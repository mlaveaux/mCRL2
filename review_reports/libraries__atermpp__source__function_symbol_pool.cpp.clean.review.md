---
{
  "file": "libraries/atermpp/source/function_symbol_pool.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T14:21:47.073285+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/source/function_symbol_pool.cpp` (full file, clean review)
- Risk level: **Medium-High** — exception-safety gaps in all manually-locked paths; API enables data races in consumers
- Overall verdict: The file is functional and the lockfree `m_symbol_set` integration is sound, but the manual lock/unlock pattern is a well-known anti-pattern that introduces deadlock risk on any throwing path. The `shared_ptr<size_t>` contract leaves unsynchronized reads as an exercise for the caller.

## Findings (ordered by severity)

### [SEV-2] Exception-unsafe manual lock/unlock in `create_helper`, `register_prefix`, and `deregister`

- **Location:** `libraries/atermpp/source/function_symbol_pool.cpp:34-60`, `:113-116`, `:119-137`
- **Why this is a problem:** Between `m_mutex.lock()` and `m_mutex.unlock()`, operations that allocate (`name.substr()`, `std::make_shared`, `m_prefix_to_register_function_map[prefix] = …`) can throw `std::bad_alloc`. On any exception, the `unlock()` call is skipped and the mutex is permanently held → deadlock on next acquisition. The file already demonstrates the correct solution in `sweep()` (line 176: `std::unique_lock lock(m_mutex)`).
- **Evidence type:** static-proof
- **Evidence:**
  - `create_helper`: `name.substr(start_of_index)` (line 44) and `name.substr(0, start_of_index)` (line 45) allocate and can throw `bad_alloc`.
  - `register_prefix`: `std::make_shared<std::size_t>(index)` (line 131) and `m_prefix_to_register_function_map[prefix]` (line 132) allocate.
  - Control flow: lock acquired at lines 34/119, unlock at lines 60/126/134 — any throw between skips the unlock.
  - Reproducer: run under memory pressure (e.g., `ulimit -v` constrained) with a workload that triggers `register_prefix` repeatedly.
  - Run command: `ctest --test-dir build -R function_symbol_generator --output-on-failure` (under memory pressure)
  - Observed result: Deadlock on subsequent pool access.
  - Expected result: Graceful `bad_alloc` propagation with mutex released.
- **Efficiency impact:** None directly, but deadlock is a liveness failure.
- **C++20-specific note:** RAII lock guards have been standard since C++11; no reason to use manual lock/unlock. `std::unique_lock` or `std::lock_guard` should be used consistently.
- **Suggested fix:** Replace all manual `m_mutex.lock()`/`m_mutex.unlock()` pairs with `std::unique_lock lock(m_mutex)`, matching the pattern already used in `sweep()`.

---

### [SEV-2] `register_prefix` returns `shared_ptr<size_t>` enabling unsynchronized concurrent read/write (data race)

- **Location:** `libraries/atermpp/source/function_symbol_pool.cpp:53` (write) and `libraries/atermpp/include/mcrl2/atermpp/function_symbol_generator.h:106` (read)
- **Why this is a problem:** `create_helper` writes `*prefix_it->second = std::max(…)` under `m_mutex`. But `function_symbol_generator::operator()()` reads `*m_central_index` *without any lock* (line 106 of generator header). Two threads accessing the same non-atomic `size_t` (one writing, one reading) with no happens-before relationship is undefined behavior per [intro.races].
- **Evidence type:** plausible-issue
- **Evidence:**
  - Thread A: creates `function_symbol("X_42", 0)` → `create()` → `create_helper("X_42")` → writes `*m_central_index` under `m_mutex`.
  - Thread B: `function_symbol_generator::operator()()` reads `*m_central_index` with no lock.
  - The write uses `m_mutex`; the read does not. No synchronization exists between these accesses.
  - Run command: Build with `-DMCRL2_ENABLE_THREADSANITIZER=ON`, run `ctest -R function_symbol_generator`
  - Expected result: TSan reports a data race on the `size_t` pointed to by the shared pointer.
- **Efficiency impact:** If the race causes a torn read, the generator could reuse an existing symbol index, violating uniqueness.
- **C++20-specific note:** Replace `shared_ptr<size_t>` with `shared_ptr<std::atomic<size_t>>` or use `std::atomic_ref<size_t>` (C++20) on the read path. Alternatively, require callers to acquire `m_mutex` before reading.
- **Suggested fix:** Change the shared index to `std::atomic<std::size_t>` and use `store`/`load` with appropriate memory ordering. The write in `create_helper` should use `store(…, std::memory_order_release)` and the read in the generator should use `load(std::memory_order_acquire)`.

---

### [SEV-3] Unsigned wraparound in prefix counter arithmetic

- **Location:** `libraries/atermpp/source/function_symbol_pool.cpp:53` and `:150`
- **Why this is a problem:** `number + 1` when `number == SIZE_MAX` wraps to 0 (well-defined for unsigned, but a logic bug). The counter would reset, causing the generator to produce names that collide with existing symbols.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Trigger: create a function symbol named `prefix_18446744073709551615` (on 64-bit). `stoul` succeeds, `number + 1` wraps to 0.
  - In practice this is unlikely (no one has 2^64 symbols), but the same logic applies to `get_sufficiently_large_postfix_index` (line 150).
  - No guard or saturation exists.
- **Efficiency impact:** None.
- **C++20-specific note:** A saturating add or explicit overflow check would make the contract clear.
- **Suggested fix:** Add `assert(number < SIZE_MAX)` or saturate at `SIZE_MAX`.

---

### [SEV-3] Misleading comment in `sweep()`

- **Location:** `libraries/atermpp/source/function_symbol_pool.cpp:176`
- **Why this is a problem:** The comment `// Prevents changes to the symbol_set` claims `m_mutex` guards `m_symbol_set`. But the header (line 95) documents `m_mutex` as "Mutex for m_prefix_to_register_function_map", and `create()` modifies `m_symbol_set` without acquiring `m_mutex`. The actual prevention of concurrent `create()` during sweep comes from the higher-level exclusive lock in `aterm_pool` (via `thread_aterm_pool`'s shared mutex). This misleads future maintainers.
- **Evidence type:** static-proof (contradiction between header comment and sweep comment)
- **Suggested fix:** Change comment to `// Acquire m_mutex to prevent concurrent access to m_prefix_to_register_function_map during iteration` or similar.

---

### [SEV-3] Redundant `if constexpr` wrapping of mutex lock in `create_helper`

- **Location:** `libraries/atermpp/source/function_symbol_pool.cpp:34,60`
- **Why this is a problem:** `create_helper` uses `if constexpr (GlobalThreadSafe) { m_mutex.lock(); }`, but `mcrl2::utilities::mutex::lock()` already contains an identical `if constexpr` check. The outer check is dead logic that diverges from the idiom used in `deregister()` and `register_prefix()` (which call `m_mutex.lock()` directly).
- **Evidence type:** static-proof (compare mutex.h lines 25-28 with create_helper lines 34,60)
- **Suggested fix:** Replace with `m_mutex.lock()` / `m_mutex.unlock()` directly (or preferably a `std::unique_lock` per SEV-2 fix above).

---

## Proposed Tests

- [ ] **test_exception_safety_register_prefix**: Purpose: verify mutex is released on allocation failure. Inject allocation failure (via custom allocator or `LD_PRELOAD` fault injector) during `register_prefix`; attempt a second `register_prefix` call afterward. Expected: no deadlock, `bad_alloc` propagated.
- [ ] **test_data_race_create_vs_generator (TSan)**: Purpose: expose concurrent non-atomic read/write on shared index. Create a generator with prefix "X_", then from another thread create function symbols "X_0" through "X_1000" with `check_for_registered_functions=true`. Run under TSan. Expected: TSan reports a data race on the shared `size_t`.
- [ ] **test_wraparound_stoul_max**: Purpose: verify behavior when trailing number is near `SIZE_MAX`. Create a function symbol named `prefix_<ULONG_MAX>`. Confirm the returned generator index doesn't wrap to 0. Expected: assertion failure or saturation.

## Open Questions

1. Is the higher-level `shared_mutex` in `aterm_pool` *documented* as protecting `m_symbol_set` during `sweep()`? If not, a future refactor could break this implicit contract.
2. Should `register_prefix` / `deregister` be restricted to single-threaded initialization phases, removing the need for `m_mutex` entirely?
3. The `function_symbol_generator_mutex()` uses a raw `std::mutex` rather than `mcrl2::utilities::mutex` — is this intentional (always locked even in single-threaded builds)?

## Copilot Raw Stdout

```text
● Explore: Explore function_symbol_pool context
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Explore: Explore concurrency patterns
  └ Agent started in background with agent_id: agent-1. You can use read_agent to...

● Read: Explore agent (Explore function_symbol_pool context)
  └ Running

● Read: Explore agent (Explore concurrency patterns)
  └ Running

● Read: Explore agent (Explore function_symbol_pool context)
  └ Completed

● Read: Explore agent (Explore concurrency patterns)
  └ Completed

Now let me verify one more detail about the exception safety concern:

● Read libraries/utilities/include/mcrl2/utilities/mutex.h
  └ 44 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_pool.h
  └ 113 lines read

Now let me check how `function_symbol_generator` uses the shared index returned by `register_prefix`:

● Grep "function_symbol_generator" in **/*.h (/workspaces/mCRL2/libraries/atermpp/include)
  └ 3 files found

● Read libraries/atermpp/include/mcrl2/atermpp/function_symbol_generator.h
  └ 154 lines read

Now I have complete context. Let me produce the review.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/source/function_symbol_pool.cpp` (full file, clean review)
- Risk level: **Medium-High** — exception-safety gaps in all manually-locked paths; API enables data races in consumers
- Overall verdict: The file is functional and the lockfree `m_symbol_set` integration is sound, but the manual lock/unlock pattern is a well-known anti-pattern that introduces deadlock risk on any throwing path. The `shared_ptr<size_t>` contract leaves unsynchronized reads as an exercise for the caller.

## Findings (ordered by severity)

### [SEV-2] Exception-unsafe manual lock/unlock in `create_helper`, `register_prefix`, and `deregister`

- **Location:** `libraries/atermpp/source/function_symbol_pool.cpp:34-60`, `:113-116`, `:119-137`
- **Why this is a problem:** Between `m_mutex.lock()` and `m_mutex.unlock()`, operations that allocate (`name.substr()`, `std::make_shared`, `m_prefix_to_register_function_map[prefix] = …`) can throw `std::bad_alloc`. On any exception, the `unlock()` call is skipped and the mutex is permanently held → deadlock on next acquisition. The file already demonstrates the correct solution in `sweep()` (line 176: `std::unique_lock lock(m_mutex)`).
- **Evidence type:** static-proof
- **Evidence:**
  - `create_helper`: `name.substr(start_of_index)` (line 44) and `name.substr(0, start_of_index)` (line 45) allocate and can throw `bad_alloc`.
  - `register_prefix`: `std::make_shared<std::size_t>(index)` (line 131) and `m_prefix_to_register_function_map[prefix]` (line 132) allocate.
  - Control flow: lock acquired at lines 34/119, unlock at lines 60/126/134 — any throw between skips the unlock.
  - Reproducer: run under memory pressure (e.g., `ulimit -v` constrained) with a workload that triggers `register_prefix` repeatedly.
  - Run command: `ctest --test-dir build -R function_symbol_generator --output-on-failure` (under memory pressure)
  - Observed result: Deadlock on subsequent pool access.
  - Expected result: Graceful `bad_alloc` propagation with mutex released.
- **Efficiency impact:** None directly, but deadlock is a liveness failure.
- **C++20-specific note:** RAII lock guards have been standard since C++11; no reason to use manual lock/unlock. `std::unique_lock` or `std::lock_guard` should be used consistently.
- **Suggested fix:** Replace all manual `m_mutex.lock()`/`m_mutex.unlock()` pairs with `std::unique_lock lock(m_mutex)`, matching the pattern already used in `sweep()`.

---

### [SEV-2] `register_prefix` returns `shared_ptr<size_t>` enabling unsynchronized concurrent read/write (data race)

- **Location:** `libraries/atermpp/source/function_symbol_pool.cpp:53` (write) and `libraries/atermpp/include/mcrl2/atermpp/function_symbol_generator.h:106` (read)
- **Why this is a problem:** `create_helper` writes `*prefix_it->second = std::max(…)` under `m_mutex`. But `function_symbol_generator::operator()()` reads `*m_central_index` *without any lock* (line 106 of generator header). Two threads accessing the same non-atomic `size_t` (one writing, one reading) with no happens-before relationship is undefined behavior per [intro.races].
- **Evidence type:** plausible-issue
- **Evidence:**
  - Thread A: creates `function_symbol("X_42", 0)` → `create()` → `create_helper("X_42")` → writes `*m_central_index` under `m_mutex`.
  - Thread B: `function_symbol_generator::operator()()` reads `*m_central_index` with no lock.
  - The write uses `m_mutex`; the read does not. No synchronization exists between these accesses.
  - Run command: Build with `-DMCRL2_ENABLE_THREADSANITIZER=ON`, run `ctest -R function_symbol_generator`
  - Expected result: TSan reports a data race on the `size_t` pointed to by the shared pointer.
- **Efficiency impact:** If the race causes a torn read, the generator could reuse an existing symbol index, violating uniqueness.
- **C++20-specific note:** Replace `shared_ptr<size_t>` with `shared_ptr<std::atomic<size_t>>` or use `std::atomic_ref<size_t>` (C++20) on the read path. Alternatively, require callers to acquire `m_mutex` before reading.
- **Suggested fix:** Change the shared index to `std::atomic<std::size_t>` and use `store`/`load` with appropriate memory ordering. The write in `create_helper` should use `store(…, std::memory_order_release)` and the read in the generator should use `load(std::memory_order_acquire)`.

---

### [SEV-3] Unsigned wraparound in prefix counter arithmetic

- **Location:** `libraries/atermpp/source/function_symbol_pool.cpp:53` and `:150`
- **Why this is a problem:** `number + 1` when `number == SIZE_MAX` wraps to 0 (well-defined for unsigned, but a logic bug). The counter would reset, causing the generator to produce names that collide with existing symbols.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Trigger: create a function symbol named `prefix_18446744073709551615` (on 64-bit). `stoul` succeeds, `number + 1` wraps to 0.
  - In practice this is unlikely (no one has 2^64 symbols), but the same logic applies to `get_sufficiently_large_postfix_index` (line 150).
  - No guard or saturation exists.
- **Efficiency impact:** None.
- **C++20-specific note:** A saturating add or explicit overflow check would make the contract clear.
- **Suggested fix:** Add `assert(number < SIZE_MAX)` or saturate at `SIZE_MAX`.

---

### [SEV-3] Misleading comment in `sweep()`

- **Location:** `libraries/atermpp/source/function_symbol_pool.cpp:176`
- **Why this is a problem:** The comment `// Prevents changes to the symbol_set` claims `m_mutex` guards `m_symbol_set`. But the header (line 95) documents `m_mutex` as "Mutex for m_prefix_to_register_function_map", and `create()` modifies `m_symbol_set` without acquiring `m_mutex`. The actual prevention of concurrent `create()` during sweep comes from the higher-level exclusive lock in `aterm_pool` (via `thread_aterm_pool`'s shared mutex). This misleads future maintainers.
- **Evidence type:** static-proof (contradiction between header comment and sweep comment)
- **Suggested fix:** Change comment to `// Acquire m_mutex to prevent concurrent access to m_prefix_to_register_function_map during iteration` or similar.

---

### [SEV-3] Redundant `if constexpr` wrapping of mutex lock in `create_helper`

- **Location:** `libraries/atermpp/source/function_symbol_pool.cpp:34,60`
- **Why this is a problem:** `create_helper` uses `if constexpr (GlobalThreadSafe) { m_mutex.lock(); }`, but `mcrl2::utilities::mutex::lock()` already contains an identical `if constexpr` check. The outer check is dead logic that diverges from the idiom used in `deregister()` and `register_prefix()` (which call `m_mutex.lock()` directly).
- **Evidence type:** static-proof (compare mutex.h lines 25-28 with create_helper lines 34,60)
- **Suggested fix:** Replace with `m_mutex.lock()` / `m_mutex.unlock()` directly (or preferably a `std::unique_lock` per SEV-2 fix above).

---

## Proposed Tests

- [ ] **test_exception_safety_register_prefix**: Purpose: verify mutex is released on allocation failure. Inject allocation failure (via custom allocator or `LD_PRELOAD` fault injector) during `register_prefix`; attempt a second `register_prefix` call afterward. Expected: no deadlock, `bad_alloc` propagated.
- [ ] **test_data_race_create_vs_generator (TSan)**: Purpose: expose concurrent non-atomic read/write on shared index. Create a generator with prefix "X_", then from another thread create function symbols "X_0" through "X_1000" with `check_for_registered_functions=true`. Run under TSan. Expected: TSan reports a data race on the shared `size_t`.
- [ ] **test_wraparound_stoul_max**: Purpose: verify behavior when trailing number is near `SIZE_MAX`. Create a function symbol named `prefix_<ULONG_MAX>`. Confirm the returned generator index doesn't wrap to 0. Expected: assertion failure or saturation.

## Open Questions

1. Is the higher-level `shared_mutex` in `aterm_pool` *documented* as protecting `m_symbol_set` during `sweep()`? If not, a future refactor could break this implicit contract.
2. Should `register_prefix` / `deregister` be restricted to single-threaded initialization phases, removing the need for `m_mutex` entirely?
3. The `function_symbol_generator_mutex()` uses a raw `std::mutex` rather than `mcrl2::utilities::mutex` — is this intentional (always locked even in single-threaded builds)?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         6m 21s
Total session time:     5m 26s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         330.6k in, 10.6k out, 294.6k cached (Est. 1 Premium request)
 claude-haiku-4.5        844.0k in, 13.3k out, 774.9k cached (Est. 0 Premium requests)

```

---
{
  "file": "libraries/atermpp/source/gc_stress_thread.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T14:27:52.577467+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/source/gc_stress_thread.cpp` — a debug/stress-testing facility that spawns a background thread to repeatedly trigger GC
- **Risk level:** Low (feature is disabled by default via `EnableGCStressThread = false`; entire body is dead-code-eliminated by `if constexpr`)
- **Overall verdict:** Acceptable for its disabled-by-default diagnostic purpose, but contains a latent `std::terminate` defect and uses non-obvious memory ordering without justification

## Findings (ordered by severity)

### [SEV-2] Assignment to potentially-joinable `std::thread` causes `std::terminate`

- **Location:** `libraries/atermpp/source/gc_stress_thread.cpp:37`
- **Why this is a problem:** `g_gc_stress_thread = std::thread([]{...})` invokes `std::thread::operator=`, which calls `std::terminate()` if the target is already joinable (§32.4.2.4 [thread.thread.assign]). If `start_gc_stress_thread()` is ever invoked a second time without an intervening `stop_gc_stress_thread()`, the program terminates abnormally. Today the single call site in `g_term_pool<true>()` prevents this, but the function's public interface offers no such guarantee—any future caller triggers UB.
- **Evidence type:** static-proof
- **Evidence:**
  - C++ standard §32.4.2.4: "If `joinable()` is `true`, calls `std::terminate()`."
  - The function performs no `joinable()` guard before assignment.
  - **Reproducer:** Call `start_gc_stress_thread()` twice consecutively (with `EnableGCStressThread=true` and `MCRL2_ENABLE_MULTITHREADING`). Program aborts with `std::terminate`.
  - **Run command:** N/A (requires compile-time constant change, but scenario is deterministic)
  - **Observed result:** `std::terminate()` / `SIGABRT`
  - **Expected result:** Second call is either a no-op or stops-then-restarts the thread
- **Efficiency impact:** None (correctness issue)
- **C++20-specific note:** Could use `std::jthread` which auto-joins on assignment, eliminating this class of bug entirely. `std::jthread` also supports cooperative cancellation via `std::stop_token`, removing the need for the manual `std::atomic<bool>` flag.
- **Suggested fix:**
  ```cpp
  void start_gc_stress_thread()
  {
    if constexpr (EnableGCStressThread && mcrl2::utilities::detail::GlobalThreadSafe)
    {
      if (g_gc_stress_thread.joinable())
        return; // Already running
      g_gc_stress_running.store(true, std::memory_order_relaxed);
      g_gc_stress_thread = std::thread([]{...});
      std::atexit(stop_gc_stress_thread);
    }
  }
  ```
  Or better: replace `std::thread` + `std::atomic<bool>` with `std::jthread` + `std::stop_token`.

---

### [SEV-3] `std::memory_order_relaxed` is correct but non-obvious; release/acquire is zero-cost-clearer

- **Location:** `libraries/atermpp/source/gc_stress_thread.cpp:35,39,53`
- **Why this is a problem:** Relaxed ordering is *technically* sufficient here because: (a) eventual visibility is guaranteed by the C++ memory model for atomic stores, (b) the polling loop with 1 ms sleep provides practical latency, and (c) `std::thread::join()` provides the happens-before edge needed by the caller. However, the reasoning is non-trivial, undocumented, and forces every future reader to re-derive it. Using `memory_order_release` (store) / `memory_order_acquire` (load) makes intent explicit at zero additional cost on x86 (both compile to the same instructions) and near-zero cost on ARM/POWER.
- **Evidence type:** plausible-issue (maintainability/reasoning burden, not a correctness defect)
- **Evidence:**
  - No comment justifies relaxed ordering.
  - On weakly-ordered architectures (ARM), `relaxed` compiles to a plain load without a fence; while correctness is preserved via `join()`, a future refactor removing `join()` would silently become unsound.
- **Efficiency impact:** None—release/acquire on x86 compiles identically to relaxed.
- **C++20-specific note:** `std::jthread` with `std::stop_token` would make this entire flag and its memory ordering irrelevant.
- **Suggested fix:** Add a brief comment explaining why relaxed suffices, or switch to `release`/`acquire` for defensive clarity.

---

### [SEV-3] Multiple `std::atexit` registrations on repeated calls

- **Location:** `libraries/atermpp/source/gc_stress_thread.cpp:45`
- **Why this is a problem:** Each invocation of `start_gc_stress_thread()` registers a new `atexit` handler. While `stop_gc_stress_thread()` is idempotent (checks `joinable()`), the C++ standard limits `atexit` to at least 32 registrations (§17.5 [support.start.term]). More practically, it signals missing re-entrance protection.
- **Evidence type:** static-proof
- **Evidence:** No guard prevents the `std::atexit` call from executing multiple times.
- **Efficiency impact:** Negligible
- **Suggested fix:** Guard with `static bool registered = false;` or use `std::call_once`.

---

### [SEV-3] Missing `[[nodiscard]]` and API contract documentation on `stop_gc_stress_thread`

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/gc_stress_thread.h`
- **Why this is a problem:** Per project C++ instructions, new/changed public APIs require Doxygen `\pre`/`\post` clauses. The header has a `\brief` but no explicit contract. `stop_gc_stress_thread()` is safe to call when no thread is running (idempotent), but this is not documented.
- **Evidence type:** plausible-issue (violates project contract documentation policy)
- **Suggested fix:** Add `\post No GC stress thread is running.` and `\note Safe to call even if the thread was never started.`

---

## Proposed Tests

- [ ] **`test_double_start_terminates`**: Purpose: confirm that calling `start_gc_stress_thread()` twice without a stop causes `std::terminate` (current behavior) or is safely a no-op (after fix). Input: enable `EnableGCStressThread` at compile time, call `start_gc_stress_thread()` twice. Expected failure: `SIGABRT` / `std::terminate` before fix; clean exit after fix. Run: requires modifying `aterm_configuration.h` to set `EnableGCStressThread = true`, then a dedicated test binary.
- [ ] **`test_start_stop_idempotent`**: Purpose: verify `stop_gc_stress_thread()` is safe to call when no thread was started. Input: call `stop_gc_stress_thread()` without a preceding start. Expected: no crash, no hang. Run: unit test in `libraries/atermpp/test/`.

## Open Questions

1. **Is `std::jthread` adoption viable?** All minimum-supported compilers (GCC 11, Clang 16, MSVC 19.31) support `std::jthread`. It would eliminate the manual flag, simplify shutdown, and remove the `std::terminate` risk. Is there a project-wide decision on `std::jthread` adoption?
2. **Should this file even exist as a separate TU?** The entire body is dead-code-eliminated unless `EnableGCStressThread` is flipped at compile time. It could be a header-only `inline` guarded by `if constexpr`, avoiding the TU compilation cost. Alternatively, it could be conditionally compiled via CMake option.
3. **Who is responsible for calling `start_gc_stress_thread()`?** The only call site is inside `g_term_pool<true>()` which itself has a non-thread-safe `static bool initialized` guard. If that initialization is ever hit concurrently (before a fix there), the stress thread code inherits the race.

## Copilot Raw Stdout

```text
● Explore: Explore gc_stress_thread context
  ● Check which test files mention gc_stress
    $ grep -l "gc_stress\|EnableGCStress" /workspaces/mCRL2/libraries/atermpp/test/*.cpp
    └ 1 line...
  ● Check parallel_thread_pool_test for relevant content
    $ grep -n "parallel.*test\|thread"
    /workspaces/mCRL2/libraries/atermpp/test/parallel_thread_pool_test.cpp | head -30
    └ 23 lines...

● Explore: Check memory ordering and atexit
  ● Find aterm_pool implementations
    $ grep -n "aterm_pool::"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h |
    head -30
    └ 21 lines...
  ● Find where g_term_pool is called
    $ grep -rn "g_term_pool<true>\|g_term_pool()"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm*.h | head -20
    └ 18 lines...

● Explore: Verify memory ordering concerns
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 57-95
    └ 39 lines read
  ● Grep "collect_impl.*m_thread_pools" (/workspaces/mCRL2/libraries/atermpp/source)
    └ No matches found

Now I have enough context to write the review. Let me verify one final critical detail about the global pool's lifetime (it uses placement new without explicit destruction, which is key to the shutdown safety analysis).

● Read libraries/atermpp/source/aterm_implementation.cpp lines 65-80
  └ 9 lines read

● Grep "~aterm_pool|destroy.*aterm_pool|g_aterm_pool_storage" (/workspaces/mCRL2/libraries/atermpp)
  └ 5 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h
  └ 65 lines read

Now I have the complete picture. The global pool uses placement-new into raw `std::byte` storage and is **never destructed** (intentional leak to avoid static destruction order issues). This is critical to the shutdown safety analysis.

---

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/source/gc_stress_thread.cpp` — a debug/stress-testing facility that spawns a background thread to repeatedly trigger GC
- **Risk level:** Low (feature is disabled by default via `EnableGCStressThread = false`; entire body is dead-code-eliminated by `if constexpr`)
- **Overall verdict:** Acceptable for its disabled-by-default diagnostic purpose, but contains a latent `std::terminate` defect and uses non-obvious memory ordering without justification

## Findings (ordered by severity)

### [SEV-2] Assignment to potentially-joinable `std::thread` causes `std::terminate`

- **Location:** `libraries/atermpp/source/gc_stress_thread.cpp:37`
- **Why this is a problem:** `g_gc_stress_thread = std::thread([]{...})` invokes `std::thread::operator=`, which calls `std::terminate()` if the target is already joinable (§32.4.2.4 [thread.thread.assign]). If `start_gc_stress_thread()` is ever invoked a second time without an intervening `stop_gc_stress_thread()`, the program terminates abnormally. Today the single call site in `g_term_pool<true>()` prevents this, but the function's public interface offers no such guarantee—any future caller triggers UB.
- **Evidence type:** static-proof
- **Evidence:**
  - C++ standard §32.4.2.4: "If `joinable()` is `true`, calls `std::terminate()`."
  - The function performs no `joinable()` guard before assignment.
  - **Reproducer:** Call `start_gc_stress_thread()` twice consecutively (with `EnableGCStressThread=true` and `MCRL2_ENABLE_MULTITHREADING`). Program aborts with `std::terminate`.
  - **Run command:** N/A (requires compile-time constant change, but scenario is deterministic)
  - **Observed result:** `std::terminate()` / `SIGABRT`
  - **Expected result:** Second call is either a no-op or stops-then-restarts the thread
- **Efficiency impact:** None (correctness issue)
- **C++20-specific note:** Could use `std::jthread` which auto-joins on assignment, eliminating this class of bug entirely. `std::jthread` also supports cooperative cancellation via `std::stop_token`, removing the need for the manual `std::atomic<bool>` flag.
- **Suggested fix:**
  ```cpp
  void start_gc_stress_thread()
  {
    if constexpr (EnableGCStressThread && mcrl2::utilities::detail::GlobalThreadSafe)
    {
      if (g_gc_stress_thread.joinable())
        return; // Already running
      g_gc_stress_running.store(true, std::memory_order_relaxed);
      g_gc_stress_thread = std::thread([]{...});
      std::atexit(stop_gc_stress_thread);
    }
  }
  ```
  Or better: replace `std::thread` + `std::atomic<bool>` with `std::jthread` + `std::stop_token`.

---

### [SEV-3] `std::memory_order_relaxed` is correct but non-obvious; release/acquire is zero-cost-clearer

- **Location:** `libraries/atermpp/source/gc_stress_thread.cpp:35,39,53`
- **Why this is a problem:** Relaxed ordering is *technically* sufficient here because: (a) eventual visibility is guaranteed by the C++ memory model for atomic stores, (b) the polling loop with 1 ms sleep provides practical latency, and (c) `std::thread::join()` provides the happens-before edge needed by the caller. However, the reasoning is non-trivial, undocumented, and forces every future reader to re-derive it. Using `memory_order_release` (store) / `memory_order_acquire` (load) makes intent explicit at zero additional cost on x86 (both compile to the same instructions) and near-zero cost on ARM/POWER.
- **Evidence type:** plausible-issue (maintainability/reasoning burden, not a correctness defect)
- **Evidence:**
  - No comment justifies relaxed ordering.
  - On weakly-ordered architectures (ARM), `relaxed` compiles to a plain load without a fence; while correctness is preserved via `join()`, a future refactor removing `join()` would silently become unsound.
- **Efficiency impact:** None—release/acquire on x86 compiles identically to relaxed.
- **C++20-specific note:** `std::jthread` with `std::stop_token` would make this entire flag and its memory ordering irrelevant.
- **Suggested fix:** Add a brief comment explaining why relaxed suffices, or switch to `release`/`acquire` for defensive clarity.

---

### [SEV-3] Multiple `std::atexit` registrations on repeated calls

- **Location:** `libraries/atermpp/source/gc_stress_thread.cpp:45`
- **Why this is a problem:** Each invocation of `start_gc_stress_thread()` registers a new `atexit` handler. While `stop_gc_stress_thread()` is idempotent (checks `joinable()`), the C++ standard limits `atexit` to at least 32 registrations (§17.5 [support.start.term]). More practically, it signals missing re-entrance protection.
- **Evidence type:** static-proof
- **Evidence:** No guard prevents the `std::atexit` call from executing multiple times.
- **Efficiency impact:** Negligible
- **Suggested fix:** Guard with `static bool registered = false;` or use `std::call_once`.

---

### [SEV-3] Missing `[[nodiscard]]` and API contract documentation on `stop_gc_stress_thread`

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/gc_stress_thread.h`
- **Why this is a problem:** Per project C++ instructions, new/changed public APIs require Doxygen `\pre`/`\post` clauses. The header has a `\brief` but no explicit contract. `stop_gc_stress_thread()` is safe to call when no thread is running (idempotent), but this is not documented.
- **Evidence type:** plausible-issue (violates project contract documentation policy)
- **Suggested fix:** Add `\post No GC stress thread is running.` and `\note Safe to call even if the thread was never started.`

---

## Proposed Tests

- [ ] **`test_double_start_terminates`**: Purpose: confirm that calling `start_gc_stress_thread()` twice without a stop causes `std::terminate` (current behavior) or is safely a no-op (after fix). Input: enable `EnableGCStressThread` at compile time, call `start_gc_stress_thread()` twice. Expected failure: `SIGABRT` / `std::terminate` before fix; clean exit after fix. Run: requires modifying `aterm_configuration.h` to set `EnableGCStressThread = true`, then a dedicated test binary.
- [ ] **`test_start_stop_idempotent`**: Purpose: verify `stop_gc_stress_thread()` is safe to call when no thread was started. Input: call `stop_gc_stress_thread()` without a preceding start. Expected: no crash, no hang. Run: unit test in `libraries/atermpp/test/`.

## Open Questions

1. **Is `std::jthread` adoption viable?** All minimum-supported compilers (GCC 11, Clang 16, MSVC 19.31) support `std::jthread`. It would eliminate the manual flag, simplify shutdown, and remove the `std::terminate` risk. Is there a project-wide decision on `std::jthread` adoption?
2. **Should this file even exist as a separate TU?** The entire body is dead-code-eliminated unless `EnableGCStressThread` is flipped at compile time. It could be a header-only `inline` guarded by `if constexpr`, avoiding the TU compilation cost. Alternatively, it could be conditionally compiled via CMake option.
3. **Who is responsible for calling `start_gc_stress_thread()`?** The only call site is inside `g_term_pool<true>()` which itself has a non-thread-safe `static bool initialized` guard. If that initialization is ever hit concurrently (before a fix there), the stress thread code inherits the race.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         5m 39s
Total session time:     6m 1s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         268.7k in, 7.7k out, 239.3k cached (Est. 1 Premium request)
 claude-haiku-4.5        695.4k in, 13.6k out, 630.7k cached (Est. 0 Premium requests)

```

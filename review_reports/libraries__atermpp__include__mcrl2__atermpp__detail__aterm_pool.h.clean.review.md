---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T09:48:54.848875+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h` — the global aterm pool header declaring `aterm_pool` and `thread_aterm_pool_interface`.
- **Risk level:** Medium — contains a strict-aliasing violation acknowledged by a TODO, potential signed-overflow UB, and a data race on a diagnostic path.
- **Overall verdict:** The design is sound and the concurrency model for the hot path (term creation/GC) is correctly serialized. However, there are several formally-undefined behaviors and missing modern C++ safety annotations that should be addressed.

## Findings (ordered by severity)

### [SEV-1] `reinterpret_cast<aterm&>(m_empty_list)` is undefined behavior (strict aliasing)

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h:118` and `aterm_pool_implementation.h:43`
- **Why this is a problem:** `m_empty_list` is declared as `aterm_core` (line 227). `aterm` is *derived* from `aterm_core`. Accessing an object of type `aterm_core` through a glvalue of type `aterm&` violates [basic.lval] — you cannot treat a base-class object as a derived-class object. While the memory layout is identical (both are a single pointer), this is formally UB and compilers are permitted to optimize based on the assumption it does not happen. Both the member access at line 118 and the initialization at `aterm_pool_implementation.h:43` trigger this.
- **Evidence type:** static-proof
- **Evidence:**
  - Per C++20 [basic.lval]/11, accessing the stored value of an object through a glvalue requires the type to be the dynamic type, a cv-qualified version thereof, or a base class type. `aterm` is a *derived* type of `aterm_core`, not a base, so the access is UB.
  - The TODO comment on line 118 acknowledges this issue explicitly.
  - Reproducer: Compile with `-fstrict-aliasing -Wstrict-aliasing=1` (GCC) or run under UBSan with `-fsanitize=undefined` and use an optimization level that enables type-based alias analysis (`-O2`).
  - Run command: `cmake --build build-asan -j && ctest --test-dir build-asan -R atermpp --output-on-failure`
  - Expected result: No UBSan violation (currently this is unlikely to trigger UBSan as it targets type punning via memcpy/unions, not inheritance casts — but TBAA optimizations could miscompile).
- **Efficiency impact:** None (same machine code either way).
- **C++20-specific note:** The fix is trivial: change `aterm_core m_empty_list;` to `aterm m_empty_list;` — since `aterm`'s default constructor just calls `aterm_core()`, there is no semantic difference and the cast can be removed entirely.
- **Suggested fix:**
  ```cpp
  // Declaration (line 227):
  aterm m_empty_list;
  
  // Accessor (line 118):
  aterm& empty_list() noexcept { return m_empty_list; }
  
  // Initialization (aterm_pool_implementation.h:43):
  create_appl(m_empty_list, m_function_symbol_pool.as_empty_list());
  ```

---

### [SEV-2] Signed integer overflow in GC threshold calculation

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h:276`
- **Why this is a problem:** `size()` and `protection_set_size()` both return `std::size_t`. Their sum is `std::size_t`, which is then cast to `long` via `static_cast<long>(size() + protection_set_size())`. If the total exceeds `LONG_MAX` (~2.1 billion on LP64, ~9.2×10¹⁸ on LLP64/ILP64), this is implementation-defined behavior (C++20 [conv.integral]/3 for unsigned-to-signed narrowing). More critically, the `fetch_sub` on line 178 subtracts a `long` (`count_until_check`) from `m_count_until_collection` — if the stored value is close to `LONG_MIN`, the subtraction can overflow `long`, which *is* UB for signed types stored in an `std::atomic<long>`.
- **Evidence type:** plausible-issue
- **Evidence:**
  - For models with >2 billion live terms (feasible in large state-space explorations with many workers), `size() + protection_set_size()` can exceed `LONG_MAX` on 64-bit Linux (where `long` is 64-bit, this requires ~9.2×10¹⁸ — unlikely). On Windows/LLP64 where `long` is 32-bit, >2 billion terms is plausible.
  - **Mechanism:** On platforms where `sizeof(long) == 4` (Windows), a pool of ~2.5 billion terms causes `static_cast<long>` to produce a negative value, immediately triggering GC on the next `created_term` call, creating a GC storm.
  - **Test plan:** Build on Windows (MSVC) or cross-compile for ILP32, create a test allocating >2^31 integer terms, observe whether GC threshold becomes negative.
- **Efficiency impact:** On affected platforms, perpetual GC storms would cripple performance.
- **C++20-specific note:** Consider `std::atomic<std::ptrdiff_t>` or `std::atomic<std::int64_t>` to guarantee sufficient range on all platforms.
- **Suggested fix:** Change `std::atomic<long> m_count_until_collection` to `std::atomic<std::int64_t>` and update all related `long` types to `std::int64_t`.

---

### [SEV-2] Data race in `print_performance_statistics()` and `protection_set_size()`

- **Location:** `aterm_pool_implementation.h:129` and `aterm_pool_implementation.h:445`
- **Why this is a problem:** Both functions iterate `m_thread_pools` without holding any lock. They are called from `collect_impl` (which holds the exclusive lock — safe) but are also reachable via the public `print_performance_statistics()` and private `protection_set_size()`. If `print_performance_statistics()` is ever called from outside the GC path (e.g., a diagnostic signal handler or on-demand stats dump) while another thread calls `register_thread_aterm_pool` or `remove_thread_aterm_pool`, this is a data race on the vector.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `print_performance_statistics()` is declared `public` (aterm_pool.h:115).
  - `register_thread_aterm_pool`/`remove_thread_aterm_pool` mutate the vector under exclusive lock.
  - Concurrent read (iteration in `print_performance_statistics`) + write (push/erase in register/remove) on `std::vector` is a data race.
  - **Test plan:** Add a test that spawns threads registering/unregistering while the main thread calls `print_performance_statistics()` in a loop. Run under TSan (`-DMCRL2_ENABLE_THREADSANITIZER=ON`).
  - Run command: `cmake --build build-tsan -j && ctest --test-dir build-tsan -R <new_test_name> --output-on-failure`
- **Efficiency impact:** Diagnostic-only path; unlikely to affect normal performance.
- **Suggested fix:** Acquire a shared lock at the start of `print_performance_statistics()`:
  ```cpp
  void aterm_pool::print_performance_statistics() const
  {
    // Ensure thread-safe iteration of m_thread_pools.
    mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
    // ... existing iteration ...
  }
  ```
  Similarly for `protection_set_size()` if it can be reached outside `collect_impl`.

---

### [SEV-3] `std::function` type-erasure overhead in `thread_aterm_pool_interface`

- **Location:** `aterm_pool.h:66-68`
- **Why this is a problem:** Three `std::function` members store trivial forwarding lambdas (each captures only a `this` pointer). `std::function` carries type-erasure overhead: typically 32–48 bytes per instance (on libstdc++/libc++) plus an indirect call. During GC, `mark()` is called for every registered thread — the indirect call prevents inlining of what is ultimately a simple method forwarding. The circular-dependency comment motivating this design could be resolved with a forward-declared virtual interface or CRTP instead.
- **Evidence type:** static-proof
- **Evidence:**
  - `thread_aterm_pool_interface` is constructed in `thread_aterm_pool.h` with `[this] { mark(); }` lambdas.
  - These are called in the GC hot path (`collect_impl` line 219: `pool->mark()`).
  - `sizeof(std::function<void()>)` is 32 bytes on libstdc++ (GCC 11+), times 3 = 96 bytes of overhead per thread pool interface, plus virtual dispatch cost per call.
- **Efficiency impact:** Low-to-moderate. GC frequency depends on workload, but `mark()` is on the critical GC path. For workloads with many threads and frequent GC, the indirect call overhead accumulates.
- **C++20-specific note:** A cleaner C++20 approach would be to make `thread_aterm_pool_interface` a virtual base class (with `mark()`, `print_local_performance_statistics()`, `protection_set_size()` as pure virtual methods) and have `thread_aterm_pool` inherit from it. This replaces 3 function objects with a single vtable pointer.
- **Suggested fix:** Refactor to virtual interface or store a `thread_aterm_pool*` directly (resolving the circular dependency with a forward declaration + out-of-line definitions).

---

### [SEV-3] Missing `[[nodiscard]]` on value-returning public methods

- **Location:** `aterm_pool.h:109-128` (multiple methods)
- **Why this is a problem:** `capacity()`, `size()`, `empty_list()`, `as_int()`, `as_list()`, `as_empty_list()`, and `get_symbol_pool()` return values that callers must use. Discarding them is always a bug. Without `[[nodiscard]]`, compilers won't warn on accidental discard.
- **Evidence type:** static-proof
- **Evidence:** The C++ Core Guidelines (F.21) and mCRL2's own `.github/instructions/cpp.instructions.md` both require `[[nodiscard]]` on pure queries and must-use results.
- **Suggested fix:** Add `[[nodiscard]]` to all accessors returning non-void.

---

### [SEV-3] Missing Doxygen `\pre`/`\post` contracts on public API

- **Location:** `aterm_pool.h:88-128` (public section)
- **Why this is a problem:** Per the mCRL2 coding instructions, new/changed public APIs require explicit precondition/postcondition documentation. `create_function_symbol` has no documented preconditions (e.g., is empty `name` valid? is arity 0 valid?). `register_thread_aterm_pool` documents "threadsafe" but not pre/postconditions. `empty_list()` has no `\post` specifying the returned term is defined.
- **Evidence type:** static-proof (guideline violation)
- **Suggested fix:** Add `\pre` and `\post` clauses to all public methods.

---

### [SEV-3] `m_shared_mutex` declared mutable would allow `const`-correct locking

- **Location:** `aterm_pool.h:224`
- **Why this is a problem:** `print_performance_statistics()` is `const`-qualified but cannot acquire a lock on `m_shared_mutex` without `mutable`. If the SEV-2 fix is applied (adding shared lock to `print_performance_statistics`), `m_shared_mutex` must be declared `mutable` to allow locking in const methods. This is the standard pattern for synchronization primitives.
- **Evidence type:** static-proof
- **Suggested fix:** `mutable mcrl2::utilities::shared_mutex m_shared_mutex;`

## Proposed Tests

- [ ] `test_empty_list_type_safety`: Verify `empty_list()` returns a valid `aterm` with correct function symbol `as_empty_list()`, arity 0, and is `defined()`. Cast correctness test. Run existing atermpp unit tests under UBSan.
- [ ] `test_gc_threshold_overflow_32bit` (plausible): On a 32-bit `long` platform, create >2^31 terms and verify GC threshold remains positive and GC does not storm. Input: loop creating integer terms. Expected failure: GC triggered every term creation after threshold wraps negative.
- [ ] `test_print_stats_concurrent` (plausible): Spawn N threads that register/unregister `thread_aterm_pool_interface` objects while main thread calls `print_performance_statistics()` in a loop. Run under TSan. Expected failure: TSan reports data race on `m_thread_pools`.

## Open Questions

1. Is `print_performance_statistics()` ever called outside the GC exclusive lock in practice? If so, the SEV-2 data race is confirmed rather than plausible.
2. What platforms does mCRL2 target where `sizeof(long) == 4`? If Windows/MSVC is actively supported, the SEV-2 overflow becomes a confirmed issue for large models.
3. The circular dependency between `aterm_pool` and `thread_aterm_pool` motivates the `std::function` design — has a virtual interface approach been considered and rejected for a specific reason (e.g., preventing header inclusion loops)?
4. The TODO on line 118 suggests changing `m_empty_list` to `aterm` — is there a known obstacle preventing this fix?

## Copilot Raw Stdout

```text
● Explore: Explore aterm_pool context
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Explore: Explore thread safety details
  └ Agent started in background with agent_id: agent-1. You can use read_agent to...

● Explore: Explore aterm type hierarchy
  └ Agent started in background with agent_id: agent-2. You can use read_agent to...

● Read: Explore agent (Explore aterm_pool context)
  └ Completed

● Read: Explore agent (Explore thread safety details)
  └ Completed

● Read: Explore agent (Explore aterm type hierarchy)
  └ Completed

Now let me verify one specific detail about the `print_performance_statistics` locking:

● Grep "print_performance_statistics"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ 4 lines found

● Grep "protection_set_size"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h)
  └ 21 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h lines 114-165
  └ 52 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h lines 25-50
  └ 26 lines read

Now I have all the context needed. Let me produce the review.

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h` — the global aterm pool header declaring `aterm_pool` and `thread_aterm_pool_interface`.
- **Risk level:** Medium — contains a strict-aliasing violation acknowledged by a TODO, potential signed-overflow UB, and a data race on a diagnostic path.
- **Overall verdict:** The design is sound and the concurrency model for the hot path (term creation/GC) is correctly serialized. However, there are several formally-undefined behaviors and missing modern C++ safety annotations that should be addressed.

## Findings (ordered by severity)

### [SEV-1] `reinterpret_cast<aterm&>(m_empty_list)` is undefined behavior (strict aliasing)

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h:118` and `aterm_pool_implementation.h:43`
- **Why this is a problem:** `m_empty_list` is declared as `aterm_core` (line 227). `aterm` is *derived* from `aterm_core`. Accessing an object of type `aterm_core` through a glvalue of type `aterm&` violates [basic.lval] — you cannot treat a base-class object as a derived-class object. While the memory layout is identical (both are a single pointer), this is formally UB and compilers are permitted to optimize based on the assumption it does not happen. Both the member access at line 118 and the initialization at `aterm_pool_implementation.h:43` trigger this.
- **Evidence type:** static-proof
- **Evidence:**
  - Per C++20 [basic.lval]/11, accessing the stored value of an object through a glvalue requires the type to be the dynamic type, a cv-qualified version thereof, or a base class type. `aterm` is a *derived* type of `aterm_core`, not a base, so the access is UB.
  - The TODO comment on line 118 acknowledges this issue explicitly.
  - Reproducer: Compile with `-fstrict-aliasing -Wstrict-aliasing=1` (GCC) or run under UBSan with `-fsanitize=undefined` and use an optimization level that enables type-based alias analysis (`-O2`).
  - Run command: `cmake --build build-asan -j && ctest --test-dir build-asan -R atermpp --output-on-failure`
  - Expected result: No UBSan violation (currently this is unlikely to trigger UBSan as it targets type punning via memcpy/unions, not inheritance casts — but TBAA optimizations could miscompile).
- **Efficiency impact:** None (same machine code either way).
- **C++20-specific note:** The fix is trivial: change `aterm_core m_empty_list;` to `aterm m_empty_list;` — since `aterm`'s default constructor just calls `aterm_core()`, there is no semantic difference and the cast can be removed entirely.
- **Suggested fix:**
  ```cpp
  // Declaration (line 227):
  aterm m_empty_list;
  
  // Accessor (line 118):
  aterm& empty_list() noexcept { return m_empty_list; }
  
  // Initialization (aterm_pool_implementation.h:43):
  create_appl(m_empty_list, m_function_symbol_pool.as_empty_list());
  ```

---

### [SEV-2] Signed integer overflow in GC threshold calculation

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h:276`
- **Why this is a problem:** `size()` and `protection_set_size()` both return `std::size_t`. Their sum is `std::size_t`, which is then cast to `long` via `static_cast<long>(size() + protection_set_size())`. If the total exceeds `LONG_MAX` (~2.1 billion on LP64, ~9.2×10¹⁸ on LLP64/ILP64), this is implementation-defined behavior (C++20 [conv.integral]/3 for unsigned-to-signed narrowing). More critically, the `fetch_sub` on line 178 subtracts a `long` (`count_until_check`) from `m_count_until_collection` — if the stored value is close to `LONG_MIN`, the subtraction can overflow `long`, which *is* UB for signed types stored in an `std::atomic<long>`.
- **Evidence type:** plausible-issue
- **Evidence:**
  - For models with >2 billion live terms (feasible in large state-space explorations with many workers), `size() + protection_set_size()` can exceed `LONG_MAX` on 64-bit Linux (where `long` is 64-bit, this requires ~9.2×10¹⁸ — unlikely). On Windows/LLP64 where `long` is 32-bit, >2 billion terms is plausible.
  - **Mechanism:** On platforms where `sizeof(long) == 4` (Windows), a pool of ~2.5 billion terms causes `static_cast<long>` to produce a negative value, immediately triggering GC on the next `created_term` call, creating a GC storm.
  - **Test plan:** Build on Windows (MSVC) or cross-compile for ILP32, create a test allocating >2^31 integer terms, observe whether GC threshold becomes negative.
- **Efficiency impact:** On affected platforms, perpetual GC storms would cripple performance.
- **C++20-specific note:** Consider `std::atomic<std::ptrdiff_t>` or `std::atomic<std::int64_t>` to guarantee sufficient range on all platforms.
- **Suggested fix:** Change `std::atomic<long> m_count_until_collection` to `std::atomic<std::int64_t>` and update all related `long` types to `std::int64_t`.

---

### [SEV-2] Data race in `print_performance_statistics()` and `protection_set_size()`

- **Location:** `aterm_pool_implementation.h:129` and `aterm_pool_implementation.h:445`
- **Why this is a problem:** Both functions iterate `m_thread_pools` without holding any lock. They are called from `collect_impl` (which holds the exclusive lock — safe) but are also reachable via the public `print_performance_statistics()` and private `protection_set_size()`. If `print_performance_statistics()` is ever called from outside the GC path (e.g., a diagnostic signal handler or on-demand stats dump) while another thread calls `register_thread_aterm_pool` or `remove_thread_aterm_pool`, this is a data race on the vector.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `print_performance_statistics()` is declared `public` (aterm_pool.h:115).
  - `register_thread_aterm_pool`/`remove_thread_aterm_pool` mutate the vector under exclusive lock.
  - Concurrent read (iteration in `print_performance_statistics`) + write (push/erase in register/remove) on `std::vector` is a data race.
  - **Test plan:** Add a test that spawns threads registering/unregistering while the main thread calls `print_performance_statistics()` in a loop. Run under TSan (`-DMCRL2_ENABLE_THREADSANITIZER=ON`).
  - Run command: `cmake --build build-tsan -j && ctest --test-dir build-tsan -R <new_test_name> --output-on-failure`
- **Efficiency impact:** Diagnostic-only path; unlikely to affect normal performance.
- **Suggested fix:** Acquire a shared lock at the start of `print_performance_statistics()`:
  ```cpp
  void aterm_pool::print_performance_statistics() const
  {
    // Ensure thread-safe iteration of m_thread_pools.
    mcrl2::utilities::shared_guard guard = m_shared_mutex.lock_shared();
    // ... existing iteration ...
  }
  ```
  Similarly for `protection_set_size()` if it can be reached outside `collect_impl`.

---

### [SEV-3] `std::function` type-erasure overhead in `thread_aterm_pool_interface`

- **Location:** `aterm_pool.h:66-68`
- **Why this is a problem:** Three `std::function` members store trivial forwarding lambdas (each captures only a `this` pointer). `std::function` carries type-erasure overhead: typically 32–48 bytes per instance (on libstdc++/libc++) plus an indirect call. During GC, `mark()` is called for every registered thread — the indirect call prevents inlining of what is ultimately a simple method forwarding. The circular-dependency comment motivating this design could be resolved with a forward-declared virtual interface or CRTP instead.
- **Evidence type:** static-proof
- **Evidence:**
  - `thread_aterm_pool_interface` is constructed in `thread_aterm_pool.h` with `[this] { mark(); }` lambdas.
  - These are called in the GC hot path (`collect_impl` line 219: `pool->mark()`).
  - `sizeof(std::function<void()>)` is 32 bytes on libstdc++ (GCC 11+), times 3 = 96 bytes of overhead per thread pool interface, plus virtual dispatch cost per call.
- **Efficiency impact:** Low-to-moderate. GC frequency depends on workload, but `mark()` is on the critical GC path. For workloads with many threads and frequent GC, the indirect call overhead accumulates.
- **C++20-specific note:** A cleaner C++20 approach would be to make `thread_aterm_pool_interface` a virtual base class (with `mark()`, `print_local_performance_statistics()`, `protection_set_size()` as pure virtual methods) and have `thread_aterm_pool` inherit from it. This replaces 3 function objects with a single vtable pointer.
- **Suggested fix:** Refactor to virtual interface or store a `thread_aterm_pool*` directly (resolving the circular dependency with a forward declaration + out-of-line definitions).

---

### [SEV-3] Missing `[[nodiscard]]` on value-returning public methods

- **Location:** `aterm_pool.h:109-128` (multiple methods)
- **Why this is a problem:** `capacity()`, `size()`, `empty_list()`, `as_int()`, `as_list()`, `as_empty_list()`, and `get_symbol_pool()` return values that callers must use. Discarding them is always a bug. Without `[[nodiscard]]`, compilers won't warn on accidental discard.
- **Evidence type:** static-proof
- **Evidence:** The C++ Core Guidelines (F.21) and mCRL2's own `.github/instructions/cpp.instructions.md` both require `[[nodiscard]]` on pure queries and must-use results.
- **Suggested fix:** Add `[[nodiscard]]` to all accessors returning non-void.

---

### [SEV-3] Missing Doxygen `\pre`/`\post` contracts on public API

- **Location:** `aterm_pool.h:88-128` (public section)
- **Why this is a problem:** Per the mCRL2 coding instructions, new/changed public APIs require explicit precondition/postcondition documentation. `create_function_symbol` has no documented preconditions (e.g., is empty `name` valid? is arity 0 valid?). `register_thread_aterm_pool` documents "threadsafe" but not pre/postconditions. `empty_list()` has no `\post` specifying the returned term is defined.
- **Evidence type:** static-proof (guideline violation)
- **Suggested fix:** Add `\pre` and `\post` clauses to all public methods.

---

### [SEV-3] `m_shared_mutex` declared mutable would allow `const`-correct locking

- **Location:** `aterm_pool.h:224`
- **Why this is a problem:** `print_performance_statistics()` is `const`-qualified but cannot acquire a lock on `m_shared_mutex` without `mutable`. If the SEV-2 fix is applied (adding shared lock to `print_performance_statistics`), `m_shared_mutex` must be declared `mutable` to allow locking in const methods. This is the standard pattern for synchronization primitives.
- **Evidence type:** static-proof
- **Suggested fix:** `mutable mcrl2::utilities::shared_mutex m_shared_mutex;`

## Proposed Tests

- [ ] `test_empty_list_type_safety`: Verify `empty_list()` returns a valid `aterm` with correct function symbol `as_empty_list()`, arity 0, and is `defined()`. Cast correctness test. Run existing atermpp unit tests under UBSan.
- [ ] `test_gc_threshold_overflow_32bit` (plausible): On a 32-bit `long` platform, create >2^31 terms and verify GC threshold remains positive and GC does not storm. Input: loop creating integer terms. Expected failure: GC triggered every term creation after threshold wraps negative.
- [ ] `test_print_stats_concurrent` (plausible): Spawn N threads that register/unregister `thread_aterm_pool_interface` objects while main thread calls `print_performance_statistics()` in a loop. Run under TSan. Expected failure: TSan reports data race on `m_thread_pools`.

## Open Questions

1. Is `print_performance_statistics()` ever called outside the GC exclusive lock in practice? If so, the SEV-2 data race is confirmed rather than plausible.
2. What platforms does mCRL2 target where `sizeof(long) == 4`? If Windows/MSVC is actively supported, the SEV-2 overflow becomes a confirmed issue for large models.
3. The circular dependency between `aterm_pool` and `thread_aterm_pool` motivates the `std::function` design — has a virtual interface approach been considered and rejected for a specific reason (e.g., preventing header inclusion loops)?
4. The TODO on line 118 suggests changing `m_empty_list` to `aterm` — is there a known obstacle preventing this fix?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         6m 24s
Total session time:     5m 0s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         376.9k in, 10.5k out, 337.2k cached (Est. 1 Premium request)
 claude-haiku-4.5        519.6k in, 16.1k out, 451.2k cached (Est. 0 Premium requests)

```

---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T12:15:33.828540+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h` (full file, clean review)
- Risk level: **Medium-High** (one confirmed UB, one exception-safety leak, missing safety annotations)
- Overall verdict: The design is sound for its intended protocol (main-thread-first initialization), but contains a confirmed undefined behavior bug, an exception-safety defect, and missing modern C++ safety annotations.

## Findings (ordered by severity)

### [SEV-1] `m_count_until_check` is uninitialized — undefined behavior on first increment
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:181`
- Why this is a problem: `long m_count_until_check;` has no initializer. The first use is `count_until_check++` inside `aterm_pool::created_term()` (aterm_pool_implementation.h:171). Incrementing an indeterminate value is undefined behavior per [dcl.init]/12 and [expr.pre.incr].
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Any term creation on a freshly constructed `thread_aterm_pool` when `m_enable_garbage_collection || m_enable_resize` is true.
  - Run command: `cmake --build build && ctest --test-dir build -R aterm -j1 --output-on-failure` under MSan (`-DMCRL2_ENABLE_MEMORYSANITIZER=ON`)
  - Observed result: MSan will report "use of uninitialized value" on the `count_until_check++` line.
  - Expected result: Deterministic counter starting at 0.
- Efficiency impact: None beyond the UB itself; practically, the garbage in the counter can delay or prematurely trigger GC checks (non-deterministic behavior).
- C++20-specific note: Use a default member initializer: `long m_count_until_check = 0;`
- Suggested fix:
  ```cpp
  long m_count_until_check = 0;
  ```

### [SEV-2] Exception-unsafe raw `new` in member initializer list
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:40-46`
- Why this is a problem: `m_variables` and `m_containers` are allocated with raw `new`. `m_thread_interface` (constructed last, per declaration order) calls `register_thread_aterm_pool` which does a `std::vector::insert` that can throw `std::bad_alloc`. If that throws, the destructor does NOT run, and both heap-allocated hashtables leak. Similarly, if the `new` for `m_containers` throws, the already-allocated `m_variables` leaks.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Allocate under memory pressure such that the `vector::insert` inside `register_thread_aterm_pool` throws.
  - Run command: Use a custom allocator or fault-injection framework to force `std::bad_alloc` during construction.
  - Observed result: Valgrind/LSan reports leaked hashtable allocations.
  - Expected result: No leak; RAII ensures cleanup.
- Efficiency impact: None at runtime; design robustness issue.
- C++20-specific note: Replace raw pointers with `std::unique_ptr<mcrl2::utilities::hashtable<...>>`. This eliminates the leak on throw, the manual `delete` calls in the destructor, and matches Rust-level ownership clarity.
- Suggested fix:
  ```cpp
  std::unique_ptr<mcrl2::utilities::hashtable<aterm_core*>> m_variables;
  std::unique_ptr<mcrl2::utilities::hashtable<detail::aterm_container*>> m_containers;
  ```
  Update destructor to use `.reset()` or let the unique_ptr destructor handle cleanup.

### [SEV-2] Missing `[[nodiscard]]` on `lock_shared()` and `lock()`
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:152-155`
- Why this is a problem: Both return RAII lock guards. If a caller discards the return value (e.g., `pool.lock_shared();` as a statement), the lock is acquired and immediately released — a silent logic error that would eliminate all mutual-exclusion guarantees.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Write `g_thread_term_pool().lock_shared();` without storing the result. The lock is released before the next statement.
  - Run command: Compile with `-Wunused-result` (implicit with `[[nodiscard]]`).
  - Observed result: No compiler warning without the attribute.
  - Expected result: Compiler error/warning on discarded lock guard.
- Efficiency impact: N/A (API safety)
- C++20-specific note: `[[nodiscard]]` is standard since C++17 and strongly expected on lock-guard factory functions.
- Suggested fix:
  ```cpp
  [[nodiscard]] inline mcrl2::utilities::shared_guard lock_shared() { ... }
  [[nodiscard]] inline mcrl2::utilities::lock_guard lock() { ... }
  ```

### [SEV-3] `static bool is_main_thread` relies on undocumented sequencing protocol
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:49-54`
- Why this is a problem: The non-atomic `static bool` is read and written without synchronization. Correctness depends on the implicit protocol "main thread constructs first." If a library consumer spawns threads before calling `g_thread_term_pool()`, two threads could simultaneously read `true` and both set `m_is_main_thread = true`. The design also duplicates the main-thread detection in `g_thread_term_pool()` (via `g_main_thread_pool == nullptr` check), making the static-local redundant and potentially inconsistent.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Launch two threads that race to call `g_thread_term_pool()` before any other mCRL2 initialization.
  - Run command: Build with `-DMCRL2_ENABLE_THREADSANITIZER=ON`, run the race test.
  - Observed result: TSan reports data race on `is_main_thread`.
  - Expected result: No race; at most one thread identified as main.
- Efficiency impact: Negligible.
- C++20-specific note: If the protocol is intentional, document it with `assert` or `static_assert` in `g_thread_term_pool()`. Consider using `std::atomic<bool>` with `exchange` for defense-in-depth, or remove the static and rely solely on the `g_main_thread_pool` mechanism.
- Suggested fix: Replace with `std::atomic<bool>` + `exchange`, or remove and use `g_main_thread_pool` as the single source of truth for main-thread identity.

### [SEV-3] Header guard `#endif` comment mismatch
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:197`
- Why this is a problem: The `#endif` comment says `ATERMPP_DETAIL_ATERM_POOL_H` but the actual guard is `ATERMPP_DETAIL_THREAD_ATERM_POOL_H`. Misleading during maintenance.
- Evidence type: static-proof
- Evidence:
  - Line 11: `#ifndef ATERMPP_DETAIL_THREAD_ATERM_POOL_H`
  - Line 197: `#endif // ATERMPP_DETAIL_ATERM_POOL_H`
- Suggested fix: `#endif // ATERMPP_DETAIL_THREAD_ATERM_POOL_H`

### [SEV-3] Deliberate resource leak for main thread pool — undocumented
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:60-62`
- Why this is a problem: When `m_is_main_thread` is true, neither `m_variables` nor `m_containers` are freed. The comment explains the rationale (global destruction order), but this causes persistent LeakSanitizer reports unless suppressed. Valgrind/LSan noise reduces developer trust in sanitizer output.
- Evidence type: plausible-issue
- C++20-specific note: Consider documenting this in a `LSAN_DISABLE` annotation or suppression file, or use `std::unique_ptr` with a conditional release in the destructor.
- Suggested fix: Add an LSan suppression or annotate with `__lsan_ignore_object`.

## Proposed Tests
- [ ] `test_uninitialized_count_until_check`: Build with MSan, create a single term on a fresh `thread_aterm_pool`. MSan should flag "use of uninitialized value" on `m_count_until_check++`. After fix, no report. Command: `cmake --build build-msan && ctest --test-dir build-msan -R aterm --output-on-failure`
- [ ] `test_lock_nodiscard_compilation`: Add a static analysis check (or `-Werror=unused-result` with `[[nodiscard]]`) that `pool.lock_shared()` without storing the result produces a compile error.
- [ ] `test_constructor_exception_safety`: Use a custom allocator / mock that forces OOM during `register_thread_aterm_pool`. Verify with LSan that no hashtables leak.

## Open Questions
- Is there a reason `m_variables` and `m_containers` are heap-allocated behind raw pointers rather than direct members? The `delete` in the destructor suggests the intent is to allow the main-thread pool to leak them, but `std::unique_ptr` with conditional `release()` would be safer.
- Should the `is_main_thread` static local be removed entirely in favor of deriving main-thread identity from `g_main_thread_pool == &instance` (already computed in `g_thread_term_pool()`)?
- Is there an existing LSan suppression file for the intentional main-thread leak?

## Copilot Raw Stdout

```text
● Explore: Explore thread_aterm_pool context
  ● Check the static bool is_main_thread logic in constructor
    $ grep -B 5 -A 5 "static bool is_main_thread"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h
    └ 12 lines...
  ● Read libraries/utilities/include/mcrl2/utilities/detail/hashtable.h lines 1-100
    └ 100 lines read

● Explore: Explore initialization and usage
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h lines 1-150
    └ 150 lines read
  ● Find created_term declarations
    $ grep -n "created_term"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h
    └ 2 lines...

● Explore: Check register_variable locking
  ● Find documentation about thread lifecycle
    $ grep -rn "thread.*destructor\|destructor.*thread\|join\|running\|active"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h
    └ 3 lines...
  ● Generate comprehensive analysis
    $ cat > /tmp/thread_safety_analysis.md << 'EOF'
    # Thread Safety Analysis of thread_aterm_pool

    ## 1. IMPLEMENTATION OF register_variable AND register_container...
    └ 556 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool_implementation.h lines
  127-170
  └ 44 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h lines 167-195
  └ 29 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool_implementation.h lines 1-50
  └ 49 lines read

● Read libraries/atermpp/source/aterm_implementation.cpp lines 30-60
  └ 31 lines read

Now I have enough context to produce the review. Let me verify one more detail about exception safety:

● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 36-56
  └ 20 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 170-190
  └ 21 lines read

● Check all uses of m_count_until_check
  $ cd /workspaces/mCRL2 && grep -n "m_count_until_check"
  libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool_implementation.h | head -20
  └ 7 lines...

I now have all the evidence needed. Here is the review:

---

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h` (full file, clean review)
- Risk level: **Medium-High** (one confirmed UB, one exception-safety leak, missing safety annotations)
- Overall verdict: The design is sound for its intended protocol (main-thread-first initialization), but contains a confirmed undefined behavior bug, an exception-safety defect, and missing modern C++ safety annotations.

## Findings (ordered by severity)

### [SEV-1] `m_count_until_check` is uninitialized — undefined behavior on first increment
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:181`
- Why this is a problem: `long m_count_until_check;` has no initializer. The first use is `count_until_check++` inside `aterm_pool::created_term()` (aterm_pool_implementation.h:171). Incrementing an indeterminate value is undefined behavior per [dcl.init]/12 and [expr.pre.incr].
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Any term creation on a freshly constructed `thread_aterm_pool` when `m_enable_garbage_collection || m_enable_resize` is true.
  - Run command: `cmake --build build && ctest --test-dir build -R aterm -j1 --output-on-failure` under MSan (`-DMCRL2_ENABLE_MEMORYSANITIZER=ON`)
  - Observed result: MSan will report "use of uninitialized value" on the `count_until_check++` line.
  - Expected result: Deterministic counter starting at 0.
- Efficiency impact: None beyond the UB itself; practically, the garbage in the counter can delay or prematurely trigger GC checks (non-deterministic behavior).
- C++20-specific note: Use a default member initializer: `long m_count_until_check = 0;`
- Suggested fix:
  ```cpp
  long m_count_until_check = 0;
  ```

### [SEV-2] Exception-unsafe raw `new` in member initializer list
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:40-46`
- Why this is a problem: `m_variables` and `m_containers` are allocated with raw `new`. `m_thread_interface` (constructed last, per declaration order) calls `register_thread_aterm_pool` which does a `std::vector::insert` that can throw `std::bad_alloc`. If that throws, the destructor does NOT run, and both heap-allocated hashtables leak. Similarly, if the `new` for `m_containers` throws, the already-allocated `m_variables` leaks.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Allocate under memory pressure such that the `vector::insert` inside `register_thread_aterm_pool` throws.
  - Run command: Use a custom allocator or fault-injection framework to force `std::bad_alloc` during construction.
  - Observed result: Valgrind/LSan reports leaked hashtable allocations.
  - Expected result: No leak; RAII ensures cleanup.
- Efficiency impact: None at runtime; design robustness issue.
- C++20-specific note: Replace raw pointers with `std::unique_ptr<mcrl2::utilities::hashtable<...>>`. This eliminates the leak on throw, the manual `delete` calls in the destructor, and matches Rust-level ownership clarity.
- Suggested fix:
  ```cpp
  std::unique_ptr<mcrl2::utilities::hashtable<aterm_core*>> m_variables;
  std::unique_ptr<mcrl2::utilities::hashtable<detail::aterm_container*>> m_containers;
  ```
  Update destructor to use `.reset()` or let the unique_ptr destructor handle cleanup.

### [SEV-2] Missing `[[nodiscard]]` on `lock_shared()` and `lock()`
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:152-155`
- Why this is a problem: Both return RAII lock guards. If a caller discards the return value (e.g., `pool.lock_shared();` as a statement), the lock is acquired and immediately released — a silent logic error that would eliminate all mutual-exclusion guarantees.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Write `g_thread_term_pool().lock_shared();` without storing the result. The lock is released before the next statement.
  - Run command: Compile with `-Wunused-result` (implicit with `[[nodiscard]]`).
  - Observed result: No compiler warning without the attribute.
  - Expected result: Compiler error/warning on discarded lock guard.
- Efficiency impact: N/A (API safety)
- C++20-specific note: `[[nodiscard]]` is standard since C++17 and strongly expected on lock-guard factory functions.
- Suggested fix:
  ```cpp
  [[nodiscard]] inline mcrl2::utilities::shared_guard lock_shared() { ... }
  [[nodiscard]] inline mcrl2::utilities::lock_guard lock() { ... }
  ```

### [SEV-3] `static bool is_main_thread` relies on undocumented sequencing protocol
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:49-54`
- Why this is a problem: The non-atomic `static bool` is read and written without synchronization. Correctness depends on the implicit protocol "main thread constructs first." If a library consumer spawns threads before calling `g_thread_term_pool()`, two threads could simultaneously read `true` and both set `m_is_main_thread = true`. The design also duplicates the main-thread detection in `g_thread_term_pool()` (via `g_main_thread_pool == nullptr` check), making the static-local redundant and potentially inconsistent.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Launch two threads that race to call `g_thread_term_pool()` before any other mCRL2 initialization.
  - Run command: Build with `-DMCRL2_ENABLE_THREADSANITIZER=ON`, run the race test.
  - Observed result: TSan reports data race on `is_main_thread`.
  - Expected result: No race; at most one thread identified as main.
- Efficiency impact: Negligible.
- C++20-specific note: If the protocol is intentional, document it with `assert` or `static_assert` in `g_thread_term_pool()`. Consider using `std::atomic<bool>` with `exchange` for defense-in-depth, or remove the static and rely solely on the `g_main_thread_pool` mechanism.
- Suggested fix: Replace with `std::atomic<bool>` + `exchange`, or remove and use `g_main_thread_pool` as the single source of truth for main-thread identity.

### [SEV-3] Header guard `#endif` comment mismatch
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:197`
- Why this is a problem: The `#endif` comment says `ATERMPP_DETAIL_ATERM_POOL_H` but the actual guard is `ATERMPP_DETAIL_THREAD_ATERM_POOL_H`. Misleading during maintenance.
- Evidence type: static-proof
- Evidence:
  - Line 11: `#ifndef ATERMPP_DETAIL_THREAD_ATERM_POOL_H`
  - Line 197: `#endif // ATERMPP_DETAIL_ATERM_POOL_H`
- Suggested fix: `#endif // ATERMPP_DETAIL_THREAD_ATERM_POOL_H`

### [SEV-3] Deliberate resource leak for main thread pool — undocumented
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h:60-62`
- Why this is a problem: When `m_is_main_thread` is true, neither `m_variables` nor `m_containers` are freed. The comment explains the rationale (global destruction order), but this causes persistent LeakSanitizer reports unless suppressed. Valgrind/LSan noise reduces developer trust in sanitizer output.
- Evidence type: plausible-issue
- C++20-specific note: Consider documenting this in a `LSAN_DISABLE` annotation or suppression file, or use `std::unique_ptr` with a conditional release in the destructor.
- Suggested fix: Add an LSan suppression or annotate with `__lsan_ignore_object`.

## Proposed Tests
- [ ] `test_uninitialized_count_until_check`: Build with MSan, create a single term on a fresh `thread_aterm_pool`. MSan should flag "use of uninitialized value" on `m_count_until_check++`. After fix, no report. Command: `cmake --build build-msan && ctest --test-dir build-msan -R aterm --output-on-failure`
- [ ] `test_lock_nodiscard_compilation`: Add a static analysis check (or `-Werror=unused-result` with `[[nodiscard]]`) that `pool.lock_shared()` without storing the result produces a compile error.
- [ ] `test_constructor_exception_safety`: Use a custom allocator / mock that forces OOM during `register_thread_aterm_pool`. Verify with LSan that no hashtables leak.

## Open Questions
- Is there a reason `m_variables` and `m_containers` are heap-allocated behind raw pointers rather than direct members? The `delete` in the destructor suggests the intent is to allow the main-thread pool to leak them, but `std::unique_ptr` with conditional `release()` would be safer.
- Should the `is_main_thread` static local be removed entirely in favor of deriving main-thread identity from `g_main_thread_pool == &instance` (already computed in `g_thread_term_pool()`)?
- Is there an existing LSan suppression file for the intentional main-thread leak?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         7m 9s
Total session time:     7m 47s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         387.1k in, 8.7k out, 352.4k cached (Est. 1 Premium request)
 claude-haiku-4.5        1.1m in, 18.4k out, 988.4k cached (Est. 0 Premium requests)

```

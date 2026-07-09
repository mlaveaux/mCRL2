---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T09:57:51.880826+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h` — the central term pool implementation including garbage collection, term creation routing, and hash table resizing.
- **Risk level:** Medium — one confirmed UB (uninitialized read via caller), one acknowledged strict aliasing violation, and fragile synchronization patterns that work but are maintenance-hazardous.
- **Overall verdict:** The code is functional and its concurrency design (exclusive-lock GC, shared-lock creation) is sound. However, it has one real UB reachable from this file, relies on non-obvious temporary-lifetime semantics for synchronization, and uses repetitive boilerplate that invites copy-paste errors.

## Findings (ordered by severity)

### [SEV-2] Uninitialized `count_until_check` read in `created_term()`

- **Location:** `aterm_pool_implementation.h:171` (via `thread_aterm_pool.h:181`)
- **Why this is a problem:** `created_term()` receives `count_until_check` by reference from `thread_aterm_pool::m_count_until_check`, which is declared as `long m_count_until_check;` without initialization and never assigned before first use. The `count_until_check++` at line 171 reads an indeterminate value — undefined behavior per [basic.indet] in C++20.
- **Evidence type:** static-proof
- **Evidence:**
  - `thread_aterm_pool.h:181` — declaration without initializer
  - No assignment to `m_count_until_check` found anywhere before first call to `created_term()`
  - Other members at lines 177-178 are properly initialized (`= 0`)
  - Run command: `cmake -B build-asan -DMCRL2_ENABLE_ADDRESSSANITIZER=ON && cmake --build build-asan && ctest --test-dir build-asan -R atermpp -j4` (MSan would be ideal here: `-DMCRL2_ENABLE_MEMORYSANITIZER=ON`)
  - Expected result: MSan reports "use of uninitialized value" on first term creation per thread
- **Efficiency impact:** May cause spurious early GC (if garbage > 10000) or delayed first GC (if garbage < 0).
- **C++20-specific note:** In C++20, reading an indeterminate value of type `long` is UB (not merely unspecified).
- **Suggested fix:** In `thread_aterm_pool.h:181`, change to `long m_count_until_check = 0;`

---

### [SEV-2] `reinterpret_cast<aterm&>(m_empty_list)` strict aliasing violation

- **Location:** `aterm_pool_implementation.h:39` (constructor), `aterm_pool.h:118`
- **Why this is a problem:** `m_empty_list` is declared as `aterm_core` but accessed through an `aterm&` reference. While `aterm` derives from `aterm_core` and adds no data members, the actual dynamic type of the object is `aterm_core`, not `aterm`. Accessing through a derived-type reference when the object is not actually a derived instance violates [basic.lval] p11 (strict aliasing). The compiler is permitted to assume no aliasing occurs between `aterm` and `aterm_core` objects.
- **Evidence type:** static-proof
- **Evidence:**
  - `aterm_pool.h:227` declares `aterm_core m_empty_list;`
  - `aterm_pool.h:118` has `TODO remove this reinterpret cast` acknowledging the issue
  - Under `-fstrict-aliasing` (default at -O2+), the compiler may misoptimize
  - Run command: Build with `-O2 -fstrict-aliasing -Wstrict-aliasing=1` and run aterm-related tests
  - Expected result: Potential miscompilation or warning; concrete failure depends on optimizer decisions
- **Efficiency impact:** None in practice.
- **C++20-specific note:** Should use placement-new to construct an `aterm` in-place, or change the member type to `aterm`.
- **Suggested fix:** Change `aterm_core m_empty_list;` to `aterm m_empty_list;` as the TODO suggests.

---

### [SEV-3] Discarded `shared_guard` temporary as synchronization barrier

- **Location:** `aterm_pool_implementation.h:208` and `aterm_pool_implementation.h:395`
- **Why this is a problem:** After `try_lock()` fails, the code calls `shared_mutex.lock_shared();` without capturing the returned `shared_guard`. This creates a temporary that is immediately destroyed, functioning as a "wait for exclusive lock release" barrier. While correct (blocks in `lock_shared_impl()` until the forbidden flag clears, then immediately releases in destructor), the pattern is fragile — any refactoring that moves the expression into a variable or changes temporary lifetime rules would alter semantics.
- **Evidence type:** plausible-issue
- **Evidence:**
  - The pattern relies on: (1) `lock_shared()` blocking until exclusive lock is free, (2) the returned temporary being destroyed at end-of-expression, releasing the lock
  - A `[[nodiscard]]` annotation on `lock_shared()` would flag this as suspicious
  - No test specifically validates this barrier behavior
  - Test plan: Under TSan, have Thread A hold exclusive lock doing slow work, Thread B call `collect_impl()`, verify B does not proceed until A releases
- **Efficiency impact:** Minimal — brief lock/unlock cycle.
- **C++20-specific note:** Consider adding an explicit `wait_for_exclusive_release()` method to `shared_mutex` to make intent clear.
- **Suggested fix:** Add a named method, or at minimum a `[[maybe_unused]]` capture:
  ```cpp
  [[maybe_unused]] auto wait = shared_mutex.lock_shared();
  ```

---

### [SEV-3] Mixed memory orderings on `m_count_until_collection`

- **Location:** `aterm_pool_implementation.h:178` (`relaxed`), lines 90, 197, 276 (implicit `seq_cst`)
- **Why this is a problem:** `fetch_sub` uses `memory_order_relaxed` while the direct loads (`m_count_until_collection > 0` at line 197) and stores (line 90, 276) use the default `seq_cst`. The relaxed subtraction means multiple threads can observe the counter as ≤ 0 simultaneously, causing redundant `collect_impl()` calls. While the `try_lock()` serializes actual GC, it adds unnecessary contention.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Under high thread counts (>8), multiple threads may race through `fetch_sub` and all attempt `try_lock()`, causing lock contention
  - Test plan: Run a multi-threaded term-creation stress test under TSan with 16+ threads, measure GC contention via `EnableGarbageCollectionMetrics`
- **Efficiency impact:** Potential unnecessary lock contention under high parallelism.
- **C++20-specific note:** Consider `memory_order_acq_rel` for `fetch_sub` to ensure the subsequent `collect_impl` reads observe prior writes.
- **Suggested fix:** Use `memory_order_acq_rel` for the `fetch_sub`, or document the intentional relaxed ordering with a rationale comment.

---

### [SEV-3] Repetitive arity-dispatch boilerplate

- **Location:** Lines 50-85, 99-130, 231-270, 285-320, 337-370
- **Why this is a problem:** Five separate functions contain near-identical `switch(arity)` or `std::get<N>` patterns dispatching over arities 0-7. This is error-prone (any new arity-specific method requires updating all switches) and adds ~200 lines of mechanical code.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Maintenance risk: a future change to arity handling that misses one switch case would cause a bug
  - No test specifically validates all arity paths are consistent
- **Efficiency impact:** None at runtime; compile-time and maintenance cost.
- **C++20-specific note:** C++20 `std::apply` with a lambda or a fold-expression over `std::index_sequence` could replace all switch statements. Example:
  ```cpp
  template<typename F>
  auto dispatch_by_arity(std::size_t arity, F&& f) {
    // Use index_sequence to generate switch at compile time
  }
  ```
- **Suggested fix:** Extract a helper template that dispatches by arity index, reducing each call site to one line.

---

### [SEV-3] Missing `m_int_storage.verify_mark()` assertion

- **Location:** `aterm_pool_implementation.h:222-230` (mark verification block)
- **Why this is a problem:** After marking, all `m_appl_storage` entries and `m_appl_dynamic_storage` have `verify_mark()` assertions, but `m_int_storage` does not. While `verify_mark()` only checks terms with arity > 0 (making it vacuous for integer storage with arity 0), the asymmetry is confusing and could mask a future bug if the check's semantics change.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `aterm_pool_storage_implementation.h:365-380`: `verify_mark()` exists for all storage types (including `integer_term_storage`)
  - `m_int_storage.verify_sweep()` IS called at line 242, showing the asymmetry is unintentional
- **Efficiency impact:** None.
- **Suggested fix:** Add `assert(m_int_storage.verify_mark());` at line 222 for consistency.

## Proposed Tests

- [ ] `test_uninitialized_count`: Build with MSan (`-DMCRL2_ENABLE_MEMORYSANITIZER=ON`), create a new thread that constructs terms. MSan should flag the first `count_until_check++` as use-of-uninitialized-value.
- [ ] `test_concurrent_gc_barrier`: With 4+ threads creating terms rapidly, verify under TSan that the `lock_shared()` barrier correctly synchronizes. No data races should be reported.
- [ ] `test_all_arity_paths`: Create terms with arities 0 through 8+ and verify correct routing through the dispatch switches. Ensures no arity path is dead or misrouted.
- [ ] `test_empty_list_aliasing`: Under `-O2 -fstrict-aliasing`, create a term using `empty_list()`, verify it compares equal to a separately constructed empty list. Exercises the `reinterpret_cast` path.

## Open Questions

1. Is the `reinterpret_cast` for `m_empty_list` blocked on a larger refactor, or could `m_empty_list` simply be declared as `aterm` today?
2. Why was `memory_order_relaxed` chosen for `fetch_sub` — is there a measured performance benefit vs. `acq_rel`?
3. Should `m_count_until_check` initialization be added in the declaration (zero-init) or in the constructor's member-init-list with a specific starting value?

## Copilot Raw Stdout

```text
● Explore: Explore aterm_pool class definition
  ● Grep "EnableGarbageCollection|EnableAggressiveGarbageCollection|EnableGarbageCollectionMetrics"
  in h files (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
    └ 9 lines found
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_configuration.h
    └ 47 lines read

● Explore: Investigate shared_mutex lock patterns
  ● Check sweep and print_performance_statistics context in function_symbol_pool.cpp
    $ grep -r "\.sweep\(\)" /workspaces/mCRL2/libraries/atermpp/source/function_symbol_pool.cpp -B5
    -n | head -20
    └ 1 line...
  ● Read libraries/atermpp/source/function_symbol_pool.cpp lines 150-220
    └ 71 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h lines 148-175
  └ 28 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h lines 82-230
  └ 149 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h lines 167-215
  └ 48 lines read

● Grep "protection_set_size" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ 13 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h lines 438-455
  └ 16 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool_implementation.h lines 30-65
  └ 36 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 34-70
  └ 37 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 170-195
  └ 24 lines read

● Explore: Check aterm inheritance and verify_mark
  ● Check verify_mark implementation
    $ grep -A3 "bool verify_mark()" /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail
    /aterm_pool_storage_implementation.h | head -20
    └ 1 line...
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
    365-400
    └ 36 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 49-56
  └ 7 lines read

● Grep "verify_mark" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_st
orage_implementation.h)
  └ 21 lines found

Now I have all the information needed to compile my review. Let me verify one more concern about the `lock_shared()` temporary pattern:

● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 204-215
  └ 11 lines read

● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 31-60
  └ 30 lines read

● Grep "m_count_until_check"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h)
  └ 1 line found

● Grep "m_count_until_check\s*=" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ No matches found

● Grep "count_until_check" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_
aterm_pool_implementation.h)
  └ 6 lines found

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h` — the central term pool implementation including garbage collection, term creation routing, and hash table resizing.
- **Risk level:** Medium — one confirmed UB (uninitialized read via caller), one acknowledged strict aliasing violation, and fragile synchronization patterns that work but are maintenance-hazardous.
- **Overall verdict:** The code is functional and its concurrency design (exclusive-lock GC, shared-lock creation) is sound. However, it has one real UB reachable from this file, relies on non-obvious temporary-lifetime semantics for synchronization, and uses repetitive boilerplate that invites copy-paste errors.

## Findings (ordered by severity)

### [SEV-2] Uninitialized `count_until_check` read in `created_term()`

- **Location:** `aterm_pool_implementation.h:171` (via `thread_aterm_pool.h:181`)
- **Why this is a problem:** `created_term()` receives `count_until_check` by reference from `thread_aterm_pool::m_count_until_check`, which is declared as `long m_count_until_check;` without initialization and never assigned before first use. The `count_until_check++` at line 171 reads an indeterminate value — undefined behavior per [basic.indet] in C++20.
- **Evidence type:** static-proof
- **Evidence:**
  - `thread_aterm_pool.h:181` — declaration without initializer
  - No assignment to `m_count_until_check` found anywhere before first call to `created_term()`
  - Other members at lines 177-178 are properly initialized (`= 0`)
  - Run command: `cmake -B build-asan -DMCRL2_ENABLE_ADDRESSSANITIZER=ON && cmake --build build-asan && ctest --test-dir build-asan -R atermpp -j4` (MSan would be ideal here: `-DMCRL2_ENABLE_MEMORYSANITIZER=ON`)
  - Expected result: MSan reports "use of uninitialized value" on first term creation per thread
- **Efficiency impact:** May cause spurious early GC (if garbage > 10000) or delayed first GC (if garbage < 0).
- **C++20-specific note:** In C++20, reading an indeterminate value of type `long` is UB (not merely unspecified).
- **Suggested fix:** In `thread_aterm_pool.h:181`, change to `long m_count_until_check = 0;`

---

### [SEV-2] `reinterpret_cast<aterm&>(m_empty_list)` strict aliasing violation

- **Location:** `aterm_pool_implementation.h:39` (constructor), `aterm_pool.h:118`
- **Why this is a problem:** `m_empty_list` is declared as `aterm_core` but accessed through an `aterm&` reference. While `aterm` derives from `aterm_core` and adds no data members, the actual dynamic type of the object is `aterm_core`, not `aterm`. Accessing through a derived-type reference when the object is not actually a derived instance violates [basic.lval] p11 (strict aliasing). The compiler is permitted to assume no aliasing occurs between `aterm` and `aterm_core` objects.
- **Evidence type:** static-proof
- **Evidence:**
  - `aterm_pool.h:227` declares `aterm_core m_empty_list;`
  - `aterm_pool.h:118` has `TODO remove this reinterpret cast` acknowledging the issue
  - Under `-fstrict-aliasing` (default at -O2+), the compiler may misoptimize
  - Run command: Build with `-O2 -fstrict-aliasing -Wstrict-aliasing=1` and run aterm-related tests
  - Expected result: Potential miscompilation or warning; concrete failure depends on optimizer decisions
- **Efficiency impact:** None in practice.
- **C++20-specific note:** Should use placement-new to construct an `aterm` in-place, or change the member type to `aterm`.
- **Suggested fix:** Change `aterm_core m_empty_list;` to `aterm m_empty_list;` as the TODO suggests.

---

### [SEV-3] Discarded `shared_guard` temporary as synchronization barrier

- **Location:** `aterm_pool_implementation.h:208` and `aterm_pool_implementation.h:395`
- **Why this is a problem:** After `try_lock()` fails, the code calls `shared_mutex.lock_shared();` without capturing the returned `shared_guard`. This creates a temporary that is immediately destroyed, functioning as a "wait for exclusive lock release" barrier. While correct (blocks in `lock_shared_impl()` until the forbidden flag clears, then immediately releases in destructor), the pattern is fragile — any refactoring that moves the expression into a variable or changes temporary lifetime rules would alter semantics.
- **Evidence type:** plausible-issue
- **Evidence:**
  - The pattern relies on: (1) `lock_shared()` blocking until exclusive lock is free, (2) the returned temporary being destroyed at end-of-expression, releasing the lock
  - A `[[nodiscard]]` annotation on `lock_shared()` would flag this as suspicious
  - No test specifically validates this barrier behavior
  - Test plan: Under TSan, have Thread A hold exclusive lock doing slow work, Thread B call `collect_impl()`, verify B does not proceed until A releases
- **Efficiency impact:** Minimal — brief lock/unlock cycle.
- **C++20-specific note:** Consider adding an explicit `wait_for_exclusive_release()` method to `shared_mutex` to make intent clear.
- **Suggested fix:** Add a named method, or at minimum a `[[maybe_unused]]` capture:
  ```cpp
  [[maybe_unused]] auto wait = shared_mutex.lock_shared();
  ```

---

### [SEV-3] Mixed memory orderings on `m_count_until_collection`

- **Location:** `aterm_pool_implementation.h:178` (`relaxed`), lines 90, 197, 276 (implicit `seq_cst`)
- **Why this is a problem:** `fetch_sub` uses `memory_order_relaxed` while the direct loads (`m_count_until_collection > 0` at line 197) and stores (line 90, 276) use the default `seq_cst`. The relaxed subtraction means multiple threads can observe the counter as ≤ 0 simultaneously, causing redundant `collect_impl()` calls. While the `try_lock()` serializes actual GC, it adds unnecessary contention.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Under high thread counts (>8), multiple threads may race through `fetch_sub` and all attempt `try_lock()`, causing lock contention
  - Test plan: Run a multi-threaded term-creation stress test under TSan with 16+ threads, measure GC contention via `EnableGarbageCollectionMetrics`
- **Efficiency impact:** Potential unnecessary lock contention under high parallelism.
- **C++20-specific note:** Consider `memory_order_acq_rel` for `fetch_sub` to ensure the subsequent `collect_impl` reads observe prior writes.
- **Suggested fix:** Use `memory_order_acq_rel` for the `fetch_sub`, or document the intentional relaxed ordering with a rationale comment.

---

### [SEV-3] Repetitive arity-dispatch boilerplate

- **Location:** Lines 50-85, 99-130, 231-270, 285-320, 337-370
- **Why this is a problem:** Five separate functions contain near-identical `switch(arity)` or `std::get<N>` patterns dispatching over arities 0-7. This is error-prone (any new arity-specific method requires updating all switches) and adds ~200 lines of mechanical code.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Maintenance risk: a future change to arity handling that misses one switch case would cause a bug
  - No test specifically validates all arity paths are consistent
- **Efficiency impact:** None at runtime; compile-time and maintenance cost.
- **C++20-specific note:** C++20 `std::apply` with a lambda or a fold-expression over `std::index_sequence` could replace all switch statements. Example:
  ```cpp
  template<typename F>
  auto dispatch_by_arity(std::size_t arity, F&& f) {
    // Use index_sequence to generate switch at compile time
  }
  ```
- **Suggested fix:** Extract a helper template that dispatches by arity index, reducing each call site to one line.

---

### [SEV-3] Missing `m_int_storage.verify_mark()` assertion

- **Location:** `aterm_pool_implementation.h:222-230` (mark verification block)
- **Why this is a problem:** After marking, all `m_appl_storage` entries and `m_appl_dynamic_storage` have `verify_mark()` assertions, but `m_int_storage` does not. While `verify_mark()` only checks terms with arity > 0 (making it vacuous for integer storage with arity 0), the asymmetry is confusing and could mask a future bug if the check's semantics change.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `aterm_pool_storage_implementation.h:365-380`: `verify_mark()` exists for all storage types (including `integer_term_storage`)
  - `m_int_storage.verify_sweep()` IS called at line 242, showing the asymmetry is unintentional
- **Efficiency impact:** None.
- **Suggested fix:** Add `assert(m_int_storage.verify_mark());` at line 222 for consistency.

## Proposed Tests

- [ ] `test_uninitialized_count`: Build with MSan (`-DMCRL2_ENABLE_MEMORYSANITIZER=ON`), create a new thread that constructs terms. MSan should flag the first `count_until_check++` as use-of-uninitialized-value.
- [ ] `test_concurrent_gc_barrier`: With 4+ threads creating terms rapidly, verify under TSan that the `lock_shared()` barrier correctly synchronizes. No data races should be reported.
- [ ] `test_all_arity_paths`: Create terms with arities 0 through 8+ and verify correct routing through the dispatch switches. Ensures no arity path is dead or misrouted.
- [ ] `test_empty_list_aliasing`: Under `-O2 -fstrict-aliasing`, create a term using `empty_list()`, verify it compares equal to a separately constructed empty list. Exercises the `reinterpret_cast` path.

## Open Questions

1. Is the `reinterpret_cast` for `m_empty_list` blocked on a larger refactor, or could `m_empty_list` simply be declared as `aterm` today?
2. Why was `memory_order_relaxed` chosen for `fetch_sub` — is there a measured performance benefit vs. `acq_rel`?
3. Should `m_count_until_check` initialization be added in the declaration (zero-init) or in the constructor's member-init-list with a specific starting value?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         8m 26s
Total session time:     8m 53s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         723.9k in, 12.7k out, 675.3k cached (Est. 1 Premium request)
 claude-haiku-4.5        905.5k in, 12.9k out, 825.9k cached (Est. 0 Premium requests)

```

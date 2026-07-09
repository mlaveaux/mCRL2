---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T12:03:02.004060+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h` (65 lines) — the global aterm pool singleton accessor, placement-new lifecycle, and lazy initialization logic.
- **Risk level:** Medium — correctness issues are mitigated by current usage patterns but represent latent UB that could manifest under different initialization ordering or future refactoring.
- **Overall verdict:** The design is intentional and works in practice, but it contains two instances of technically undefined behavior (missing `std::launder`, unsynchronized flag access) and lacks contracts/documentation for its safety invariants.

## Findings (ordered by severity)

### [SEV-2] Missing `std::launder` after placement-new — undefined behavior per C++20
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:34`
- **Why this is a problem:** After `new (&g_aterm_pool_instance()) aterm_pool()` constructs an object in byte storage, accessing that object through a pointer obtained from `reinterpret_cast<aterm_pool*>(g_aterm_pool_storage.data())` is undefined behavior per \[basic.life\]/8 and \[ptr.launder\]. The pointer from `data()` originally designates `std::byte` subobjects; it does not automatically point to the newly-created `aterm_pool`. `aterm_pool` is not an implicit-lifetime type (user-declared constructor), so C++20 P0593R6 does not rescue this. The codebase has **zero** uses of `std::launder` anywhere.
- **Evidence type:** static-proof
- **Evidence:**
  - C++20 \[ptr.launder\]/3 requires `std::launder` to obtain a valid pointer to an object created via placement-new when accessed through the original storage address.
  - `g_aterm_pool_instance()` is called hundreds of times across the codebase (every `g_term_pool()` call), each time dereferencing a formally invalid pointer.
  - Run command: Build with `-fsanitize=undefined` (UBSan doesn't currently flag this, but it's a standard violation that compilers may exploit in future optimization passes).
  - Observed: All compilers generate correct code today.
  - Expected: Formally correct code should use `std::launder`.
- **Efficiency impact:** None — `std::launder` is a compiler barrier only.
- **C++20-specific note:** `std::launder` was introduced in C++17; all mCRL2 target compilers support it. This is the canonical use case.
- **Suggested fix:**
  ```cpp
  #include <new> // for std::launder
  
  inline aterm_pool& g_aterm_pool_instance()
  {
    return *std::launder(reinterpret_cast<aterm_pool*>(g_aterm_pool_storage.data()));
  }
  ```

---

### [SEV-3] Unsynchronized read/write of `static bool initialized` — data race (UB)
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:45–51`
- **Why this is a problem:** `static bool initialized` is a block-scope static whose *initialization* to `false` is thread-safe (C++11 magic statics guarantee). However, subsequent reads (`if (!initialized)`) and the write (`initialized = true`) are ordinary, unsynchronized memory accesses. If two threads call `g_term_pool<true>()` concurrently before `initialized` is `true`, both read `false`, both construct the pool (double placement-new over a live object = UB), and both call `start_gc_stress_thread()` (which assigns to a non-joinable-checked `std::thread`, causing `std::terminate` if the stress thread is enabled).
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: `global_function_symbol` constructor calls `g_term_pool<true>()`. `global_function_symbol` instances are used as function-local statics in `function_symbols.h`. If two threads trigger different function-local statics before the file-scope statics in `function_symbol.cpp` have run, concurrent calls to `g_term_pool<true>()` occur.
  - Mitigation in practice: File-scope statics in `function_symbol.cpp` (lines 15–17) call `g_term_pool<true>()` during dynamic initialization of that TU, which runs single-threaded before `main()`. After that, `initialized == true` and the branch is never taken. Thread-creation after `main()` provides a happens-before edge making the `true` value visible.
  - This relies on an **undocumented invariant**: no thread spawns before `function_symbol.cpp`'s statics initialize. If a future change spawns threads earlier, or if link ordering changes cause another TU's statics to run first and trigger `function_symbols.h` inline functions, the race materializes.
  - Test plan: Compile with `EnableGCStressThread = true` and TSan, spawn two threads that both call `g_term_pool<true>()` simultaneously before the file-scope statics have run.
  - Run command: `cmake -S . -B build-tsan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_THREADSANITIZER=ON -DMCRL2_ENABLE_TESTS=ON && cmake --build build-tsan && ctest --test-dir build-tsan -R <new_test> --output-on-failure`
- **Efficiency impact:** None.
- **C++20-specific note:** The idiomatic C++11+ fix is `std::call_once` or relying on magic statics for the entire initialization block (not just the flag). Alternatively, `std::atomic<bool>` with double-checked locking or `std::once_flag`.
- **Suggested fix:**
  ```cpp
  #include <mutex>
  
  template<bool lazy = false>
  inline aterm_pool& g_term_pool()
  {
    if constexpr (lazy)
    {
      static std::once_flag flag;
      std::call_once(flag, [] {
        new (&g_aterm_pool_instance()) aterm_pool();
        start_gc_stress_thread();
      });
    }
    return g_aterm_pool_instance();
  }
  ```

---

### [SEV-3] `g_aterm_pool_instance()` is callable before placement-new — no precondition enforcement
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:32–35`
- **Why this is a problem:** `g_aterm_pool_instance()` is an `inline` function in a public header with no precondition check. Any code calling `g_term_pool<false>()` (the default, used hundreds of times) or `g_aterm_pool_instance()` directly before `g_term_pool<true>()` has been called operates on zero-initialized byte storage cast to `aterm_pool&` — UB (accessing a non-existent object). The correctness of the entire codebase depends on an implicit ordering guarantee that is nowhere documented or asserted.
- **Evidence type:** static-proof
- **Evidence:**
  - `g_aterm_pool_storage` is zero-initialized at link time (`= {}`).
  - Before `g_term_pool<true>()` runs, no `aterm_pool` object exists at that address.
  - `g_term_pool<false>()` unconditionally returns `g_aterm_pool_instance()` — no guard.
  - Any tool or library that triggers `g_term_pool()` before `function_symbol.cpp`'s dynamic init runs will operate on non-object memory.
  - Observed: Works in practice because linkers initialize `function_symbol.cpp` statics early. But this is not guaranteed across TUs by the standard (SIOF territory).
- **Efficiency impact:** An `assert` on a flag in debug builds costs nothing meaningful.
- **C++20-specific note:** `[[assume(initialized)]]` (C++23) or `assert` are appropriate tools.
- **Suggested fix:** Add a debug assertion to `g_aterm_pool_instance()`:
  ```cpp
  inline aterm_pool& g_aterm_pool_instance()
  {
    assert(/* pool has been initialized */);
    return *std::launder(reinterpret_cast<aterm_pool*>(g_aterm_pool_storage.data()));
  }
  ```
  This requires exposing the initialization state (e.g., an `std::atomic<bool>` set after construction).

---

### [SEV-3] No documented lifetime contract — deliberate leak without rationale
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:48`
- **Why this is a problem:** The `aterm_pool` is constructed via placement-new but never destroyed. This is a valid pattern to avoid the static destruction order fiasco, but it is not documented. The Doxygen comment on `g_term_pool()` says "provides lazy initialization" but does not mention:
  - The object is intentionally leaked.
  - No destructor is ever called.
  - Resources held by `aterm_pool` (mutexes, allocations) are abandoned at exit.
  - The design decision rationale.
- **Evidence type:** static-proof (absence of documentation).
- **Efficiency impact:** Memory is not reclaimed at exit (acceptable for process-lifetime singletons, but confuses LeakSanitizer unless suppressed).
- **C++20-specific note:** N/A.
- **Suggested fix:** Add `\post` documentation and a LeakSanitizer suppression comment:
  ```cpp
  /// \brief obtain a reference to the global aterm pool.
  /// \details Provides lazy initialization on first call with lazy=true.
  ///          The pool is intentionally never destroyed (leaky singleton pattern)
  ///          to avoid static destruction order issues. LeakSanitizer findings
  ///          for this allocation are expected and benign.
  /// \pre For lazy=false: the pool must have been previously initialized
  ///      via a call with lazy=true.
  ```

---

## Proposed Tests
- [ ] `test_concurrent_lazy_pool_init`: Spawn two threads that both call `g_term_pool<true>()` simultaneously (requires making a test-only entry point or using `global_function_symbol` constructors from two threads). Run under TSan. Signal: TSan data-race report on `initialized`. Command: `ctest --test-dir build-tsan -R test_concurrent_lazy_pool_init --output-on-failure`
- [ ] `test_pool_access_before_init`: In an isolated TU with controlled init ordering, call `g_term_pool<false>()` before `g_term_pool<true>()`. Signal: assertion failure (after fix), or ASan/UBSan report on current code. Command: `ctest --test-dir build-asan -R test_pool_access_before_init --output-on-failure`

## Open Questions
1. Is the undocumented invariant "no thread spawns before `function_symbol.cpp` dynamic init completes" intentional and stable? If so, should it be enforced (e.g., `std::call_once`) or at least asserted?
2. Has the interaction with LeakSanitizer been considered? The intentionally-leaked `aterm_pool` will appear in LSan reports unless suppressed. Is there an existing suppression file?
3. Given that `global_function_symbol` objects are function-local statics (thread-safe init guaranteed by the compiler), could the lazy initialization of the global pool itself be restructured as a Meyers singleton (`static aterm_pool& instance() { static aterm_pool p; return p; }`) to eliminate both the manual flag and the `std::launder` requirement? The current byte-array approach was presumably chosen to avoid destruction, but `[[gnu::no_destroy]]` / `[[msvc::no_destroy]]` or `std::optional` with no reset could achieve the same.

## Copilot Raw Stdout

```text
● Explore: Explore aterm pool usage
  ● Read libraries/atermpp/source/function_symbol.cpp
    └ 29 lines read
  ● Read libraries/atermpp/source/gc_stress_thread.cpp
    └ 58 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h
  └ 65 lines read

● Read libraries/atermpp/source/aterm_implementation.cpp lines 60-80
  └ 14 lines read

● Read libraries/atermpp/source/gc_stress_thread.cpp
  └ 58 lines read

● Grep "g_term_pool<true>"
  └ 75 lines found

Now I have all the information needed to produce the review. Let me verify one more detail about whether `global_function_symbol` can be called after static init:

● Grep "global_function_symbol" (/workspaces/mCRL2/libraries)
  └ 4 files found

● Grep "std::launder" (/workspaces/mCRL2/libraries)
  └ No matches found

● Grep "class global_function_symbol"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/function_symbol.h)
  └ 20 lines found

● Grep "global_function_symbol"
(/workspaces/mCRL2/libraries/core/include/mcrl2/core/detail/function_symbols.h)
  └ 20 lines found

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h` (65 lines) — the global aterm pool singleton accessor, placement-new lifecycle, and lazy initialization logic.
- **Risk level:** Medium — correctness issues are mitigated by current usage patterns but represent latent UB that could manifest under different initialization ordering or future refactoring.
- **Overall verdict:** The design is intentional and works in practice, but it contains two instances of technically undefined behavior (missing `std::launder`, unsynchronized flag access) and lacks contracts/documentation for its safety invariants.

## Findings (ordered by severity)

### [SEV-2] Missing `std::launder` after placement-new — undefined behavior per C++20
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:34`
- **Why this is a problem:** After `new (&g_aterm_pool_instance()) aterm_pool()` constructs an object in byte storage, accessing that object through a pointer obtained from `reinterpret_cast<aterm_pool*>(g_aterm_pool_storage.data())` is undefined behavior per \[basic.life\]/8 and \[ptr.launder\]. The pointer from `data()` originally designates `std::byte` subobjects; it does not automatically point to the newly-created `aterm_pool`. `aterm_pool` is not an implicit-lifetime type (user-declared constructor), so C++20 P0593R6 does not rescue this. The codebase has **zero** uses of `std::launder` anywhere.
- **Evidence type:** static-proof
- **Evidence:**
  - C++20 \[ptr.launder\]/3 requires `std::launder` to obtain a valid pointer to an object created via placement-new when accessed through the original storage address.
  - `g_aterm_pool_instance()` is called hundreds of times across the codebase (every `g_term_pool()` call), each time dereferencing a formally invalid pointer.
  - Run command: Build with `-fsanitize=undefined` (UBSan doesn't currently flag this, but it's a standard violation that compilers may exploit in future optimization passes).
  - Observed: All compilers generate correct code today.
  - Expected: Formally correct code should use `std::launder`.
- **Efficiency impact:** None — `std::launder` is a compiler barrier only.
- **C++20-specific note:** `std::launder` was introduced in C++17; all mCRL2 target compilers support it. This is the canonical use case.
- **Suggested fix:**
  ```cpp
  #include <new> // for std::launder
  
  inline aterm_pool& g_aterm_pool_instance()
  {
    return *std::launder(reinterpret_cast<aterm_pool*>(g_aterm_pool_storage.data()));
  }
  ```

---

### [SEV-3] Unsynchronized read/write of `static bool initialized` — data race (UB)
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:45–51`
- **Why this is a problem:** `static bool initialized` is a block-scope static whose *initialization* to `false` is thread-safe (C++11 magic statics guarantee). However, subsequent reads (`if (!initialized)`) and the write (`initialized = true`) are ordinary, unsynchronized memory accesses. If two threads call `g_term_pool<true>()` concurrently before `initialized` is `true`, both read `false`, both construct the pool (double placement-new over a live object = UB), and both call `start_gc_stress_thread()` (which assigns to a non-joinable-checked `std::thread`, causing `std::terminate` if the stress thread is enabled).
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: `global_function_symbol` constructor calls `g_term_pool<true>()`. `global_function_symbol` instances are used as function-local statics in `function_symbols.h`. If two threads trigger different function-local statics before the file-scope statics in `function_symbol.cpp` have run, concurrent calls to `g_term_pool<true>()` occur.
  - Mitigation in practice: File-scope statics in `function_symbol.cpp` (lines 15–17) call `g_term_pool<true>()` during dynamic initialization of that TU, which runs single-threaded before `main()`. After that, `initialized == true` and the branch is never taken. Thread-creation after `main()` provides a happens-before edge making the `true` value visible.
  - This relies on an **undocumented invariant**: no thread spawns before `function_symbol.cpp`'s statics initialize. If a future change spawns threads earlier, or if link ordering changes cause another TU's statics to run first and trigger `function_symbols.h` inline functions, the race materializes.
  - Test plan: Compile with `EnableGCStressThread = true` and TSan, spawn two threads that both call `g_term_pool<true>()` simultaneously before the file-scope statics have run.
  - Run command: `cmake -S . -B build-tsan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_THREADSANITIZER=ON -DMCRL2_ENABLE_TESTS=ON && cmake --build build-tsan && ctest --test-dir build-tsan -R <new_test> --output-on-failure`
- **Efficiency impact:** None.
- **C++20-specific note:** The idiomatic C++11+ fix is `std::call_once` or relying on magic statics for the entire initialization block (not just the flag). Alternatively, `std::atomic<bool>` with double-checked locking or `std::once_flag`.
- **Suggested fix:**
  ```cpp
  #include <mutex>
  
  template<bool lazy = false>
  inline aterm_pool& g_term_pool()
  {
    if constexpr (lazy)
    {
      static std::once_flag flag;
      std::call_once(flag, [] {
        new (&g_aterm_pool_instance()) aterm_pool();
        start_gc_stress_thread();
      });
    }
    return g_aterm_pool_instance();
  }
  ```

---

### [SEV-3] `g_aterm_pool_instance()` is callable before placement-new — no precondition enforcement
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:32–35`
- **Why this is a problem:** `g_aterm_pool_instance()` is an `inline` function in a public header with no precondition check. Any code calling `g_term_pool<false>()` (the default, used hundreds of times) or `g_aterm_pool_instance()` directly before `g_term_pool<true>()` has been called operates on zero-initialized byte storage cast to `aterm_pool&` — UB (accessing a non-existent object). The correctness of the entire codebase depends on an implicit ordering guarantee that is nowhere documented or asserted.
- **Evidence type:** static-proof
- **Evidence:**
  - `g_aterm_pool_storage` is zero-initialized at link time (`= {}`).
  - Before `g_term_pool<true>()` runs, no `aterm_pool` object exists at that address.
  - `g_term_pool<false>()` unconditionally returns `g_aterm_pool_instance()` — no guard.
  - Any tool or library that triggers `g_term_pool()` before `function_symbol.cpp`'s dynamic init runs will operate on non-object memory.
  - Observed: Works in practice because linkers initialize `function_symbol.cpp` statics early. But this is not guaranteed across TUs by the standard (SIOF territory).
- **Efficiency impact:** An `assert` on a flag in debug builds costs nothing meaningful.
- **C++20-specific note:** `[[assume(initialized)]]` (C++23) or `assert` are appropriate tools.
- **Suggested fix:** Add a debug assertion to `g_aterm_pool_instance()`:
  ```cpp
  inline aterm_pool& g_aterm_pool_instance()
  {
    assert(/* pool has been initialized */);
    return *std::launder(reinterpret_cast<aterm_pool*>(g_aterm_pool_storage.data()));
  }
  ```
  This requires exposing the initialization state (e.g., an `std::atomic<bool>` set after construction).

---

### [SEV-3] No documented lifetime contract — deliberate leak without rationale
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:48`
- **Why this is a problem:** The `aterm_pool` is constructed via placement-new but never destroyed. This is a valid pattern to avoid the static destruction order fiasco, but it is not documented. The Doxygen comment on `g_term_pool()` says "provides lazy initialization" but does not mention:
  - The object is intentionally leaked.
  - No destructor is ever called.
  - Resources held by `aterm_pool` (mutexes, allocations) are abandoned at exit.
  - The design decision rationale.
- **Evidence type:** static-proof (absence of documentation).
- **Efficiency impact:** Memory is not reclaimed at exit (acceptable for process-lifetime singletons, but confuses LeakSanitizer unless suppressed).
- **C++20-specific note:** N/A.
- **Suggested fix:** Add `\post` documentation and a LeakSanitizer suppression comment:
  ```cpp
  /// \brief obtain a reference to the global aterm pool.
  /// \details Provides lazy initialization on first call with lazy=true.
  ///          The pool is intentionally never destroyed (leaky singleton pattern)
  ///          to avoid static destruction order issues. LeakSanitizer findings
  ///          for this allocation are expected and benign.
  /// \pre For lazy=false: the pool must have been previously initialized
  ///      via a call with lazy=true.
  ```

---

## Proposed Tests
- [ ] `test_concurrent_lazy_pool_init`: Spawn two threads that both call `g_term_pool<true>()` simultaneously (requires making a test-only entry point or using `global_function_symbol` constructors from two threads). Run under TSan. Signal: TSan data-race report on `initialized`. Command: `ctest --test-dir build-tsan -R test_concurrent_lazy_pool_init --output-on-failure`
- [ ] `test_pool_access_before_init`: In an isolated TU with controlled init ordering, call `g_term_pool<false>()` before `g_term_pool<true>()`. Signal: assertion failure (after fix), or ASan/UBSan report on current code. Command: `ctest --test-dir build-asan -R test_pool_access_before_init --output-on-failure`

## Open Questions
1. Is the undocumented invariant "no thread spawns before `function_symbol.cpp` dynamic init completes" intentional and stable? If so, should it be enforced (e.g., `std::call_once`) or at least asserted?
2. Has the interaction with LeakSanitizer been considered? The intentionally-leaked `aterm_pool` will appear in LSan reports unless suppressed. Is there an existing suppression file?
3. Given that `global_function_symbol` objects are function-local statics (thread-safe init guaranteed by the compiler), could the lazy initialization of the global pool itself be restructured as a Meyers singleton (`static aterm_pool& instance() { static aterm_pool p; return p; }`) to eliminate both the manual flag and the `std::launder` requirement? The current byte-array approach was presumably chosen to avoid destruction, but `[[gnu::no_destroy]]` / `[[msvc::no_destroy]]` or `std::optional` with no reset could achieve the same.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         4m 55s
Total session time:     5m 15s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         273.0k in, 9.0k out, 241.9k cached (Est. 1 Premium request)
 claude-haiku-4.5        638.2k in, 9.4k out, 600.7k cached (Est. 0 Premium requests)

```

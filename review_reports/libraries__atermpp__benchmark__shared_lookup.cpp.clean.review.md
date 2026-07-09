---
{
  "file": "libraries/atermpp/benchmark/shared_lookup.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-sonnet-4.6",
  "thinking": "deep",
  "effort": "high",
  "generated_at_utc": "2026-07-07T12:20:36.535897+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/benchmark/shared_lookup.cpp` plus transitive shared context in `benchmark_shared.h`, `aterm_pool.h`, and `aterm_configuration.h`
- **Risk level:** High — benchmark is broken and uncompilable when benchmarks are enabled; behavioral correctness also compromised
- **Overall verdict:** Do not ship. One regression introduced by commit `31c2b9f6ed` (Jun 30 2026) makes every benchmark in this directory fail to compile. Secondary findings affect measurement correctness and code safety.

---

## Findings (ordered by severity)

---

### [SEV-1] `enable_garbage_collection` method no longer exists — all benchmarks fail to compile

- **Location:** `libraries/atermpp/benchmark/shared_lookup.cpp:16` (and all 5 sibling benchmarks)
- **Why this is a problem:**
  Commit `31c2b9f6ed` ("Relaxed/removed some atomic counters", Jun 30 2026) changed `m_enable_garbage_collection` from `std::atomic<bool>` to `static constexpr bool` and **deleted** the public method `enable_garbage_collection(bool)` from `aterm_pool`. Every benchmark in the directory still calls `detail::g_term_pool().enable_garbage_collection(false)`. When `MCRL2_ENABLE_BENCHMARKS=ON` this is a hard compilation error.

- **Evidence type:** static-proof

- **Evidence:**
  - **Reproducer:** Configure with `-DMCRL2_ENABLE_BENCHMARKS=ON` and build benchmark targets.
  - **Run command:**
    ```bash
    cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug \
          -DMCRL2_ENABLE_GUI_TOOLS=OFF -DMCRL2_ENABLE_BENCHMARKS=ON
    cmake --build build --target benchmark_target_atermpp_shared_lookup
    ```
  - **Observed result:** `error: 'class atermpp::detail::aterm_pool' has no member named 'enable_garbage_collection'`
  - **Expected result:** Successful compilation

- **Behavioral impact:** GC is now compile-time-controlled via `constexpr static bool EnableGarbageCollection = true` (`aterm_configuration.h:19`). There is no runtime API to disable it. With GC enabled, the "lookup" benchmark silently degrades into a mixed creation-and-GC benchmark: terms built during one `create_nested_function` iteration can be reclaimed between iterations (the outer `aterm f` only roots the final term, not the intermediates), altering what is being measured.

- **C++20-specific note:** `static constexpr` is correct C++20 for compile-time constants; the original `std::atomic<bool>` member for a flag that was only written once at startup was itself unnecessary overhead.

- **Suggested fix:** Remove the `enable_garbage_collection(false)` call from all 6 benchmark files. If runtime GC control is required for benchmark isolation, add a `enable_garbage_collection(bool)` accessor back to `aterm_pool` that writes `EnableGarbageCollection` through a non-`constexpr` flag, or add a new `aterm_configuration.h` override mechanism. At minimum, add a comment explaining that GC remains active and how that affects measurement.

---

### [SEV-2] `create_nested_function<N>` hardcodes arity 2 regardless of template parameter `N`

- **Location:** `libraries/atermpp/benchmark/benchmark_shared.h` (template used in `shared_lookup.cpp:30` and `shared_lookup.cpp:38`)

- **Why this is a problem:**
  ```cpp
  template<std::size_t N>
  aterm create_nested_function(...)
  {
    function_symbol f(function_name, N);   // arity = N
    aterm f_term(f, c_term, c_term);       // ALWAYS 2 children
    for (...) {
      make_term_appl(f_term, f, f_term, f_term);  // ALWAYS 2 children
    }
  }
  ```
  `f` is created with arity `N` but the term is always constructed with 2 children. For `N ≠ 2`, the arity of the function symbol and the number of arguments do not match. The atermpp pool stores the term using the function symbol's arity — a mismatch will cause an assertion failure or silent data corruption depending on build flags.

- **Evidence type:** static-proof (arity mismatch is unconditional for N ≠ 2)

- **Evidence:**
  - **Run command:**
    ```bash
    # add a call: create_nested_function<3>("f", "c", 5); to a test binary
    cmake --build build --target benchmark_target_atermpp_shared_lookup && \
    ./build/benchmarks/benchmark_target_atermpp_shared_lookup 1
    ```
  - **Observed result:** assertion failure or mismatched term storage (N=3 symbol stored in 2-arg pool slot)
  - **Expected result:** `N` children used consistently

- **Currently masked because:** Every caller in the suite uses `N=2` explicitly. The latent bug will surface if any new benchmark or test uses a different arity template instantiation.

- **Suggested fix:**
  ```cpp
  // Replace hardcoded 2-arg construction with N-arg construction:
  std::array<aterm, N> args;
  args.fill(c_term);
  make_term_appl(f_term, f, args.begin(), args.end());
  // and in the loop:
  args.fill(f_term);
  make_term_appl(f_term, f, args.begin(), args.end());
  ```

---

### [SEV-2] `std::stoi` → `static_cast<std::size_t>` allows negative thread counts to wrap

- **Location:** `libraries/atermpp/benchmark/shared_lookup.cpp:21–23`

- **Why this is a problem:**
  ```cpp
  number_of_threads = static_cast<std::size_t>(std::stoi(argv[1]));
  ```
  `std::stoi` returns `int`. A negative string argument (e.g., `"-1"`) returns `-1`, which wraps to `SIZE_MAX` (~1.8×10¹⁹ on 64-bit). `benchmark_threads` then allocates `SIZE_MAX - 1` `std::thread` objects, exhausting all memory before hitting OS thread limits.

  A string exceeding `INT_MAX` throws `std::out_of_range` uncaught, producing a non-descriptive `terminate`.

- **Evidence type:** static-proof

- **Evidence:**
  - **Run command:** `./build/.../benchmark_target_atermpp_shared_lookup -- -1`
  - **Observed result:** `std::bad_alloc` or OOM kill
  - **Expected result:** Validated input, at most one error message and graceful exit

- **Suggested fix:** Use `std::stoull`, validate the result is ≥ 1 and within a sane bound, or emit an error on non-positive input. This pattern is identical across all 6 benchmarks and should be fixed in a shared argument-parsing helper.

---

### [SEV-3] Integer truncation in `iterations / number_of_threads` silently produces zero-iteration benchmark

- **Location:** `libraries/atermpp/benchmark/shared_lookup.cpp:34`

- **Why this is a problem:**
  `iterations = 1000`. If `number_of_threads > 1000`, unsigned integer division produces `0`. The benchmark runs, reports a time, and the user assumes it measured something — it measured zero work. There is no assertion or warning.

- **Evidence type:** static-proof

- **Evidence:**
  - **Run command:** `./build/.../benchmark_target_atermpp_shared_lookup 1001`
  - **Observed result:** near-zero reported time with no error, silently invalid result
  - **Expected result:** diagnostic or iteration count clamped to at least 1

- **Suggested fix:** Assert or warn when `iterations / number_of_threads == 0`. Consider rounding up: `(iterations + number_of_threads - 1) / number_of_threads`.

---

### [SEV-3] Unnamed lambda parameter — thread ID is silently discarded; inconsistency with sibling benchmarks

- **Location:** `libraries/atermpp/benchmark/shared_lookup.cpp:31`

- **Why this is a problem:**
  ```cpp
  auto nested_function = [&](std::size_t) -> void   // unnamed parameter
  ```
  vs. `unique_lookup.cpp` and `unique_creation.cpp` which use `[&](int id)`. The unnamed parameter signals intent that thread ID is irrelevant (all threads work on the same terms), but the type is `std::size_t` while `benchmark_threads` passes `int`. The implicit narrowing-to-widening conversion is safe but unintentional-looking.

  More critically, if someone adds per-thread logic to this benchmark (as in `unique_lookup.cpp`), the unnamed parameter invites a silent bug where the wrong ID is used or ID is ignored.

- **Evidence type:** plausible-issue (no current failure, correctness risk on future change)

- **Suggested fix:** Either name the parameter `/*thread_id*/` to document the deliberate discard, or align the type to `int` to match `benchmark_threads`'s actual argument type. Mark with a comment: `// all threads contend on the same shared term — id unused`.

---

## Proposed Tests

- [ ] **benchmark_compile_test**: Configure with `MCRL2_ENABLE_BENCHMARKS=ON` and build all targets; expected failure on `enable_garbage_collection` call, confirms SEV-1. Command: `cmake --build build --target benchmark_target_atermpp_shared_lookup 2>&1 | grep "has no member"`
- [ ] **negative_thread_count_test**: Run any benchmark with `-- -1`; should reject input gracefully rather than OOM. Command: `./benchmark_target_atermpp_shared_lookup -- -1`
- [ ] **arity_mismatch_test**: Instantiate `create_nested_function<3>(...)` in a unit test; expected: assertion failure or sanitizer error when arity 3 function symbol is stored in a 2-argument term. Run with `-DCMAKE_BUILD_TYPE=Debug` to enable assertions.
- [ ] **zero_iteration_test**: Run with thread count > `iterations` (e.g., 2000 threads, 1000 iterations); verify that elapsed time is not near-zero without a warning. Command: `./benchmark_target_atermpp_shared_lookup 2000`

---

## Open Questions

1. **Intent of `enable_garbage_collection(false)` removal:** Was the intent to make GC behavior purely compile-time, or should a runtime override remain available for benchmarking isolation? If the former, all benchmark files need to be updated to document that GC is active during measurement.
2. **Behavioral validity with active GC:** With GC enabled, is the "shared lookup" benchmark still measuring what it intends (hash-table lookup contention), or does GC interference introduce noise that invalidates comparisons against `unique_lookup`?
3. **`create_nested_function<N>` ownership:** The template lives in `benchmark_shared.h`, which is not under test. Should this header be promoted to a proper tested utility, or kept as benchmark-only infrastructure? The N≠2 latent bug argues for at least a `static_assert(N == 2)` if 2 is the only supported arity.
4. **Timer placement in `benchmark_threads`:** The `stopwatch` starts before threads are launched and the first thread only reaches `f(0)` after the spawn loop. For a large `number_of_threads`, thread-launch overhead is included in the reported benchmark time. Is this intentional?

## Copilot Raw Stdout

```text
● Glob "/workspaces/mCRL2/libraries/atermpp/**/*.h"
  └ No matches found

● Find benchmark files and read shared header
  $ find /workspaces/mCRL2/libraries/atermpp/benchmark -type f | sort && cat
  /workspaces/mCRL2/libraries/atermpp/benchmark/benchmark_shared.h
  └ 109 lines...

● Read CMakeLists for benchmark targets
  $ cat /workspaces/mCRL2/libraries/atermpp/benchmark/CMakeLists.txt
  └ 30 lines...

● Look at sister benchmarks for comparison
  $ find /workspaces/mCRL2/libraries/atermpp/include -type f | sort | head -40 && cat
  /workspaces/mCRL2/libraries/atermpp/benchmark/shared_creation.cpp && cat
  /workspaces/mCRL2/libraries/atermpp/benchmark/unique_lookup.cpp
  └ 123 lines...

● Read remaining benchmark files for comparison
  $ cat /workspaces/mCRL2/libraries/atermpp/benchmark/shared_inspect.cpp && cat
  /workspaces/mCRL2/libraries/atermpp/benchmark/unique_creation.cpp && cat
  /workspaces/mCRL2/libraries/atermpp/benchmark/unique_inspect.cpp
  └ 170 lines...

● Check global term pool interface
  $ grep -n "enable_garbage_collection\|g_term_pool\|thread_aterm_pool\|make_term_appl\|create_neste
  d_function" /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h |
   head -40
  └ 4 lines...

● Read global aterm pool header
  $ cat /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h && cat
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h | head -80
  └ 146 lines...

● Look for enable_garbage_collection signature
  $ grep -n "enable_garbage_collection\|iterations\|number_of_threads\|std::stoi"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h 2>/dev/null | head
  -20
  └ 2 lines...

● Find enable_garbage_collection definition
  $ grep -rn "enable_garbage_collection" /workspaces/mCRL2/libraries/atermpp/include/ | head -20
  └ 4 lines...

● Find enable_garbage_collection method
  $ grep -n "enable_garbage_collection"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h | head -20 && grep
  -n "enable_garbage_collection"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h
  2>/dev/null | head -20
  └ 4 lines...

● Find all usages of enable_garbage_collection
  $ grep -rn "enable_garbage_collection" /workspaces/mCRL2/libraries/atermpp/ | grep -v ".pyc" |
  head -30
  └ 10 lines...

● Look for enable_garbage_collection as a method
  $ grep -n "enable_garbage_collection"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h && grep -n
  "enable_garbage_collection" /workspaces/mCRL2/libraries/atermpp/source/*.cpp 2>/dev/null | head
  -20
  └ 2 lines...

● Check if enable_garbage_collection is a method in aterm_pool
  $ grep -n "void enable_garbage\|enable_garbage_collection\|m_enable_garbage"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h | head -30
  └ 2 lines...

● Find aterm_pool class extent
  $ grep -n "class aterm_pool"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h && wc -l
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h
  └ 4 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h lines 82-259
  └ 178 lines read

● Search for enable_garbage_collection method definition
  $ grep -rn "void enable_garbage_collection\|enable_garbage_collection("
  /workspaces/mCRL2/libraries/atermpp/ | grep -v "\.pyc" | head -30
  └ 7 lines...

● Check EnableGarbageCollection usage in headers
  $ grep -rn "EnableGarbageCollection\|enable_garbage" /workspaces/mCRL2/libraries/atermpp/include/
  | head -20
  └ 9 lines...

● Read aterm configuration header
  $ cat /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_configuration.h
  └ 48 lines...

● Check for existing build directory
  $ find /workspaces/mCRL2 -name "build" -maxdepth 2 -type d 2>/dev/null && ls
  /workspaces/mCRL2/build 2>/dev/null | head -20
  └ 1 line...

✗ Find CMakeLists mentioning benchmarks
  $ find /workspaces/mCRL2 -name "CMakeLists.txt" | xargs grep -l "benchmark\|atermpp" 2>/dev/null |
   head -10
  Permission denied and could not request permission from user

● Check atermpp top-level CMakeLists for benchmark inclusion
  $ grep -rn "benchmark" /workspaces/mCRL2/libraries/atermpp/CMakeLists.txt 2>/dev/null | head -20
  && cat /workspaces/mCRL2/libraries/atermpp/CMakeLists.txt
  └ 20 lines...

● Search all non-benchmark source files for enable_garbage_collection
  $ grep -rn "enable_garbage_collection" /workspaces/mCRL2/ --include="*.h" --include="*.cpp"
  --include="*.hpp" | grep -v benchmark | head -20
  └ 4 lines...

● Check when enable_garbage_collection was removed
  $ cd /workspaces/mCRL2 && git log --oneline --follow --
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h | head -15
  └ 16 lines...

● Check git history for the relevant files
  $ cd /workspaces/mCRL2 && git log --oneline -10 --
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h && git log --oneline -10 --
  libraries/atermpp/benchmark/shared_lookup.cpp
  └ 15 lines...

● Check what was changed when disable automatic resizing was added
  $ cd /workspaces/mCRL2 && git show 7b43922d32 --stat | head -30
  └ 11 lines...

● Find commit that removed enable_garbage_collection method
  $ cd /workspaces/mCRL2 && git log --all --oneline --follow --
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h | grep -i
  "garbage\|gc\|collect\|enable" | head -10
  └ 4 lines...

● Check the commit that added test resizing and garbage collection functions
  $ cd /workspaces/mCRL2 && git show e4231250f9 --stat | head -20 && git show e4231250f9 --
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h | head -60
  └ 42 lines...

● Check the most recent commit that changed aterm_pool
  $ cd /workspaces/mCRL2 && git show 31c2b9f6ed --stat | head -20
  └ 21 lines...

● Check thread_aterm_pool for GC control methods
  $ grep -n "enable_garbage_collection\|test_garbage\|collect"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h | head -20
  └ 6 lines...

● See what changed in the most recent aterm_pool commit
  $ cd /workspaces/mCRL2 && git show 31c2b9f6ed --
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h | head -80
  └ 63 lines...

Now I have full context. Let me write the review.

---

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/benchmark/shared_lookup.cpp` plus transitive shared context in `benchmark_shared.h`, `aterm_pool.h`, and `aterm_configuration.h`
- **Risk level:** High — benchmark is broken and uncompilable when benchmarks are enabled; behavioral correctness also compromised
- **Overall verdict:** Do not ship. One regression introduced by commit `31c2b9f6ed` (Jun 30 2026) makes every benchmark in this directory fail to compile. Secondary findings affect measurement correctness and code safety.

---

## Findings (ordered by severity)

---

### [SEV-1] `enable_garbage_collection` method no longer exists — all benchmarks fail to compile

- **Location:** `libraries/atermpp/benchmark/shared_lookup.cpp:16` (and all 5 sibling benchmarks)
- **Why this is a problem:**
  Commit `31c2b9f6ed` ("Relaxed/removed some atomic counters", Jun 30 2026) changed `m_enable_garbage_collection` from `std::atomic<bool>` to `static constexpr bool` and **deleted** the public method `enable_garbage_collection(bool)` from `aterm_pool`. Every benchmark in the directory still calls `detail::g_term_pool().enable_garbage_collection(false)`. When `MCRL2_ENABLE_BENCHMARKS=ON` this is a hard compilation error.

- **Evidence type:** static-proof

- **Evidence:**
  - **Reproducer:** Configure with `-DMCRL2_ENABLE_BENCHMARKS=ON` and build benchmark targets.
  - **Run command:**
    ```bash
    cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug \
          -DMCRL2_ENABLE_GUI_TOOLS=OFF -DMCRL2_ENABLE_BENCHMARKS=ON
    cmake --build build --target benchmark_target_atermpp_shared_lookup
    ```
  - **Observed result:** `error: 'class atermpp::detail::aterm_pool' has no member named 'enable_garbage_collection'`
  - **Expected result:** Successful compilation

- **Behavioral impact:** GC is now compile-time-controlled via `constexpr static bool EnableGarbageCollection = true` (`aterm_configuration.h:19`). There is no runtime API to disable it. With GC enabled, the "lookup" benchmark silently degrades into a mixed creation-and-GC benchmark: terms built during one `create_nested_function` iteration can be reclaimed between iterations (the outer `aterm f` only roots the final term, not the intermediates), altering what is being measured.

- **C++20-specific note:** `static constexpr` is correct C++20 for compile-time constants; the original `std::atomic<bool>` member for a flag that was only written once at startup was itself unnecessary overhead.

- **Suggested fix:** Remove the `enable_garbage_collection(false)` call from all 6 benchmark files. If runtime GC control is required for benchmark isolation, add a `enable_garbage_collection(bool)` accessor back to `aterm_pool` that writes `EnableGarbageCollection` through a non-`constexpr` flag, or add a new `aterm_configuration.h` override mechanism. At minimum, add a comment explaining that GC remains active and how that affects measurement.

---

### [SEV-2] `create_nested_function<N>` hardcodes arity 2 regardless of template parameter `N`

- **Location:** `libraries/atermpp/benchmark/benchmark_shared.h` (template used in `shared_lookup.cpp:30` and `shared_lookup.cpp:38`)

- **Why this is a problem:**
  ```cpp
  template<std::size_t N>
  aterm create_nested_function(...)
  {
    function_symbol f(function_name, N);   // arity = N
    aterm f_term(f, c_term, c_term);       // ALWAYS 2 children
    for (...) {
      make_term_appl(f_term, f, f_term, f_term);  // ALWAYS 2 children
    }
  }
  ```
  `f` is created with arity `N` but the term is always constructed with 2 children. For `N ≠ 2`, the arity of the function symbol and the number of arguments do not match. The atermpp pool stores the term using the function symbol's arity — a mismatch will cause an assertion failure or silent data corruption depending on build flags.

- **Evidence type:** static-proof (arity mismatch is unconditional for N ≠ 2)

- **Evidence:**
  - **Run command:**
    ```bash
    # add a call: create_nested_function<3>("f", "c", 5); to a test binary
    cmake --build build --target benchmark_target_atermpp_shared_lookup && \
    ./build/benchmarks/benchmark_target_atermpp_shared_lookup 1
    ```
  - **Observed result:** assertion failure or mismatched term storage (N=3 symbol stored in 2-arg pool slot)
  - **Expected result:** `N` children used consistently

- **Currently masked because:** Every caller in the suite uses `N=2` explicitly. The latent bug will surface if any new benchmark or test uses a different arity template instantiation.

- **Suggested fix:**
  ```cpp
  // Replace hardcoded 2-arg construction with N-arg construction:
  std::array<aterm, N> args;
  args.fill(c_term);
  make_term_appl(f_term, f, args.begin(), args.end());
  // and in the loop:
  args.fill(f_term);
  make_term_appl(f_term, f, args.begin(), args.end());
  ```

---

### [SEV-2] `std::stoi` → `static_cast<std::size_t>` allows negative thread counts to wrap

- **Location:** `libraries/atermpp/benchmark/shared_lookup.cpp:21–23`

- **Why this is a problem:**
  ```cpp
  number_of_threads = static_cast<std::size_t>(std::stoi(argv[1]));
  ```
  `std::stoi` returns `int`. A negative string argument (e.g., `"-1"`) returns `-1`, which wraps to `SIZE_MAX` (~1.8×10¹⁹ on 64-bit). `benchmark_threads` then allocates `SIZE_MAX - 1` `std::thread` objects, exhausting all memory before hitting OS thread limits.

  A string exceeding `INT_MAX` throws `std::out_of_range` uncaught, producing a non-descriptive `terminate`.

- **Evidence type:** static-proof

- **Evidence:**
  - **Run command:** `./build/.../benchmark_target_atermpp_shared_lookup -- -1`
  - **Observed result:** `std::bad_alloc` or OOM kill
  - **Expected result:** Validated input, at most one error message and graceful exit

- **Suggested fix:** Use `std::stoull`, validate the result is ≥ 1 and within a sane bound, or emit an error on non-positive input. This pattern is identical across all 6 benchmarks and should be fixed in a shared argument-parsing helper.

---

### [SEV-3] Integer truncation in `iterations / number_of_threads` silently produces zero-iteration benchmark

- **Location:** `libraries/atermpp/benchmark/shared_lookup.cpp:34`

- **Why this is a problem:**
  `iterations = 1000`. If `number_of_threads > 1000`, unsigned integer division produces `0`. The benchmark runs, reports a time, and the user assumes it measured something — it measured zero work. There is no assertion or warning.

- **Evidence type:** static-proof

- **Evidence:**
  - **Run command:** `./build/.../benchmark_target_atermpp_shared_lookup 1001`
  - **Observed result:** near-zero reported time with no error, silently invalid result
  - **Expected result:** diagnostic or iteration count clamped to at least 1

- **Suggested fix:** Assert or warn when `iterations / number_of_threads == 0`. Consider rounding up: `(iterations + number_of_threads - 1) / number_of_threads`.

---

### [SEV-3] Unnamed lambda parameter — thread ID is silently discarded; inconsistency with sibling benchmarks

- **Location:** `libraries/atermpp/benchmark/shared_lookup.cpp:31`

- **Why this is a problem:**
  ```cpp
  auto nested_function = [&](std::size_t) -> void   // unnamed parameter
  ```
  vs. `unique_lookup.cpp` and `unique_creation.cpp` which use `[&](int id)`. The unnamed parameter signals intent that thread ID is irrelevant (all threads work on the same terms), but the type is `std::size_t` while `benchmark_threads` passes `int`. The implicit narrowing-to-widening conversion is safe but unintentional-looking.

  More critically, if someone adds per-thread logic to this benchmark (as in `unique_lookup.cpp`), the unnamed parameter invites a silent bug where the wrong ID is used or ID is ignored.

- **Evidence type:** plausible-issue (no current failure, correctness risk on future change)

- **Suggested fix:** Either name the parameter `/*thread_id*/` to document the deliberate discard, or align the type to `int` to match `benchmark_threads`'s actual argument type. Mark with a comment: `// all threads contend on the same shared term — id unused`.

---

## Proposed Tests

- [ ] **benchmark_compile_test**: Configure with `MCRL2_ENABLE_BENCHMARKS=ON` and build all targets; expected failure on `enable_garbage_collection` call, confirms SEV-1. Command: `cmake --build build --target benchmark_target_atermpp_shared_lookup 2>&1 | grep "has no member"`
- [ ] **negative_thread_count_test**: Run any benchmark with `-- -1`; should reject input gracefully rather than OOM. Command: `./benchmark_target_atermpp_shared_lookup -- -1`
- [ ] **arity_mismatch_test**: Instantiate `create_nested_function<3>(...)` in a unit test; expected: assertion failure or sanitizer error when arity 3 function symbol is stored in a 2-argument term. Run with `-DCMAKE_BUILD_TYPE=Debug` to enable assertions.
- [ ] **zero_iteration_test**: Run with thread count > `iterations` (e.g., 2000 threads, 1000 iterations); verify that elapsed time is not near-zero without a warning. Command: `./benchmark_target_atermpp_shared_lookup 2000`

---

## Open Questions

1. **Intent of `enable_garbage_collection(false)` removal:** Was the intent to make GC behavior purely compile-time, or should a runtime override remain available for benchmarking isolation? If the former, all benchmark files need to be updated to document that GC is active during measurement.
2. **Behavioral validity with active GC:** With GC enabled, is the "shared lookup" benchmark still measuring what it intends (hash-table lookup contention), or does GC interference introduce noise that invalidates comparisons against `unique_lookup`?
3. **`create_nested_function<N>` ownership:** The template lives in `benchmark_shared.h`, which is not under test. Should this header be promoted to a proper tested utility, or kept as benchmark-only infrastructure? The N≠2 latent bug argues for at least a `static_assert(N == 2)` if 2 is the only supported arity.
4. **Timer placement in `benchmark_threads`:** The `stopwatch` starts before threads are launched and the first thread only reaches `f(0)` after the spawn loop. For a large `number_of_threads`, thread-launch overhead is included in the reported benchmark time. Is this intentional?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         4m 5s
Total session time:     4m 24s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-sonnet-4.6       674.0k in, 13.5k out, 636.2k cached (Est. 1 Premium request)

```

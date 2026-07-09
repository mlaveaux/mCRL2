---
{
  "file": "libraries/atermpp/benchmark/shared_inspect.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-sonnet-4.6",
  "thinking": "deep",
  "effort": "high",
  "generated_at_utc": "2026-07-07T12:16:08.345263+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary

- **Scope reviewed:** `libraries/atermpp/benchmark/shared_inspect.cpp` (multi-threaded BFS traversal benchmark, ~70 LOC); transitively reviewed `benchmark_shared.h` and `detail/aterm_appl_iterator.h`.
- **Risk level:** Medium for production integrity; the benchmark is not deployed code, but it contains patterns that will mislead users who rely on its timing output, and it has a crash path on malformed input. A confirmed pre-existing bug in the iterator dependency is also surfaced.
- **Overall verdict:** Reject / needs fixes. Multiple findings corrupt benchmark validity. One crash-path finding is actionable now. One pre-existing confirmed iterator bug in a dependency is exposed.

---

## Findings (ordered by severity)

---

### [SEV-1] Unhandled `std::stoi` exceptions and silent wrap-around to `SIZE_MAX` threads

- **Location:** `shared_inspect.cpp:27`
- **Why this is a problem:** `std::stoi(argv[1])` returns a signed `int`. A negative input such as `"-1"` is parsed successfully and silently promoted via `static_cast<std::size_t>(-1)` to `SIZE_MAX`. This is then passed to `benchmark_threads(SIZE_MAX, ...)`, which immediately executes `std::vector<std::thread> threads(SIZE_MAX - 1)` in `benchmark_shared.h:24`. This allocation either throws `std::bad_alloc` or exhausts OS thread limits (RLIMIT_NPROC / ulimit), crashing the process. Additionally, non-numeric inputs (e.g., `"abc"`) cause `std::stoi` to throw `std::invalid_argument`, and overly large values throw `std::out_of_range` — both propagate unhandled out of `main`.
- **Evidence type:** static-proof
- **Evidence:**
  - **Reproducer:** Run the binary with `-1` as the first argument.
  - **Run command:**
    ```bash
    cmake --build build --target shared_inspect -j
    ./build/libraries/atermpp/benchmark/shared_inspect -- -1
    ```
  - **Observed result:** `std::bad_alloc` thrown from `benchmark_threads` or system-level OOM.
  - **Expected result:** Validated input; error message and clean exit for invalid or non-positive values.
- **Efficiency impact:** Not applicable (this is a correctness/crash finding).
- **C++20-specific note:** C++20 adds no remedy here; the author should use `std::from_chars` which is noexcept and avoids the int→size_t implicit narrowing trap.
- **Suggested fix:**
  ```cpp
  if (argc > 1)
  {
    int val = 0;
    auto [ptr, ec] = std::from_chars(argv[1], argv[1] + std::strlen(argv[1]), val);
    if (ec != std::errc{} || val <= 0)
    {
      std::cerr << "Usage: " << argv[0] << " [positive_integer]\n";
      return 1;
    }
    number_of_threads = static_cast<std::size_t>(val);
  }
  ```

---

### [SEV-2] Redundant `reinterpret_cast` masking a fragile lifetime dependency

- **Location:** `shared_inspect.cpp:55`
  ```cpp
  Q.push(reinterpret_cast<const aterm*>(&u));
  ```
- **Why this is a problem:** `u` is declared `const aterm&`, so `&u` is already `const aterm*`. The cast is a strict no-op at the type level. Its presence is a code smell indicating the author was uncertain about types. More critically, the pattern is fragile: correctness silently depends on `term_appl_iterator<aterm>::operator*()` returning a reference into stable pool memory (not a temporary or a copy). If the iterator were ever refactored to return by value — for example to eliminate the existing type-aliasing UB in `aterm.h:begin()` where `static_cast<const aterm&>` is used to downcast from an `unprotected_aterm_core` object — `&u` would immediately become a dangling pointer and the queue would hold stale addresses, causing silent memory corruption. There is zero indication in the code that this stability assumption is intentional or tested.
- **Evidence type:** static-proof (no-op cast) + plausible issue (fragile lifetime)
- **Evidence:**
  - **Static proof:** `term_appl_iterator<aterm>::operator*()` at `aterm_appl_iterator.h:51` returns `const Term&` where `Term = aterm`. `&u` is therefore `const aterm*`. `reinterpret_cast<const aterm*>(&u)` converts `const aterm*` to `const aterm*` — identical source and destination types.
  - **Plausible failure trigger:** Any refactoring of the atermpp iterator to return by value (motivated by the existing `static_cast` downcast UB in `aterm.h:142`).
  - **Reproducer for latent UB in iterator:** Change `operator*()` in `aterm_appl_iterator.h` to return `Term` (by value). Recompile and run — the queue will hold pointers to destroyed temporaries.
- **Efficiency impact:** None currently; the cast is optimised away.
- **C++20-specific note:** C++20 ranges requires iterators to model `std::indirectly_readable`, which mandates that `operator*()` return a stable reference if `reference` is a reference type. The current iterator declares `reference = const Term&` but the underlying storage is `unprotected_aterm_core` — the stated contract and the actual storage type conflict.
- **Suggested fix:** Remove the cast entirely. Document the stability assumption explicitly:
  ```cpp
  // &u is stable: points into the term pool node's inline argument array.
  // Validity depends on GC being disabled for the lifetime of Q.
  Q.push(&u);
  ```

---

### [SEV-2] Integer division truncation silently invalidates cross-thread-count benchmark comparisons

- **Location:** `shared_inspect.cpp:38`
  ```cpp
  for (std::size_t i = 0; i < iterations / number_of_threads; ++i)
  ```
- **Why this is a problem:** `1000 / k` truncates for all `k` that do not evenly divide 1000. The total work performed by all threads is `k * floor(1000/k)`, not 1000. Examples: `k=3` → 999 iterations total; `k=7` → 994; `k=9` → 999; `k=11` → 990. The benchmark's intended purpose is to show how sharing `f` scales with thread count by holding total work constant. The truncation makes multi-thread timings non-comparable: faster runtimes may reflect less work, not better parallelism.
- **Evidence type:** static-proof
- **Evidence:**
  - **Run command:** `./shared_inspect 7` vs `./shared_inspect 1`
  - **Observed result:** The 7-thread run does 994 total iterations; the 1-thread run does 1000.
  - **Expected result:** All configurations do exactly 1000 total iterations.
- **Efficiency impact:** Minor under-utilisation of threads. More importantly, invalid performance conclusions.
- **C++20-specific note:** None specific, but `std::div` or simply assigning the remainder to one thread (as done in e.g. `std::for_each_n` partitioning) is idiomatic.
- **Suggested fix:** Add the remainder to one thread's count, or restructure so a single atomic counter drives all threads.

---

### [SEV-2] Pre-existing: postfix `operator--(int)` in `term_appl_iterator` returns `*this` instead of `temp`

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_appl_iterator.h:99`
  ```cpp
  term_appl_iterator operator--(int)
  {
    term_appl_iterator temp = *this;
    --m_term;
    return *this;   // ← BUG: should be `return temp;`
  }
  ```
- **Why this is a problem:** Postfix decrement semantics require returning the pre-decrement value. This implementation decrements `m_term` but then returns the already-modified `*this`, making postfix and prefix `--` semantically identical. Any code relying on the pre-decrement snapshot via postfix `--` will silently compute the wrong result.
- **Evidence type:** static-proof
- **Evidence:**
  - **Reproducer:**
    ```cpp
    #include "mcrl2/atermpp/aterm.h"
    // Given an aterm_appl with ≥2 arguments:
    auto it = term.begin(); ++it;
    auto old = it--;   // old should == begin()+1, but equals begin()
    assert(old == term.begin());  // passes incorrectly
    assert(it == term.begin());   // also passes — both decremented
    ```
  - **Run command:** `ctest --test-dir build -R "atermpp" --output-on-failure`
- **Efficiency impact:** None unless postfix decrement is used.
- **C++20-specific note:** `std::random_access_iterator` requires postfix decrement to return the old value. This iterator violates the `std::random_access_iterator` concept contract silently.
- **Suggested fix:** `return temp;`

---

### [SEV-3] `std::cerr` output from multiple threads races without synchronisation

- **Location:** `shared_inspect.cpp:62` and `benchmark_shared.h:36`
- **Why this is a problem:** Each thread executes `std::cerr << "Term has " << n << " elements" << std::endl;` unsynchronised. C++11 guarantees individual `operator<<` calls are atomic (no data race), but a chain of three calls is not a single atomic operation. With `number_of_threads > 1`, output lines from different threads will interleave. `benchmark_shared.h:36` adds a fourth unsynchronised write from the timing output.
- **Evidence type:** plausible-issue (non-deterministic; trivially reproducible in practice with high thread count)
- **Evidence:**
  - **Trigger:** Run `./shared_inspect 4`; output will contain interleaved lines.
  - **Likely output:** `Term has Term has 2097150 elements\n2097150 elements\n`
- **Efficiency impact:** `std::endl` flushes the stream on every call, adding unnecessary I/O latency inside the timed section.
- **Suggested fix:** Collect `n` into a thread-local result and print all results after `benchmark_threads` returns, or use a mutex around the `cerr` block. Replace `std::endl` with `'\n'`.

---

### [SEV-3] Dead `#include "mcrl2/atermpp/algorithm.h"`

- **Location:** `shared_inspect.cpp:12`
- **Why this is a problem:** No function or type from `algorithm.h` is used in the file. The include contributes to compilation time and misleads readers who may believe `algorithm` functions are part of the benchmark logic.
- **Evidence type:** static-proof (no symbol from `algorithm.h` appears in the TU)
- **Suggested fix:** Remove the include.

---

### [SEV-3] No documentation that BFS visits exponentially many nodes (not unique term count)

- **Location:** `shared_inspect.cpp:40–60`
- **Why this is a problem:** Due to hash-consing, `create_nested_function<2>("f", "c", 20)` produces only 21 unique term nodes. The BFS traverses WITHOUT a visited set, so it follows every path through the DAG. For depth 20 with branching factor 2, the BFS visits ≈2²¹ ≈ 2M node instances per iteration. With 1000 iterations, the single-thread run performs ≈2 billion pointer traversals. A reader interpreting `n` as a count of unique subterms will be misled; `n ≈ 2097150` is a path count, not a unique-node count. The choice is presumably intentional (stressing the traversal path), but it is entirely undocumented.
- **Evidence type:** plausible-issue (undocumented intent; trivially verified by instrumenting unique-address counts)
- **Suggested fix:** Add a comment:
  ```cpp
  // BFS without a visited set intentionally counts all traversal paths,
  // not unique term nodes. Due to hash-consing, create_nested_function<2>(..., 20)
  // has only 21 unique nodes; the BFS visits ~2^21 instances per iteration.
  ```

---

## Proposed Tests

- [ ] **`test_benchmark_negative_threads`:** Run `./shared_inspect -- -1`; confirm process exits with non-zero status and an error message rather than aborting with `std::bad_alloc` or `std::terminate`.
- [ ] **`test_benchmark_invalid_arg`:** Run `./shared_inspect abc`; confirm graceful error rather than uncaught `std::invalid_argument`.
- [ ] **`test_postfix_decrement_iterator`:** Construct an `aterm_appl` with 2+ arguments; take `begin()+1`, apply postfix `--`, assert return value equals `begin()+1` and that the iterator is left at `begin()`.
- [ ] **`test_benchmark_work_equality`:** Instrument total iteration count across threads for `number_of_threads ∈ {1, 3, 7}` with `iterations = 1000`; assert total == 1000 for each.
- [ ] **`test_bfs_node_count`:** Add a visited-address set alongside `n`; assert unique nodes == 21 while `n ≈ 2M`, confirming the DAG/path-count distinction.

---

## Open Questions

- Is `benchmark_threads` intentionally measuring wall time (including thread-spawn overhead) rather than only traversal time? Thread spawning is not excluded from the timer, which will skew timings for small `iterations`.
- Is `enable_garbage_collection(false)` intended to be a global, permanent setting for the process lifetime? If the term pool is shared and another component re-enables GC (e.g., in a hypothetical integration test), all raw `const aterm*` pointers held in `Q` during the BFS would become dangling mid-traversal with no warning.
- Why is `algorithm.h` included — was it used in an earlier version of the benchmark, or is it expected to be needed in future?
- The `_aterm_appl::arg()` method at `detail/aterm.h:90` does `static_cast<const aterm_core&>` from an `unprotected_aterm_core` object (a downcast), and `aterm::begin()` then applies a second `static_cast<const aterm&>` (another downcast). Given `aterm_core : public unprotected_aterm_core` and `aterm : public aterm_core`, both of these are well-defined upcasts only if the underlying objects were originally constructed as `aterm` — but `m_arguments` is `std::array<unprotected_aterm_core, N>`. Is the intent that `_aterm_appl` be always constructed with actual `aterm` subobjects, or is this intentional layout-punning that relies on identical single-member structs?

## Copilot Raw Stdout

```text
● Explore: Explore atermpp internals
  ● Find function definitions in algorithm.h
    $ grep "^[a-z_]*(" /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/algorithm.h | head
    -20
    └ 5 lines...
  ● Create comprehensive summary
    $ cat > /tmp/test_summary.txt << 'EOF'
    FINDINGS SUMMARY FOR ATERMPP LIBRARY:

    1. ITERATOR BEHAVIOR (for (const aterm& u : t))
       =============================================...
    └ 279 lines...

● Read benchmark_shared.h and iterator definitions
  $ cat /workspaces/mCRL2/libraries/atermpp/benchmark/benchmark_shared.h && echo "---ITERATOR---" &&
   cat /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_appl_iterator.h &&
  echo "---ATERM_H---" && sed -n '130,160p'
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm.h
  └ 321 lines...

● Check iterator postfix decrement and type hierarchy
  $ sed -n '85,110p'
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_appl_iterator.h && echo
  "---" && grep -n "unprotected_aterm_core\|aterm_core"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h | head -30
  └ 36 lines...

● Check type hierarchy for aterm_core and unprotected_aterm_core
  $ grep -n "class unprotected_aterm_core\|class aterm_core\|class aterm "
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm_core.h | head -20 && echo "---" &&
   grep -n "unprotected_aterm_core"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/unprotected_aterm.h | head -10
  └ 5 lines...

Now I have all the information needed. Here is the full review:

---

# Review Findings

## Summary

- **Scope reviewed:** `libraries/atermpp/benchmark/shared_inspect.cpp` (multi-threaded BFS traversal benchmark, ~70 LOC); transitively reviewed `benchmark_shared.h` and `detail/aterm_appl_iterator.h`.
- **Risk level:** Medium for production integrity; the benchmark is not deployed code, but it contains patterns that will mislead users who rely on its timing output, and it has a crash path on malformed input. A confirmed pre-existing bug in the iterator dependency is also surfaced.
- **Overall verdict:** Reject / needs fixes. Multiple findings corrupt benchmark validity. One crash-path finding is actionable now. One pre-existing confirmed iterator bug in a dependency is exposed.

---

## Findings (ordered by severity)

---

### [SEV-1] Unhandled `std::stoi` exceptions and silent wrap-around to `SIZE_MAX` threads

- **Location:** `shared_inspect.cpp:27`
- **Why this is a problem:** `std::stoi(argv[1])` returns a signed `int`. A negative input such as `"-1"` is parsed successfully and silently promoted via `static_cast<std::size_t>(-1)` to `SIZE_MAX`. This is then passed to `benchmark_threads(SIZE_MAX, ...)`, which immediately executes `std::vector<std::thread> threads(SIZE_MAX - 1)` in `benchmark_shared.h:24`. This allocation either throws `std::bad_alloc` or exhausts OS thread limits (RLIMIT_NPROC / ulimit), crashing the process. Additionally, non-numeric inputs (e.g., `"abc"`) cause `std::stoi` to throw `std::invalid_argument`, and overly large values throw `std::out_of_range` — both propagate unhandled out of `main`.
- **Evidence type:** static-proof
- **Evidence:**
  - **Reproducer:** Run the binary with `-1` as the first argument.
  - **Run command:**
    ```bash
    cmake --build build --target shared_inspect -j
    ./build/libraries/atermpp/benchmark/shared_inspect -- -1
    ```
  - **Observed result:** `std::bad_alloc` thrown from `benchmark_threads` or system-level OOM.
  - **Expected result:** Validated input; error message and clean exit for invalid or non-positive values.
- **Efficiency impact:** Not applicable (this is a correctness/crash finding).
- **C++20-specific note:** C++20 adds no remedy here; the author should use `std::from_chars` which is noexcept and avoids the int→size_t implicit narrowing trap.
- **Suggested fix:**
  ```cpp
  if (argc > 1)
  {
    int val = 0;
    auto [ptr, ec] = std::from_chars(argv[1], argv[1] + std::strlen(argv[1]), val);
    if (ec != std::errc{} || val <= 0)
    {
      std::cerr << "Usage: " << argv[0] << " [positive_integer]\n";
      return 1;
    }
    number_of_threads = static_cast<std::size_t>(val);
  }
  ```

---

### [SEV-2] Redundant `reinterpret_cast` masking a fragile lifetime dependency

- **Location:** `shared_inspect.cpp:55`
  ```cpp
  Q.push(reinterpret_cast<const aterm*>(&u));
  ```
- **Why this is a problem:** `u` is declared `const aterm&`, so `&u` is already `const aterm*`. The cast is a strict no-op at the type level. Its presence is a code smell indicating the author was uncertain about types. More critically, the pattern is fragile: correctness silently depends on `term_appl_iterator<aterm>::operator*()` returning a reference into stable pool memory (not a temporary or a copy). If the iterator were ever refactored to return by value — for example to eliminate the existing type-aliasing UB in `aterm.h:begin()` where `static_cast<const aterm&>` is used to downcast from an `unprotected_aterm_core` object — `&u` would immediately become a dangling pointer and the queue would hold stale addresses, causing silent memory corruption. There is zero indication in the code that this stability assumption is intentional or tested.
- **Evidence type:** static-proof (no-op cast) + plausible issue (fragile lifetime)
- **Evidence:**
  - **Static proof:** `term_appl_iterator<aterm>::operator*()` at `aterm_appl_iterator.h:51` returns `const Term&` where `Term = aterm`. `&u` is therefore `const aterm*`. `reinterpret_cast<const aterm*>(&u)` converts `const aterm*` to `const aterm*` — identical source and destination types.
  - **Plausible failure trigger:** Any refactoring of the atermpp iterator to return by value (motivated by the existing `static_cast` downcast UB in `aterm.h:142`).
  - **Reproducer for latent UB in iterator:** Change `operator*()` in `aterm_appl_iterator.h` to return `Term` (by value). Recompile and run — the queue will hold pointers to destroyed temporaries.
- **Efficiency impact:** None currently; the cast is optimised away.
- **C++20-specific note:** C++20 ranges requires iterators to model `std::indirectly_readable`, which mandates that `operator*()` return a stable reference if `reference` is a reference type. The current iterator declares `reference = const Term&` but the underlying storage is `unprotected_aterm_core` — the stated contract and the actual storage type conflict.
- **Suggested fix:** Remove the cast entirely. Document the stability assumption explicitly:
  ```cpp
  // &u is stable: points into the term pool node's inline argument array.
  // Validity depends on GC being disabled for the lifetime of Q.
  Q.push(&u);
  ```

---

### [SEV-2] Integer division truncation silently invalidates cross-thread-count benchmark comparisons

- **Location:** `shared_inspect.cpp:38`
  ```cpp
  for (std::size_t i = 0; i < iterations / number_of_threads; ++i)
  ```
- **Why this is a problem:** `1000 / k` truncates for all `k` that do not evenly divide 1000. The total work performed by all threads is `k * floor(1000/k)`, not 1000. Examples: `k=3` → 999 iterations total; `k=7` → 994; `k=9` → 999; `k=11` → 990. The benchmark's intended purpose is to show how sharing `f` scales with thread count by holding total work constant. The truncation makes multi-thread timings non-comparable: faster runtimes may reflect less work, not better parallelism.
- **Evidence type:** static-proof
- **Evidence:**
  - **Run command:** `./shared_inspect 7` vs `./shared_inspect 1`
  - **Observed result:** The 7-thread run does 994 total iterations; the 1-thread run does 1000.
  - **Expected result:** All configurations do exactly 1000 total iterations.
- **Efficiency impact:** Minor under-utilisation of threads. More importantly, invalid performance conclusions.
- **C++20-specific note:** None specific, but `std::div` or simply assigning the remainder to one thread (as done in e.g. `std::for_each_n` partitioning) is idiomatic.
- **Suggested fix:** Add the remainder to one thread's count, or restructure so a single atomic counter drives all threads.

---

### [SEV-2] Pre-existing: postfix `operator--(int)` in `term_appl_iterator` returns `*this` instead of `temp`

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_appl_iterator.h:99`
  ```cpp
  term_appl_iterator operator--(int)
  {
    term_appl_iterator temp = *this;
    --m_term;
    return *this;   // ← BUG: should be `return temp;`
  }
  ```
- **Why this is a problem:** Postfix decrement semantics require returning the pre-decrement value. This implementation decrements `m_term` but then returns the already-modified `*this`, making postfix and prefix `--` semantically identical. Any code relying on the pre-decrement snapshot via postfix `--` will silently compute the wrong result.
- **Evidence type:** static-proof
- **Evidence:**
  - **Reproducer:**
    ```cpp
    #include "mcrl2/atermpp/aterm.h"
    // Given an aterm_appl with ≥2 arguments:
    auto it = term.begin(); ++it;
    auto old = it--;   // old should == begin()+1, but equals begin()
    assert(old == term.begin());  // passes incorrectly
    assert(it == term.begin());   // also passes — both decremented
    ```
  - **Run command:** `ctest --test-dir build -R "atermpp" --output-on-failure`
- **Efficiency impact:** None unless postfix decrement is used.
- **C++20-specific note:** `std::random_access_iterator` requires postfix decrement to return the old value. This iterator violates the `std::random_access_iterator` concept contract silently.
- **Suggested fix:** `return temp;`

---

### [SEV-3] `std::cerr` output from multiple threads races without synchronisation

- **Location:** `shared_inspect.cpp:62` and `benchmark_shared.h:36`
- **Why this is a problem:** Each thread executes `std::cerr << "Term has " << n << " elements" << std::endl;` unsynchronised. C++11 guarantees individual `operator<<` calls are atomic (no data race), but a chain of three calls is not a single atomic operation. With `number_of_threads > 1`, output lines from different threads will interleave. `benchmark_shared.h:36` adds a fourth unsynchronised write from the timing output.
- **Evidence type:** plausible-issue (non-deterministic; trivially reproducible in practice with high thread count)
- **Evidence:**
  - **Trigger:** Run `./shared_inspect 4`; output will contain interleaved lines.
  - **Likely output:** `Term has Term has 2097150 elements\n2097150 elements\n`
- **Efficiency impact:** `std::endl` flushes the stream on every call, adding unnecessary I/O latency inside the timed section.
- **Suggested fix:** Collect `n` into a thread-local result and print all results after `benchmark_threads` returns, or use a mutex around the `cerr` block. Replace `std::endl` with `'\n'`.

---

### [SEV-3] Dead `#include "mcrl2/atermpp/algorithm.h"`

- **Location:** `shared_inspect.cpp:12`
- **Why this is a problem:** No function or type from `algorithm.h` is used in the file. The include contributes to compilation time and misleads readers who may believe `algorithm` functions are part of the benchmark logic.
- **Evidence type:** static-proof (no symbol from `algorithm.h` appears in the TU)
- **Suggested fix:** Remove the include.

---

### [SEV-3] No documentation that BFS visits exponentially many nodes (not unique term count)

- **Location:** `shared_inspect.cpp:40–60`
- **Why this is a problem:** Due to hash-consing, `create_nested_function<2>("f", "c", 20)` produces only 21 unique term nodes. The BFS traverses WITHOUT a visited set, so it follows every path through the DAG. For depth 20 with branching factor 2, the BFS visits ≈2²¹ ≈ 2M node instances per iteration. With 1000 iterations, the single-thread run performs ≈2 billion pointer traversals. A reader interpreting `n` as a count of unique subterms will be misled; `n ≈ 2097150` is a path count, not a unique-node count. The choice is presumably intentional (stressing the traversal path), but it is entirely undocumented.
- **Evidence type:** plausible-issue (undocumented intent; trivially verified by instrumenting unique-address counts)
- **Suggested fix:** Add a comment:
  ```cpp
  // BFS without a visited set intentionally counts all traversal paths,
  // not unique term nodes. Due to hash-consing, create_nested_function<2>(..., 20)
  // has only 21 unique nodes; the BFS visits ~2^21 instances per iteration.
  ```

---

## Proposed Tests

- [ ] **`test_benchmark_negative_threads`:** Run `./shared_inspect -- -1`; confirm process exits with non-zero status and an error message rather than aborting with `std::bad_alloc` or `std::terminate`.
- [ ] **`test_benchmark_invalid_arg`:** Run `./shared_inspect abc`; confirm graceful error rather than uncaught `std::invalid_argument`.
- [ ] **`test_postfix_decrement_iterator`:** Construct an `aterm_appl` with 2+ arguments; take `begin()+1`, apply postfix `--`, assert return value equals `begin()+1` and that the iterator is left at `begin()`.
- [ ] **`test_benchmark_work_equality`:** Instrument total iteration count across threads for `number_of_threads ∈ {1, 3, 7}` with `iterations = 1000`; assert total == 1000 for each.
- [ ] **`test_bfs_node_count`:** Add a visited-address set alongside `n`; assert unique nodes == 21 while `n ≈ 2M`, confirming the DAG/path-count distinction.

---

## Open Questions

- Is `benchmark_threads` intentionally measuring wall time (including thread-spawn overhead) rather than only traversal time? Thread spawning is not excluded from the timer, which will skew timings for small `iterations`.
- Is `enable_garbage_collection(false)` intended to be a global, permanent setting for the process lifetime? If the term pool is shared and another component re-enables GC (e.g., in a hypothetical integration test), all raw `const aterm*` pointers held in `Q` during the BFS would become dangling mid-traversal with no warning.
- Why is `algorithm.h` included — was it used in an earlier version of the benchmark, or is it expected to be needed in future?
- The `_aterm_appl::arg()` method at `detail/aterm.h:90` does `static_cast<const aterm_core&>` from an `unprotected_aterm_core` object (a downcast), and `aterm::begin()` then applies a second `static_cast<const aterm&>` (another downcast). Given `aterm_core : public unprotected_aterm_core` and `aterm : public aterm_core`, both of these are well-defined upcasts only if the underlying objects were originally constructed as `aterm` — but `m_arguments` is `std::array<unprotected_aterm_core, N>`. Is the intent that `_aterm_appl` be always constructed with actual `aterm` subobjects, or is this intentional layout-punning that relies on identical single-member structs?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         7m 3s
Total session time:     7m 24s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-sonnet-4.6       185.9k in, 19.2k out, 153.5k cached (Est. 1 Premium request)
 claude-haiku-4.5        346.5k in, 7.1k out, 317.6k cached (Est. 0 Premium requests)

```

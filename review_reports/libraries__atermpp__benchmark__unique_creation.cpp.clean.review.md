---
{
  "file": "libraries/atermpp/benchmark/unique_creation.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-sonnet-4.6",
  "thinking": "deep",
  "effort": "high",
  "generated_at_utc": "2026-07-07T12:23:56.065850+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/benchmark/unique_creation.cpp` with full transitive context from `benchmark_shared.h`, `shared_creation.cpp`, `unique_lookup.cpp`, `unique_inspect.cpp`
- **Risk level:** Medium (benchmark-only code, not production, but contains confirmed UB)
- **Overall verdict:** Reject / Revise. Three confirmed unsafe paths, one misleading design discrepancy. The code would crash or exhibit undefined behaviour under trivially reachable inputs.

---

## Findings (ordered by severity)

---

### [SEV-1] Division-by-zero UB and `size_t` underflow when thread count is zero

- **Location:** `unique_creation.cpp:22,30` and `benchmark_shared.h:23`
- **Why this is a problem:**  
  `std::stoi(argv[1])` legitimately returns `0` (input `"0"`).  
  `static_cast<std::size_t>(0)` sets `number_of_threads = 0`.  
  Two separate UB paths then fire:
  1. `size / number_of_threads` (line 30) — unsigned integer division by zero is **undefined behaviour** in C++ (standard [expr.mul]/4 states the behaviour is undefined when the right operand is zero regardless of signedness).
  2. Inside `benchmark_threads`, `std::vector<std::thread> threads(number_of_threads - 1)` (benchmark_shared.h:23) — `0 - 1` on a `std::size_t` wraps to `SIZE_MAX` (well-defined wrap, but then the vector constructor attempts to allocate ~1.8×10¹⁹ `std::thread` objects, which will throw `std::bad_alloc` or `std::length_error` before even reaching the UB division).

- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer: run the built binary with argument `"0"`
  - Run command: `./unique_creation 0`
  - Observed result: UB or `std::bad_alloc`/`std::length_error` (platform-dependent crash or termination)
  - Expected result: meaningful error message and clean exit
- **Efficiency impact:** N/A (crash path)
- **C++20-specific note:** C++20 does not change UB rules for integer division by zero. Using `[[likely]]`/`[[unlikely]]` on the guard would be appropriate once the guard is added.
- **Suggested fix:**
  ```cpp
  if (argc > 1)
  {
    const int n = std::stoi(argv[1]);
    if (n <= 0) { std::cerr << "Thread count must be >= 1\n"; return 1; }
    number_of_threads = static_cast<std::size_t>(n);
  }
  ```

---

### [SEV-1] Negative thread count silently wraps to `SIZE_MAX`, spawning catastrophic thread count

- **Location:** `unique_creation.cpp:22`
- **Why this is a problem:**  
  `std::stoi(argv[1])` can return any value in `[INT_MIN, INT_MAX]`.  
  For negative inputs (e.g., `"-1"`), `static_cast<std::size_t>(-1)` yields `SIZE_MAX` (≈1.8×10¹⁹ on 64-bit).  
  This is then passed to `benchmark_threads`, which attempts `std::vector<std::thread>(SIZE_MAX - 1)` — guaranteed `std::bad_alloc` or OS resource exhaustion before any graceful recovery. No diagnostic message is printed.

- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer: `./unique_creation -1`
  - Observed result: `std::bad_alloc` or OS-level resource exhaustion
  - Expected result: validated error message and graceful exit
- **Efficiency impact:** N/A (crash path)
- **C++20-specific note:** None specific; this is a long-standing pitfall with signed-to-unsigned cast.
- **Suggested fix:** Same guard as SEV-1 above (`n <= 0` check before cast).

---

### [SEV-2] `std::stoi` throws on malformed or out-of-range input; no exception handler

- **Location:** `unique_creation.cpp:22`
- **Why this is a problem:**  
  `std::stoi` throws `std::invalid_argument` on non-numeric input and `std::out_of_range` when the value exceeds `INT_MAX`. Both are uncaught, terminating the process with a noisy `terminate()` call and no user-facing message.  
  Every other benchmark in the directory (`shared_creation.cpp:22`, `unique_lookup.cpp:22`, etc.) has the same flaw — but this file is the reviewed target.

- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer: `./unique_creation abc` or `./unique_creation 99999999999999`
  - Observed result: `terminate called after throwing an instance of 'std::invalid_argument'`
  - Expected result: clean diagnostic and non-zero exit
- **Efficiency impact:** N/A
- **C++20-specific note:** `std::from_chars` (C++17/20) avoids exceptions entirely and is preferable for CLI parsing:
  ```cpp
  std::size_t n = 1;
  if (argc > 1) {
    auto [ptr, ec] = std::from_chars(argv[1], argv[1] + std::strlen(argv[1]), n);
    if (ec != std::errc{} || n == 0) { std::cerr << "Invalid thread count\n"; return 1; }
    number_of_threads = n;
  }
  ```
- **Suggested fix:** Replace `std::stoi` with `std::from_chars` and add range validation.

---

### [SEV-2] `benchmark_threads` internal `id` counter uses `int`; overflows for large thread counts

- **Location:** `benchmark_shared.h:24`
- **Why this is a problem:**  
  ```cpp
  int id = 1;
  for (auto& thread : threads) { thread = std::thread(f, id); ++id; }
  ```
  `id` is `int` but `threads.size()` is `std::size_t`. For `number_of_threads > INT_MAX`, `++id` overflows, which is signed integer overflow — **undefined behaviour** in C++. Even at practical scales, `id` silently wraps at 2^31 on 32-bit or is platform-dependent on 64-bit.  
  The lambda in `unique_creation.cpp` uses the `id` as a leaf name discriminator: collisions produce incorrect pool-sharing behaviour that invalidates the "unique" benchmark semantics.

- **Evidence type:** plausible-issue (thread counts this large are impractical, but the type mismatch propagates to correctness of thread identity)
- **Evidence:**
  - Mechanism: signed integer overflow is UB; id collision corrupts benchmark uniqueness.
  - Triggering scenario: theoretically `number_of_threads > INT_MAX`; practically, any `number_of_threads > 2^31`.
  - Test plan: Pass `number_of_threads = INT_MAX + 2` on a 64-bit system, observe whether id wraps and two threads get the same leaf name.
- **Efficiency impact:** Incorrect benchmark measurement if id collision occurs.
- **C++20-specific note:** `benchmark_threads` should use `std::size_t id` consistently.
- **Suggested fix:** Change `int id = 1;` to `std::size_t id = 1;` in `benchmark_shared.h`; update all lambda signatures to `std::size_t`.

---

### [SEV-3] Integer truncation in `size / number_of_threads` makes total work non-reproducible

- **Location:** `unique_creation.cpp:30`
- **Why this is a problem:**  
  When `size` (400000) is not evenly divisible by `number_of_threads`, integer truncation silently reduces depth: e.g., `400000 / 3 = 133333` (loses 1 layer). The total number of unique pool entries created across all threads is `N × ⌊400000/N⌋`, not `400000`. For a benchmark tool this means successive runs with different thread counts are not measuring equivalent amounts of work, making scaling measurements unreliable.

- **Evidence type:** static-proof
- **Evidence:**
  - `400000 / 3 * 3 = 399999 ≠ 400000`
  - The benchmark's `std::cerr` line only prints elapsed time — no work-item count — so the discrepancy is silent.
- **Efficiency impact:** Benchmark measurements are not directly comparable across thread counts.
- **C++20-specific note:** None.
- **Suggested fix:**  
  ```cpp
  // Distribute work evenly, ensuring total depth ≥ size
  const std::size_t depth_per_thread = (size + number_of_threads - 1) / number_of_threads;
  ```
  Or document that the division is intentional (with a comment).

---

### [SEV-3] Lambda captures `size` and `number_of_threads` by mutable reference; no `const` protection

- **Location:** `unique_creation.cpp:28–31`
- **Why this is a problem:**  
  The lambda `[&](int id) -> void` captures both `size` and `number_of_threads` by mutable reference. These variables are only read inside the lambda, but there is nothing preventing accidental mutation if this benchmark is extended. When multiple threads simultaneously read these through references, the absence of `const` makes the intent unclear and the code fragile.

- **Evidence type:** plausible-issue  
- **Efficiency impact:** None in current code, but data-race risk if a future edit writes to either variable after thread creation.
- **C++20-specific note:** C++20 lambda `[size, number_of_threads](int id) -> void` (value capture) is the correct idiom here since both variables are of value type (`std::size_t`) and are not modified.  
  Alternatively: `[size = size, not = number_of_threads]` or capture the computed `depth` directly:
  ```cpp
  const std::size_t depth = size / number_of_threads;
  auto nested_function = [depth](int id) -> void { ... };
  ```
- **Suggested fix:** Capture `size` and `number_of_threads` by value, or precompute `depth` as `const` and capture only that.

---

### [SEV-3] "Unique" per-thread semantics are undocumented; design intent is non-obvious

- **Location:** `unique_creation.cpp:28–31`
- **Why this is a problem:**  
  The file is named `unique_creation.cpp` and uses `std::to_string(id)` as the leaf name — making each thread's term tree disjoint from other threads' trees (no shared pool nodes). This is the deliberate contrast to `shared_creation.cpp` (leaf `"c"` is shared across threads). However, there is no comment explaining this design, so future maintainers may silently break the benchmark's distinctiveness by changing the leaf name.

- **Evidence type:** plausible-issue
- **Suggested fix:** Add a comment:
  ```cpp
  // Each thread uses a distinct leaf name so no pool entries are shared between
  // threads. This measures creation throughput without inter-thread contention on
  // the term pool.
  ```

---

## Proposed Tests

- [ ] `test_zero_threads`: Run `./unique_creation 0`; expected: graceful error message, exit code != 0, no crash/UB.
- [ ] `test_negative_threads`: Run `./unique_creation -1`; expected: graceful error message, exit code != 0, no resource exhaustion.
- [ ] `test_invalid_threads_arg`: Run `./unique_creation abc`; expected: no unhandled exception, graceful exit with message.
- [ ] `test_work_invariance`: Run `./unique_creation 1`, `./unique_creation 2`, `./unique_creation 4`; verify total unique pool terms created (e.g., via pool size query before/after) equals `size` in all cases — catches the truncation bug.
- [ ] `test_unique_vs_shared_no_collision`: After running with 2 threads, verify terms from thread 0 (`leaf="0"`) and thread 1 (`leaf="1"`) are provably disjoint in the pool; expected: no node shared between the two subtrees.

---

## Open Questions

1. **Why does `unique_creation.cpp` use `size / number_of_threads` (total work constant) while `shared_creation.cpp` uses `size` directly (total work scales with threads)?** The inconsistency suggests different authors or different unspoken intents — but there is no documentation. If the goal is weak-scaling measurement, that should be stated explicitly.
2. **Why is `benchmark_threads` internal `id` typed as `int` when `number_of_threads` is `std::size_t`?** This is an API inconsistency across all benchmark files that should be resolved globally in `benchmark_shared.h`.
3. **Is GC being disabled (`enable_garbage_collection(false)`) intentional for benchmarking creation cost specifically, or was it copied from another benchmark without consideration?** Disabling GC means pool memory grows unboundedly for the duration of the benchmark, which could invalidate memory-pressure performance conclusions.
4. **Is the timer in `benchmark_threads` intended to include thread-spawn overhead?** It starts before `std::thread` construction — this may inflate results for large `number_of_threads`, conflating thread creation latency with term creation latency.

## Copilot Raw Stdout

```text
● List directory libraries/atermpp/benchmark
  └ 8 files found

● Glob "libraries/atermpp/benchmark/*.{h,hpp,cpp}"
  └ 7 files found

● Read libraries/atermpp/benchmark/shared_creation.cpp
  └ 36 lines read

● Read libraries/atermpp/benchmark/unique_lookup.cpp
  └ 46 lines read

● Read libraries/atermpp/benchmark/benchmark_shared.h
  └ 100 lines read

● Read libraries/atermpp/benchmark/unique_inspect.cpp
  └ 66 lines read

● Grep "g_term_pool|enable_garbage_collection" in **/*.{h,hpp,cpp}
(/workspaces/mCRL2/libraries/atermpp)
  └ 16 files found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h
  └ 65 lines read

● Read libraries/atermpp/benchmark/shared_inspect.cpp
  └ 67 lines read

● Read libraries/atermpp/benchmark/shared_lookup.cpp
  └ 42 lines read

● Grep "stoi|number_of_threads|size_t"
(/workspaces/mCRL2/libraries/atermpp/benchmark/unique_creation.cpp)
  └ 5 lines found

Now I have full context. Let me write the review.

---

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/benchmark/unique_creation.cpp` with full transitive context from `benchmark_shared.h`, `shared_creation.cpp`, `unique_lookup.cpp`, `unique_inspect.cpp`
- **Risk level:** Medium (benchmark-only code, not production, but contains confirmed UB)
- **Overall verdict:** Reject / Revise. Three confirmed unsafe paths, one misleading design discrepancy. The code would crash or exhibit undefined behaviour under trivially reachable inputs.

---

## Findings (ordered by severity)

---

### [SEV-1] Division-by-zero UB and `size_t` underflow when thread count is zero

- **Location:** `unique_creation.cpp:22,30` and `benchmark_shared.h:23`
- **Why this is a problem:**  
  `std::stoi(argv[1])` legitimately returns `0` (input `"0"`).  
  `static_cast<std::size_t>(0)` sets `number_of_threads = 0`.  
  Two separate UB paths then fire:
  1. `size / number_of_threads` (line 30) — unsigned integer division by zero is **undefined behaviour** in C++ (standard [expr.mul]/4 states the behaviour is undefined when the right operand is zero regardless of signedness).
  2. Inside `benchmark_threads`, `std::vector<std::thread> threads(number_of_threads - 1)` (benchmark_shared.h:23) — `0 - 1` on a `std::size_t` wraps to `SIZE_MAX` (well-defined wrap, but then the vector constructor attempts to allocate ~1.8×10¹⁹ `std::thread` objects, which will throw `std::bad_alloc` or `std::length_error` before even reaching the UB division).

- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer: run the built binary with argument `"0"`
  - Run command: `./unique_creation 0`
  - Observed result: UB or `std::bad_alloc`/`std::length_error` (platform-dependent crash or termination)
  - Expected result: meaningful error message and clean exit
- **Efficiency impact:** N/A (crash path)
- **C++20-specific note:** C++20 does not change UB rules for integer division by zero. Using `[[likely]]`/`[[unlikely]]` on the guard would be appropriate once the guard is added.
- **Suggested fix:**
  ```cpp
  if (argc > 1)
  {
    const int n = std::stoi(argv[1]);
    if (n <= 0) { std::cerr << "Thread count must be >= 1\n"; return 1; }
    number_of_threads = static_cast<std::size_t>(n);
  }
  ```

---

### [SEV-1] Negative thread count silently wraps to `SIZE_MAX`, spawning catastrophic thread count

- **Location:** `unique_creation.cpp:22`
- **Why this is a problem:**  
  `std::stoi(argv[1])` can return any value in `[INT_MIN, INT_MAX]`.  
  For negative inputs (e.g., `"-1"`), `static_cast<std::size_t>(-1)` yields `SIZE_MAX` (≈1.8×10¹⁹ on 64-bit).  
  This is then passed to `benchmark_threads`, which attempts `std::vector<std::thread>(SIZE_MAX - 1)` — guaranteed `std::bad_alloc` or OS resource exhaustion before any graceful recovery. No diagnostic message is printed.

- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer: `./unique_creation -1`
  - Observed result: `std::bad_alloc` or OS-level resource exhaustion
  - Expected result: validated error message and graceful exit
- **Efficiency impact:** N/A (crash path)
- **C++20-specific note:** None specific; this is a long-standing pitfall with signed-to-unsigned cast.
- **Suggested fix:** Same guard as SEV-1 above (`n <= 0` check before cast).

---

### [SEV-2] `std::stoi` throws on malformed or out-of-range input; no exception handler

- **Location:** `unique_creation.cpp:22`
- **Why this is a problem:**  
  `std::stoi` throws `std::invalid_argument` on non-numeric input and `std::out_of_range` when the value exceeds `INT_MAX`. Both are uncaught, terminating the process with a noisy `terminate()` call and no user-facing message.  
  Every other benchmark in the directory (`shared_creation.cpp:22`, `unique_lookup.cpp:22`, etc.) has the same flaw — but this file is the reviewed target.

- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer: `./unique_creation abc` or `./unique_creation 99999999999999`
  - Observed result: `terminate called after throwing an instance of 'std::invalid_argument'`
  - Expected result: clean diagnostic and non-zero exit
- **Efficiency impact:** N/A
- **C++20-specific note:** `std::from_chars` (C++17/20) avoids exceptions entirely and is preferable for CLI parsing:
  ```cpp
  std::size_t n = 1;
  if (argc > 1) {
    auto [ptr, ec] = std::from_chars(argv[1], argv[1] + std::strlen(argv[1]), n);
    if (ec != std::errc{} || n == 0) { std::cerr << "Invalid thread count\n"; return 1; }
    number_of_threads = n;
  }
  ```
- **Suggested fix:** Replace `std::stoi` with `std::from_chars` and add range validation.

---

### [SEV-2] `benchmark_threads` internal `id` counter uses `int`; overflows for large thread counts

- **Location:** `benchmark_shared.h:24`
- **Why this is a problem:**  
  ```cpp
  int id = 1;
  for (auto& thread : threads) { thread = std::thread(f, id); ++id; }
  ```
  `id` is `int` but `threads.size()` is `std::size_t`. For `number_of_threads > INT_MAX`, `++id` overflows, which is signed integer overflow — **undefined behaviour** in C++. Even at practical scales, `id` silently wraps at 2^31 on 32-bit or is platform-dependent on 64-bit.  
  The lambda in `unique_creation.cpp` uses the `id` as a leaf name discriminator: collisions produce incorrect pool-sharing behaviour that invalidates the "unique" benchmark semantics.

- **Evidence type:** plausible-issue (thread counts this large are impractical, but the type mismatch propagates to correctness of thread identity)
- **Evidence:**
  - Mechanism: signed integer overflow is UB; id collision corrupts benchmark uniqueness.
  - Triggering scenario: theoretically `number_of_threads > INT_MAX`; practically, any `number_of_threads > 2^31`.
  - Test plan: Pass `number_of_threads = INT_MAX + 2` on a 64-bit system, observe whether id wraps and two threads get the same leaf name.
- **Efficiency impact:** Incorrect benchmark measurement if id collision occurs.
- **C++20-specific note:** `benchmark_threads` should use `std::size_t id` consistently.
- **Suggested fix:** Change `int id = 1;` to `std::size_t id = 1;` in `benchmark_shared.h`; update all lambda signatures to `std::size_t`.

---

### [SEV-3] Integer truncation in `size / number_of_threads` makes total work non-reproducible

- **Location:** `unique_creation.cpp:30`
- **Why this is a problem:**  
  When `size` (400000) is not evenly divisible by `number_of_threads`, integer truncation silently reduces depth: e.g., `400000 / 3 = 133333` (loses 1 layer). The total number of unique pool entries created across all threads is `N × ⌊400000/N⌋`, not `400000`. For a benchmark tool this means successive runs with different thread counts are not measuring equivalent amounts of work, making scaling measurements unreliable.

- **Evidence type:** static-proof
- **Evidence:**
  - `400000 / 3 * 3 = 399999 ≠ 400000`
  - The benchmark's `std::cerr` line only prints elapsed time — no work-item count — so the discrepancy is silent.
- **Efficiency impact:** Benchmark measurements are not directly comparable across thread counts.
- **C++20-specific note:** None.
- **Suggested fix:**  
  ```cpp
  // Distribute work evenly, ensuring total depth ≥ size
  const std::size_t depth_per_thread = (size + number_of_threads - 1) / number_of_threads;
  ```
  Or document that the division is intentional (with a comment).

---

### [SEV-3] Lambda captures `size` and `number_of_threads` by mutable reference; no `const` protection

- **Location:** `unique_creation.cpp:28–31`
- **Why this is a problem:**  
  The lambda `[&](int id) -> void` captures both `size` and `number_of_threads` by mutable reference. These variables are only read inside the lambda, but there is nothing preventing accidental mutation if this benchmark is extended. When multiple threads simultaneously read these through references, the absence of `const` makes the intent unclear and the code fragile.

- **Evidence type:** plausible-issue  
- **Efficiency impact:** None in current code, but data-race risk if a future edit writes to either variable after thread creation.
- **C++20-specific note:** C++20 lambda `[size, number_of_threads](int id) -> void` (value capture) is the correct idiom here since both variables are of value type (`std::size_t`) and are not modified.  
  Alternatively: `[size = size, not = number_of_threads]` or capture the computed `depth` directly:
  ```cpp
  const std::size_t depth = size / number_of_threads;
  auto nested_function = [depth](int id) -> void { ... };
  ```
- **Suggested fix:** Capture `size` and `number_of_threads` by value, or precompute `depth` as `const` and capture only that.

---

### [SEV-3] "Unique" per-thread semantics are undocumented; design intent is non-obvious

- **Location:** `unique_creation.cpp:28–31`
- **Why this is a problem:**  
  The file is named `unique_creation.cpp` and uses `std::to_string(id)` as the leaf name — making each thread's term tree disjoint from other threads' trees (no shared pool nodes). This is the deliberate contrast to `shared_creation.cpp` (leaf `"c"` is shared across threads). However, there is no comment explaining this design, so future maintainers may silently break the benchmark's distinctiveness by changing the leaf name.

- **Evidence type:** plausible-issue
- **Suggested fix:** Add a comment:
  ```cpp
  // Each thread uses a distinct leaf name so no pool entries are shared between
  // threads. This measures creation throughput without inter-thread contention on
  // the term pool.
  ```

---

## Proposed Tests

- [ ] `test_zero_threads`: Run `./unique_creation 0`; expected: graceful error message, exit code != 0, no crash/UB.
- [ ] `test_negative_threads`: Run `./unique_creation -1`; expected: graceful error message, exit code != 0, no resource exhaustion.
- [ ] `test_invalid_threads_arg`: Run `./unique_creation abc`; expected: no unhandled exception, graceful exit with message.
- [ ] `test_work_invariance`: Run `./unique_creation 1`, `./unique_creation 2`, `./unique_creation 4`; verify total unique pool terms created (e.g., via pool size query before/after) equals `size` in all cases — catches the truncation bug.
- [ ] `test_unique_vs_shared_no_collision`: After running with 2 threads, verify terms from thread 0 (`leaf="0"`) and thread 1 (`leaf="1"`) are provably disjoint in the pool; expected: no node shared between the two subtrees.

---

## Open Questions

1. **Why does `unique_creation.cpp` use `size / number_of_threads` (total work constant) while `shared_creation.cpp` uses `size` directly (total work scales with threads)?** The inconsistency suggests different authors or different unspoken intents — but there is no documentation. If the goal is weak-scaling measurement, that should be stated explicitly.
2. **Why is `benchmark_threads` internal `id` typed as `int` when `number_of_threads` is `std::size_t`?** This is an API inconsistency across all benchmark files that should be resolved globally in `benchmark_shared.h`.
3. **Is GC being disabled (`enable_garbage_collection(false)`) intentional for benchmarking creation cost specifically, or was it copied from another benchmark without consideration?** Disabling GC means pool memory grows unboundedly for the duration of the benchmark, which could invalidate memory-pressure performance conclusions.
4. **Is the timer in `benchmark_threads` intended to include thread-spawn overhead?** It starts before `std::thread` construction — this may inflate results for large `number_of_threads`, conflating thread creation latency with term creation latency.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         3m 7s
Total session time:     3m 15s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-sonnet-4.6       176.0k in, 10.8k out, 151.0k cached (Est. 1 Premium request)

```

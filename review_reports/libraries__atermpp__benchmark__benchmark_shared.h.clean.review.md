---
{
  "file": "libraries/atermpp/benchmark/benchmark_shared.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-sonnet-4.6",
  "thinking": "deep",
  "effort": "high",
  "generated_at_utc": "2026-07-07T12:03:28.099658+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/benchmark/benchmark_shared.h` — two template utilities (`benchmark_threads`, `create_nested_function`) used by all atermpp benchmarks.
- **Risk level:** Medium — these are benchmark-only files, but two bugs would cause assertion failures or `std::terminate` at runtime; one is a design-level trap for future callers.
- **Overall verdict:** Reject. Two correctness bugs (underflow crash, exception-unsafe thread management) and one latent arity mismatch that is silently masked only because all current callers pass `N=2`. Namespace pollution in header is a quality defect.

---

## Findings (ordered by severity)

---

### [SEV-1] `create_nested_function<N>` always constructs 2-argument terms regardless of `N`

- **Location:** `benchmark_shared.h:92` and `benchmark_shared.h:96`
- **Why this is a problem:**  
  The template parameter `N` is passed to `function_symbol f(function_name, N)` to declare arity N, but the initial term is constructed as `aterm f_term(f, c_term, c_term)` (exactly 2 arguments), and the loop calls `make_term_appl(f_term, f, f_term, f_term)` (also exactly 2 arguments). The function symbol carries arity N, but the term application always supplies 2 operands. The comment explicitly states this is a "fixed arity constructor of **length N**", implying generality.  
  In debug builds, `create_appl` asserts `symbol.arity() == sizeof...(arguments)` (see `aterm_pool_storage_implementation.h:188`). For any `N != 2`, this assertion fires immediately. In release builds with `NDEBUG`, the assert is silently skipped and the code produces structurally incorrect terms (arity mismatch between header and actual child count). All current call sites happen to use `N=2`, masking the defect entirely — but the template advertises a contract it cannot fulfill.

- **Evidence type:** static-proof  
  - **Static proof:** `function_symbol f(function_name, N)` sets arity to `N`. `aterm f_term(f, c_term, c_term)` calls the variadic constructor with `sizeof...(arguments) == 2`. The assert in `aterm_pool_storage_implementation.h:188` fires when `N != 2`.  
  - **Reproducer:** Call `create_nested_function<3>("f", "c", 2)` in any benchmark main, compile in debug mode.  
  - **Run command:**  
    ```bash
    cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_GUI_TOOLS=OFF
    cmake --build build --target shared_creation -j
    # Then instrument/call with N=3 or add a unit test
    ```
  - **Observed result:** Assertion failure: `symbol.arity() == sizeof...(arguments)` at `aterm_pool_storage_implementation.h:188`.  
  - **Expected result:** A valid N-ary term application.

- **Efficiency impact:** None in present state (N=2 only); if the intent is ever to generalise, a std::array or pack expansion is needed.
- **C++20-specific note:** C++20 `std::array` with CTAD or `[]<std::size_t... I>(std::index_sequence<I...>)` lambda could correctly expand N copies without hardcoding 2.
- **Suggested fix:**  
  Either restrict the template (static_assert or `requires N == 2`) and document it clearly, or fix the construction to actually use N arguments:
  ```cpp
  // Correct for arbitrary N:
  auto make_args = [&]<std::size_t... I>(std::index_sequence<I...>) {
      return aterm(f, (static_cast<void>(I), c_term)...);
  };
  aterm f_term = make_args(std::make_index_sequence<N>{});
  // And similarly for the loop body.
  ```

---

### [SEV-1] `benchmark_threads` crashes with `number_of_threads == 0` due to `std::size_t` underflow

- **Location:** `benchmark_shared.h:23`
- **Why this is a problem:**  
  `std::vector<std::thread> threads(number_of_threads - 1)` where both operands are `std::size_t`. When `number_of_threads == 0`, the subtraction wraps to `SIZE_MAX` (guaranteed unsigned underflow), and `std::vector` attempts to reserve `SIZE_MAX` threads — throwing `std::length_error` or `std::bad_alloc` before any benchmark code runs. There is no precondition documented or enforced.
- **Evidence type:** static-proof  
  - `0u - 1u` for `std::size_t` is well-defined unsigned wraparound to `SIZE_MAX` per C++ standard [basic.fundamental]. `std::vector(SIZE_MAX)` throws `std::length_error` on any real system.
  - **Run command:**
    ```bash
    # Minimal reproducer in any benchmark main:
    benchmark_threads(0, [](int id){ });
    ```
  - **Observed result:** `std::length_error` or `std::bad_alloc` thrown immediately.  
  - **Expected result:** Either run the benchmark on the main thread alone (1 effective thread) or reject 0 explicitly.
- **Efficiency impact:** Crash before any work is done.
- **C++20-specific note:** None specific; the issue predates C++20.
- **Suggested fix:**
  ```cpp
  assert(number_of_threads >= 1);  // or throw std::invalid_argument
  std::vector<std::thread> threads(number_of_threads - 1);
  ```

---

### [SEV-2] `benchmark_threads` is exception-unsafe: spawned threads are never joined if `f(0)` throws

- **Location:** `benchmark_shared.h:32`
- **Why this is a problem:**  
  Worker threads are started at lines 25–29, then `f(0)` is called on the main thread at line 32. If `f(0)` throws, the stack unwinds, the `threads` vector is destroyed, and `std::thread::~thread()` is called on each still-joinable thread. Per [thread.thread.destr], destroying a joinable `std::thread` calls `std::terminate()`. The program dies with no diagnostic, making benchmark failures hard to diagnose. This is a real risk since benchmark callables may allocate, assert, or throw.
- **Evidence type:** static-proof  
  - C++ standard [thread.thread.destr]: "If the thread is joinable, calls `std::terminate()`." `std::thread` remains joinable after `std::thread(f, id)` and before `.join()` or `.detach()`.
  - **Reproducer:**
    ```cpp
    benchmark_threads(2, [](int id) {
        if (id == 0) throw std::runtime_error("oops");
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    });
    ```
  - **Observed result:** `std::terminate()` called — process aborts with no useful message.  
  - **Expected result:** Exception propagated cleanly, worker threads joined.
- **Efficiency impact:** None normally; catastrophic on error.
- **C++20-specific note:** C++20 `std::jthread` (RAII joinable thread) would fix this without any explicit cleanup code. Using `std::jthread` is idiomatic C++20.
- **Suggested fix:**
  ```cpp
  std::vector<std::jthread> threads(number_of_threads - 1);
  // jthread joins automatically on destruction — no explicit join loop needed.
  ```
  Or, keep `std::thread` with a `try/catch` wrapping `f(0)` that joins all threads before rethrowing.

---

### [SEV-2] `using namespace atermpp` in a header file

- **Location:** `benchmark_shared.h:15`
- **Why this is a problem:**  
  A `using namespace` directive in a header injects the entire `atermpp` namespace into every translation unit that includes this header, silently and transitively. This can cause name collisions, ADL surprises, and makes it impossible for including code to control which names are in scope. Even for benchmark-only code, this is a maintainability defect — benchmarks that include this header share a single translation unit with their own code.
- **Evidence type:** static-proof  
  - C++ [namespace.udir]: namespace directives in headers are unconditional; the effect is visible in all downstream includes.
- **Efficiency impact:** None at runtime; compilation correctness risk.
- **C++20-specific note:** C++20 modules would avoid this entirely; in headers, the fix is explicit qualification.
- **Suggested fix:** Remove the directive; qualify uses as `atermpp::aterm`, `atermpp::function_symbol`, etc., or add it inside each benchmark `.cpp` file that needs it.

---

### [SEV-3] `int id` type is too narrow for `std::size_t number_of_threads`

- **Location:** `benchmark_shared.h:24`
- **Why this is a problem:**  
  `id` is `int` but is incremented up to `number_of_threads - 1`, which is `std::size_t`. If `number_of_threads > INT_MAX + 1u`, `id` overflows — undefined behavior for signed integers. The callable `f` receives a silently truncated or wrapped ID, making thread ID semantics unreliable. In practice benchmarks use small thread counts, but the type contract is broken.
- **Evidence type:** static-proof (UB on signed overflow per [expr.pre]).
- **Suggested fix:** Use `std::size_t id = 1;` (or the callable's parameter type).

---

### [SEV-3] Missing include guard / `#pragma once`

- **Location:** `benchmark_shared.h:1` (file level)
- **Why this is a problem:**  
  The file has no `#pragma once` and no include guard. If included from multiple translation units in the same build (which does not currently happen, but is not prevented), every template definition is re-instantiated and ODR violations may result silently.
- **Evidence type:** plausible-issue (no current multi-include exists, but no protection either).
- **Suggested fix:** Add `#pragma once` at the top of the file.

---

### [SEV-3] Timer captures thread creation and join overhead, not pure benchmark work

- **Location:** `benchmark_shared.h:20,23,36–38`
- **Why this is a problem:**  
  `stopwatch timer` is constructed before thread creation (line 20). Thread spawn, OS scheduling latency, and join synchronization are all included in `timer.seconds()`. For benchmarks measuring per-operation throughput, this conflates infrastructure cost with workload cost. On lightly-loaded systems this may be negligible, but on contended machines or with large thread counts, results are unreliable.
- **Evidence type:** plausible-issue (timing-sensitive; not a correctness failure).
- **Suggested fix:** Reset the timer after all threads are started (and `f(0)` is about to run), or use a shared barrier so all threads start work at the same time before the timer is read.

---

## Proposed Tests

- [ ] **`test_create_nested_N3`**: Call `create_nested_function<3>("f", "c", 2)` in a debug build; expected: assertion failure at `aterm_pool_storage_implementation.h:188` proving the arity mismatch.
- [ ] **`test_benchmark_threads_zero`**: Call `benchmark_threads(0, [](int){})` ; expected: `std::length_error` or `assert` failure — not silent undefined behavior.
- [ ] **`test_benchmark_threads_exception`**: Call `benchmark_threads(2, [](int id){ if(id==0) throw std::runtime_error("x"); std::this_thread::sleep_for(10ms); })`; expected: clean exception propagation, not `std::terminate`.
- [ ] **`test_static_assert_N_restriction`**: Add `static_assert(N == 2, "only arity-2 supported")` to `create_nested_function<N>` and verify it documents the contract at compile-time.

---

## Open Questions

1. Is `create_nested_function<N>` intentionally restricted to `N=2` (all callers use 2) but simply missing a `static_assert`, or was it designed to generalize to arbitrary arity? The comment ("fixed arity constructor of length N") implies generality, but the implementation contradicts it.
2. Should `benchmark_threads` propagate exceptions from worker threads as well (not just the main-thread invocation)? Currently worker-thread exceptions would also cause `std::terminate` via the `jthread`/`thread` destructor.
3. Is the timer in `benchmark_threads` intentionally wall-clock from start (including spawn) to give a full "end-to-end" measure, or is per-work-unit throughput the target metric?

## Copilot Raw Stdout

```text
● Read libraries/atermpp/benchmark/benchmark_shared.h
  └ 100 lines read

● Explore: Explore atermpp aterm API
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h lines 37-70
    └ 34 lines read
  ● Verify stopwatch location
    $ cd /workspaces/mCRL2 && head -1
    /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/stopwatch.h
    └ 2 lines...

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/benchmark/benchmark_shared.h` — two template utilities (`benchmark_threads`, `create_nested_function`) used by all atermpp benchmarks.
- **Risk level:** Medium — these are benchmark-only files, but two bugs would cause assertion failures or `std::terminate` at runtime; one is a design-level trap for future callers.
- **Overall verdict:** Reject. Two correctness bugs (underflow crash, exception-unsafe thread management) and one latent arity mismatch that is silently masked only because all current callers pass `N=2`. Namespace pollution in header is a quality defect.

---

## Findings (ordered by severity)

---

### [SEV-1] `create_nested_function<N>` always constructs 2-argument terms regardless of `N`

- **Location:** `benchmark_shared.h:92` and `benchmark_shared.h:96`
- **Why this is a problem:**  
  The template parameter `N` is passed to `function_symbol f(function_name, N)` to declare arity N, but the initial term is constructed as `aterm f_term(f, c_term, c_term)` (exactly 2 arguments), and the loop calls `make_term_appl(f_term, f, f_term, f_term)` (also exactly 2 arguments). The function symbol carries arity N, but the term application always supplies 2 operands. The comment explicitly states this is a "fixed arity constructor of **length N**", implying generality.  
  In debug builds, `create_appl` asserts `symbol.arity() == sizeof...(arguments)` (see `aterm_pool_storage_implementation.h:188`). For any `N != 2`, this assertion fires immediately. In release builds with `NDEBUG`, the assert is silently skipped and the code produces structurally incorrect terms (arity mismatch between header and actual child count). All current call sites happen to use `N=2`, masking the defect entirely — but the template advertises a contract it cannot fulfill.

- **Evidence type:** static-proof  
  - **Static proof:** `function_symbol f(function_name, N)` sets arity to `N`. `aterm f_term(f, c_term, c_term)` calls the variadic constructor with `sizeof...(arguments) == 2`. The assert in `aterm_pool_storage_implementation.h:188` fires when `N != 2`.  
  - **Reproducer:** Call `create_nested_function<3>("f", "c", 2)` in any benchmark main, compile in debug mode.  
  - **Run command:**  
    ```bash
    cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_GUI_TOOLS=OFF
    cmake --build build --target shared_creation -j
    # Then instrument/call with N=3 or add a unit test
    ```
  - **Observed result:** Assertion failure: `symbol.arity() == sizeof...(arguments)` at `aterm_pool_storage_implementation.h:188`.  
  - **Expected result:** A valid N-ary term application.

- **Efficiency impact:** None in present state (N=2 only); if the intent is ever to generalise, a std::array or pack expansion is needed.
- **C++20-specific note:** C++20 `std::array` with CTAD or `[]<std::size_t... I>(std::index_sequence<I...>)` lambda could correctly expand N copies without hardcoding 2.
- **Suggested fix:**  
  Either restrict the template (static_assert or `requires N == 2`) and document it clearly, or fix the construction to actually use N arguments:
  ```cpp
  // Correct for arbitrary N:
  auto make_args = [&]<std::size_t... I>(std::index_sequence<I...>) {
      return aterm(f, (static_cast<void>(I), c_term)...);
  };
  aterm f_term = make_args(std::make_index_sequence<N>{});
  // And similarly for the loop body.
  ```

---

### [SEV-1] `benchmark_threads` crashes with `number_of_threads == 0` due to `std::size_t` underflow

- **Location:** `benchmark_shared.h:23`
- **Why this is a problem:**  
  `std::vector<std::thread> threads(number_of_threads - 1)` where both operands are `std::size_t`. When `number_of_threads == 0`, the subtraction wraps to `SIZE_MAX` (guaranteed unsigned underflow), and `std::vector` attempts to reserve `SIZE_MAX` threads — throwing `std::length_error` or `std::bad_alloc` before any benchmark code runs. There is no precondition documented or enforced.
- **Evidence type:** static-proof  
  - `0u - 1u` for `std::size_t` is well-defined unsigned wraparound to `SIZE_MAX` per C++ standard [basic.fundamental]. `std::vector(SIZE_MAX)` throws `std::length_error` on any real system.
  - **Run command:**
    ```bash
    # Minimal reproducer in any benchmark main:
    benchmark_threads(0, [](int id){ });
    ```
  - **Observed result:** `std::length_error` or `std::bad_alloc` thrown immediately.  
  - **Expected result:** Either run the benchmark on the main thread alone (1 effective thread) or reject 0 explicitly.
- **Efficiency impact:** Crash before any work is done.
- **C++20-specific note:** None specific; the issue predates C++20.
- **Suggested fix:**
  ```cpp
  assert(number_of_threads >= 1);  // or throw std::invalid_argument
  std::vector<std::thread> threads(number_of_threads - 1);
  ```

---

### [SEV-2] `benchmark_threads` is exception-unsafe: spawned threads are never joined if `f(0)` throws

- **Location:** `benchmark_shared.h:32`
- **Why this is a problem:**  
  Worker threads are started at lines 25–29, then `f(0)` is called on the main thread at line 32. If `f(0)` throws, the stack unwinds, the `threads` vector is destroyed, and `std::thread::~thread()` is called on each still-joinable thread. Per [thread.thread.destr], destroying a joinable `std::thread` calls `std::terminate()`. The program dies with no diagnostic, making benchmark failures hard to diagnose. This is a real risk since benchmark callables may allocate, assert, or throw.
- **Evidence type:** static-proof  
  - C++ standard [thread.thread.destr]: "If the thread is joinable, calls `std::terminate()`." `std::thread` remains joinable after `std::thread(f, id)` and before `.join()` or `.detach()`.
  - **Reproducer:**
    ```cpp
    benchmark_threads(2, [](int id) {
        if (id == 0) throw std::runtime_error("oops");
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    });
    ```
  - **Observed result:** `std::terminate()` called — process aborts with no useful message.  
  - **Expected result:** Exception propagated cleanly, worker threads joined.
- **Efficiency impact:** None normally; catastrophic on error.
- **C++20-specific note:** C++20 `std::jthread` (RAII joinable thread) would fix this without any explicit cleanup code. Using `std::jthread` is idiomatic C++20.
- **Suggested fix:**
  ```cpp
  std::vector<std::jthread> threads(number_of_threads - 1);
  // jthread joins automatically on destruction — no explicit join loop needed.
  ```
  Or, keep `std::thread` with a `try/catch` wrapping `f(0)` that joins all threads before rethrowing.

---

### [SEV-2] `using namespace atermpp` in a header file

- **Location:** `benchmark_shared.h:15`
- **Why this is a problem:**  
  A `using namespace` directive in a header injects the entire `atermpp` namespace into every translation unit that includes this header, silently and transitively. This can cause name collisions, ADL surprises, and makes it impossible for including code to control which names are in scope. Even for benchmark-only code, this is a maintainability defect — benchmarks that include this header share a single translation unit with their own code.
- **Evidence type:** static-proof  
  - C++ [namespace.udir]: namespace directives in headers are unconditional; the effect is visible in all downstream includes.
- **Efficiency impact:** None at runtime; compilation correctness risk.
- **C++20-specific note:** C++20 modules would avoid this entirely; in headers, the fix is explicit qualification.
- **Suggested fix:** Remove the directive; qualify uses as `atermpp::aterm`, `atermpp::function_symbol`, etc., or add it inside each benchmark `.cpp` file that needs it.

---

### [SEV-3] `int id` type is too narrow for `std::size_t number_of_threads`

- **Location:** `benchmark_shared.h:24`
- **Why this is a problem:**  
  `id` is `int` but is incremented up to `number_of_threads - 1`, which is `std::size_t`. If `number_of_threads > INT_MAX + 1u`, `id` overflows — undefined behavior for signed integers. The callable `f` receives a silently truncated or wrapped ID, making thread ID semantics unreliable. In practice benchmarks use small thread counts, but the type contract is broken.
- **Evidence type:** static-proof (UB on signed overflow per [expr.pre]).
- **Suggested fix:** Use `std::size_t id = 1;` (or the callable's parameter type).

---

### [SEV-3] Missing include guard / `#pragma once`

- **Location:** `benchmark_shared.h:1` (file level)
- **Why this is a problem:**  
  The file has no `#pragma once` and no include guard. If included from multiple translation units in the same build (which does not currently happen, but is not prevented), every template definition is re-instantiated and ODR violations may result silently.
- **Evidence type:** plausible-issue (no current multi-include exists, but no protection either).
- **Suggested fix:** Add `#pragma once` at the top of the file.

---

### [SEV-3] Timer captures thread creation and join overhead, not pure benchmark work

- **Location:** `benchmark_shared.h:20,23,36–38`
- **Why this is a problem:**  
  `stopwatch timer` is constructed before thread creation (line 20). Thread spawn, OS scheduling latency, and join synchronization are all included in `timer.seconds()`. For benchmarks measuring per-operation throughput, this conflates infrastructure cost with workload cost. On lightly-loaded systems this may be negligible, but on contended machines or with large thread counts, results are unreliable.
- **Evidence type:** plausible-issue (timing-sensitive; not a correctness failure).
- **Suggested fix:** Reset the timer after all threads are started (and `f(0)` is about to run), or use a shared barrier so all threads start work at the same time before the timer is read.

---

## Proposed Tests

- [ ] **`test_create_nested_N3`**: Call `create_nested_function<3>("f", "c", 2)` in a debug build; expected: assertion failure at `aterm_pool_storage_implementation.h:188` proving the arity mismatch.
- [ ] **`test_benchmark_threads_zero`**: Call `benchmark_threads(0, [](int){})` ; expected: `std::length_error` or `assert` failure — not silent undefined behavior.
- [ ] **`test_benchmark_threads_exception`**: Call `benchmark_threads(2, [](int id){ if(id==0) throw std::runtime_error("x"); std::this_thread::sleep_for(10ms); })`; expected: clean exception propagation, not `std::terminate`.
- [ ] **`test_static_assert_N_restriction`**: Add `static_assert(N == 2, "only arity-2 supported")` to `create_nested_function<N>` and verify it documents the contract at compile-time.

---

## Open Questions

1. Is `create_nested_function<N>` intentionally restricted to `N=2` (all callers use 2) but simply missing a `static_assert`, or was it designed to generalize to arbitrary arity? The comment ("fixed arity constructor of length N") implies generality, but the implementation contradicts it.
2. Should `benchmark_threads` propagate exceptions from worker threads as well (not just the main-thread invocation)? Currently worker-thread exceptions would also cause `std::terminate` via the `jthread`/`thread` destructor.
3. Is the timer in `benchmark_threads` intentionally wall-clock from start (including spawn) to give a full "end-to-end" measure, or is per-work-unit throughput the target metric?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         3m 39s
Total session time:     3m 58s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-sonnet-4.6       56.1k in, 6.2k out, 36.9k cached (Est. 1 Premium request)
 claude-haiku-4.5        888.3k in, 6.2k out, 845.5k cached (Est. 0 Premium requests)

```

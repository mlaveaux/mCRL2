---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/function_symbol_generator.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T12:32:50.394522+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/function_symbol_generator.h` — clean/full-file review
- **Risk level:** HIGH — multiple concurrency defects, an ODR violation, and exception-safety holes in a header included from at least 4 distinct TUs
- **Overall verdict:** This file has several correctness defects ranging from technically-UB (but mitigated by other mechanisms) to genuinely exploitable data races. The most critical issues are the `static inline` ODR violation, exception-unsafe locking, and the unsynchronized read of `*m_central_index` in `operator()`.

---

## Findings (ordered by severity)

### [SEV-1] `static inline` free functions create per-TU singletons — ODR violation

- **Location:** `function_symbol_generator.h:21-31`
- **Why this is a problem:** `static inline` at namespace scope gives internal linkage. Each translation unit that includes this header gets its **own** `std::mutex` and its **own** `size_t` sequence counter. The inline member functions of `function_symbol_generator` (defined in-class, thus implicitly inline with external linkage) reference these internal-linkage entities — this is an ODR violation \[basic.def.odr\]/13.2, which is undefined behavior. In practice, different TUs use different mutexes (no cross-TU synchronization) and independent counters (non-unique prefixes across TUs).
- **Evidence type:** static-proof
- **Evidence:**
  - The C++20 standard states: an inline function defined in multiple TUs shall not use a name with internal linkage unless it names the same entity in all TUs. `static` functions have internal linkage; each TU names a different entity.
  - The header is included transitively from at least: `enumerator_test.cpp`, `jitty.cpp`, `function_symbol_generator_test.cpp`, and any TU including `pres2res.h`.
  - Two TUs constructing a generator with prefix `"x_"` both start `generator_sequence_number()` at 0, producing identical prefix `"x_"` in both TUs. While `register_prefix` partially compensates via `get_sufficiently_large_postfix_index`, the mutex provides no cross-TU protection.
- **Efficiency impact:** None directly; the real cost is correctness.
- **C++20-specific note:** Should use plain `inline` (external linkage) to guarantee a single static-local instance across all TUs, which C++17+ guarantees for inline functions.
- **Suggested fix:**
  ```cpp
  inline std::mutex& function_symbol_generator_mutex()
  {
    static std::mutex m;
    return m;
  }

  inline std::size_t& generator_sequence_number()
  {
    static std::size_t n = 0;
    return n;
  }
  ```

---

### [SEV-1] Manual lock/unlock without RAII — exception-unsafe

- **Location:** `function_symbol_generator.h:60-78` (constructor), `function_symbol_generator.h:90-98` (destructor), `function_symbol_generator.h:127-146` (operator with prefix)
- **Why this is a problem:** Between `lock()` and `unlock()`, the constructor performs string concatenation (`operator+`, `std::to_string`) and calls `register_prefix`, all of which may throw `std::bad_alloc`. If an exception propagates, the mutex is never unlocked, causing permanent deadlock for all subsequent generator constructions.
- **Evidence type:** static-proof
- **Evidence:**
  - `m_prefix = prefix + (...)` (line 65): `std::string::operator+` can throw `std::bad_alloc`.
  - `register_prefix(m_prefix)` (line 70): `std::make_shared` inside this function can throw.
  - No `try`/`catch`, no `std::lock_guard`, no RAII wrapper is present.
- **Efficiency impact:** N/A (correctness issue).
- **C++20-specific note:** Use `std::lock_guard` or the project's `mcrl2::utilities::mutex` wrapper which provides RAII.
- **Suggested fix:** Replace manual lock/unlock with `std::lock_guard`:
  ```cpp
  if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
  {
    std::lock_guard lock(function_symbol_generator_mutex());
    // ... body ...
  }
  ```
  Or restructure to use `mcrl2::utilities::mutex` which already handles the `GlobalThreadSafe` conditionality internally.

---

### [SEV-2] Data race on `*m_central_index` in `operator()(std::size_t)`

- **Location:** `function_symbol_generator.h:106`
- **Why this is a problem:** `operator()(std::size_t)` reads `*m_central_index` without holding any lock. This `shared_ptr<std::size_t>` points to an object that `register_prefix` and `get_sufficiently_large_postfix_index` in the symbol pool can modify (the pool holds its own copy of the `shared_ptr`). If another thread concurrently calls `register_prefix` or constructs a `function_symbol` that triggers `create_helper`, a concurrent write to `*m_central_index` races with the unsynchronized read here.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `register_prefix` returns the same `shared_ptr` to all generators sharing a prefix. The symbol pool's `create_helper` may increment `*shared_index` when new symbols matching the prefix are created.
  - Two threads: T1 calls `generator()` reading `*m_central_index`; T2 concurrently creates a function symbol triggering `create_helper` which writes to the same pointee.
  - **Test plan:** Build with `MCRL2_ENABLE_THREADSANITIZER=ON`. Create two threads: one continuously generating symbols, the other continuously creating function symbols with the same prefix via the pool directly. TSan should report a data race on the `size_t` pointed to by the shared index.
  - **Run command:** `cmake --build build-tsan && ctest --test-dir build-tsan -R function_symbol_generator`
- **Efficiency impact:** None; the fix is a lock or atomic.
- **C++20-specific note:** Consider `std::atomic<std::size_t>` for the central index, or hold a lock during the read.
- **Suggested fix:** Either make the central index `std::shared_ptr<std::atomic<std::size_t>>` or hold `function_symbol_generator_mutex()` during the check.

---

### [SEV-2] `operator()(const std::string&)` corrupts shared generator state

- **Location:** `function_symbol_generator.h:125-148`
- **Why this is a problem:** This overload:
  1. Overwrites `m_central_index` (line 136) with a new prefix's index — silently breaking `operator()(std::size_t)` which relies on `m_central_index` corresponding to `m_prefix`.
  2. Sets `m_initial_index = m_index` (line 139) which has no clear semantic purpose — it doesn't advance the central index for the new prefix, so subsequent calls to this overload with the same prefix may generate duplicates.
  3. Does not update `m_prefix` or `m_string_buffer`, leaving the object in an inconsistent state where member variables reference different prefixes.
- **Evidence type:** static-proof
- **Evidence:**
  - After `gen("foo")` is called, `m_central_index` points to the index for `"foo"`, but `m_prefix` is still the original prefix. A subsequent call to `gen(0)` will use the wrong `m_central_index` to validate `m_index`.
  - The central index `*m_central_index` is never incremented, so two calls to `gen("bar")` both read the same value and generate the same function symbol name.
- **C++20-specific note:** N/A — logic bug independent of language version.
- **Suggested fix:** This overload appears fundamentally broken. Either make it a free function that creates a temporary generator, or properly update all relevant state and increment the central index.

---

### [SEV-2] Class is implicitly copyable — double deregister on destruction

- **Location:** `function_symbol_generator.h:35` (commented-out `noncopyable`)
- **Why this is a problem:** Without the `noncopyable` base, the class has compiler-generated copy constructor and copy-assignment operator. Copying a generator produces two objects that both hold the same `m_prefix`. Both destructors call `deregister(m_prefix)` — the second call erases a non-existent entry (benign on its own), but any third generator sharing that prefix loses its registration.
- **Evidence type:** static-proof
- **Evidence:**
  - Copy construct: both original and copy have identical `m_prefix`. First destructor calls `deregister` successfully. Second destructor calls `deregister` on an already-removed entry — `std::map::erase` with a missing key is a no-op, so no crash, but the design contract is violated.
  - If a third generator registered the same prefix between the two destructions, the second `deregister` would remove its registration.
- **C++20-specific note:** Use `= delete` directly instead of inheritance for clarity in C++20:
  ```cpp
  function_symbol_generator(const function_symbol_generator&) = delete;
  function_symbol_generator& operator=(const function_symbol_generator&) = delete;
  ```
- **Suggested fix:** Uncomment the `noncopyable` base or explicitly delete copy operations.

---

### [SEV-3] Uses raw `std::mutex` instead of project wrappers

- **Location:** `function_symbol_generator.h:21-25`
- **Why this is a problem:** The project convention (per `mcrl2/utilities/mutex.h`) wraps mutexes so that single-threaded builds compile away all synchronization. The manual `if constexpr (GlobalThreadSafe) { mutex.lock(); }` pattern is error-prone and duplicates logic that the wrappers already encapsulate.
- **Evidence type:** plausible-issue (style/convention violation)
- **Suggested fix:** Replace with `mcrl2::utilities::mutex` and use `std::lock_guard` over it.

---

### [SEV-3] Missing `[[nodiscard]]` on `operator()` overloads

- **Location:** `function_symbol_generator.h:102, 125`
- **Why this is a problem:** The entire purpose of calling `operator()` is to obtain the generated `function_symbol`. Discarding the return value is certainly a bug. `[[nodiscard]]` would catch this at compile time.
- **C++20-specific note:** `[[nodiscard]]` is standard since C++17 and should be applied to any function whose return value must not be ignored.
- **Suggested fix:**
  ```cpp
  [[nodiscard]] function_symbol operator()(std::size_t arity = 0)
  [[nodiscard]] function_symbol operator()(const std::string& prefix)
  ```

---

### [SEV-3] Incomplete Doxygen — missing `\post` on all mutating operations

- **Location:** `function_symbol_generator.h:53-55, 81, 101, 122`
- **Why this is a problem:** Per project conventions, public APIs must document postconditions. The constructor documents `\pre` but no `\post`. `clear()` has no contract documentation. `operator()` doesn't specify the uniqueness guarantee or its conditions.
- **Suggested fix:** Add `\post` clauses specifying:
  - Constructor: generator is ready to produce symbols starting at the central index.
  - `clear()`: next generated symbol will be the same as the first one produced after construction.
  - `operator()`: returned symbol has a unique name not previously generated by any generator sharing this prefix (conditional on `clear()` not being called).

---

## Proposed Tests

- [ ] **`test_cross_tu_uniqueness`**: Link two TUs that each construct a `function_symbol_generator("x_")` and generate 1000 symbols. Collect all names into a shared set and verify no duplicates. Expected failure signal: duplicate name detected (assertion failure).
- [ ] **`test_exception_safety`**: Mock/intercept allocation to throw `std::bad_alloc` during constructor. Verify that subsequent generator constructions do not deadlock (timeout = failure signal). Run with: `ctest --test-dir build -R function_symbol_generator --timeout 5`
- [ ] **`test_tsan_concurrent_generation`**: Two threads sharing a generator prefix, one calling `operator()(0)` and the other calling `function_symbol(prefix + "5", 0)` to trigger `create_helper`. Build with `-DMCRL2_ENABLE_THREADSANITIZER=ON`. Expected failure: TSan race report.
- [ ] **`test_operator_prefix_state_corruption`**: Call `gen("foo")` then `gen(0)`. Verify the name of the second symbol still uses the original prefix, not "foo". Expected failure: symbol name mismatch.
- [ ] **`test_copy_double_deregister`**: Copy a generator, destroy the copy, then verify the original still generates valid unique symbols. Expected failure: prefix lookup fails or duplicates are generated.

---

## Open Questions

1. Is the `operator()(const std::string& prefix)` overload still actively used? It appears fundamentally broken (corrupts object state, doesn't increment central index). If unused, it should be removed.
2. Was the `noncopyable` base deliberately commented out, or was it accidentally lost during a refactor? The class is clearly not safely copyable.
3. Should `*m_central_index` be `std::atomic<std::size_t>` given that it's shared across generators and the symbol pool can modify it concurrently?
4. Is the `clear()` functionality (intentional reuse of previously generated symbol names) actually safe in all call sites, or does it silently create name collisions?

## Copilot Raw Stdout

```text
● Explore: Gather context for review
  ● Find all .cpp files with direct includes of function_symbol_generator
    $ grep -r "#include.*function_symbol_generator" /workspaces/mCRL2 --include="*.cpp" 2>/dev/null
    └ 2 lines...
  ● Find all .h files with direct includes of function_symbol_generator
    $ grep -r "#include.*function_symbol_generator" /workspaces/mCRL2 --include="*.h" 2>/dev/null
    └ 3 lines...

● Read libraries/atermpp/source/function_symbol_pool.cpp lines 100-145
  └ 46 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/function_symbol_generator.h
  └ 154 lines read

● Read libraries/atermpp/test/function_symbol_generator_test.cpp
  └ 53 lines read

Now let me check one more detail about the ODR implications:

● Grep "function_symbol_generator" in **/*.cpp (/workspaces/mCRL2/libraries)
  └ 1 file found

● Grep "enumerator_identifier_generator" in **/*.cpp (/workspaces/mCRL2/libraries)
  └ 3 files found

● Read libraries/data/include/mcrl2/data/detail/enumerator_identifier_generator.h lines 1-40
  └ 40 lines read

● Grep "#include.*pres2res"" in **/*.cpp (/workspaces/mCRL2/libraries)
  └ No matches found

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/function_symbol_generator.h` — clean/full-file review
- **Risk level:** HIGH — multiple concurrency defects, an ODR violation, and exception-safety holes in a header included from at least 4 distinct TUs
- **Overall verdict:** This file has several correctness defects ranging from technically-UB (but mitigated by other mechanisms) to genuinely exploitable data races. The most critical issues are the `static inline` ODR violation, exception-unsafe locking, and the unsynchronized read of `*m_central_index` in `operator()`.

---

## Findings (ordered by severity)

### [SEV-1] `static inline` free functions create per-TU singletons — ODR violation

- **Location:** `function_symbol_generator.h:21-31`
- **Why this is a problem:** `static inline` at namespace scope gives internal linkage. Each translation unit that includes this header gets its **own** `std::mutex` and its **own** `size_t` sequence counter. The inline member functions of `function_symbol_generator` (defined in-class, thus implicitly inline with external linkage) reference these internal-linkage entities — this is an ODR violation \[basic.def.odr\]/13.2, which is undefined behavior. In practice, different TUs use different mutexes (no cross-TU synchronization) and independent counters (non-unique prefixes across TUs).
- **Evidence type:** static-proof
- **Evidence:**
  - The C++20 standard states: an inline function defined in multiple TUs shall not use a name with internal linkage unless it names the same entity in all TUs. `static` functions have internal linkage; each TU names a different entity.
  - The header is included transitively from at least: `enumerator_test.cpp`, `jitty.cpp`, `function_symbol_generator_test.cpp`, and any TU including `pres2res.h`.
  - Two TUs constructing a generator with prefix `"x_"` both start `generator_sequence_number()` at 0, producing identical prefix `"x_"` in both TUs. While `register_prefix` partially compensates via `get_sufficiently_large_postfix_index`, the mutex provides no cross-TU protection.
- **Efficiency impact:** None directly; the real cost is correctness.
- **C++20-specific note:** Should use plain `inline` (external linkage) to guarantee a single static-local instance across all TUs, which C++17+ guarantees for inline functions.
- **Suggested fix:**
  ```cpp
  inline std::mutex& function_symbol_generator_mutex()
  {
    static std::mutex m;
    return m;
  }

  inline std::size_t& generator_sequence_number()
  {
    static std::size_t n = 0;
    return n;
  }
  ```

---

### [SEV-1] Manual lock/unlock without RAII — exception-unsafe

- **Location:** `function_symbol_generator.h:60-78` (constructor), `function_symbol_generator.h:90-98` (destructor), `function_symbol_generator.h:127-146` (operator with prefix)
- **Why this is a problem:** Between `lock()` and `unlock()`, the constructor performs string concatenation (`operator+`, `std::to_string`) and calls `register_prefix`, all of which may throw `std::bad_alloc`. If an exception propagates, the mutex is never unlocked, causing permanent deadlock for all subsequent generator constructions.
- **Evidence type:** static-proof
- **Evidence:**
  - `m_prefix = prefix + (...)` (line 65): `std::string::operator+` can throw `std::bad_alloc`.
  - `register_prefix(m_prefix)` (line 70): `std::make_shared` inside this function can throw.
  - No `try`/`catch`, no `std::lock_guard`, no RAII wrapper is present.
- **Efficiency impact:** N/A (correctness issue).
- **C++20-specific note:** Use `std::lock_guard` or the project's `mcrl2::utilities::mutex` wrapper which provides RAII.
- **Suggested fix:** Replace manual lock/unlock with `std::lock_guard`:
  ```cpp
  if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
  {
    std::lock_guard lock(function_symbol_generator_mutex());
    // ... body ...
  }
  ```
  Or restructure to use `mcrl2::utilities::mutex` which already handles the `GlobalThreadSafe` conditionality internally.

---

### [SEV-2] Data race on `*m_central_index` in `operator()(std::size_t)`

- **Location:** `function_symbol_generator.h:106`
- **Why this is a problem:** `operator()(std::size_t)` reads `*m_central_index` without holding any lock. This `shared_ptr<std::size_t>` points to an object that `register_prefix` and `get_sufficiently_large_postfix_index` in the symbol pool can modify (the pool holds its own copy of the `shared_ptr`). If another thread concurrently calls `register_prefix` or constructs a `function_symbol` that triggers `create_helper`, a concurrent write to `*m_central_index` races with the unsynchronized read here.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `register_prefix` returns the same `shared_ptr` to all generators sharing a prefix. The symbol pool's `create_helper` may increment `*shared_index` when new symbols matching the prefix are created.
  - Two threads: T1 calls `generator()` reading `*m_central_index`; T2 concurrently creates a function symbol triggering `create_helper` which writes to the same pointee.
  - **Test plan:** Build with `MCRL2_ENABLE_THREADSANITIZER=ON`. Create two threads: one continuously generating symbols, the other continuously creating function symbols with the same prefix via the pool directly. TSan should report a data race on the `size_t` pointed to by the shared index.
  - **Run command:** `cmake --build build-tsan && ctest --test-dir build-tsan -R function_symbol_generator`
- **Efficiency impact:** None; the fix is a lock or atomic.
- **C++20-specific note:** Consider `std::atomic<std::size_t>` for the central index, or hold a lock during the read.
- **Suggested fix:** Either make the central index `std::shared_ptr<std::atomic<std::size_t>>` or hold `function_symbol_generator_mutex()` during the check.

---

### [SEV-2] `operator()(const std::string&)` corrupts shared generator state

- **Location:** `function_symbol_generator.h:125-148`
- **Why this is a problem:** This overload:
  1. Overwrites `m_central_index` (line 136) with a new prefix's index — silently breaking `operator()(std::size_t)` which relies on `m_central_index` corresponding to `m_prefix`.
  2. Sets `m_initial_index = m_index` (line 139) which has no clear semantic purpose — it doesn't advance the central index for the new prefix, so subsequent calls to this overload with the same prefix may generate duplicates.
  3. Does not update `m_prefix` or `m_string_buffer`, leaving the object in an inconsistent state where member variables reference different prefixes.
- **Evidence type:** static-proof
- **Evidence:**
  - After `gen("foo")` is called, `m_central_index` points to the index for `"foo"`, but `m_prefix` is still the original prefix. A subsequent call to `gen(0)` will use the wrong `m_central_index` to validate `m_index`.
  - The central index `*m_central_index` is never incremented, so two calls to `gen("bar")` both read the same value and generate the same function symbol name.
- **C++20-specific note:** N/A — logic bug independent of language version.
- **Suggested fix:** This overload appears fundamentally broken. Either make it a free function that creates a temporary generator, or properly update all relevant state and increment the central index.

---

### [SEV-2] Class is implicitly copyable — double deregister on destruction

- **Location:** `function_symbol_generator.h:35` (commented-out `noncopyable`)
- **Why this is a problem:** Without the `noncopyable` base, the class has compiler-generated copy constructor and copy-assignment operator. Copying a generator produces two objects that both hold the same `m_prefix`. Both destructors call `deregister(m_prefix)` — the second call erases a non-existent entry (benign on its own), but any third generator sharing that prefix loses its registration.
- **Evidence type:** static-proof
- **Evidence:**
  - Copy construct: both original and copy have identical `m_prefix`. First destructor calls `deregister` successfully. Second destructor calls `deregister` on an already-removed entry — `std::map::erase` with a missing key is a no-op, so no crash, but the design contract is violated.
  - If a third generator registered the same prefix between the two destructions, the second `deregister` would remove its registration.
- **C++20-specific note:** Use `= delete` directly instead of inheritance for clarity in C++20:
  ```cpp
  function_symbol_generator(const function_symbol_generator&) = delete;
  function_symbol_generator& operator=(const function_symbol_generator&) = delete;
  ```
- **Suggested fix:** Uncomment the `noncopyable` base or explicitly delete copy operations.

---

### [SEV-3] Uses raw `std::mutex` instead of project wrappers

- **Location:** `function_symbol_generator.h:21-25`
- **Why this is a problem:** The project convention (per `mcrl2/utilities/mutex.h`) wraps mutexes so that single-threaded builds compile away all synchronization. The manual `if constexpr (GlobalThreadSafe) { mutex.lock(); }` pattern is error-prone and duplicates logic that the wrappers already encapsulate.
- **Evidence type:** plausible-issue (style/convention violation)
- **Suggested fix:** Replace with `mcrl2::utilities::mutex` and use `std::lock_guard` over it.

---

### [SEV-3] Missing `[[nodiscard]]` on `operator()` overloads

- **Location:** `function_symbol_generator.h:102, 125`
- **Why this is a problem:** The entire purpose of calling `operator()` is to obtain the generated `function_symbol`. Discarding the return value is certainly a bug. `[[nodiscard]]` would catch this at compile time.
- **C++20-specific note:** `[[nodiscard]]` is standard since C++17 and should be applied to any function whose return value must not be ignored.
- **Suggested fix:**
  ```cpp
  [[nodiscard]] function_symbol operator()(std::size_t arity = 0)
  [[nodiscard]] function_symbol operator()(const std::string& prefix)
  ```

---

### [SEV-3] Incomplete Doxygen — missing `\post` on all mutating operations

- **Location:** `function_symbol_generator.h:53-55, 81, 101, 122`
- **Why this is a problem:** Per project conventions, public APIs must document postconditions. The constructor documents `\pre` but no `\post`. `clear()` has no contract documentation. `operator()` doesn't specify the uniqueness guarantee or its conditions.
- **Suggested fix:** Add `\post` clauses specifying:
  - Constructor: generator is ready to produce symbols starting at the central index.
  - `clear()`: next generated symbol will be the same as the first one produced after construction.
  - `operator()`: returned symbol has a unique name not previously generated by any generator sharing this prefix (conditional on `clear()` not being called).

---

## Proposed Tests

- [ ] **`test_cross_tu_uniqueness`**: Link two TUs that each construct a `function_symbol_generator("x_")` and generate 1000 symbols. Collect all names into a shared set and verify no duplicates. Expected failure signal: duplicate name detected (assertion failure).
- [ ] **`test_exception_safety`**: Mock/intercept allocation to throw `std::bad_alloc` during constructor. Verify that subsequent generator constructions do not deadlock (timeout = failure signal). Run with: `ctest --test-dir build -R function_symbol_generator --timeout 5`
- [ ] **`test_tsan_concurrent_generation`**: Two threads sharing a generator prefix, one calling `operator()(0)` and the other calling `function_symbol(prefix + "5", 0)` to trigger `create_helper`. Build with `-DMCRL2_ENABLE_THREADSANITIZER=ON`. Expected failure: TSan race report.
- [ ] **`test_operator_prefix_state_corruption`**: Call `gen("foo")` then `gen(0)`. Verify the name of the second symbol still uses the original prefix, not "foo". Expected failure: symbol name mismatch.
- [ ] **`test_copy_double_deregister`**: Copy a generator, destroy the copy, then verify the original still generates valid unique symbols. Expected failure: prefix lookup fails or duplicates are generated.

---

## Open Questions

1. Is the `operator()(const std::string& prefix)` overload still actively used? It appears fundamentally broken (corrupts object state, doesn't increment central index). If unused, it should be removed.
2. Was the `noncopyable` base deliberately commented out, or was it accidentally lost during a refactor? The class is clearly not safely copyable.
3. Should `*m_central_index` be `std::atomic<std::size_t>` given that it's shared across generators and the symbol pool can modify it concurrently?
4. Is the `clear()` functionality (intentional reuse of previously generated symbol names) actually safe in all call sites, or does it silently create name collisions?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         3m 55s
Total session time:     4m 8s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         223.5k in, 8.5k out, 195.6k cached (Est. 1 Premium request)
 claude-haiku-4.5        221.2k in, 4.6k out, 191.9k cached (Est. 0 Premium requests)

```

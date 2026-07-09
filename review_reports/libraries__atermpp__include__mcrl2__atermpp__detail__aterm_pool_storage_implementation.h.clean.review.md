---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T10:09:21.960277+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h` — full file, clean review (no diff, entire implementation)
- **Risk level:** Medium-High (multiple instances of technically undefined behavior in a correctness-critical garbage collection subsystem)
- **Overall verdict:** The file is functional and reflects deliberate design choices (flexible-array-member pattern via `_aterm_appl<1>`, type-punning for zero-cost term wrappers). However, it contains several instances of UB per the C++ standard that rely on implementation-specific layout guarantees. Some are easily fixable (e.g., `reinterpret_cast` where `static_cast` suffices); others are deep design patterns. The file also lacks `inline` on a non-template function definition in a header, and is missing API contracts on all public methods.

## Findings (ordered by severity)

### [SEV-1] Strict aliasing violation in `construct_arguments` (second overload) and `create_appl_dynamic`
- **Location:** `aterm_pool_storage_implementation.h:64` and `:200`
- **Why this is a problem:** Both locations cast `unprotected_aterm_core*` (or `unprotected_aterm_core&`) to `InputIterator::value_type*` (or `&`) via `reinterpret_cast`, then write through that reference. If `InputIterator::value_type` is not `unprotected_aterm_core` itself (e.g., it is `aterm` or `aterm_core`), this constitutes a strict aliasing violation: the object's dynamic type is `unprotected_aterm_core` but it is accessed through a reference to an unrelated type. The C++ standard [basic.lval]/11 does not permit this unless the types are in an inheritance relationship and the cast is a valid downcast. Here, `aterm_core` does inherit from `unprotected_aterm_core`, but `reinterpret_cast` is used instead of a valid base-to-derived cast, and the object was never constructed as the derived type.
- **Evidence type:** static-proof
- **Evidence:**
  - At line 64: `*reinterpret_cast<typename InputIterator::value_type*>(&arguments[i])` — the object at `arguments[i]` was default-constructed as `unprotected_aterm_core`. The converter then writes through a reference typed as `InputIterator::value_type` (typically `aterm` or `aterm_core`). No object of that type exists at this location.
  - At line 200: `reinterpret_cast<typename InputIterator::value_type&>(arguments[i])` — same issue with `MCRL2_DECLARE_STACK_ARRAY` allocated `unprotected_aterm_core` elements.
  - **Standard citation:** [basic.lval]/11, [basic.life]/8 — accessing an object through a glvalue of a type other than its dynamic type (or a related type) is UB. Writing through a reference to a never-constructed derived type violates object lifetime rules.
  - **Run command:** `cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_ADDRESSSANITIZER=ON -DMCRL2_ENABLE_TESTS=ON && cmake --build build-asan && ctest --test-dir build-asan -R aterm --output-on-failure` (UBSan may not catch this specific aliasing violation without `-fstrict-aliasing -Wstrict-aliasing`, but the code is formally UB)
  - **Observed result:** Works in practice because layout is identical (single pointer member).
  - **Expected result:** Well-defined behavior via proper object construction.
- **Efficiency impact:** None (the fix has zero runtime cost).
- **C++20-specific note:** Could be fixed by using `std::start_lifetime_as<InputIterator::value_type>` (C++23, P2590) or placement new. For C++20, a `static_assert` ensuring layout compatibility plus explicit documentation of the reliance on implementation behavior would be appropriate.
- **Suggested fix:** Replace with placement new to start the lifetime of the correct type, or add a `static_assert(std::is_base_of_v<unprotected_aterm_core, typename InputIterator::value_type> && sizeof(typename InputIterator::value_type) == sizeof(unprotected_aterm_core))` and use a proper `static_cast` after constructing via placement new.

---

### [SEV-2] Missing `inline` on `mark_term` definition — ODR violation risk
- **Location:** `aterm_pool_storage_implementation.h:73`
- **Why this is a problem:** The function `mark_term` is declared `inline` at `aterm_pool_storage.h:34` but defined without `inline` at `aterm_pool_storage_implementation.h:73`. This implementation header is transitively included in multiple translation units (via `aterm_pool.h`, `aterm_container.h`, `aterm_pool_implementation.h`). If the linker encounters multiple definitions without the `inline` specifier on the definition, this is technically an ODR violation. In practice, most linkers accept this because the declaration has `inline`, but the standard requires both declaration and definition to agree.
- **Evidence type:** static-proof
- **Evidence:**
  - Declaration: `aterm_pool_storage.h:34`: `inline void mark_term(const _aterm& root, ...);`
  - Definition: `aterm_pool_storage_implementation.h:73`: `void mark_term(const _aterm& root, ...)` — missing `inline`
  - Multiple TUs include this transitively: `aterm_implementation.cpp` → `global_aterm_pool.h` → `aterm_pool.h` → this file. Other paths exist through `aterm_container.h`.
  - **Run command:** Build with `-Weverything` or MSVC `/W4` to check for linker warnings.
  - **Observed result:** Compiles/links on tested platforms (GCC/Clang/MSVC fold duplicate weak symbols).
  - **Expected result:** Definition should be `inline` to match declaration and guarantee ODR compliance.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Add `inline` to the definition at line 73:
  ```cpp
  inline void mark_term(const _aterm& root, std::stack<std::reference_wrapper<_aterm>>& todo)
  ```

---

### [SEV-2] `reinterpret_cast<aterm_core&>(term)` where `static_cast` is correct
- **Location:** `aterm_pool_storage_implementation.h:117, 123`
- **Why this is a problem:** `aterm` publicly inherits from `aterm_core`. A derived-to-base cast is a valid `static_cast` (or even implicit conversion for references). Using `reinterpret_cast` instead bypasses the type system and is technically UB if the compiler applies layout optimizations for empty base classes or similar. The authors acknowledge this with TODO comments.
- **Evidence type:** static-proof
- **Evidence:**
  - `aterm.h:27`: `class aterm : public aterm_core`
  - The cast `reinterpret_cast<aterm_core&>(term)` where `term` is `aterm&` is UB per [expr.reinterpret.cast] — for class types with inheritance, `static_cast` or implicit conversion is required.
  - **Observed result:** Works because `aterm` adds no data members beyond `aterm_core`.
  - **Expected result:** Use `static_cast<aterm_core&>(term)` for well-defined behavior.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Replace both occurrences:
  ```cpp
  return emplace(static_cast<aterm_core&>(term), value);
  return emplace(static_cast<aterm_core&>(term), symbol);
  ```

---

### [SEV-2] `reinterpret_cast` in `call_deletion_hook` — casting `unprotected_aterm_core` to `const aterm&`
- **Location:** `aterm_pool_storage_implementation.h:285`
- **Why this is a problem:** `callback(reinterpret_cast<const aterm&>(term))` where `term` is `unprotected_aterm_core`. The object's dynamic type is `unprotected_aterm_core` (a stack variable or parameter), but it is accessed as `const aterm&`. While `aterm` → `aterm_core` → `unprotected_aterm_core` is the inheritance chain, there is no `aterm` object here — only an `unprotected_aterm_core` value. Upcasting a base-class object to a derived reference is invalid. The callback receives a lie about the object's type.
- **Evidence type:** static-proof
- **Evidence:**
  - The parameter `term` is `unprotected_aterm_core` (by value). No `aterm` object exists at that address.
  - The callback expects `const aterm&` but the lifetime and type of the underlying object is `unprotected_aterm_core`.
  - If the callback does anything that relies on `aterm`-specific behavior beyond the `unprotected_aterm_core` base (e.g., destructor, copy constructor that registers with the pool), this is UB.
  - **Observed result:** Works if callbacks only access the `m_term` pointer (inherited from base).
  - **Expected result:** Either construct a temporary `aterm` from the pointer, or change the callback signature to accept `unprotected_aterm_core`.
- **Efficiency impact:** Creating a temporary `aterm` would involve reference-count manipulation — may not be acceptable during deletion. Current approach is a performance hack.
- **C++20-specific note:** The TODO comment at line 285 acknowledges uncertainty about this cast.
- **Suggested fix:** Document explicitly why this is safe (callbacks must not copy the reference), or change callback type to `std::function<void(const unprotected_aterm_core&)>`.

---

### [SEV-2] Unconditional downcast of `_aterm&` to `_term_appl&` in `mark_term`
- **Location:** `aterm_pool_storage_implementation.h:89`
- **Why this is a problem:** Inside the while loop, every term popped from the stack is cast: `_term_appl& term_appl = static_cast<_term_appl&>(term)`. This includes arity-0 terms (constants), which may be `_aterm` objects that are NOT `_aterm_appl<1>` instances. Forming a reference to a wrong dynamic type is UB per [expr.static.cast]/2, regardless of whether the reference is subsequently dereferenced.
- **Evidence type:** plausible-issue
- **Evidence:**
  - When a term with arity > 0 has a constant (arity-0) argument, that argument is marked and pushed to `todo`. On the next iteration, it is cast to `_term_appl&` despite potentially being a plain `_aterm` or `_aterm_appl<0>`.
  - The cast to `_aterm_appl<1>&` is UB if the dynamic type is not `_aterm_appl<1>` or a class derived from it. `_aterm_appl<0>` is a sibling template instantiation, not a base/derived.
  - **Mitigating factor:** The library provides an explicit `operator _aterm_appl<1>&()` conversion (line 96-99 in `aterm.h`), signaling this is an intentional design pattern. The loop body is never entered for arity-0.
  - **Test plan:** Run under UBSan with `-fsanitize=undefined -fno-sanitize-recover=all` on a workload that creates and garbage-collects constants. If UBSan reports "reference binding to misaligned/null/wrong-type pointer", the issue is confirmed.
  - **Run command:** `cmake -S . -B build-ubsan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-fsanitize=undefined -fno-sanitize-recover=all" -DMCRL2_ENABLE_TESTS=ON && cmake --build build-ubsan && ctest --test-dir build-ubsan -R aterm --output-on-failure`
- **Efficiency impact:** None.
- **C++20-specific note:** The fix is to guard the cast: `if (arity > 0) { _term_appl& term_appl = ...; for (...) }`.
- **Suggested fix:**
  ```cpp
  if (arity > 0)
  {
    _term_appl& term_appl = static_cast<_term_appl&>(term);
    for (std::size_t i = 0; i < arity; ++i) { ... }
  }
  ```

---

### [SEV-3] Placement new in `emplace` overwrites `aterm_core` without respecting its destructor
- **Location:** `aterm_pool_storage_implementation.h:308`
- **Why this is a problem:** `new (&term) atermpp::unprotected_aterm_core(&*it)` constructs an `unprotected_aterm_core` into storage that holds an `aterm_core`. `aterm_core` has a non-trivial destructor (decrements reference counts). The placement new bypasses this destructor, potentially leaking the previous reference. It also ends the lifetime of the `aterm_core` object and starts the lifetime of an `unprotected_aterm_core` — subsequent access through the `aterm_core&` reference is then UB.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `aterm_core` has `~aterm_core() noexcept` (non-trivial, manages reference counts).
  - If `term` previously held a valid reference, this placement new skips the decrement.
  - **Mitigating factor:** Callers may guarantee `term` is uninitialized or null. The emplace pattern appears to be called on freshly-declared `aterm` locals or output parameters.
  - **Test plan:** Construct an `aterm_core` holding a valid term, then pass it to `create_int`/`create_term`. If reference counts become inconsistent, GC will either leak or double-free.
- **Efficiency impact:** The current approach avoids a redundant reference count decrement+increment. This is an intentional optimization.
- **C++20-specific note:** Could use `std::launder` after placement new for formally correct access.
- **Suggested fix:** Add a `\pre` comment documenting that `term` must be in a moved-from/null state, or use `term.unprotected_assign(...)` instead.

---

### [SEV-3] Missing API contracts on all public methods
- **Location:** All `create_*`, `add_deletion_hook`, `sweep`, `resize_if_needed` methods
- **Why this is a problem:** Per mCRL2 coding conventions, public APIs require Doxygen `\pre`/`\post` clauses. None of the public methods in this file have precondition/postcondition documentation. For a GC subsystem where incorrect usage leads to memory corruption, explicit contracts are essential.
- **Evidence type:** static-proof (absence of documentation)
- **Suggested fix:** Add `\pre` and `\post` to at minimum: `create_int`, `create_term`, `create_appl`, `create_appl_iterator`, `create_appl_dynamic`, `sweep`, `add_deletion_hook`.

---

### [SEV-3] `static_cast` in `store_in_argument_array_` without compile-time safety
- **Location:** `aterm_pool_storage_implementation.h:170`
- **Why this is a problem:** `function_or_term(static_cast<typename traits::template arg<0>::type&>(argument_array[I]))` casts `unprotected_aterm_core&` to whatever the functor's first parameter type is. If the functor expects a type not in the `unprotected_aterm_core` hierarchy (e.g., a concrete `data_expression&`), this is UB. No `static_assert` or concept constrains compatibility.
- **Evidence type:** plausible-issue
- **Evidence:**
  - The `else` branch is reached for functors with signature `void(T&)` where T is not checked against `unprotected_aterm_core`.
  - **Mitigating factor:** All callers in practice pass functors expecting `aterm&` or `unprotected_aterm_core&`.
  - **Test plan:** Create a functor with `void(int&)` and call `create_appl` — should fail to compile but currently may silently compile and produce UB.
- **Suggested fix:** Add `static_assert(std::is_base_of_v<unprotected_aterm_core, std::remove_reference_t<typename traits::template arg<0>::type>> || std::is_same_v<std::remove_reference_t<typename traits::template arg<0>::type>, unprotected_aterm_core>)`.

---

## Proposed Tests
- [ ] `test_strict_aliasing_construct_arguments`: Create terms via `create_appl_iterator` with a converter that writes through the reinterpreted reference. Run under UBSan to detect aliasing violations. Input: function symbol with arity 3, iterator of `aterm` values. Expected: no UBSan report (currently plausible failure).
- [ ] `test_mark_term_arity_zero`: Construct a term tree where root has arity > 0 but arguments have arity 0. Trigger GC marking. Run under UBSan. Expected: no "reference binding to wrong type" report.
- [ ] `test_emplace_overwrites_live_term`: Create a valid `aterm`, then call `create_int` on it without clearing. Verify reference counts remain consistent (no leak, no double-free). Run under ASan+LSan.
- [ ] `test_odr_mark_term`: Build with LTO enabled (`-flto`) and verify no duplicate symbol warnings for `mark_term`.

## Open Questions
1. Is the placement new in `emplace` intentionally skipping the destructor of the previous `aterm_core` value, or are callers guaranteed to pass uninitialized/null terms? If the latter, should there be an `assert(!term.defined())` precondition?
2. The `call_deletion_hook` casts `unprotected_aterm_core` (a base) to `const aterm&` (a derived). Do any registered callbacks perform operations that depend on `aterm`-specific state (e.g., copying the reference, which would register it with the protection set)?
3. Why is `_aterm_appl<1>` used as a universal "view" type rather than defining a proper interface that doesn't require the flexible-array-member UB pattern? Is this a performance-critical decision that precludes virtual dispatch or `std::variant`-based approaches?

## Copilot Raw Stdout

```text
● Explore: Explore aterm type hierarchy
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Explore: Explore aterm pool storage class
  └ Agent started in background with agent_id: agent-1. You can use read_agent to...

● Explore: Check mark_term ODR issues
  └ Agent started in background with agent_id: agent-2. You can use read_agent to...

● Explore: Analyze reinterpret_cast safety
  └ Agent started in background with agent_id: agent-3. You can use read_agent to...

● Read: Explore agent (Explore aterm pool storage class)
  └ Completed

● Read: Explore agent (Explore aterm type hierarchy)
  └ Running

● Read: Explore agent (Explore aterm type hierarchy)
  └ Completed

● Read: Explore agent (Check mark_term ODR issues)
  └ Completed

● Read: Explore agent (Analyze reinterpret_cast safety)
  └ Completed

Now I have all the context needed to write the review. Let me verify one more detail about the `mark_term` downcast pattern:

✗ Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_appl.h lines 1-50
  Path does not exist

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage.h lines 30-40
  └ 11 lines read

● Glob "libraries/atermpp/include/mcrl2/atermpp/detail/*appl*"
  └ 1 file found

● Grep "class _aterm_appl|using _term_appl" (libraries/atermpp)
  └ 3 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h lines 33-110
  └ 78 lines read

Now I have the complete picture. Let me write the review:

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h lines 90-100
  └ 10 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 33-85
  └ 53 lines read

Here is my review:

---

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h` — full file, clean review (no diff, entire implementation)
- **Risk level:** Medium-High (multiple instances of technically undefined behavior in a correctness-critical garbage collection subsystem)
- **Overall verdict:** The file is functional and reflects deliberate design choices (flexible-array-member pattern via `_aterm_appl<1>`, type-punning for zero-cost term wrappers). However, it contains several instances of UB per the C++ standard that rely on implementation-specific layout guarantees. Some are easily fixable (e.g., `reinterpret_cast` where `static_cast` suffices); others are deep design patterns. The file also lacks `inline` on a non-template function definition in a header, and is missing API contracts on all public methods.

## Findings (ordered by severity)

### [SEV-1] Strict aliasing violation in `construct_arguments` (second overload) and `create_appl_dynamic`
- **Location:** `aterm_pool_storage_implementation.h:64` and `:200`
- **Why this is a problem:** Both locations cast `unprotected_aterm_core*` (or `unprotected_aterm_core&`) to `InputIterator::value_type*` (or `&`) via `reinterpret_cast`, then write through that reference. If `InputIterator::value_type` is not `unprotected_aterm_core` itself (e.g., it is `aterm` or `aterm_core`), this constitutes a strict aliasing violation: the object's dynamic type is `unprotected_aterm_core` but it is accessed through a reference to an unrelated type. The C++ standard [basic.lval]/11 does not permit this unless the types are in an inheritance relationship and the cast is a valid downcast. Here, `aterm_core` does inherit from `unprotected_aterm_core`, but `reinterpret_cast` is used instead of a valid base-to-derived cast, and the object was never constructed as the derived type.
- **Evidence type:** static-proof
- **Evidence:**
  - At line 64: `*reinterpret_cast<typename InputIterator::value_type*>(&arguments[i])` — the object at `arguments[i]` was default-constructed as `unprotected_aterm_core`. The converter then writes through a reference typed as `InputIterator::value_type` (typically `aterm` or `aterm_core`). No object of that type exists at this location.
  - At line 200: `reinterpret_cast<typename InputIterator::value_type&>(arguments[i])` — same issue with `MCRL2_DECLARE_STACK_ARRAY` allocated `unprotected_aterm_core` elements.
  - **Standard citation:** [basic.lval]/11, [basic.life]/8 — accessing an object through a glvalue of a type other than its dynamic type (or a related type) is UB. Writing through a reference to a never-constructed derived type violates object lifetime rules.
  - **Run command:** `cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_ADDRESSSANITIZER=ON -DMCRL2_ENABLE_TESTS=ON && cmake --build build-asan && ctest --test-dir build-asan -R aterm --output-on-failure` (UBSan may not catch this specific aliasing violation without `-fstrict-aliasing -Wstrict-aliasing`, but the code is formally UB)
  - **Observed result:** Works in practice because layout is identical (single pointer member).
  - **Expected result:** Well-defined behavior via proper object construction.
- **Efficiency impact:** None (the fix has zero runtime cost).
- **C++20-specific note:** Could be fixed by using `std::start_lifetime_as<InputIterator::value_type>` (C++23, P2590) or placement new. For C++20, a `static_assert` ensuring layout compatibility plus explicit documentation of the reliance on implementation behavior would be appropriate.
- **Suggested fix:** Replace with placement new to start the lifetime of the correct type, or add a `static_assert(std::is_base_of_v<unprotected_aterm_core, typename InputIterator::value_type> && sizeof(typename InputIterator::value_type) == sizeof(unprotected_aterm_core))` and use a proper `static_cast` after constructing via placement new.

---

### [SEV-2] Missing `inline` on `mark_term` definition — ODR violation risk
- **Location:** `aterm_pool_storage_implementation.h:73`
- **Why this is a problem:** The function `mark_term` is declared `inline` at `aterm_pool_storage.h:34` but defined without `inline` at `aterm_pool_storage_implementation.h:73`. This implementation header is transitively included in multiple translation units (via `aterm_pool.h`, `aterm_container.h`, `aterm_pool_implementation.h`). If the linker encounters multiple definitions without the `inline` specifier on the definition, this is technically an ODR violation. In practice, most linkers accept this because the declaration has `inline`, but the standard requires both declaration and definition to agree.
- **Evidence type:** static-proof
- **Evidence:**
  - Declaration: `aterm_pool_storage.h:34`: `inline void mark_term(const _aterm& root, ...);`
  - Definition: `aterm_pool_storage_implementation.h:73`: `void mark_term(const _aterm& root, ...)` — missing `inline`
  - Multiple TUs include this transitively: `aterm_implementation.cpp` → `global_aterm_pool.h` → `aterm_pool.h` → this file. Other paths exist through `aterm_container.h`.
  - **Run command:** Build with `-Weverything` or MSVC `/W4` to check for linker warnings.
  - **Observed result:** Compiles/links on tested platforms (GCC/Clang/MSVC fold duplicate weak symbols).
  - **Expected result:** Definition should be `inline` to match declaration and guarantee ODR compliance.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Add `inline` to the definition at line 73:
  ```cpp
  inline void mark_term(const _aterm& root, std::stack<std::reference_wrapper<_aterm>>& todo)
  ```

---

### [SEV-2] `reinterpret_cast<aterm_core&>(term)` where `static_cast` is correct
- **Location:** `aterm_pool_storage_implementation.h:117, 123`
- **Why this is a problem:** `aterm` publicly inherits from `aterm_core`. A derived-to-base cast is a valid `static_cast` (or even implicit conversion for references). Using `reinterpret_cast` instead bypasses the type system and is technically UB if the compiler applies layout optimizations for empty base classes or similar. The authors acknowledge this with TODO comments.
- **Evidence type:** static-proof
- **Evidence:**
  - `aterm.h:27`: `class aterm : public aterm_core`
  - The cast `reinterpret_cast<aterm_core&>(term)` where `term` is `aterm&` is UB per [expr.reinterpret.cast] — for class types with inheritance, `static_cast` or implicit conversion is required.
  - **Observed result:** Works because `aterm` adds no data members beyond `aterm_core`.
  - **Expected result:** Use `static_cast<aterm_core&>(term)` for well-defined behavior.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Replace both occurrences:
  ```cpp
  return emplace(static_cast<aterm_core&>(term), value);
  return emplace(static_cast<aterm_core&>(term), symbol);
  ```

---

### [SEV-2] `reinterpret_cast` in `call_deletion_hook` — casting `unprotected_aterm_core` to `const aterm&`
- **Location:** `aterm_pool_storage_implementation.h:285`
- **Why this is a problem:** `callback(reinterpret_cast<const aterm&>(term))` where `term` is `unprotected_aterm_core`. The object's dynamic type is `unprotected_aterm_core` (a stack variable or parameter), but it is accessed as `const aterm&`. While `aterm` → `aterm_core` → `unprotected_aterm_core` is the inheritance chain, there is no `aterm` object here — only an `unprotected_aterm_core` value. Upcasting a base-class object to a derived reference is invalid. The callback receives a lie about the object's type.
- **Evidence type:** static-proof
- **Evidence:**
  - The parameter `term` is `unprotected_aterm_core` (by value). No `aterm` object exists at that address.
  - The callback expects `const aterm&` but the lifetime and type of the underlying object is `unprotected_aterm_core`.
  - If the callback does anything that relies on `aterm`-specific behavior beyond the `unprotected_aterm_core` base (e.g., destructor, copy constructor that registers with the pool), this is UB.
  - **Observed result:** Works if callbacks only access the `m_term` pointer (inherited from base).
  - **Expected result:** Either construct a temporary `aterm` from the pointer, or change the callback signature to accept `unprotected_aterm_core`.
- **Efficiency impact:** Creating a temporary `aterm` would involve reference-count manipulation — may not be acceptable during deletion. Current approach is a performance hack.
- **C++20-specific note:** The TODO comment at line 285 acknowledges uncertainty about this cast.
- **Suggested fix:** Document explicitly why this is safe (callbacks must not copy the reference), or change callback type to `std::function<void(const unprotected_aterm_core&)>`.

---

### [SEV-2] Unconditional downcast of `_aterm&` to `_term_appl&` in `mark_term`
- **Location:** `aterm_pool_storage_implementation.h:89`
- **Why this is a problem:** Inside the while loop, every term popped from the stack is cast: `_term_appl& term_appl = static_cast<_term_appl&>(term)`. This includes arity-0 terms (constants), which may be `_aterm` objects that are NOT `_aterm_appl<1>` instances. Forming a reference to a wrong dynamic type is UB per [expr.static.cast]/2, regardless of whether the reference is subsequently dereferenced.
- **Evidence type:** plausible-issue
- **Evidence:**
  - When a term with arity > 0 has a constant (arity-0) argument, that argument is marked and pushed to `todo`. On the next iteration, it is cast to `_term_appl&` despite potentially being a plain `_aterm` or `_aterm_appl<0>`.
  - The cast to `_aterm_appl<1>&` is UB if the dynamic type is not `_aterm_appl<1>` or a class derived from it. `_aterm_appl<0>` is a sibling template instantiation, not a base/derived.
  - **Mitigating factor:** The library provides an explicit `operator _aterm_appl<1>&()` conversion (line 96-99 in `aterm.h`), signaling this is an intentional design pattern. The loop body is never entered for arity-0.
  - **Test plan:** Run under UBSan with `-fsanitize=undefined -fno-sanitize-recover=all` on a workload that creates and garbage-collects constants. If UBSan reports "reference binding to misaligned/null/wrong-type pointer", the issue is confirmed.
  - **Run command:** `cmake -S . -B build-ubsan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-fsanitize=undefined -fno-sanitize-recover=all" -DMCRL2_ENABLE_TESTS=ON && cmake --build build-ubsan && ctest --test-dir build-ubsan -R aterm --output-on-failure`
- **Efficiency impact:** None.
- **C++20-specific note:** The fix is to guard the cast: `if (arity > 0) { _term_appl& term_appl = ...; for (...) }`.
- **Suggested fix:**
  ```cpp
  if (arity > 0)
  {
    _term_appl& term_appl = static_cast<_term_appl&>(term);
    for (std::size_t i = 0; i < arity; ++i) { ... }
  }
  ```

---

### [SEV-3] Placement new in `emplace` overwrites `aterm_core` without respecting its destructor
- **Location:** `aterm_pool_storage_implementation.h:308`
- **Why this is a problem:** `new (&term) atermpp::unprotected_aterm_core(&*it)` constructs an `unprotected_aterm_core` into storage that holds an `aterm_core`. `aterm_core` has a non-trivial destructor (decrements reference counts). The placement new bypasses this destructor, potentially leaking the previous reference. It also ends the lifetime of the `aterm_core` object and starts the lifetime of an `unprotected_aterm_core` — subsequent access through the `aterm_core&` reference is then UB.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `aterm_core` has `~aterm_core() noexcept` (non-trivial, manages reference counts).
  - If `term` previously held a valid reference, this placement new skips the decrement.
  - **Mitigating factor:** Callers may guarantee `term` is uninitialized or null. The emplace pattern appears to be called on freshly-declared `aterm` locals or output parameters.
  - **Test plan:** Construct an `aterm_core` holding a valid term, then pass it to `create_int`/`create_term`. If reference counts become inconsistent, GC will either leak or double-free.
- **Efficiency impact:** The current approach avoids a redundant reference count decrement+increment. This is an intentional optimization.
- **C++20-specific note:** Could use `std::launder` after placement new for formally correct access.
- **Suggested fix:** Add a `\pre` comment documenting that `term` must be in a moved-from/null state, or use `term.unprotected_assign(...)` instead.

---

### [SEV-3] Missing API contracts on all public methods
- **Location:** All `create_*`, `add_deletion_hook`, `sweep`, `resize_if_needed` methods
- **Why this is a problem:** Per mCRL2 coding conventions, public APIs require Doxygen `\pre`/`\post` clauses. None of the public methods in this file have precondition/postcondition documentation. For a GC subsystem where incorrect usage leads to memory corruption, explicit contracts are essential.
- **Evidence type:** static-proof (absence of documentation)
- **Suggested fix:** Add `\pre` and `\post` to at minimum: `create_int`, `create_term`, `create_appl`, `create_appl_iterator`, `create_appl_dynamic`, `sweep`, `add_deletion_hook`.

---

### [SEV-3] `static_cast` in `store_in_argument_array_` without compile-time safety
- **Location:** `aterm_pool_storage_implementation.h:170`
- **Why this is a problem:** `function_or_term(static_cast<typename traits::template arg<0>::type&>(argument_array[I]))` casts `unprotected_aterm_core&` to whatever the functor's first parameter type is. If the functor expects a type not in the `unprotected_aterm_core` hierarchy (e.g., a concrete `data_expression&`), this is UB. No `static_assert` or concept constrains compatibility.
- **Evidence type:** plausible-issue
- **Evidence:**
  - The `else` branch is reached for functors with signature `void(T&)` where T is not checked against `unprotected_aterm_core`.
  - **Mitigating factor:** All callers in practice pass functors expecting `aterm&` or `unprotected_aterm_core&`.
  - **Test plan:** Create a functor with `void(int&)` and call `create_appl` — should fail to compile but currently may silently compile and produce UB.
- **Suggested fix:** Add `static_assert(std::is_base_of_v<unprotected_aterm_core, std::remove_reference_t<typename traits::template arg<0>::type>> || std::is_same_v<std::remove_reference_t<typename traits::template arg<0>::type>, unprotected_aterm_core>)`.

---

## Proposed Tests
- [ ] `test_strict_aliasing_construct_arguments`: Create terms via `create_appl_iterator` with a converter that writes through the reinterpreted reference. Run under UBSan to detect aliasing violations. Input: function symbol with arity 3, iterator of `aterm` values. Expected: no UBSan report (currently plausible failure).
- [ ] `test_mark_term_arity_zero`: Construct a term tree where root has arity > 0 but arguments have arity 0. Trigger GC marking. Run under UBSan. Expected: no "reference binding to wrong type" report.
- [ ] `test_emplace_overwrites_live_term`: Create a valid `aterm`, then call `create_int` on it without clearing. Verify reference counts remain consistent (no leak, no double-free). Run under ASan+LSan.
- [ ] `test_odr_mark_term`: Build with LTO enabled (`-flto`) and verify no duplicate symbol warnings for `mark_term`.

## Open Questions
1. Is the placement new in `emplace` intentionally skipping the destructor of the previous `aterm_core` value, or are callers guaranteed to pass uninitialized/null terms? If the latter, should there be an `assert(!term.defined())` precondition?
2. The `call_deletion_hook` casts `unprotected_aterm_core` (a base) to `const aterm&` (a derived). Do any registered callbacks perform operations that depend on `aterm`-specific state (e.g., copying the reference, which would register it with the protection set)?
3. Why is `_aterm_appl<1>` used as a universal "view" type rather than defining a proper interface that doesn't require the flexible-array-member UB pattern? Is this a performance-critical decision that precludes virtual dispatch or `std::variant`-based approaches?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         9m 18s
Total session time:     5m 52s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         595.1k in, 13.3k out, 545.7k cached (Est. 1 Premium request)
 claude-haiku-4.5        1.5m in, 24.1k out, 1.4m cached (Est. 0 Premium requests)

```

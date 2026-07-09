---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool_implementation.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T12:23:06.788013+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool_implementation.h` — the per-thread aterm creation and GC protection set management layer.
- **Risk level:** HIGH — contains `reinterpret_cast`-based type punning that constitutes undefined behavior under C++ strict aliasing rules, in a project that mandates "No undefined behavior."
- **Overall verdict:** The synchronization design is sound (shared/exclusive lock pattern correctly prevents data races), but the `create_appl_index` function contains formally undefined behavior via `reinterpret_cast` aliasing violations that are not guarded by any `static_assert` for the `std::pair` case. The output-parameter aliasing pattern is fragile. Missing API contracts violate project conventions.

## Findings (ordered by severity)

### [SEV-1] Strict aliasing violation via `reinterpret_cast` of `std::pair` over `unprotected_aterm_core` array

- **Location:** `thread_aterm_pool_implementation.h:91-92`
- **Why this is a problem:**  
  The code reinterprets two consecutive `unprotected_aterm_core` objects (each containing a single `const _aterm*`) as a `std::pair<core::identifier_string, sort_expression>`. Per [basic.lval]/11, accessing an object through a glvalue of an incompatible type is undefined behavior. `std::pair<T1, T2>` is entirely unrelated to `unprotected_aterm_core` in the type hierarchy.  
  Additionally, no `static_assert` validates that `sizeof(std::pair<INDEX_TYPE_first, INDEX_TYPE_second>) == 2 * sizeof(unprotected_aterm_core)` or that the pair has no inter-member padding. While this likely holds on all target platforms (both pair members are pointer-sized standard-layout types), the assumption is unverified at compile time.
- **Evidence type:** static-proof
- **Evidence:**
  - Per C++20 [basic.lval]/11, the permitted access types for an object of dynamic type `unprotected_aterm_core` do NOT include `std::pair<identifier_string, sort_expression>`.
  - The `IsATerm` concept (concepts.h:20-30) enforces standard-layout and size constraints on individual aterm types but has no coverage for `std::pair` composites.
  - UBSan's `-fsanitize=undefined` does not catch strict aliasing violations at runtime; this requires `-fstrict-aliasing` compiler optimizations to manifest as miscompilation (reordered/eliminated reads).
  - Run command: Build with `cmake -DCMAKE_CXX_FLAGS="-O2 -Wall -Wstrict-aliasing=1"` and inspect diagnostics.
  - Observed result: Formal UB; compiler may assume `std::pair*` and `unprotected_aterm_core*` don't alias.
  - Expected result: Well-defined behavior via `std::bit_cast`, `memcpy`, or constructing the actual pair.
- **Efficiency impact:** The reinterpret_cast avoids constructing a `std::pair` copy (the commented-out "more elegant" code), saving one pair construction per indexed term creation. This is a micro-optimization.
- **C++20-specific note:** C++20 `std::bit_cast` would be type-safe for trivially-copyable types of same size, but requires a value (not reference) result. A `memcpy` into a local `INDEX_TYPE` is the correct fix that any optimizer will elide.
- **Suggested fix:**
  ```cpp
  INDEX_TYPE key;
  static_assert(sizeof(key) == sizeof...(arguments) * sizeof(unprotected_aterm_core),
                "INDEX_TYPE must match argument array layout");
  std::memcpy(&key, &argument_array[0], sizeof(key));
  m_pool.create_int(term, atermpp::detail::index_traits<Term, INDEX_TYPE, 2>::insert(key));
  ```

---

### [SEV-2] Strict aliasing violation via `reinterpret_cast` of single `unprotected_aterm_core` as derived type

- **Location:** `thread_aterm_pool_implementation.h:79-80`
- **Why this is a problem:**  
  `*reinterpret_cast<INDEX_TYPE*>(&(argument_array[0]))` accesses an object of dynamic type `unprotected_aterm_core` through a pointer to `INDEX_TYPE` (a derived class like `identifier_string`/`aterm_string`). While the `IsATerm` concept guarantees size and layout compatibility, the C++ object model requires that a derived-type access is only valid when a derived object actually exists at that address. No `INDEX_TYPE` object was constructed in `argument_array[0]`.
- **Evidence type:** static-proof
- **Evidence:**
  - [basic.lval]/11: accessing a base-class object through a derived-class pointer is not in the permitted list (the inverse — accessing through a base-class type — is permitted).
  - The `IsATerm` concept's comment at concepts.h:28 ("A pointer to a standard-layout class may be converted to a pointer to its first non-static data member") references [class.mem] pointer-interconvertibility, which governs struct↔member conversion, not base↔derived object identity.
  - Reproducer: A compiler exploiting this via devirtualization or reordering is theoretically possible under LTO.
  - Run command: `cmake --build build --target mcrl2_test_atermpp` under aggressive optimization (`-O3 -flto`).
- **Efficiency impact:** Saves one copy-construction of an aterm type (pointer copy + potential refcount bump).
- **C++20-specific note:** Same fix as SEV-1 applies: use `memcpy` or construct the derived type properly via `static_cast` from a properly-formed base reference.
- **Suggested fix:**
  ```cpp
  INDEX_TYPE key;
  std::memcpy(&key, &argument_array[0], sizeof(key));
  m_pool.create_int(term, atermpp::detail::index_traits<Term, INDEX_TYPE, 1>::insert(key));
  ```

---

### [SEV-2] Output parameter `term` aliased as both intermediate result and final argument to `create_appl`

- **Location:** `thread_aterm_pool_implementation.h:80-81` and `91-93`
- **Why this is a problem:**  
  ```cpp
  m_pool.create_int(term, ...);  // writes int-aterm into term
  added = m_pool.create_appl(term, sym, argument_array[0], term);  // term is BOTH output AND input
  ```
  `term` is passed by reference as the output parameter AND as one of the variadic arguments. The correctness depends on `emplace` fully consuming all argument values (reading `term.m_term`) before the placement `new (&term)` overwrites the output. This ordering dependency is implicit and undocumented.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Current `emplace` implementation (aterm_pool_storage_implementation.h:411-427) first calls `m_term_set.emplace(...)` which hashes/copies all arguments, THEN writes `term`. The ordering is currently safe.
  - However, any refactoring of `emplace` or the hash set that defers argument evaluation (e.g., lazy construction, move semantics on arguments) would silently corrupt `term`.
  - No assertion or comment documents this ordering requirement.
  - Test plan: Instrument `emplace` to write `term` before reading arguments; observe corruption.
- **Efficiency impact:** None; this is a correctness/maintainability concern.
- **Suggested fix:** Store the int-aterm result in a separate local variable:
  ```cpp
  aterm index_term;
  m_pool.create_int(index_term, atermpp::detail::index_traits<Term, INDEX_TYPE, 1>::insert(...));
  added = m_pool.create_appl(term, sym, argument_array[0], index_term);
  ```

---

### [SEV-3] Missing precondition assertion and documentation on `mark()` and `protection_set_size()`

- **Location:** `thread_aterm_pool_implementation.h:159` and `210`
- **Why this is a problem:**  
  Both `mark()` and `protection_set_size()` iterate per-thread protection sets (`m_variables`, `m_containers`) WITHOUT acquiring any lock. Their correctness requires that an exclusive lock is held by the caller (confirmed: called from `collect_impl` which holds the exclusive lock). This invariant is neither documented with Doxygen `\pre` nor enforced with an assertion.
- **Evidence type:** static-proof (violation of project API contract conventions)
- **Evidence:**
  - Project convention (cpp.instructions.md): "Document every non-trivial public function with Doxygen comments" and "Enforce cheap preconditions with assert(...)."
  - A future caller invoking `mark()` without the exclusive lock would cause a data race with concurrent `register_variable`/`deregister_variable` calls.
- **Suggested fix:**
  ```cpp
  void thread_aterm_pool::mark()
  {
    // \pre The global shared_mutex must be exclusively locked (GC phase).
    assert(!m_shared_mutex.is_shared_locked()); // Not in shared section; exclusive must be held
    // ... existing code ...
  }
  ```

---

### [SEV-3] Missing `[[nodiscard]]` on value-returning functions

- **Location:** `thread_aterm_pool_implementation.h:24,33,210`
- **Why this is a problem:**  
  `create_function_symbol` (returns `function_symbol`) and `protection_set_size` (returns `std::size_t`) lack `[[nodiscard]]`. Discarding a function_symbol result would be a bug (wasted allocation + lost handle). Project conventions require `[[nodiscard]]` on "pure queries and must-use results."
- **Evidence type:** static-proof (convention violation)
- **Suggested fix:** Add `[[nodiscard]]` to declarations in `thread_aterm_pool.h`.

---

### [SEV-3] No Doxygen `\pre`/`\post` contracts on any public methods

- **Location:** entire file
- **Why this is a problem:**  
  None of the public methods (`create_int`, `create_term`, `create_appl`, `register_variable`, etc.) have Doxygen documentation specifying preconditions (e.g., "arguments must be valid terms"), postconditions (e.g., "`term` is a valid, GC-protected aterm"), or thread-safety guarantees. Per project conventions, "Document every non-trivial public function with Doxygen comments: `\brief`, `\param`, `\returns`, and explicit `\pre` and `\post` clauses."
- **Evidence type:** static-proof (convention violation)
- **Suggested fix:** Add minimal Doxygen to the declarations in `thread_aterm_pool.h`.

---

## Proposed Tests

- [ ] **`test_create_appl_index_pair_layout`:** Compile-time `static_assert` that `sizeof(std::pair<identifier_string, sort_expression>) == 2 * sizeof(unprotected_aterm_core)` and that `offsetof(pair, second) == sizeof(unprotected_aterm_core)`. Place in the existing atermpp test suite. Expected failure signal: compilation error if pair layout assumption is violated on any platform.
- [ ] **`test_create_appl_index_aliasing`:** Unit test that creates an indexed term where the `term` output parameter is pre-populated with a specific int-aterm, then verifies the resulting application's last argument matches the expected index value. Run: `ctest --test-dir build -R atermpp`. Expected failure signal: incorrect term structure if output/input ordering is violated.
- [ ] **`test_strict_aliasing_memcpy_equivalence`:** After applying the `memcpy` fix, verify that `index_traits::insert` receives the same key value as the `reinterpret_cast` path would produce, using a randomized set of aterm arguments. Run under `-O2 -fstrict-aliasing`. Expected failure: index mismatch.

## Open Questions

1. Has the `reinterpret_cast` path ever been tested under `-fstrict-aliasing -Wstrict-aliasing=1` with GCC/Clang at `-O2`/`-O3`? If yes, are there suppression pragmas?
2. Is the `std::pair` layout guarantee relied upon elsewhere in the codebase, or is this the only instance? If widespread, a project-wide `static_assert` or wrapper type would be warranted.
3. The `emplace` function uses `new (&term) unprotected_aterm_core(...)` which constructs a base-class object in memory that holds a derived `aterm` — is there a deliberate policy about implicit object creation / `std::launder` requirements here, or is this additional UB that the `// TODO remove reinterpret_cast` comment alludes to?

## Copilot Raw Stdout

```text
● Explore: Explore aterm pool types
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Explore: Explore reinterpret_cast safety
  └ Agent started in background with agent_id: agent-1. You can use read_agent to...

● Explore: Find mark() call context
  └ Agent started in background with agent_id: agent-2. You can use read_agent to...

● Read: Explore agent (Explore aterm pool types)
  └ Running

● Read: Explore agent (Explore reinterpret_cast safety)
  └ Running

● Read: Explore agent (Explore aterm pool types)
  └ Completed

● Read: Explore agent (Explore reinterpret_cast safety)
  └ Completed

● Read: Explore agent (Find mark() call context)
  └ Running

● Read: Explore agent (Find mark() call context)
  └ Completed

Now I have comprehensive context. Let me verify one more detail about the `std::pair` layout assumption:

● Find create_appl_index instantiations
  $ cd /workspaces/mCRL2 && grep -rn "create_appl_index\|make_term_appl_with_index" libraries/
  --include="*.h" --include="*.cpp" | grep -v "^Binary" | head -30
  └ 9 lines...

● Check for static_asserts on pair layout
  $ cd /workspaces/mCRL2 && grep -rn "static_assert.*pair\|static_assert.*INDEX" libraries/atermpp/
  --include="*.h" | head -20
  └ 1 line...

● All INDEX_TYPE usages
  $ cd /workspaces/mCRL2 && grep -rn "make_term_appl_with_index" libraries/ --include="*.h"
  --include="*.cpp" | grep -v "^Binary" | grep -v "detail/thread_aterm_pool" | grep -v
  "detail/index_traits"
  └ 5 lines...

● Check INDEX_TYPE component types
  $ cd /workspaces/mCRL2 && grep -rn "class identifier_string" libraries/core/ --include="*.h" |
  head -5 && echo "---" && grep -rn "class sort_expression" libraries/data/ --include="*.h" | head
  -5
  └ 8 lines...

● Read libraries/core/include/mcrl2/core/identifier_string.h lines 25-52
  └ 28 lines read

● Check if variable type also uses index
  $ cd /workspaces/mCRL2 && grep -n "class variable " libraries/data/include/mcrl2/data/variable.h |
   head -3 && echo "---" && grep -rn
  "make_term_appl_with_index.*variable\|make_term_appl_with_index<variable" libraries/
  --include="*.h" | head -5
  └ 4 lines...

● Find identifier_string class definition
  $ cd /workspaces/mCRL2 && grep -B5 "class identifier_string"
  libraries/core/include/mcrl2/core/identifier_string.h | head -20
  └ 13 lines...

● Check IsATerm concept
  $ cd /workspaces/mCRL2 && cat libraries/atermpp/include/mcrl2/atermpp/detail/aterm_configuration.h
   2>/dev/null | head -40 && echo "---" && grep -n "IsATerm"
  libraries/atermpp/include/mcrl2/atermpp/ -r
  └ 53 lines...

● Check aterm_string definition
  $ cd /workspaces/mCRL2 && grep -n "class aterm_string" libraries/atermpp/include/mcrl2/atermpp/ -r
   | head -5 && echo "---" && grep -A10 "class aterm_string"
  libraries/atermpp/include/mcrl2/atermpp/aterm_string.h | head -15
  └ 14 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/concepts.h lines 1-45
  └ 44 lines read

Now I have all the information needed. Let me compile the review.

● Check create_appl signature for aliasing concern
  $ cd /workspaces/mCRL2 && grep -n "create_appl"
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h | head -10
  └ 7 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
  186-205
  └ 20 lines read

● Find emplace implementation
  $ cd /workspaces/mCRL2 && grep -n "bool.*emplace"
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h | head -5
  └ 2 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
  411-460
  └ 49 lines read

✗ Quick check on pair layout assumption
  $ cd /workspaces/mCRL2 && python3 -c "
  import struct, sys
  # Check if std::pair<T,T> for pointer-sized T has padding on typical systems
  # For standard layout types with same alignment, no padding expected
  # But we should verify the assumption with a compile test...
  Permission denied and could not request permission from user

✗ Verify pair layout assumption compiles
  $ cd /workspaces/mCRL2 && cat > /tmp/check_pair_layout.cpp << 'EOF'
  #include <utility>
  #include <cstddef>
  #include <type_traits>...
  Permission denied and could not request permission from user

● Check are_terms trait
  $ cd /workspaces/mCRL2 && grep -rn "are_terms" libraries/atermpp/include/mcrl2/atermpp/detail/ |
  head -5
  └ 5 lines...

I now have all the context needed to produce the review. Here is my critical analysis:

---

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool_implementation.h` — the per-thread aterm creation and GC protection set management layer.
- **Risk level:** HIGH — contains `reinterpret_cast`-based type punning that constitutes undefined behavior under C++ strict aliasing rules, in a project that mandates "No undefined behavior."
- **Overall verdict:** The synchronization design is sound (shared/exclusive lock pattern correctly prevents data races), but the `create_appl_index` function contains formally undefined behavior via `reinterpret_cast` aliasing violations that are not guarded by any `static_assert` for the `std::pair` case. The output-parameter aliasing pattern is fragile. Missing API contracts violate project conventions.

## Findings (ordered by severity)

### [SEV-1] Strict aliasing violation via `reinterpret_cast` of `std::pair` over `unprotected_aterm_core` array

- **Location:** `thread_aterm_pool_implementation.h:91-92`
- **Why this is a problem:**  
  The code reinterprets two consecutive `unprotected_aterm_core` objects (each containing a single `const _aterm*`) as a `std::pair<core::identifier_string, sort_expression>`. Per [basic.lval]/11, accessing an object through a glvalue of an incompatible type is undefined behavior. `std::pair<T1, T2>` is entirely unrelated to `unprotected_aterm_core` in the type hierarchy.  
  Additionally, no `static_assert` validates that `sizeof(std::pair<INDEX_TYPE_first, INDEX_TYPE_second>) == 2 * sizeof(unprotected_aterm_core)` or that the pair has no inter-member padding. While this likely holds on all target platforms (both pair members are pointer-sized standard-layout types), the assumption is unverified at compile time.
- **Evidence type:** static-proof
- **Evidence:**
  - Per C++20 [basic.lval]/11, the permitted access types for an object of dynamic type `unprotected_aterm_core` do NOT include `std::pair<identifier_string, sort_expression>`.
  - The `IsATerm` concept (concepts.h:20-30) enforces standard-layout and size constraints on individual aterm types but has no coverage for `std::pair` composites.
  - UBSan's `-fsanitize=undefined` does not catch strict aliasing violations at runtime; this requires `-fstrict-aliasing` compiler optimizations to manifest as miscompilation (reordered/eliminated reads).
  - Run command: Build with `cmake -DCMAKE_CXX_FLAGS="-O2 -Wall -Wstrict-aliasing=1"` and inspect diagnostics.
  - Observed result: Formal UB; compiler may assume `std::pair*` and `unprotected_aterm_core*` don't alias.
  - Expected result: Well-defined behavior via `std::bit_cast`, `memcpy`, or constructing the actual pair.
- **Efficiency impact:** The reinterpret_cast avoids constructing a `std::pair` copy (the commented-out "more elegant" code), saving one pair construction per indexed term creation. This is a micro-optimization.
- **C++20-specific note:** C++20 `std::bit_cast` would be type-safe for trivially-copyable types of same size, but requires a value (not reference) result. A `memcpy` into a local `INDEX_TYPE` is the correct fix that any optimizer will elide.
- **Suggested fix:**
  ```cpp
  INDEX_TYPE key;
  static_assert(sizeof(key) == sizeof...(arguments) * sizeof(unprotected_aterm_core),
                "INDEX_TYPE must match argument array layout");
  std::memcpy(&key, &argument_array[0], sizeof(key));
  m_pool.create_int(term, atermpp::detail::index_traits<Term, INDEX_TYPE, 2>::insert(key));
  ```

---

### [SEV-2] Strict aliasing violation via `reinterpret_cast` of single `unprotected_aterm_core` as derived type

- **Location:** `thread_aterm_pool_implementation.h:79-80`
- **Why this is a problem:**  
  `*reinterpret_cast<INDEX_TYPE*>(&(argument_array[0]))` accesses an object of dynamic type `unprotected_aterm_core` through a pointer to `INDEX_TYPE` (a derived class like `identifier_string`/`aterm_string`). While the `IsATerm` concept guarantees size and layout compatibility, the C++ object model requires that a derived-type access is only valid when a derived object actually exists at that address. No `INDEX_TYPE` object was constructed in `argument_array[0]`.
- **Evidence type:** static-proof
- **Evidence:**
  - [basic.lval]/11: accessing a base-class object through a derived-class pointer is not in the permitted list (the inverse — accessing through a base-class type — is permitted).
  - The `IsATerm` concept's comment at concepts.h:28 ("A pointer to a standard-layout class may be converted to a pointer to its first non-static data member") references [class.mem] pointer-interconvertibility, which governs struct↔member conversion, not base↔derived object identity.
  - Reproducer: A compiler exploiting this via devirtualization or reordering is theoretically possible under LTO.
  - Run command: `cmake --build build --target mcrl2_test_atermpp` under aggressive optimization (`-O3 -flto`).
- **Efficiency impact:** Saves one copy-construction of an aterm type (pointer copy + potential refcount bump).
- **C++20-specific note:** Same fix as SEV-1 applies: use `memcpy` or construct the derived type properly via `static_cast` from a properly-formed base reference.
- **Suggested fix:**
  ```cpp
  INDEX_TYPE key;
  std::memcpy(&key, &argument_array[0], sizeof(key));
  m_pool.create_int(term, atermpp::detail::index_traits<Term, INDEX_TYPE, 1>::insert(key));
  ```

---

### [SEV-2] Output parameter `term` aliased as both intermediate result and final argument to `create_appl`

- **Location:** `thread_aterm_pool_implementation.h:80-81` and `91-93`
- **Why this is a problem:**  
  ```cpp
  m_pool.create_int(term, ...);  // writes int-aterm into term
  added = m_pool.create_appl(term, sym, argument_array[0], term);  // term is BOTH output AND input
  ```
  `term` is passed by reference as the output parameter AND as one of the variadic arguments. The correctness depends on `emplace` fully consuming all argument values (reading `term.m_term`) before the placement `new (&term)` overwrites the output. This ordering dependency is implicit and undocumented.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Current `emplace` implementation (aterm_pool_storage_implementation.h:411-427) first calls `m_term_set.emplace(...)` which hashes/copies all arguments, THEN writes `term`. The ordering is currently safe.
  - However, any refactoring of `emplace` or the hash set that defers argument evaluation (e.g., lazy construction, move semantics on arguments) would silently corrupt `term`.
  - No assertion or comment documents this ordering requirement.
  - Test plan: Instrument `emplace` to write `term` before reading arguments; observe corruption.
- **Efficiency impact:** None; this is a correctness/maintainability concern.
- **Suggested fix:** Store the int-aterm result in a separate local variable:
  ```cpp
  aterm index_term;
  m_pool.create_int(index_term, atermpp::detail::index_traits<Term, INDEX_TYPE, 1>::insert(...));
  added = m_pool.create_appl(term, sym, argument_array[0], index_term);
  ```

---

### [SEV-3] Missing precondition assertion and documentation on `mark()` and `protection_set_size()`

- **Location:** `thread_aterm_pool_implementation.h:159` and `210`
- **Why this is a problem:**  
  Both `mark()` and `protection_set_size()` iterate per-thread protection sets (`m_variables`, `m_containers`) WITHOUT acquiring any lock. Their correctness requires that an exclusive lock is held by the caller (confirmed: called from `collect_impl` which holds the exclusive lock). This invariant is neither documented with Doxygen `\pre` nor enforced with an assertion.
- **Evidence type:** static-proof (violation of project API contract conventions)
- **Evidence:**
  - Project convention (cpp.instructions.md): "Document every non-trivial public function with Doxygen comments" and "Enforce cheap preconditions with assert(...)."
  - A future caller invoking `mark()` without the exclusive lock would cause a data race with concurrent `register_variable`/`deregister_variable` calls.
- **Suggested fix:**
  ```cpp
  void thread_aterm_pool::mark()
  {
    // \pre The global shared_mutex must be exclusively locked (GC phase).
    assert(!m_shared_mutex.is_shared_locked()); // Not in shared section; exclusive must be held
    // ... existing code ...
  }
  ```

---

### [SEV-3] Missing `[[nodiscard]]` on value-returning functions

- **Location:** `thread_aterm_pool_implementation.h:24,33,210`
- **Why this is a problem:**  
  `create_function_symbol` (returns `function_symbol`) and `protection_set_size` (returns `std::size_t`) lack `[[nodiscard]]`. Discarding a function_symbol result would be a bug (wasted allocation + lost handle). Project conventions require `[[nodiscard]]` on "pure queries and must-use results."
- **Evidence type:** static-proof (convention violation)
- **Suggested fix:** Add `[[nodiscard]]` to declarations in `thread_aterm_pool.h`.

---

### [SEV-3] No Doxygen `\pre`/`\post` contracts on any public methods

- **Location:** entire file
- **Why this is a problem:**  
  None of the public methods (`create_int`, `create_term`, `create_appl`, `register_variable`, etc.) have Doxygen documentation specifying preconditions (e.g., "arguments must be valid terms"), postconditions (e.g., "`term` is a valid, GC-protected aterm"), or thread-safety guarantees. Per project conventions, "Document every non-trivial public function with Doxygen comments: `\brief`, `\param`, `\returns`, and explicit `\pre` and `\post` clauses."
- **Evidence type:** static-proof (convention violation)
- **Suggested fix:** Add minimal Doxygen to the declarations in `thread_aterm_pool.h`.

---

## Proposed Tests

- [ ] **`test_create_appl_index_pair_layout`:** Compile-time `static_assert` that `sizeof(std::pair<identifier_string, sort_expression>) == 2 * sizeof(unprotected_aterm_core)` and that `offsetof(pair, second) == sizeof(unprotected_aterm_core)`. Place in the existing atermpp test suite. Expected failure signal: compilation error if pair layout assumption is violated on any platform.
- [ ] **`test_create_appl_index_aliasing`:** Unit test that creates an indexed term where the `term` output parameter is pre-populated with a specific int-aterm, then verifies the resulting application's last argument matches the expected index value. Run: `ctest --test-dir build -R atermpp`. Expected failure signal: incorrect term structure if output/input ordering is violated.
- [ ] **`test_strict_aliasing_memcpy_equivalence`:** After applying the `memcpy` fix, verify that `index_traits::insert` receives the same key value as the `reinterpret_cast` path would produce, using a randomized set of aterm arguments. Run under `-O2 -fstrict-aliasing`. Expected failure: index mismatch.

## Open Questions

1. Has the `reinterpret_cast` path ever been tested under `-fstrict-aliasing -Wstrict-aliasing=1` with GCC/Clang at `-O2`/`-O3`? If yes, are there suppression pragmas?
2. Is the `std::pair` layout guarantee relied upon elsewhere in the codebase, or is this the only instance? If widespread, a project-wide `static_assert` or wrapper type would be warranted.
3. The `emplace` function uses `new (&term) unprotected_aterm_core(...)` which constructs a base-class object in memory that holds a derived `aterm` — is there a deliberate policy about implicit object creation / `std::launder` requirements here, or is this additional UB that the `// TODO remove reinterpret_cast` comment alludes to?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         10m 33s
Total session time:     7m 28s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         1.0m in, 15.0k out, 965.9k cached (Est. 1 Premium request)
 claude-haiku-4.5        1.1m in, 27.5k out, 1.0m cached (Est. 0 Premium requests)

```

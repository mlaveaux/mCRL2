---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-08T22:13:04.734728+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h` — the internal representation of term applications and their allocator.
- **Risk level:** Medium. The code is a low-level, performance-critical foundation of the entire aterm library and uses deliberate UB patterns (C-style flexible arrays in C++). It is well-established and architecturally protected, but multiple strict-standard violations remain that could become real bugs under aggressive optimizers or sanitizers.
- **Overall verdict:** Functional and battle-tested code, but contains several technically-undefined operations that rely on assumptions about layout and compiler behavior rather than language guarantees.

## Findings (ordered by severity)

### [SEV-2] Invalid downcast in `arg()` — UB per strict aliasing / object model

- **Location:** `aterm.h:90-93`
- **Why this is a problem:** `static_cast<const aterm_core&>(m_arguments.data()[index])` performs a downcast from `unprotected_aterm_core` to `aterm_core`. The objects stored in `m_arguments` are `unprotected_aterm_core` instances, not `aterm_core`. Per [expr.static.cast]/2 and [basic.lval], accessing an object through a reference to a derived type when the dynamic type is actually the base type is undefined behavior.
- **Evidence type:** static-proof
- **Evidence:**
  - `m_arguments` is declared as `std::array<unprotected_aterm_core, N>`. Its elements have dynamic type `unprotected_aterm_core`.
  - `aterm_core` inherits publicly from `unprotected_aterm_core` but adds non-trivial special members (copy ctor, move ctor, destructor with reference counting).
  - A downcast via `static_cast` from base to derived is only defined when the actual object is the derived type ([expr.static.cast]/2).
  - While layout-compatible today (no added data members), `aterm_core`'s non-trivial destructor means the compiler may legally assume objects accessed through `aterm_core&` have been properly constructed as `aterm_core` — enabling optimizations that break this code.
- **Efficiency impact:** None; this is a correctness issue.
- **C++20-specific note:** C++20 implicit lifetime types (P0593R6) explicitly address similar patterns but only for *implicit-lifetime types*. `aterm_core` has a user-declared destructor and is therefore NOT an implicit-lifetime type. The UB is not cured by C++20.
- **Suggested fix:** Return `const unprotected_aterm_core&` from `arg()`, or use `std::launder`/`reinterpret_cast` with a documented justification, or store `aterm_core` objects in the array (accepting the reference-counting cost).

---

### [SEV-2] Out-of-bounds pointer arithmetic in the dynamic-arity constructor

- **Location:** `aterm.h:74-87` (constructor with `bool` sentinel parameter)
- **Why this is a problem:** `m_arguments.data()[i]` for `i >= N` (when `symbol.arity() > N`) performs pointer arithmetic past the end of the `std::array` object. Per [expr.add]/4, adding an integer to a pointer that moves it beyond one-past-the-end of the array object is UB, regardless of whether extra memory exists at that address.
- **Evidence type:** static-proof
- **Evidence:**
  - `m_arguments` is `std::array<unprotected_aterm_core, N>` (N=1 for `_term_appl`).
  - The allocator allocates `sizeof(T) + (arity - 1) * sizeof(aterm_core)` bytes, placing extra memory after the array.
  - The constructor indexes `m_arguments.data()[i]` for `i` up to `symbol.arity() - 1`, which can exceed `N - 1`.
  - The standard's pointer arithmetic rules apply to the array *object*, not to the underlying storage. `data() + i` for `i > N` is UB.
  - This is the C++ "flexible array member" antipattern. It works on all current compilers, but is not portable per the standard and could break under Link-Time Optimization or compiler builtins that know the array bound.
- **Efficiency impact:** None; performance-motivated design.
- **C++20-specific note:** C++20 did not legalize this pattern. A conformant alternative would be a trailing `unprotected_aterm_core[]` flexible array member (C99 extension, supported as extension by all target compilers), or manual placement-new with `std::launder`.
- **Suggested fix:** Add a `// NOLINT` comment with explicit rationale, or refactor to use a `char[]` buffer with explicit placement-new for each argument (trading code clarity for strict conformance).

---

### [SEV-2] `reinterpret_cast` aliasing violation in conversion operator

- **Location:** `aterm.h:96-99`
- **Why this is a problem:** `reinterpret_cast<_aterm_appl<1>&>(*this)` accesses an `_aterm_appl<N>` object through a reference to `_aterm_appl<1>`, violating strict aliasing ([basic.lval]/11). These are unrelated types (different template instantiations).
- **Evidence type:** static-proof
- **Evidence:**
  - `_aterm_appl<1>` and `_aterm_appl<N>` (N≠1) are distinct types. Neither is a base of the other.
  - Accessing an object through a glvalue of a different type is UB unless it's a character type, a base/derived relationship, or the types are similar.
  - Compilers may assume `_aterm_appl<1>&` and `_aterm_appl<N>&` do not alias, enabling reordering or elimination of stores.
- **Efficiency impact:** None.
- **C++20-specific note:** No relief from C++20. `std::bit_cast` cannot help here (reference types).
- **Suggested fix:** Prefer a free function that does `reinterpret_cast` on the pointer with a `// NOLINTNEXTLINE` and explicit documentation of the standard-library-grade assumption, or use `std::launder` on the result.

---

### [SEV-3] Missing bounds assertion in `arg()`

- **Location:** `aterm.h:90-93`
- **Why this is a problem:** `arg(index)` has no bounds check, even in debug mode. An out-of-range `index` silently reads arbitrary memory. For a core accessor in a model-checking library where correctness is paramount, this is a missed diagnostic opportunity.
- **Evidence type:** plausible-issue
- **Evidence:**
  - The function has no precondition check. The arity is available through `function()` in the base class.
  - A caller that miscalculates index will get silent data corruption rather than a clear assertion failure.
  - Test plan: Add `assert(index < function().arity())` and run existing test suite; no test should fail. Any test that does reveals a latent bug.
- **Efficiency impact:** Negligible (one comparison in debug builds only).
- **C++20-specific note:** Consider `[[pre: index < function().arity()]]` when C++26 contracts become available; for now, use `assert`.
- **Suggested fix:**
  ```cpp
  const aterm_core& arg(std::size_t index) const
  {
    assert(index < function().arity());
    return static_cast<const aterm_core&>(m_arguments.data()[index]);
  }
  ```

---

### [SEV-3] `term_appl_size` has unsigned underflow for arity=0

- **Location:** `aterm.h:116-118`
- **Why this is a problem:** `term_appl_size(0)` computes `sizeof(T) + (0 - 1) * sizeof(aterm_core)`. Since the operands are `std::size_t` (unsigned), `0 - 1` wraps to `SIZE_MAX`, producing a nonsensical allocation size.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Architecturally guarded: the pool routes arity-0 to `term_storage` (which stores `_aterm`, not `_aterm_appl`). So `term_appl_size(0)` should never be called in practice.
  - However, the function is `public` within the allocator class and has no precondition guard. A future refactoring or new caller could trigger this.
  - Test plan: Call `_aterm_appl_allocator<>::allocate_args(symbol_with_arity_0, ...)` directly in a unit test. Expected: assertion failure or safe handling. Actual: attempted allocation of ~`SIZE_MAX` bytes → `std::bad_alloc` or system kill.
- **Efficiency impact:** None in normal operation.
- **Suggested fix:** Add `assert(arity >= 1)` inside `term_appl_size`, or use a saturating computation.

---

### [SEV-3] Missing `[[nodiscard]]` on `arg()` and allocator accessors

- **Location:** `aterm.h:90`, `aterm.h:177-179`
- **Why this is a problem:** `arg()` returns a const reference to critical data; discarding it is always a bug. `capacity()`, `consolidate()`, `has_free_slots()` are query functions whose results should never be ignored.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Per project conventions (`.github/instructions/cpp.instructions.md`): "Mark pure queries and must-use results `[[nodiscard]]`".
  - A caller accidentally writing `term.arg(0);` as a statement would compile silently today.
- **Suggested fix:** Add `[[nodiscard]]` to `arg()`, `capacity()`, `consolidate()`, `has_free_slots()`.

---

### [SEV-3] Missing Doxygen `\pre`/`\post` contracts on public API

- **Location:** `aterm.h:90` (`arg`), `aterm.h:116` (`term_appl_size`), `aterm.h:134` (`allocate_args`)
- **Why this is a problem:** Per project conventions, public APIs must document preconditions. `arg()` requires `index < function().arity()`. `term_appl_size()` requires `arity >= 1`. `allocate_args()` requires the symbol arity to match the iterator range. None are documented with `\pre`.
- **Evidence type:** plausible-issue (convention violation, not runtime failure)
- **Suggested fix:** Add `/// \pre index < function().arity()` etc.

---

## Proposed Tests

- [ ] `test_arg_bounds_assertion`: Call `_aterm_appl<1>::arg(1)` on a term with arity 1. Expected: assertion failure (index out of bounds). Validates the proposed `assert(index < function().arity())`.
  - Run: `ctest --test-dir build -R aterm_test --output-on-failure`

- [ ] `test_term_appl_size_arity_zero`: Call `_aterm_appl_allocator<>::term_appl_size(0)` (make it accessible via friend/test helper). Expected: assertion failure. Validates the underflow guard.
  - Run: `ctest --test-dir build -R aterm_test --output-on-failure`

- [ ] `test_arg_downcast_under_ubsan`: Create a term with arity 2, access both args via `arg(0)` and `arg(1)`, run under UBSan (`-DMCRL2_ENABLE_ADDRESSSANITIZER=ON`). Expected: no UBSan diagnostic (currently may flag the downcast under strict `-fsanitize=undefined -fno-sanitize-recover`). If it fires, confirms the SEV-2 finding.
  - Run: `cmake -B build-asan -DMCRL2_ENABLE_ADDRESSSANITIZER=ON && cmake --build build-asan && ctest --test-dir build-asan -R aterm_test --output-on-failure`

## Open Questions

1. Has this code ever triggered a UBSan diagnostic in CI under the `static_cast` downcast in `arg()`? If not, is `-fsanitize=object-size` or `-fsanitize=vptr` enabled in the sanitizer builds?
2. Is the flexible-array pattern intentionally preserved for performance, or would a `std::vector`-based approach (avoiding the out-of-bounds UB) be acceptable given the maximally-shared nature of terms?
3. The conversion operator `operator _aterm_appl<1>&()` — is it still used? If it's only for legacy compatibility, consider deprecation.

## Copilot Raw Stdout

```text
● Explore: Explore aterm context
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h
  └ 194 lines read

● List directory libraries/atermpp/include/mcrl2/atermpp/detail
  └ 25 files found

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 1-80
  └ 80 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h lines 1-80
  └ 80 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 80-160
  └ 81 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 160-220
  └ 57 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h lines 80-120
  └ 8 lines read

● Read: Explore agent (Explore aterm context)
  └ Running

● Grep "sizeof\(unprotected_aterm_core\)|sizeof\(aterm_core\)" (/workspaces/mCRL2/libraries/atermpp)
  └ 2 lines found

● Grep "term_appl_size" (/workspaces/mCRL2/libraries/atermpp)
  └ 4 lines found

● Grep "arity.*==.*0|arity\(\) == 0|arity == 0" (/workspaces/mCRL2/libraries/atermpp)
  └ 2 lines found

● Grep "allocate_args" (/workspaces/mCRL2/libraries/atermpp)
  └ 2 lines found

● Read: Explore agent (Explore aterm context)
  └ Completed

● Grep "allocate_args|term_appl_size|create_appl" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2
/atermpp/detail/aterm_pool_storage_implementation.h)
  └ 6 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
  186-270
  └ 85 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
  135-186
  └ 52 lines read

● Grep "emplace" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storag
e_implementation.h)
  └ 15 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
  405-440
  └ 36 lines read

● Grep "class.*ATERM_POOL_STORAGE|m_term_set|allocator"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage.h)
  └ 8 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage.h lines 30-70
  └ 41 lines read

● Grep "DynamicNumberOfArguments|storage|arity"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h)
  └ 30 lines found

● Grep "DynamicNumberOfArguments" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ 16 lines found

● Grep "create_appl|create_term"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h)
  └ 4 lines found

● Grep "create_appl_dynamic"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h)
  └ 88 lines found

Now I have sufficient context to produce the review. Let me verify one more detail about the `arg()` downcast and alignment:

● Check aterm_core class hierarchy
  $ cd /workspaces/mCRL2 && grep -n "class aterm_core"
  libraries/atermpp/include/mcrl2/atermpp/aterm_core.h | head -5
  └ 2 lines...

● Grep "arg\(" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h)
  └ 1 line found

Here is the critical review:

---

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h` — the internal representation of term applications and their allocator.
- **Risk level:** Medium. The code is a low-level, performance-critical foundation of the entire aterm library and uses deliberate UB patterns (C-style flexible arrays in C++). It is well-established and architecturally protected, but multiple strict-standard violations remain that could become real bugs under aggressive optimizers or sanitizers.
- **Overall verdict:** Functional and battle-tested code, but contains several technically-undefined operations that rely on assumptions about layout and compiler behavior rather than language guarantees.

## Findings (ordered by severity)

### [SEV-2] Invalid downcast in `arg()` — UB per strict aliasing / object model

- **Location:** `aterm.h:90-93`
- **Why this is a problem:** `static_cast<const aterm_core&>(m_arguments.data()[index])` performs a downcast from `unprotected_aterm_core` to `aterm_core`. The objects stored in `m_arguments` are `unprotected_aterm_core` instances, not `aterm_core`. Per [expr.static.cast]/2 and [basic.lval], accessing an object through a reference to a derived type when the dynamic type is actually the base type is undefined behavior.
- **Evidence type:** static-proof
- **Evidence:**
  - `m_arguments` is declared as `std::array<unprotected_aterm_core, N>`. Its elements have dynamic type `unprotected_aterm_core`.
  - `aterm_core` inherits publicly from `unprotected_aterm_core` but adds non-trivial special members (copy ctor, move ctor, destructor with reference counting).
  - A downcast via `static_cast` from base to derived is only defined when the actual object is the derived type ([expr.static.cast]/2).
  - While layout-compatible today (no added data members), `aterm_core`'s non-trivial destructor means the compiler may legally assume objects accessed through `aterm_core&` have been properly constructed as `aterm_core` — enabling optimizations that break this code.
- **Efficiency impact:** None; this is a correctness issue.
- **C++20-specific note:** C++20 implicit lifetime types (P0593R6) explicitly address similar patterns but only for *implicit-lifetime types*. `aterm_core` has a user-declared destructor and is therefore NOT an implicit-lifetime type. The UB is not cured by C++20.
- **Suggested fix:** Return `const unprotected_aterm_core&` from `arg()`, or use `std::launder`/`reinterpret_cast` with a documented justification, or store `aterm_core` objects in the array (accepting the reference-counting cost).

---

### [SEV-2] Out-of-bounds pointer arithmetic in the dynamic-arity constructor

- **Location:** `aterm.h:74-87` (constructor with `bool` sentinel parameter)
- **Why this is a problem:** `m_arguments.data()[i]` for `i >= N` (when `symbol.arity() > N`) performs pointer arithmetic past the end of the `std::array` object. Per [expr.add]/4, adding an integer to a pointer that moves it beyond one-past-the-end of the array object is UB, regardless of whether extra memory exists at that address.
- **Evidence type:** static-proof
- **Evidence:**
  - `m_arguments` is `std::array<unprotected_aterm_core, N>` (N=1 for `_term_appl`).
  - The allocator allocates `sizeof(T) + (arity - 1) * sizeof(aterm_core)` bytes, placing extra memory after the array.
  - The constructor indexes `m_arguments.data()[i]` for `i` up to `symbol.arity() - 1`, which can exceed `N - 1`.
  - The standard's pointer arithmetic rules apply to the array *object*, not to the underlying storage. `data() + i` for `i > N` is UB.
  - This is the C++ "flexible array member" antipattern. It works on all current compilers, but is not portable per the standard and could break under Link-Time Optimization or compiler builtins that know the array bound.
- **Efficiency impact:** None; performance-motivated design.
- **C++20-specific note:** C++20 did not legalize this pattern. A conformant alternative would be a trailing `unprotected_aterm_core[]` flexible array member (C99 extension, supported as extension by all target compilers), or manual placement-new with `std::launder`.
- **Suggested fix:** Add a `// NOLINT` comment with explicit rationale, or refactor to use a `char[]` buffer with explicit placement-new for each argument (trading code clarity for strict conformance).

---

### [SEV-2] `reinterpret_cast` aliasing violation in conversion operator

- **Location:** `aterm.h:96-99`
- **Why this is a problem:** `reinterpret_cast<_aterm_appl<1>&>(*this)` accesses an `_aterm_appl<N>` object through a reference to `_aterm_appl<1>`, violating strict aliasing ([basic.lval]/11). These are unrelated types (different template instantiations).
- **Evidence type:** static-proof
- **Evidence:**
  - `_aterm_appl<1>` and `_aterm_appl<N>` (N≠1) are distinct types. Neither is a base of the other.
  - Accessing an object through a glvalue of a different type is UB unless it's a character type, a base/derived relationship, or the types are similar.
  - Compilers may assume `_aterm_appl<1>&` and `_aterm_appl<N>&` do not alias, enabling reordering or elimination of stores.
- **Efficiency impact:** None.
- **C++20-specific note:** No relief from C++20. `std::bit_cast` cannot help here (reference types).
- **Suggested fix:** Prefer a free function that does `reinterpret_cast` on the pointer with a `// NOLINTNEXTLINE` and explicit documentation of the standard-library-grade assumption, or use `std::launder` on the result.

---

### [SEV-3] Missing bounds assertion in `arg()`

- **Location:** `aterm.h:90-93`
- **Why this is a problem:** `arg(index)` has no bounds check, even in debug mode. An out-of-range `index` silently reads arbitrary memory. For a core accessor in a model-checking library where correctness is paramount, this is a missed diagnostic opportunity.
- **Evidence type:** plausible-issue
- **Evidence:**
  - The function has no precondition check. The arity is available through `function()` in the base class.
  - A caller that miscalculates index will get silent data corruption rather than a clear assertion failure.
  - Test plan: Add `assert(index < function().arity())` and run existing test suite; no test should fail. Any test that does reveals a latent bug.
- **Efficiency impact:** Negligible (one comparison in debug builds only).
- **C++20-specific note:** Consider `[[pre: index < function().arity()]]` when C++26 contracts become available; for now, use `assert`.
- **Suggested fix:**
  ```cpp
  const aterm_core& arg(std::size_t index) const
  {
    assert(index < function().arity());
    return static_cast<const aterm_core&>(m_arguments.data()[index]);
  }
  ```

---

### [SEV-3] `term_appl_size` has unsigned underflow for arity=0

- **Location:** `aterm.h:116-118`
- **Why this is a problem:** `term_appl_size(0)` computes `sizeof(T) + (0 - 1) * sizeof(aterm_core)`. Since the operands are `std::size_t` (unsigned), `0 - 1` wraps to `SIZE_MAX`, producing a nonsensical allocation size.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Architecturally guarded: the pool routes arity-0 to `term_storage` (which stores `_aterm`, not `_aterm_appl`). So `term_appl_size(0)` should never be called in practice.
  - However, the function is `public` within the allocator class and has no precondition guard. A future refactoring or new caller could trigger this.
  - Test plan: Call `_aterm_appl_allocator<>::allocate_args(symbol_with_arity_0, ...)` directly in a unit test. Expected: assertion failure or safe handling. Actual: attempted allocation of ~`SIZE_MAX` bytes → `std::bad_alloc` or system kill.
- **Efficiency impact:** None in normal operation.
- **Suggested fix:** Add `assert(arity >= 1)` inside `term_appl_size`, or use a saturating computation.

---

### [SEV-3] Missing `[[nodiscard]]` on `arg()` and allocator accessors

- **Location:** `aterm.h:90`, `aterm.h:177-179`
- **Why this is a problem:** `arg()` returns a const reference to critical data; discarding it is always a bug. `capacity()`, `consolidate()`, `has_free_slots()` are query functions whose results should never be ignored.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Per project conventions (`.github/instructions/cpp.instructions.md`): "Mark pure queries and must-use results `[[nodiscard]]`".
  - A caller accidentally writing `term.arg(0);` as a statement would compile silently today.
- **Suggested fix:** Add `[[nodiscard]]` to `arg()`, `capacity()`, `consolidate()`, `has_free_slots()`.

---

### [SEV-3] Missing Doxygen `\pre`/`\post` contracts on public API

- **Location:** `aterm.h:90` (`arg`), `aterm.h:116` (`term_appl_size`), `aterm.h:134` (`allocate_args`)
- **Why this is a problem:** Per project conventions, public APIs must document preconditions. `arg()` requires `index < function().arity()`. `term_appl_size()` requires `arity >= 1`. `allocate_args()` requires the symbol arity to match the iterator range. None are documented with `\pre`.
- **Evidence type:** plausible-issue (convention violation, not runtime failure)
- **Suggested fix:** Add `/// \pre index < function().arity()` etc.

---

## Proposed Tests

- [ ] `test_arg_bounds_assertion`: Call `_aterm_appl<1>::arg(1)` on a term with arity 1. Expected: assertion failure (index out of bounds). Validates the proposed `assert(index < function().arity())`.
  - Run: `ctest --test-dir build -R aterm_test --output-on-failure`

- [ ] `test_term_appl_size_arity_zero`: Call `_aterm_appl_allocator<>::term_appl_size(0)` (make it accessible via friend/test helper). Expected: assertion failure. Validates the underflow guard.
  - Run: `ctest --test-dir build -R aterm_test --output-on-failure`

- [ ] `test_arg_downcast_under_ubsan`: Create a term with arity 2, access both args via `arg(0)` and `arg(1)`, run under UBSan (`-DMCRL2_ENABLE_ADDRESSSANITIZER=ON`). Expected: no UBSan diagnostic (currently may flag the downcast under strict `-fsanitize=undefined -fno-sanitize-recover`). If it fires, confirms the SEV-2 finding.
  - Run: `cmake -B build-asan -DMCRL2_ENABLE_ADDRESSSANITIZER=ON && cmake --build build-asan && ctest --test-dir build-asan -R aterm_test --output-on-failure`

## Open Questions

1. Has this code ever triggered a UBSan diagnostic in CI under the `static_cast` downcast in `arg()`? If not, is `-fsanitize=object-size` or `-fsanitize=vptr` enabled in the sanitizer builds?
2. Is the flexible-array pattern intentionally preserved for performance, or would a `std::vector`-based approach (avoiding the out-of-bounds UB) be acceptable given the maximally-shared nature of terms?
3. The conversion operator `operator _aterm_appl<1>&()` — is it still used? If it's only for legacy compatibility, consider deprecation.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         5m 26s
Total session time:     5m 6s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         929.2k in, 11.3k out, 853.9k cached (Est. 1 Premium request)
 claude-haiku-4.5        362.2k in, 7.5k out, 321.1k cached (Est. 0 Premium requests)

```

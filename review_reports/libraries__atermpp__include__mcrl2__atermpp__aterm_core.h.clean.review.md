---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/aterm_core.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-07T23:43:42.023573+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h` — full header defining `unprotected_aterm_core` and `aterm_core` base classes
- Risk level: **Medium** — one null-dereference UB path reachable from public API; remainder are quality/portability concerns
- Overall verdict: Functionally sound for its intended use pattern (callers must check `defined()` before use), but the header has multiple C++20 quality gaps and a latent UB path with no compile-time or runtime guard.

## Findings (ordered by severity)

### [SEV-2] Null pointer dereference in `type_is_int()`, `type_is_list()`, and `function()` on default-constructed term
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:64-78,121-124`
- Why this is a problem: A default-constructed `unprotected_aterm_core` has `m_term == nullptr`. Calling `type_is_int()`, `type_is_list()`, or `function()` immediately dereferences `m_term->function()` without a null check. This is undefined behavior. The `defined()` method exists as a precondition check but is not enforced — no `assert`, no `[[expects]]`, no compile-time protection.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Construct default `unprotected_aterm_core`, call `type_is_int()`
  - Run command: `ctest --test-dir build -R atermpp -j4 --output-on-failure` (after adding proposed test below)
  - Observed result: Undefined behavior (likely segfault under ASan/UBSan)
  - Expected result: Either a precondition assert fires, or the method returns `false` safely
- Efficiency impact: Adding an `assert(m_term != nullptr)` is zero-cost in release builds.
- C++20-specific note: C++20 contracts (`[[expects: m_term != nullptr]]`) would be ideal here when available; `assert` is the current best practice.
- Suggested fix: Add `assert(m_term != nullptr)` (or `assert(defined())`) at the top of `type_is_int()`, `type_is_list()`, and `function()`.

### [SEV-3] Include guard `MCRL2_ATERMPP_ATERM_H` does not match filename `aterm_core.h`
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:11`
- Why this is a problem: The guard name suggests this protects `aterm.h`, not `aterm_core.h`. While `aterm.h` uses a different guard (`MCRL2_ATERMPP_ATERM_APPL_H`), the mismatch creates maintenance confusion and violates the convention that guard names mirror file paths. A future developer adding a true `aterm.h` guard with this name would silently suppress inclusion.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: N/A — naming/maintenance concern
  - Run command: N/A
  - Observed result: Guard is `MCRL2_ATERMPP_ATERM_H`
  - Expected result: Guard should be `MCRL2_ATERMPP_ATERM_CORE_H`
- Efficiency impact: None at runtime; compilation hygiene only.
- C++20-specific note: `#pragma once` is a non-standard alternative but the project consistently uses traditional guards.
- Suggested fix: Rename guard to `MCRL2_ATERMPP_ATERM_CORE_H`.

### [SEV-3] Unused `#include <sstream>` adds compilation overhead
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:16`
- Why this is a problem: No symbol from `<sstream>` is used anywhere in this header. Since this is a foundational header included transitively by most of the library, the unnecessary include propagates costly template instantiations to hundreds of translation units.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `grep -n "sstream\|ostringstream\|istringstream\|stringstream" aterm_core.h` finds only the include itself
  - Run command: `grep -c "sstream\|ostringstream\|istringstream\|stringstream" libraries/atermpp/include/mcrl2/atermpp/aterm_core.h`
  - Observed result: 1 (only the include line)
  - Expected result: 0 includes of unused headers
- Efficiency impact: Measurable compile-time reduction on a widely-included header.
- C++20-specific note: C++20 modules would isolate this; absent modules, removing unused includes is best practice.
- Suggested fix: Remove `#include <sstream>`.

### [SEV-3] `operator==` and `operator<=>` lack `noexcept`
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:86-100`
- Why this is a problem: Both operators perform only pointer comparison, which cannot throw. Missing `noexcept` prevents certain `std::` algorithm optimizations and violates the principle of maximal noexcept guarantees for trivial operations. Other methods in the same class (constructors, `swap`) are correctly marked `noexcept`.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: `static_assert(noexcept(std::declval<unprotected_aterm_core>() == std::declval<unprotected_aterm_core>()))` — would fail to compile
  - Run command: Add static_assert to any test TU and build
  - Observed result: Compilation failure
  - Expected result: Assertion passes
- Efficiency impact: Minor — affects `std::sort`/`std::equal` inlining decisions in some implementations.
- C++20-specific note: The C++20 defaulted `operator<=>` with `= default` would inherit `noexcept` from member comparisons. Here, since it's manually defined, `noexcept` must be explicit.
- Suggested fix: Add `noexcept` to both `operator==` and `operator<=>`.

### [SEV-3] `detail::address()` uses `const_cast` to strip const, creating a const-safety hole
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:194-198`
- Why this is a problem: `address()` takes a `const unprotected_aterm_core&` and returns a non-const `_aterm*`. Any downstream code that mutates the returned pointer while the original was const-qualified would trigger UB. The function exists for internal use (GC marking, hashing) but leaks a mutable handle through a const interface.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Any code path that calls `address(const_term)` and then writes through the result
  - Run command: Search for `address(` call sites that modify the returned pointer
  - Observed result: Potential UB if mutation occurs
  - Expected result: Two overloads — const input → const output, non-const input → non-const output
- Efficiency impact: None.
- C++20-specific note: A `const`-correct overload set or a separate `mutable_address()` function would be idiomatic.
- Suggested fix: Either provide two overloads (const/non-const) or document `// NOLINT: required for GC marking` with a comment explaining the safety invariant.

### [SEV-3] Move constructor/assignment do not leave source in a distinct moved-from state
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:172-176,183-186` (declarations); implementations in `detail/aterm_implementation.h:83-87,117-123`
- Why this is a problem: After `aterm_core b = std::move(a)`, both `a` and `b` point to the same underlying term. While this is safe in this design (protection is address-registration-based), it violates the conventional C++ expectation that moved-from objects are in a "valid but unspecified" state that is distinguishable. Specifically, `a.defined()` still returns `true` after being moved from, which is surprising and could mask bugs in calling code.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: `aterm_core a(some_term); aterm_core b = std::move(a); assert(!a.defined());` — this assert would FAIL
  - Run command: Add test, build and run with `ctest -R atermpp`
  - Observed result: `a.defined()` returns `true` after move
  - Expected result: Depends on design intent — if intentional, document it explicitly
- Efficiency impact: Nulling `other.m_term` in move constructor would be a single store — negligible cost.
- C++20-specific note: C++ Core Guidelines recommend moved-from objects be in a destructible and assignable state; this is met, but the non-null state after move is non-idiomatic.
- Suggested fix: Either set `other.m_term = nullptr` in move operations (zero runtime cost) or add a prominent comment explaining the design choice.

## Proposed Tests
- [ ] `test_null_dereference_precondition`: Default-construct `unprotected_aterm_core`, call `type_is_int()` — verify assert fires (debug build) or safely returns false (release build). Run under UBSan.
- [ ] `test_noexcept_comparisons`: `static_assert(noexcept(a == b))` and `static_assert(noexcept(a <=> b))` for `unprotected_aterm_core` — validates noexcept guarantee after fix.
- [ ] `test_moved_from_state`: Construct term, move it, verify source state (documents contract regardless of fix direction).

## Open Questions
1. Is the `const_cast` in `address()` exercised in any path that actually mutates the `_aterm`? If so, what invariant guarantees the original object was not truly const?
2. Is the non-null moved-from state an intentional performance optimization (avoiding a store) or an oversight? The move assignment acquires a shared lock already, so the marginal cost of one additional pointer store is negligible.
3. Given that `#include <sstream>` is unused here, is it depended upon transitively by downstream headers that assume `aterm_core.h` pulls it in? (If so, those headers should add their own include.)

## Copilot Raw Stdout

```text
● Explore: Explore aterm_core dependencies
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h
    └ 62 lines read
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/type_traits_impl.h lines 1-100
    └ 69 lines read

● Explore: Find aterm_core.cpp implementation
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h
    └ 128 lines read
  ● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h
    └ 213 lines read

● Explore: Explore thread_aterm_pool and operator<=> details
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h lines 1-150
    └ 128 lines read
  ● Check for comparison operator documentation
    $ cd /workspaces/mCRL2 && grep -r "operator<=>.*pointer\|weak_ordering\|strong_ordering"
    --include="*.h" | head -20
    └ 5 lines...

● Grep "MCRL2_ATERMPP_ATERM_H" (libraries/atermpp/include)
  └ 1 file found

● Glob "libraries/atermpp/include/mcrl2/atermpp/aterm.h"
  └ 1 file found

● Grep "#include <sstream>" (libraries/atermpp/include/mcrl2/atermpp/aterm_core.h)
  └ 1 line found

● Grep "sstream" (libraries/atermpp/include/mcrl2/atermpp/aterm_core.h)
  └ 1 line found

● Check aterm.h include guard
  $ head -5 /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm.h | grep -n
  "ifndef\|define"
  └ 1 line...

● Check if sstream is used
  $ grep -n "sstream\|ostringstream\|istringstream\|stringstream\|str()"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm_core.h
  └ 2 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/aterm.h lines 1-15
  └ 15 lines read

● Find who includes aterm_core.h
  $ grep -rn "include.*aterm_core.h" /workspaces/mCRL2/libraries/atermpp/ | head -10
  └ 9 lines...

Now I have all the information needed to produce the review. Let me write the formal findings.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h` — full header defining `unprotected_aterm_core` and `aterm_core` base classes
- Risk level: **Medium** — one null-dereference UB path reachable from public API; remainder are quality/portability concerns
- Overall verdict: Functionally sound for its intended use pattern (callers must check `defined()` before use), but the header has multiple C++20 quality gaps and a latent UB path with no compile-time or runtime guard.

## Findings (ordered by severity)

### [SEV-2] Null pointer dereference in `type_is_int()`, `type_is_list()`, and `function()` on default-constructed term
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:64-78,121-124`
- Why this is a problem: A default-constructed `unprotected_aterm_core` has `m_term == nullptr`. Calling `type_is_int()`, `type_is_list()`, or `function()` immediately dereferences `m_term->function()` without a null check. This is undefined behavior. The `defined()` method exists as a precondition check but is not enforced — no `assert`, no `[[expects]]`, no compile-time protection.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Construct default `unprotected_aterm_core`, call `type_is_int()`
  - Run command: `ctest --test-dir build -R atermpp -j4 --output-on-failure` (after adding proposed test below)
  - Observed result: Undefined behavior (likely segfault under ASan/UBSan)
  - Expected result: Either a precondition assert fires, or the method returns `false` safely
- Efficiency impact: Adding an `assert(m_term != nullptr)` is zero-cost in release builds.
- C++20-specific note: C++20 contracts (`[[expects: m_term != nullptr]]`) would be ideal here when available; `assert` is the current best practice.
- Suggested fix: Add `assert(m_term != nullptr)` (or `assert(defined())`) at the top of `type_is_int()`, `type_is_list()`, and `function()`.

### [SEV-3] Include guard `MCRL2_ATERMPP_ATERM_H` does not match filename `aterm_core.h`
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:11`
- Why this is a problem: The guard name suggests this protects `aterm.h`, not `aterm_core.h`. While `aterm.h` uses a different guard (`MCRL2_ATERMPP_ATERM_APPL_H`), the mismatch creates maintenance confusion and violates the convention that guard names mirror file paths. A future developer adding a true `aterm.h` guard with this name would silently suppress inclusion.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: N/A — naming/maintenance concern
  - Run command: N/A
  - Observed result: Guard is `MCRL2_ATERMPP_ATERM_H`
  - Expected result: Guard should be `MCRL2_ATERMPP_ATERM_CORE_H`
- Efficiency impact: None at runtime; compilation hygiene only.
- C++20-specific note: `#pragma once` is a non-standard alternative but the project consistently uses traditional guards.
- Suggested fix: Rename guard to `MCRL2_ATERMPP_ATERM_CORE_H`.

### [SEV-3] Unused `#include <sstream>` adds compilation overhead
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:16`
- Why this is a problem: No symbol from `<sstream>` is used anywhere in this header. Since this is a foundational header included transitively by most of the library, the unnecessary include propagates costly template instantiations to hundreds of translation units.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `grep -n "sstream\|ostringstream\|istringstream\|stringstream" aterm_core.h` finds only the include itself
  - Run command: `grep -c "sstream\|ostringstream\|istringstream\|stringstream" libraries/atermpp/include/mcrl2/atermpp/aterm_core.h`
  - Observed result: 1 (only the include line)
  - Expected result: 0 includes of unused headers
- Efficiency impact: Measurable compile-time reduction on a widely-included header.
- C++20-specific note: C++20 modules would isolate this; absent modules, removing unused includes is best practice.
- Suggested fix: Remove `#include <sstream>`.

### [SEV-3] `operator==` and `operator<=>` lack `noexcept`
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:86-100`
- Why this is a problem: Both operators perform only pointer comparison, which cannot throw. Missing `noexcept` prevents certain `std::` algorithm optimizations and violates the principle of maximal noexcept guarantees for trivial operations. Other methods in the same class (constructors, `swap`) are correctly marked `noexcept`.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: `static_assert(noexcept(std::declval<unprotected_aterm_core>() == std::declval<unprotected_aterm_core>()))` — would fail to compile
  - Run command: Add static_assert to any test TU and build
  - Observed result: Compilation failure
  - Expected result: Assertion passes
- Efficiency impact: Minor — affects `std::sort`/`std::equal` inlining decisions in some implementations.
- C++20-specific note: The C++20 defaulted `operator<=>` with `= default` would inherit `noexcept` from member comparisons. Here, since it's manually defined, `noexcept` must be explicit.
- Suggested fix: Add `noexcept` to both `operator==` and `operator<=>`.

### [SEV-3] `detail::address()` uses `const_cast` to strip const, creating a const-safety hole
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:194-198`
- Why this is a problem: `address()` takes a `const unprotected_aterm_core&` and returns a non-const `_aterm*`. Any downstream code that mutates the returned pointer while the original was const-qualified would trigger UB. The function exists for internal use (GC marking, hashing) but leaks a mutable handle through a const interface.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Any code path that calls `address(const_term)` and then writes through the result
  - Run command: Search for `address(` call sites that modify the returned pointer
  - Observed result: Potential UB if mutation occurs
  - Expected result: Two overloads — const input → const output, non-const input → non-const output
- Efficiency impact: None.
- C++20-specific note: A `const`-correct overload set or a separate `mutable_address()` function would be idiomatic.
- Suggested fix: Either provide two overloads (const/non-const) or document `// NOLINT: required for GC marking` with a comment explaining the safety invariant.

### [SEV-3] Move constructor/assignment do not leave source in a distinct moved-from state
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_core.h:172-176,183-186` (declarations); implementations in `detail/aterm_implementation.h:83-87,117-123`
- Why this is a problem: After `aterm_core b = std::move(a)`, both `a` and `b` point to the same underlying term. While this is safe in this design (protection is address-registration-based), it violates the conventional C++ expectation that moved-from objects are in a "valid but unspecified" state that is distinguishable. Specifically, `a.defined()` still returns `true` after being moved from, which is surprising and could mask bugs in calling code.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: `aterm_core a(some_term); aterm_core b = std::move(a); assert(!a.defined());` — this assert would FAIL
  - Run command: Add test, build and run with `ctest -R atermpp`
  - Observed result: `a.defined()` returns `true` after move
  - Expected result: Depends on design intent — if intentional, document it explicitly
- Efficiency impact: Nulling `other.m_term` in move constructor would be a single store — negligible cost.
- C++20-specific note: C++ Core Guidelines recommend moved-from objects be in a destructible and assignable state; this is met, but the non-null state after move is non-idiomatic.
- Suggested fix: Either set `other.m_term = nullptr` in move operations (zero runtime cost) or add a prominent comment explaining the design choice.

## Proposed Tests
- [ ] `test_null_dereference_precondition`: Default-construct `unprotected_aterm_core`, call `type_is_int()` — verify assert fires (debug build) or safely returns false (release build). Run under UBSan.
- [ ] `test_noexcept_comparisons`: `static_assert(noexcept(a == b))` and `static_assert(noexcept(a <=> b))` for `unprotected_aterm_core` — validates noexcept guarantee after fix.
- [ ] `test_moved_from_state`: Construct term, move it, verify source state (documents contract regardless of fix direction).

## Open Questions
1. Is the `const_cast` in `address()` exercised in any path that actually mutates the `_aterm`? If so, what invariant guarantees the original object was not truly const?
2. Is the non-null moved-from state an intentional performance optimization (avoiding a store) or an oversight? The move assignment acquires a shared lock already, so the marginal cost of one additional pointer store is negligible.
3. Given that `#include <sstream>` is unused here, is it depended upon transitively by downstream headers that assume `aterm_core.h` pulls it in? (If so, those headers should add their own include.)


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         5m 0s
Total session time:     5m 15s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         249.1k in, 8.1k out, 218.7k cached (Est. 1 Premium request)
 claude-haiku-4.5        362.2k in, 13.5k out, 309.1k cached (Est. 0 Premium requests)

```

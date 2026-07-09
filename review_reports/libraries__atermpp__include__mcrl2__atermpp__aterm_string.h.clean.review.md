---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/aterm_string.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-08T00:08:48.404492+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h` — a thin wrapper around `aterm` representing zero-arity terms as strings
- Risk level: Low-Medium (no confirmed crash bugs; several design-level concerns with concrete failure paths in edge cases)
- Overall verdict: The class is functional and consistent with codebase conventions, but relies on runtime-only invariant checks that vanish in release builds, and exposes an implicit conversion operator with dangling-reference risk.

## Findings (ordered by severity)

### [SEV-2] Runtime-only invariant enforcement in `aterm_string(const aterm& t)` — silent corruption in release builds

- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h:28-31`
- Why this is a problem: The constructor validates `size() == 0` only via `assert()`. In release builds (`-DNDEBUG`), a caller can construct an `aterm_string` from an `aterm` with arity > 0, creating an object that violates the class invariant. Downstream code (e.g., `operator const std::string&()`) will still call `function().name()` on what it assumes is a zero-arity term, but the semantic contract is broken. Any comparison, hashing, or string extraction will silently operate on a corrupted abstraction.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Construct `aterm_string` from multi-arity aterm in release mode
  - Run command: `cmake --build build --target aterm_test && cd build && ctest -R aterm_test --output-on-failure` (after adding test below)
  - Observed result: In debug builds, `assert` fires. In release (`-DCMAKE_BUILD_TYPE=Release`), the invalid object is silently created.
  - Expected result: Either a compile-time narrowing check or a runtime exception/error that survives release builds.
- Efficiency impact: None (this is a correctness issue, not performance).
- C++20-specific note: A `requires` clause or `static_assert` cannot easily guard this at compile time since `aterm::size()` is runtime. A `[[gnu::always_inline]]` conditional throwing constructor or a factory function with `std::expected` would be more appropriate in modern C++20.
- Suggested fix: Replace `assert(size() == 0)` with a conditional that throws or calls `std::terminate` regardless of build mode:
  ```cpp
  explicit aterm_string(const aterm& t)
    : aterm(t)
  {
    if (size() != 0) { throw std::invalid_argument("aterm_string requires arity 0"); }
  }
  ```
  Alternatively, if exceptions are undesirable in this codebase, use `mcrl2_assert` or the project's own abort mechanism that is not stripped in release.

### [SEV-2] Implicit `operator const std::string&()` enables dangling references

- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h:67-70`
- Why this is a problem: The implicit conversion operator returns a `const std::string&` to an internal string owned by the function symbol in the global term pool. If an `aterm_string` temporary is implicitly converted, the resulting reference binds to data whose lifetime is managed by the term pool, not the local scope. While the term pool typically persists, the implicit nature means callers can accidentally write `const std::string& s = some_function_returning_aterm_string();` where lifetime extension does NOT apply (because the reference binds through a conversion operator, not directly to the temporary). Additionally, implicit conversion operators create overload ambiguity with `std::string` parameters.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Verify whether `const std::string& s = aterm_string("hello");` produces a dangling reference
  - Run command: Compile with `-fsanitize=address` and exercise the pattern
  - Observed result: Likely safe in practice due to term pool persistence, but the pattern is undefined behavior per the standard (temporary `aterm_string` destroyed, reference to internal data of a pool-managed object survives only coincidentally).
  - Expected result: Either mark the conversion operator `explicit` or document that the returned reference has pool-managed lifetime.
- Efficiency impact: None directly, but implicit conversions can cause unexpected temporary `std::string` copies at call sites expecting value semantics.
- C++20-specific note: In C++20, the idiomatic approach is `explicit operator const std::string&()` with a named `.str()` or `.name()` accessor for the common case.
- Suggested fix: Mark `explicit` and add a named accessor:
  ```cpp
  explicit operator const std::string&() const { return function().name(); }
  const std::string& str() const { return function().name(); }
  ```

### [SEV-3] Non-explicit converting constructors allow silent implicit conversions

- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h:35-48`
- Why this is a problem: `aterm_string(const function_symbol& f)` and `aterm_string(const std::string& s)` are both non-explicit, enabling implicit conversions. A function accepting `const aterm_string&` can silently accept a raw `std::string` or `function_symbol`, triggering allocation of a new function symbol in the global pool as a side effect of what looks like a simple function call. This is a classic "action at a distance" pattern.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Measure allocation count when passing `std::string` to a function expecting `aterm_string`
  - Run command: N/A (design review finding)
  - Observed result: The codebase deliberately relies on these implicit conversions (e.g., `pbeschain.h:112` comparisons). Changing this would be a breaking API modification.
  - Expected result: Ideally `explicit` constructors with deliberate call-site conversions.
- Efficiency impact: Each implicit conversion from `std::string` allocates or looks up a `function_symbol` in the global pool — this is hidden cost at every implicit conversion site.
- C++20-specific note: C++20 `explicit(bool)` conditional explicitness could be used to make these explicit in new code while preserving backward compatibility via a macro/policy.
- Suggested fix: This is a known design trade-off in the codebase. No immediate fix recommended without broader API discussion. Document the allocation side-effect in the constructor's doxygen comment.

### [SEV-3] `empty_string()` static destruction order concern

- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h:75-79`
- Why this is a problem: The function returns a `const` reference to a `static` local `aterm_string`. If any other static-lifetime object (in another TU) holds a reference to this and the term pool is destroyed before this static, accessing the reference after pool teardown is undefined behavior. This is the classic "static destruction order fiasco."
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Create a static object in another TU whose destructor accesses `empty_string()`
  - Run command: Requires multi-TU test with deliberate static ordering
  - Observed result: Likely safe in current usage (the term pool is typically the last thing torn down), but architecturally fragile.
  - Expected result: Use a leaked singleton (`static aterm_string& t = *new aterm_string("")`) if absolute safety is required, or document the constraint.
- Efficiency impact: None.
- C++20-specific note: `constinit` cannot be used here since `aterm_string("")` is not constexpr-constructible (it allocates in the term pool).
- Suggested fix: Low priority. Document that `empty_string()` must not be called during static destruction, or use the leak-singleton pattern if this guarantee cannot be upheld.

### [SEV-3] `std::swap` specialization in namespace `std` — deprecated practice

- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h:92-96`
- Why this is a problem: Explicit specialization of `std::swap` in namespace `std` is valid but has been superseded by hidden-friend ADL swap in modern C++. The current approach also prevents partial specialization (not applicable here, but limits future templated derivatives). More importantly, the specialization delegates to `t1.swap(t2)` which is inherited from `unprotected_aterm_core` and only swaps the internal pointer — this is correct but non-obvious.
- Evidence type: N/A (style/modernization observation)
- Efficiency impact: None.
- C++20-specific note: In C++20, the preferred pattern is a `friend void swap(aterm_string& a, aterm_string& b) noexcept` defined inline in the class body, found via ADL.
- Suggested fix: Low priority. Consider adding a `friend swap` in the class for ADL discovery while keeping the `std::swap` specialization for backward compatibility.

## Proposed Tests

- [ ] `test_aterm_string_from_nonzero_arity`: Construct `aterm_string` from an `aterm` with arity > 0 in a release build. Verify the class either throws or aborts (currently it silently succeeds). Input: `aterm(function_symbol("f", 2), aterm_string("a"), aterm_string("b"))`. Expected failure signal: exception or assertion failure.
- [ ] `test_aterm_string_implicit_conversion_allocation`: Measure that passing `std::string("x")` where `aterm_string` is expected causes exactly one function symbol pool lookup/insertion. Input: function accepting `const aterm_string&` called with `std::string`. Expected signal: pool size increases by at most 1.
- [ ] `test_empty_string_identity`: Verify `&empty_string() == &empty_string()` (same object returned on repeated calls). Input: two calls to `empty_string()`. Expected signal: pointer equality.

## Open Questions

- Is the decision to keep converting constructors non-explicit a deliberate API choice, or historical inertia? The codebase relies on it, but it hides allocation costs.
- Has the team considered marking `operator const std::string&()` as `explicit` now that C++20 is the baseline? The call sites using `std::string(x)` already perform explicit conversion.
- Are there any known cases of `aterm_string` being constructed from non-zero-arity aterms in production code paths (i.e., has the `assert` ever fired in CI)?

## Copilot Raw Stdout

```text
● Explore: Explore aterm base class
  ● Find empty_tree static variable
    $ grep -B5 -A5 "static const aterm& empty_tree"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm_balanced_tree.h
    └ 12 lines...
  ● Full test file for reference
    $ cat /workspaces/mCRL2/libraries/atermpp/test/aterm_test.cpp
    └ 81 lines...

● Explore: Check implicit conversion risks
  ● Check aterm_core class
    $ grep -rn "unprotected_aterm_core"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm_core.h | head -10
    └ 10 lines...
  ● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 30-150
    └ 120 lines read

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h` — a thin wrapper around `aterm` representing zero-arity terms as strings
- Risk level: Low-Medium (no confirmed crash bugs; several design-level concerns with concrete failure paths in edge cases)
- Overall verdict: The class is functional and consistent with codebase conventions, but relies on runtime-only invariant checks that vanish in release builds, and exposes an implicit conversion operator with dangling-reference risk.

## Findings (ordered by severity)

### [SEV-2] Runtime-only invariant enforcement in `aterm_string(const aterm& t)` — silent corruption in release builds

- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h:28-31`
- Why this is a problem: The constructor validates `size() == 0` only via `assert()`. In release builds (`-DNDEBUG`), a caller can construct an `aterm_string` from an `aterm` with arity > 0, creating an object that violates the class invariant. Downstream code (e.g., `operator const std::string&()`) will still call `function().name()` on what it assumes is a zero-arity term, but the semantic contract is broken. Any comparison, hashing, or string extraction will silently operate on a corrupted abstraction.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Construct `aterm_string` from multi-arity aterm in release mode
  - Run command: `cmake --build build --target aterm_test && cd build && ctest -R aterm_test --output-on-failure` (after adding test below)
  - Observed result: In debug builds, `assert` fires. In release (`-DCMAKE_BUILD_TYPE=Release`), the invalid object is silently created.
  - Expected result: Either a compile-time narrowing check or a runtime exception/error that survives release builds.
- Efficiency impact: None (this is a correctness issue, not performance).
- C++20-specific note: A `requires` clause or `static_assert` cannot easily guard this at compile time since `aterm::size()` is runtime. A `[[gnu::always_inline]]` conditional throwing constructor or a factory function with `std::expected` would be more appropriate in modern C++20.
- Suggested fix: Replace `assert(size() == 0)` with a conditional that throws or calls `std::terminate` regardless of build mode:
  ```cpp
  explicit aterm_string(const aterm& t)
    : aterm(t)
  {
    if (size() != 0) { throw std::invalid_argument("aterm_string requires arity 0"); }
  }
  ```
  Alternatively, if exceptions are undesirable in this codebase, use `mcrl2_assert` or the project's own abort mechanism that is not stripped in release.

### [SEV-2] Implicit `operator const std::string&()` enables dangling references

- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h:67-70`
- Why this is a problem: The implicit conversion operator returns a `const std::string&` to an internal string owned by the function symbol in the global term pool. If an `aterm_string` temporary is implicitly converted, the resulting reference binds to data whose lifetime is managed by the term pool, not the local scope. While the term pool typically persists, the implicit nature means callers can accidentally write `const std::string& s = some_function_returning_aterm_string();` where lifetime extension does NOT apply (because the reference binds through a conversion operator, not directly to the temporary). Additionally, implicit conversion operators create overload ambiguity with `std::string` parameters.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Verify whether `const std::string& s = aterm_string("hello");` produces a dangling reference
  - Run command: Compile with `-fsanitize=address` and exercise the pattern
  - Observed result: Likely safe in practice due to term pool persistence, but the pattern is undefined behavior per the standard (temporary `aterm_string` destroyed, reference to internal data of a pool-managed object survives only coincidentally).
  - Expected result: Either mark the conversion operator `explicit` or document that the returned reference has pool-managed lifetime.
- Efficiency impact: None directly, but implicit conversions can cause unexpected temporary `std::string` copies at call sites expecting value semantics.
- C++20-specific note: In C++20, the idiomatic approach is `explicit operator const std::string&()` with a named `.str()` or `.name()` accessor for the common case.
- Suggested fix: Mark `explicit` and add a named accessor:
  ```cpp
  explicit operator const std::string&() const { return function().name(); }
  const std::string& str() const { return function().name(); }
  ```

### [SEV-3] Non-explicit converting constructors allow silent implicit conversions

- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h:35-48`
- Why this is a problem: `aterm_string(const function_symbol& f)` and `aterm_string(const std::string& s)` are both non-explicit, enabling implicit conversions. A function accepting `const aterm_string&` can silently accept a raw `std::string` or `function_symbol`, triggering allocation of a new function symbol in the global pool as a side effect of what looks like a simple function call. This is a classic "action at a distance" pattern.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Measure allocation count when passing `std::string` to a function expecting `aterm_string`
  - Run command: N/A (design review finding)
  - Observed result: The codebase deliberately relies on these implicit conversions (e.g., `pbeschain.h:112` comparisons). Changing this would be a breaking API modification.
  - Expected result: Ideally `explicit` constructors with deliberate call-site conversions.
- Efficiency impact: Each implicit conversion from `std::string` allocates or looks up a `function_symbol` in the global pool — this is hidden cost at every implicit conversion site.
- C++20-specific note: C++20 `explicit(bool)` conditional explicitness could be used to make these explicit in new code while preserving backward compatibility via a macro/policy.
- Suggested fix: This is a known design trade-off in the codebase. No immediate fix recommended without broader API discussion. Document the allocation side-effect in the constructor's doxygen comment.

### [SEV-3] `empty_string()` static destruction order concern

- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h:75-79`
- Why this is a problem: The function returns a `const` reference to a `static` local `aterm_string`. If any other static-lifetime object (in another TU) holds a reference to this and the term pool is destroyed before this static, accessing the reference after pool teardown is undefined behavior. This is the classic "static destruction order fiasco."
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Create a static object in another TU whose destructor accesses `empty_string()`
  - Run command: Requires multi-TU test with deliberate static ordering
  - Observed result: Likely safe in current usage (the term pool is typically the last thing torn down), but architecturally fragile.
  - Expected result: Use a leaked singleton (`static aterm_string& t = *new aterm_string("")`) if absolute safety is required, or document the constraint.
- Efficiency impact: None.
- C++20-specific note: `constinit` cannot be used here since `aterm_string("")` is not constexpr-constructible (it allocates in the term pool).
- Suggested fix: Low priority. Document that `empty_string()` must not be called during static destruction, or use the leak-singleton pattern if this guarantee cannot be upheld.

### [SEV-3] `std::swap` specialization in namespace `std` — deprecated practice

- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_string.h:92-96`
- Why this is a problem: Explicit specialization of `std::swap` in namespace `std` is valid but has been superseded by hidden-friend ADL swap in modern C++. The current approach also prevents partial specialization (not applicable here, but limits future templated derivatives). More importantly, the specialization delegates to `t1.swap(t2)` which is inherited from `unprotected_aterm_core` and only swaps the internal pointer — this is correct but non-obvious.
- Evidence type: N/A (style/modernization observation)
- Efficiency impact: None.
- C++20-specific note: In C++20, the preferred pattern is a `friend void swap(aterm_string& a, aterm_string& b) noexcept` defined inline in the class body, found via ADL.
- Suggested fix: Low priority. Consider adding a `friend swap` in the class for ADL discovery while keeping the `std::swap` specialization for backward compatibility.

## Proposed Tests

- [ ] `test_aterm_string_from_nonzero_arity`: Construct `aterm_string` from an `aterm` with arity > 0 in a release build. Verify the class either throws or aborts (currently it silently succeeds). Input: `aterm(function_symbol("f", 2), aterm_string("a"), aterm_string("b"))`. Expected failure signal: exception or assertion failure.
- [ ] `test_aterm_string_implicit_conversion_allocation`: Measure that passing `std::string("x")` where `aterm_string` is expected causes exactly one function symbol pool lookup/insertion. Input: function accepting `const aterm_string&` called with `std::string`. Expected signal: pool size increases by at most 1.
- [ ] `test_empty_string_identity`: Verify `&empty_string() == &empty_string()` (same object returned on repeated calls). Input: two calls to `empty_string()`. Expected signal: pointer equality.

## Open Questions

- Is the decision to keep converting constructors non-explicit a deliberate API choice, or historical inertia? The codebase relies on it, but it hides allocation costs.
- Has the team considered marking `operator const std::string&()` as `explicit` now that C++20 is the baseline? The call sites using `std::string(x)` already perform explicit conversion.
- Are there any known cases of `aterm_string` being constructed from non-zero-arity aterms in production code paths (i.e., has the `assert` ever fired in CI)?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         4m 11s
Total session time:     4m 33s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         81.8k in, 4.3k out, 64.1k cached (Est. 1 Premium request)
 claude-haiku-4.5        848.6k in, 12.5k out, 781.7k cached (Est. 0 Premium requests)

```

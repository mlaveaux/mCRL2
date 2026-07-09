---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/aterm_int.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-07T23:47:11.619698+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h` — a thin wrapper around the internal `_aterm_int` storage.
- Risk level: **Medium** — contains undefined behavior in a debug-mode assertion and lacks defensive preconditions in a noexcept accessor.
- Overall verdict: Functional in the happy path, but harbors latent UB that surfaces in debug builds and provides no safety net in release builds for misuse of `value()`.

## Findings (ordered by severity)

### [SEV-1] Assert in converting constructor dereferences null pointer (UB in debug builds)

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:40`
- **Why this is a problem:** The assertion `assert(type_is_int() || !defined())` evaluates `type_is_int()` first due to left-to-right evaluation of `||`. `type_is_int()` (in `aterm_core.h:65-69`) unconditionally dereferences `m_term` via `m_term->function()`. When the source `aterm t` is default-constructed (i.e., `m_term == nullptr`), this is a null pointer dereference — undefined behavior — before `!defined()` is ever evaluated.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: `test_aterm_int_from_default_aterm`
  - Code:
    ```cpp
    atermpp::aterm default_term;        // m_term == nullptr
    atermpp::aterm_int i(default_term); // assert dereferences nullptr
    ```
  - Run command: `ctest --test-dir build -R aterm_test --output-on-failure` (after adding to `libraries/atermpp/test/aterm_test.cpp`)
  - Observed result: UB (crash/SIGSEGV in debug build, silent in release)
  - Expected result: Assert should pass for undefined terms (the `|| !defined()` branch)
- **Efficiency impact:** None (debug-only path).
- **C++20-specific note:** N/A — this is a fundamental C/C++ correctness issue.
- **Suggested fix:** Swap operand order to ensure short-circuit safety:
  ```cpp
  assert(!defined() || type_is_int());
  ```

---

### [SEV-2] `value()` has no precondition check — UB on default-constructed instance

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:47-50`
- **Why this is a problem:** Calling `value()` on a default-constructed `aterm_int` (where `m_term == nullptr`) performs `reinterpret_cast<const detail::_aterm_int*>(nullptr)->value()` — null pointer dereference, undefined behavior. There is no assert or precondition guarding this path, providing zero diagnostic even in debug builds.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: `test_aterm_int_value_on_default`
  - Code:
    ```cpp
    atermpp::aterm_int i;  // default-constructed, m_term == nullptr
    auto v = i.value();    // UB: null dereference
    ```
  - Run command: `ctest --test-dir build -R aterm_test --output-on-failure` (after adding to test suite)
  - Observed result: UB (likely crash under ASan/debug, silent corruption otherwise)
  - Expected result: Assertion failure with clear diagnostic
- **Efficiency impact:** Adding `assert(defined())` is zero-cost in release builds.
- **C++20-specific note:** A `[[pre: defined()]]` contract (C++26 contracts) would be the ideal long-term solution; `assert` is appropriate today.
- **Suggested fix:**
  ```cpp
  std::size_t value() const noexcept
  {
    assert(defined());
    return reinterpret_cast<const detail::_aterm_int*>(m_term)->value();
  }
  ```

---

### [SEV-3] `value()` accessor lacks `[[nodiscard]]`

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:47`
- **Why this is a problem:** Discarding the return value of `value()` is almost certainly a bug (the function has no side effects). `[[nodiscard]]` provides a compile-time diagnostic for this mistake.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Accidental `i.value();` statement (e.g., confused with a mutating operation) compiles silently.
  - Test plan: Compiler warning audit with `-Wunused-result` plus `[[nodiscard]]`.
- **Efficiency impact:** None.
- **C++20-specific note:** `[[nodiscard]]` has been standard since C++17 and is expected on pure accessors in modern code.
- **Suggested fix:**
  ```cpp
  [[nodiscard]] std::size_t value() const noexcept
  ```

---

### [SEV-3] `std::swap` explicit specialization — legacy pattern

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:68-76`
- **Why this is a problem:** Explicit specializations in `namespace std` are permitted but discouraged in modern C++. A hidden-friend `swap` found via ADL is the idiomatic C++20 approach and avoids potential ODR issues in complex build configurations.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: In large translation● Request failed due to a transient API error. Retrying...

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 60-70
  └ 10 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 105-120
  └ 16 lines read

Now I have all the context needed to produce the review.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h` (single header, public interface for integer terms)
- Risk level: **Low-Medium** — One confirmed UB in debug builds, minor API safety gaps.
- Overall verdict: Functional in release builds but contains a logic error in a debug assertion that causes undefined behavior on a valid input path.

## Findings (ordered by severity)

### [SEV-2] Assert operands reversed — null dereference UB in debug builds
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:40`
- Why this is a problem: The assertion `assert(type_is_int() || !defined())` evaluates `type_is_int()` first. `type_is_int()` (in `aterm_core.h:67`) unconditionally dereferences `m_term` via `m_term->function()`. When constructing `aterm_int` from a default-constructed (undefined) `aterm`, `m_term` is `nullptr`. Since `||` evaluates left-to-right, `type_is_int()` triggers a null pointer dereference **before** `!defined()` can short-circuit.
- Evidence type: **static-proof**
- Evidence:
  - Reproducer/test name: `test_aterm_int_from_undefined`
  - Run command: Add to `libraries/atermpp/test/aterm_test.cpp` and run `ctest --test-dir build -R aterm_test --output-on-failure`
  - Code:
    ```cpp
    BOOST_AUTO_TEST_CASE(test_aterm_int_from_undefined)
    {
      atermpp::aterm undefined_term;       // m_term == nullptr
      atermpp::aterm_int i(undefined_term); // triggers type_is_int() → null deref
      BOOST_CHECK(!i.defined());
    }
    ```
  - Observed result: Undefined behavior (segfault or sanitizer error) in debug builds.
  - Expected result: Construction succeeds; the resulting `aterm_int` is undefined.
- Efficiency impact: None (debug-only path).
- C++20-specific note: N/A.
- Suggested fix: Swap the operands to check `defined()` first:
  ```cpp
  assert(!defined() || type_is_int());
  ```

### [SEV-3] `value()` dereferences `m_term` without precondition guard
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:48-50`
- Why this is a problem: Calling `value()` on a default-constructed `aterm_int` (where `m_term == nullptr`) performs `reinterpret_cast` on null and dereferences it — undefined behavior. No assertion or documentation states the precondition.
- Evidence type: **plausible-issue**
- Evidence:
  - Mechanism: Any code path that calls `.value()` on a default-constructed or moved-from `aterm_int` hits UB.
  - Triggering scenario: A container of `aterm_int` is resized with default-inserted elements; iteration calls `.value()` without checking `defined()`.
  - Test plan: Call `aterm_int{}.value()` under AddressSanitizer and verify it is flagged.
  - Run command: Build with `-DCMAKE_CXX_FLAGS="-fsanitize=address,undefined"` and run `ctest --test-dir build -R aterm_test`
- Efficiency impact: None — adding an assert would be zero-cost in release.
- C++20-specific note: Could use `[[assume(m_term != nullptr)]]` (C++23) or at minimum add a debug `assert(defined())`.
- Suggested fix: Add `assert(defined());` at the start of `value()`.

### [SEV-3] Missing `[[nodiscard]]` on `value()` accessor
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:46`
- Why this is a problem: Discarding the return value of a pure accessor is almost certainly a bug at the call site. C++17/20 `[[nodiscard]]` would catch this at compile time.
- Evidence type: **plausible-issue**
- Evidence:
  - Mechanism: Caller writes `t.value();` as a no-op statement (typo or dead code).
  - Triggering scenario: Low probability but trivially preventable.
- Efficiency impact: None.
- C++20-specific note: `[[nodiscard]]` is the idiomatic C++17+ annotation for pure accessors.
- Suggested fix: `[[nodiscard]] std::size_t value() const noexcept`

## Proposed Tests
- [ ] `test_aterm_int_from_undefined`: Construct `aterm_int` from a default `aterm`; verify no crash and `!defined()` holds. Exposes the SEV-2 assert-order bug in debug builds.
- [ ] `test_aterm_int_value_precondition`: Under ASan/UBSan, call `value()` on default-constructed `aterm_int`; expect sanitizer report. Validates SEV-3 precondition gap.

## Open Questions
- Is constructing `aterm_int` from an undefined `aterm` intentionally supported, or should the constructor reject it entirely (removing the `|| !defined()` branch)?
- Should `value()` be hardened with `assert(defined())` consistently across the `aterm_*` type family, or is the precondition-on-caller pattern the project's deliberate choice?

## Copilot Raw Stdout

```text
● Explore: Explore aterm_int dependencies
  ● Check how m_term is used in aterm_int
    $ grep -n "m_term" /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm_int.h
    └ 2 lines...
  ● Read libraries/atermpp/test/constructors.cpp
    └ 39 lines read

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h` — a thin wrapper around the internal `_aterm_int` storage.
- Risk level: **Medium** — contains undefined behavior in a debug-mode assertion and lacks defensive preconditions in a noexcept accessor.
- Overall verdict: Functional in the happy path, but harbors latent UB that surfaces in debug builds and provides no safety net in release builds for misuse of `value()`.

## Findings (ordered by severity)

### [SEV-1] Assert in converting constructor dereferences null pointer (UB in debug builds)

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:40`
- **Why this is a problem:** The assertion `assert(type_is_int() || !defined())` evaluates `type_is_int()` first due to left-to-right evaluation of `||`. `type_is_int()` (in `aterm_core.h:65-69`) unconditionally dereferences `m_term` via `m_term->function()`. When the source `aterm t` is default-constructed (i.e., `m_term == nullptr`), this is a null pointer dereference — undefined behavior — before `!defined()` is ever evaluated.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: `test_aterm_int_from_default_aterm`
  - Code:
    ```cpp
    atermpp::aterm default_term;        // m_term == nullptr
    atermpp::aterm_int i(default_term); // assert dereferences nullptr
    ```
  - Run command: `ctest --test-dir build -R aterm_test --output-on-failure` (after adding to `libraries/atermpp/test/aterm_test.cpp`)
  - Observed result: UB (crash/SIGSEGV in debug build, silent in release)
  - Expected result: Assert should pass for undefined terms (the `|| !defined()` branch)
- **Efficiency impact:** None (debug-only path).
- **C++20-specific note:** N/A — this is a fundamental C/C++ correctness issue.
- **Suggested fix:** Swap operand order to ensure short-circuit safety:
  ```cpp
  assert(!defined() || type_is_int());
  ```

---

### [SEV-2] `value()` has no precondition check — UB on default-constructed instance

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:47-50`
- **Why this is a problem:** Calling `value()` on a default-constructed `aterm_int` (where `m_term == nullptr`) performs `reinterpret_cast<const detail::_aterm_int*>(nullptr)->value()` — null pointer dereference, undefined behavior. There is no assert or precondition guarding this path, providing zero diagnostic even in debug builds.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: `test_aterm_int_value_on_default`
  - Code:
    ```cpp
    atermpp::aterm_int i;  // default-constructed, m_term == nullptr
    auto v = i.value();    // UB: null dereference
    ```
  - Run command: `ctest --test-dir build -R aterm_test --output-on-failure` (after adding to test suite)
  - Observed result: UB (likely crash under ASan/debug, silent corruption otherwise)
  - Expected result: Assertion failure with clear diagnostic
- **Efficiency impact:** Adding `assert(defined())` is zero-cost in release builds.
- **C++20-specific note:** A `[[pre: defined()]]` contract (C++26 contracts) would be the ideal long-term solution; `assert` is appropriate today.
- **Suggested fix:**
  ```cpp
  std::size_t value() const noexcept
  {
    assert(defined());
    return reinterpret_cast<const detail::_aterm_int*>(m_term)->value();
  }
  ```

---

### [SEV-3] `value()` accessor lacks `[[nodiscard]]`

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:47`
- **Why this is a problem:** Discarding the return value of `value()` is almost certainly a bug (the function has no side effects). `[[nodiscard]]` provides a compile-time diagnostic for this mistake.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Accidental `i.value();` statement (e.g., confused with a mutating operation) compiles silently.
  - Test plan: Compiler warning audit with `-Wunused-result` plus `[[nodiscard]]`.
- **Efficiency impact:** None.
- **C++20-specific note:** `[[nodiscard]]` has been standard since C++17 and is expected on pure accessors in modern code.
- **Suggested fix:**
  ```cpp
  [[nodiscard]] std::size_t value() const noexcept
  ```

---

### [SEV-3] `std::swap` explicit specialization — legacy pattern

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:68-76`
- **Why this is a problem:** Explicit specializations in `namespace std` are permitted but discouraged in modern C++. A hidden-friend `swap` found via ADL is the idiomatic C++20 approach and avoids potential ODR issues in complex build configurations.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: In large translation● Request failed due to a transient API error. Retrying...

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 60-70
  └ 10 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 105-120
  └ 16 lines read

Now I have all the context needed to produce the review.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h` (single header, public interface for integer terms)
- Risk level: **Low-Medium** — One confirmed UB in debug builds, minor API safety gaps.
- Overall verdict: Functional in release builds but contains a logic error in a debug assertion that causes undefined behavior on a valid input path.

## Findings (ordered by severity)

### [SEV-2] Assert operands reversed — null dereference UB in debug builds
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:40`
- Why this is a problem: The assertion `assert(type_is_int() || !defined())` evaluates `type_is_int()` first. `type_is_int()` (in `aterm_core.h:67`) unconditionally dereferences `m_term` via `m_term->function()`. When constructing `aterm_int` from a default-constructed (undefined) `aterm`, `m_term` is `nullptr`. Since `||` evaluates left-to-right, `type_is_int()` triggers a null pointer dereference **before** `!defined()` can short-circuit.
- Evidence type: **static-proof**
- Evidence:
  - Reproducer/test name: `test_aterm_int_from_undefined`
  - Run command: Add to `libraries/atermpp/test/aterm_test.cpp` and run `ctest --test-dir build -R aterm_test --output-on-failure`
  - Code:
    ```cpp
    BOOST_AUTO_TEST_CASE(test_aterm_int_from_undefined)
    {
      atermpp::aterm undefined_term;       // m_term == nullptr
      atermpp::aterm_int i(undefined_term); // triggers type_is_int() → null deref
      BOOST_CHECK(!i.defined());
    }
    ```
  - Observed result: Undefined behavior (segfault or sanitizer error) in debug builds.
  - Expected result: Construction succeeds; the resulting `aterm_int` is undefined.
- Efficiency impact: None (debug-only path).
- C++20-specific note: N/A.
- Suggested fix: Swap the operands to check `defined()` first:
  ```cpp
  assert(!defined() || type_is_int());
  ```

### [SEV-3] `value()` dereferences `m_term` without precondition guard
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:48-50`
- Why this is a problem: Calling `value()` on a default-constructed `aterm_int` (where `m_term == nullptr`) performs `reinterpret_cast` on null and dereferences it — undefined behavior. No assertion or documentation states the precondition.
- Evidence type: **plausible-issue**
- Evidence:
  - Mechanism: Any code path that calls `.value()` on a default-constructed or moved-from `aterm_int` hits UB.
  - Triggering scenario: A container of `aterm_int` is resized with default-inserted elements; iteration calls `.value()` without checking `defined()`.
  - Test plan: Call `aterm_int{}.value()` under AddressSanitizer and verify it is flagged.
  - Run command: Build with `-DCMAKE_CXX_FLAGS="-fsanitize=address,undefined"` and run `ctest --test-dir build -R aterm_test`
- Efficiency impact: None — adding an assert would be zero-cost in release.
- C++20-specific note: Could use `[[assume(m_term != nullptr)]]` (C++23) or at minimum add a debug `assert(defined())`.
- Suggested fix: Add `assert(defined());` at the start of `value()`.

### [SEV-3] Missing `[[nodiscard]]` on `value()` accessor
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_int.h:46`
- Why this is a problem: Discarding the return value of a pure accessor is almost certainly a bug at the call site. C++17/20 `[[nodiscard]]` would catch this at compile time.
- Evidence type: **plausible-issue**
- Evidence:
  - Mechanism: Caller writes `t.value();` as a no-op statement (typo or dead code).
  - Triggering scenario: Low probability but trivially preventable.
- Efficiency impact: None.
- C++20-specific note: `[[nodiscard]]` is the idiomatic C++17+ annotation for pure accessors.
- Suggested fix: `[[nodiscard]] std::size_t value() const noexcept`

## Proposed Tests
- [ ] `test_aterm_int_from_undefined`: Construct `aterm_int` from a default `aterm`; verify no crash and `!defined()` holds. Exposes the SEV-2 assert-order bug in debug builds.
- [ ] `test_aterm_int_value_precondition`: Under ASan/UBSan, call `value()` on default-constructed `aterm_int`; expect sanitizer report. Validates SEV-3 precondition gap.

## Open Questions
- Is constructing `aterm_int` from an undefined `aterm` intentionally supported, or should the constructor reject it entirely (removing the `|| !defined()` branch)?
- Should `value()` be hardened with `assert(defined())` consistently across the `aterm_*` type family, or is the precondition-on-caller pattern the project's deliberate choice?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         2m 22s
Total session time:     3m 27s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         113.8k in, 4.1k out, 100.6k cached (Est. 1 Premium request)
 claude-haiku-4.5        386.9k in, 5.6k out, 350.4k cached (Est. 0 Premium requests)

```

---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T09:29:16.228843+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h` — a small internal header defining the integer-carrying aterm node.
- **Risk level:** Low. This is a stable, simple leaf class with minimal surface area.
- **Overall verdict:** Functionally correct. A few minor hygiene issues exist; no confirmed defects found.

## Findings (ordered by severity)

### [SEV-3] Trailing whitespace on line 33
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h:33`
- **Why this is a problem:** Violates `.clang-format` formatting expectations. The line between `value()` and `private:` contains trailing spaces.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer: `grep -nP ' +$' libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h`
  - Run command: `clang-format --dry-run -Werror libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h`
  - Observed result: Line 33 has trailing whitespace.
  - Expected result: No trailing whitespace.
- **Suggested fix:** Remove trailing spaces on line 33.

### [SEV-3] `value()` getter lacks `[[nodiscard]]`
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h:29`
- **Why this is a problem:** Calling `value()` without using the result is always a bug. `[[nodiscard]]` would catch this at compile time. The C++ instructions require `[[nodiscard]]` on pure queries.
- **Evidence type:** plausible-issue
- **Evidence:**
  - This is a style/safety improvement. No confirmed silent misuse found in the current codebase, but the project instructions mandate it for must-use results.
- **C++20-specific note:** `[[nodiscard]]` is standard since C++17 and trivially applicable here.
- **Suggested fix:** `[[nodiscard]] std::size_t value() const noexcept`

### [SEV-3] Constructor lacks `explicit` keyword
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h:24`
- **Why this is a problem:** Single-argument constructors should be `explicit` to prevent accidental implicit conversions from `std::size_t`. While `_aterm` is `noncopyable` (limiting abuse scenarios), and the existing codebase convention does not use `explicit` on sibling classes (`_aterm_appl`), this remains a deviation from modern C++ best practice.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mitigated by `noncopyable` semantics — implicit conversion would require constructing in place, which pool allocation handles explicitly. Risk is theoretical.
  - Convention: `_aterm_appl` also omits `explicit`, so this is consistent with surrounding code.
- **C++20-specific note:** Core Guidelines C.46 recommends `explicit` for single-argument constructors.
- **Suggested fix:** Add `explicit`, but note this is a broader codebase convention gap, not specific to this file.

### [SEV-3] Missing Doxygen on `value()` accessor
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h:29`
- **Why this is a problem:** Project instructions require Doxygen for public APIs. `value()` lacks `/// \brief` and `\returns`.
- **Evidence type:** static-proof (documentation requirement)
- **Suggested fix:** Add `/// \brief Returns the integer value stored in this term.`

## Non-issues noted
- **No virtual destructor in `_aterm`:** This is by design. Objects live in a pool allocator and are destroyed via `m_term_set.erase()` which knows the concrete type. No polymorphic delete occurs.
- **`static_assert` on layout:** Excellent — confirms no unexpected padding. Correct for the pool/hash-table memory model.
- **`noncopyable` inheritance:** Appropriate for pool-managed objects.
- **Constructor not `noexcept`:** The base `_aterm(const function_symbol&)` is not `noexcept`, so marking this constructor `noexcept` would be incorrect without first auditing the base. No defect here.

## Proposed Tests
- None required. This header defines a trivial data-carrying node. Existing pool allocation tests cover construction. No behavioral defects identified.

## Open Questions
- Should a codebase-wide pass add `[[nodiscard]]` to all aterm accessors (`_aterm::function()`, `_aterm_appl::arg()`, `_aterm_int::value()`) uniformly rather than piecemeal?
- Should `explicit` be adopted across all single-argument `_aterm*` constructors as a batch modernization?

## Copilot Raw Stdout

```text
● Explore: Gather context for review
  ● Check _aterm_appl constructors for explicit/noexcept
    $ grep -n "noexcept\|explicit"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h | grep
    "constructor\|_aterm_appl(" | head -20
    └ 1 line...
  ● Get _aterm_appl constructor details
    $ sed -n '37,46p' /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h
    └ 10 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h
  └ 44 lines read

✗ Check trailing whitespace on line 33
  $ cd /workspaces/mCRL2 && cat -A libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h | sed
  -n '33p'
  Permission denied and could not request permission from user

● Grep " +$" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h)
  └ 1 line found

Now I have all the context needed. Here is my review:

---

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h` — a small internal header defining the integer-carrying aterm node.
- **Risk level:** Low. This is a stable, simple leaf class with minimal surface area.
- **Overall verdict:** Functionally correct. A few minor hygiene issues exist; no confirmed defects found.

## Findings (ordered by severity)

### [SEV-3] Trailing whitespace on line 33
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h:33`
- **Why this is a problem:** Violates `.clang-format` formatting expectations. The line between `value()` and `private:` contains trailing spaces.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer: `grep -nP ' +$' libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h`
  - Run command: `clang-format --dry-run -Werror libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h`
  - Observed result: Line 33 has trailing whitespace.
  - Expected result: No trailing whitespace.
- **Suggested fix:** Remove trailing spaces on line 33.

### [SEV-3] `value()` getter lacks `[[nodiscard]]`
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h:29`
- **Why this is a problem:** Calling `value()` without using the result is always a bug. `[[nodiscard]]` would catch this at compile time. The C++ instructions require `[[nodiscard]]` on pure queries.
- **Evidence type:** plausible-issue
- **Evidence:**
  - This is a style/safety improvement. No confirmed silent misuse found in the current codebase, but the project instructions mandate it for must-use results.
- **C++20-specific note:** `[[nodiscard]]` is standard since C++17 and trivially applicable here.
- **Suggested fix:** `[[nodiscard]] std::size_t value() const noexcept`

### [SEV-3] Constructor lacks `explicit` keyword
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h:24`
- **Why this is a problem:** Single-argument constructors should be `explicit` to prevent accidental implicit conversions from `std::size_t`. While `_aterm` is `noncopyable` (limiting abuse scenarios), and the existing codebase convention does not use `explicit` on sibling classes (`_aterm_appl`), this remains a deviation from modern C++ best practice.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mitigated by `noncopyable` semantics — implicit conversion would require constructing in place, which pool allocation handles explicitly. Risk is theoretical.
  - Convention: `_aterm_appl` also omits `explicit`, so this is consistent with surrounding code.
- **C++20-specific note:** Core Guidelines C.46 recommends `explicit` for single-argument constructors.
- **Suggested fix:** Add `explicit`, but note this is a broader codebase convention gap, not specific to this file.

### [SEV-3] Missing Doxygen on `value()` accessor
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_int.h:29`
- **Why this is a problem:** Project instructions require Doxygen for public APIs. `value()` lacks `/// \brief` and `\returns`.
- **Evidence type:** static-proof (documentation requirement)
- **Suggested fix:** Add `/// \brief Returns the integer value stored in this term.`

## Non-issues noted
- **No virtual destructor in `_aterm`:** This is by design. Objects live in a pool allocator and are destroyed via `m_term_set.erase()` which knows the concrete type. No polymorphic delete occurs.
- **`static_assert` on layout:** Excellent — confirms no unexpected padding. Correct for the pool/hash-table memory model.
- **`noncopyable` inheritance:** Appropriate for pool-managed objects.
- **Constructor not `noexcept`:** The base `_aterm(const function_symbol&)` is not `noexcept`, so marking this constructor `noexcept` would be incorrect without first auditing the base. No defect here.

## Proposed Tests
- None required. This header defines a trivial data-carrying node. Existing pool allocation tests cover construction. No behavioral defects identified.

## Open Questions
- Should a codebase-wide pass add `[[nodiscard]]` to all aterm accessors (`_aterm::function()`, `_aterm_appl::arg()`, `_aterm_int::value()`) uniformly rather than piecemeal?
- Should `explicit` be adopted across all single-argument `_aterm*` constructors as a batch modernization?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         2m 10s
Total session time:     2m 26s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         196.7k in, 2.9k out, 178.0k cached (Est. 1 Premium request)
 claude-haiku-4.5        315.0k in, 6.0k out, 285.6k cached (Est. 0 Premium requests)

```

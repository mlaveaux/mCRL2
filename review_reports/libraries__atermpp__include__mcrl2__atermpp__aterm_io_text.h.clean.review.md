---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/aterm_io_text.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-08T00:00:39.228500+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_io_text.h` (header) with full implementation context from `aterm_io_text.cpp`
- Risk level: **Medium** — contains a confirmed UB path, a logic bug, and a design hazard
- Overall verdict: **Functional for well-formed inputs; several defects in edge cases and one confirmed undefined behavior**

## Findings (ordered by severity)

### [SEV-1] Buffer overwrite past array bounds in `parse_aterm_int`
- Location: `libraries/atermpp/source/aterm_io_text.cpp:220`
- Why this is a problem: `std::array<char, 32> number` is filled until `it == number.end()`. When exactly 32 digits (or 31 + sign) are consumed, `it` equals `number.end()` and then `*it = '\0'` writes one past the end of the array — undefined behavior.
- Evidence type: static-proof
- Evidence:
  - Input: `read_term_from_string("12345678901234567890123456789012")` (32 digits)
  - The loop exits with `it == number.end()` (past last valid index). Line 220 dereferences this invalid iterator.
  - Run command: `cd build && ctest -R aterm_test --output-on-failure` (after adding reproducer to `libraries/atermpp/test/aterm_test.cpp`)
  - Observed result: Undefined behavior (stack corruption, ASan error, or silent data corruption)
  - Expected result: Graceful parsing or a clear error
- Efficiency impact: N/A (correctness issue)
- C++20-specific note: Could use `std::string` with `reserve(32)` or `std::from_chars` on the accumulated string to avoid fixed-buffer risks entirely.
- Suggested fix:
```cpp
// Either: guard the null-terminator write
if (it != number.end()) {
  *it = '\0';
} else {
  throw std::runtime_error("Integer literal too long");
}
// Or: replace with std::string accumulation + std::from_chars
```

### [SEV-2] Duplicate condition in `parse_unquoted_string` — `'['` not excluded
- Location: `libraries/atermpp/source/aterm_io_text.cpp:363-364`
- Why this is a problem: The while-condition checks `character != ']'` twice (copy-paste error). The second `']'` should almost certainly be `'['`. As a result, `[` is accepted within unquoted strings, which can cause incorrect parsing when an unquoted symbol is immediately followed by a list.
- Evidence type: static-proof
- Evidence:
  - Code: `character != ']' && character != ']'` — identical conditions, one is redundant
  - Input: `read_term_from_string("f[1,2]")` — the unquoted parser consumes `f[1,2]` as a single string instead of stopping at `[`
  - Run command: Add `test_aterm_io("f[1,2]")` to test and observe it does not parse as expected
  - Observed result: Parsed as function application with name containing `[1,2]` or parse error
  - Expected result: Parser should stop at `[` since it begins a list
- Efficiency impact: None (logic bug)
- C++20-specific note: N/A
- Suggested fix: Change the second `']'` to `'['`:
```cpp
while (character != '"' && character != '(' && character != ')' && character != '['
  && character != ']' && character != ',' && ...
```

### [SEV-2] Negative integer parsed via `atol` then cast to `std::size_t` (unsigned wraparound)
- Location: `libraries/atermpp/source/aterm_io_text.cpp:206-221`
- Why this is a problem: `parse_aterm_int` accepts a leading `-` sign, calls `atol()` (returns `long`), then casts to `std::size_t`. For negative inputs, this silently wraps to a large unsigned value. `aterm_int` stores `std::size_t`, so the sign is lost.
- Evidence type: static-proof
- Evidence:
  - `atol("-1")` returns `-1L`. `static_cast<std::size_t>(-1L)` == `SIZE_MAX`.
  - `aterm_int(SIZE_MAX).value()` != `-1` in any meaningful sense.
  - Run command: `read_term_from_string("-1")` → verify `down_cast<aterm_int>(t).value()` is `SIZE_MAX`
  - Observed result: Value is `SIZE_MAX` (18446744073709551615 on 64-bit)
  - Expected result: Either reject negative integers or preserve the sign
- Efficiency impact: None (semantic correctness)
- C++20-specific note: `std::from_chars` with `std::size_t` would reject negative input at parse time.
- Suggested fix: Either remove the `-` handling (if negative ints are not supported by `aterm_int`) or use a signed type with validation.

### [SEV-3] Member variable `character` shadows method parameter `int& character`
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_io_text.h:81` vs lines 47-63
- Why this is a problem: Every parse method takes `int& character` which shadows the member `int character`. In `next_char()` (which takes no such parameter), the member is modified directly. The intent is that `get()` passes the member by reference, so modifications propagate, but this pattern is fragile — any future method that forgets the reference parameter will silently read stale state.
- Evidence type: plausible-issue
- Evidence:
  - Mechanism: A new helper method that accesses `character` without the reference parameter would use the member directly. If `next_char()` updates the member but the caller holds a local copy, state diverges.
  - Current code works because `get()` always passes `this->character` by reference into all methods.
  - Test plan: Static analysis with `-Wshadow` would flag this immediately.
  - Run command: `g++ -Wshadow -c libraries/atermpp/source/aterm_io_text.cpp -I libraries/atermpp/include/`
- Efficiency impact: None
- C++20-specific note: Modern style would use `m_character` for the member to avoid any ambiguity.
- Suggested fix: Rename member to `m_character` (consistent with `m_stream`, `m_line`, `m_column`, `m_history`).

### [SEV-3] Include guard comment mismatch
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_io_text.h:83`
- Why this is a problem: The `#endif` comment says `// ATERM_IO_TEXT_H` but the actual guard is `MCRL2_ATERMPP_ATERM_IO_TEXT_H`. Misleading during maintenance.
- Evidence type: static-proof
- Evidence: Line 11: `#ifndef MCRL2_ATERMPP_ATERM_IO_TEXT_H`; last line: `#endif // ATERM_IO_TEXT_H`
- Efficiency impact: None
- C++20-specific note: Could use `#pragma once` (non-standard but universally supported) to eliminate this class of error.
- Suggested fix: `#endif // MCRL2_ATERMPP_ATERM_IO_TEXT_H`

## Proposed Tests
- [ ] `test_parse_long_integer`: Input `"12345678901234567890123456789012"` (32+ digits), expected: either valid parse or `runtime_error`, signal: ASan/UBSan stack-buffer-overflow
- [ ] `test_parse_negative_integer`: Input `"-1"`, verify `aterm_int::value()` round-trips correctly or throws
- [ ] `test_unquoted_with_bracket`: Input `"f[1,2]"`, expected: parse error or correct decomposition, signal: wrong function name containing `[`
- [ ] `test_empty_function_symbol`: Input `""()` (quoted empty string with args), verify round-trip correctness

## Open Questions
- Is negative integer support in `parse_aterm_int` intentional given `aterm_int` is unsigned? If so, what is the expected semantics?
- Should `parse_unquoted_string` truly stop at `[`? The current test `"f([a,f(x),[]],2,[g,g(34566)])"` uses parenthesized arguments for `f`, never bare `f[...]`. Is `f[1,2]` valid syntax in the aterm text format?
- The `next_char()` conversion of trailing `\n` to `EOF` (line 310) — is this documented behavior or a legacy workaround? It means streams ending with `\n` are treated as complete, but it also silently truncates content after the last term if separated only by newlines.

## Copilot Raw Stdout

```text
● Explore: Gather context for review
  ● Find all uses of 'character' variable in implementation
    $ grep -n "int character" /workspaces/mCRL2/libraries/atermpp/source/aterm_io_text.cpp | head
    -20
    └ 1 line...
  ● Find all assignments to character
    $ grep -n "character = " /workspaces/mCRL2/libraries/atermpp/source/aterm_io_text.cpp
    └ 14 lines...

● Read libraries/atermpp/source/aterm_io_text.cpp
  └ 400 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_io.h
  └ 234 lines read

● Grep "aterm_int" in *aterm_int* (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp)
  └ 2 files found

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_int.h
  └ 86 lines read

● Grep "class aterm_int" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm_int.h)
  └ 11 lines found

Now I have comprehensive context. Let me produce the review.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_io_text.h` (header) with full implementation context from `aterm_io_text.cpp`
- Risk level: **Medium** — contains a confirmed UB path, a logic bug, and a design hazard
- Overall verdict: **Functional for well-formed inputs; several defects in edge cases and one confirmed undefined behavior**

## Findings (ordered by severity)

### [SEV-1] Buffer overwrite past array bounds in `parse_aterm_int`
- Location: `libraries/atermpp/source/aterm_io_text.cpp:220`
- Why this is a problem: `std::array<char, 32> number` is filled until `it == number.end()`. When exactly 32 digits (or 31 + sign) are consumed, `it` equals `number.end()` and then `*it = '\0'` writes one past the end of the array — undefined behavior.
- Evidence type: static-proof
- Evidence:
  - Input: `read_term_from_string("12345678901234567890123456789012")` (32 digits)
  - The loop exits with `it == number.end()` (past last valid index). Line 220 dereferences this invalid iterator.
  - Run command: `cd build && ctest -R aterm_test --output-on-failure` (after adding reproducer to `libraries/atermpp/test/aterm_test.cpp`)
  - Observed result: Undefined behavior (stack corruption, ASan error, or silent data corruption)
  - Expected result: Graceful parsing or a clear error
- Efficiency impact: N/A (correctness issue)
- C++20-specific note: Could use `std::string` with `reserve(32)` or `std::from_chars` on the accumulated string to avoid fixed-buffer risks entirely.
- Suggested fix:
```cpp
// Either: guard the null-terminator write
if (it != number.end()) {
  *it = '\0';
} else {
  throw std::runtime_error("Integer literal too long");
}
// Or: replace with std::string accumulation + std::from_chars
```

### [SEV-2] Duplicate condition in `parse_unquoted_string` — `'['` not excluded
- Location: `libraries/atermpp/source/aterm_io_text.cpp:363-364`
- Why this is a problem: The while-condition checks `character != ']'` twice (copy-paste error). The second `']'` should almost certainly be `'['`. As a result, `[` is accepted within unquoted strings, which can cause incorrect parsing when an unquoted symbol is immediately followed by a list.
- Evidence type: static-proof
- Evidence:
  - Code: `character != ']' && character != ']'` — identical conditions, one is redundant
  - Input: `read_term_from_string("f[1,2]")` — the unquoted parser consumes `f[1,2]` as a single string instead of stopping at `[`
  - Run command: Add `test_aterm_io("f[1,2]")` to test and observe it does not parse as expected
  - Observed result: Parsed as function application with name containing `[1,2]` or parse error
  - Expected result: Parser should stop at `[` since it begins a list
- Efficiency impact: None (logic bug)
- C++20-specific note: N/A
- Suggested fix: Change the second `']'` to `'['`:
```cpp
while (character != '"' && character != '(' && character != ')' && character != '['
  && character != ']' && character != ',' && ...
```

### [SEV-2] Negative integer parsed via `atol` then cast to `std::size_t` (unsigned wraparound)
- Location: `libraries/atermpp/source/aterm_io_text.cpp:206-221`
- Why this is a problem: `parse_aterm_int` accepts a leading `-` sign, calls `atol()` (returns `long`), then casts to `std::size_t`. For negative inputs, this silently wraps to a large unsigned value. `aterm_int` stores `std::size_t`, so the sign is lost.
- Evidence type: static-proof
- Evidence:
  - `atol("-1")` returns `-1L`. `static_cast<std::size_t>(-1L)` == `SIZE_MAX`.
  - `aterm_int(SIZE_MAX).value()` != `-1` in any meaningful sense.
  - Run command: `read_term_from_string("-1")` → verify `down_cast<aterm_int>(t).value()` is `SIZE_MAX`
  - Observed result: Value is `SIZE_MAX` (18446744073709551615 on 64-bit)
  - Expected result: Either reject negative integers or preserve the sign
- Efficiency impact: None (semantic correctness)
- C++20-specific note: `std::from_chars` with `std::size_t` would reject negative input at parse time.
- Suggested fix: Either remove the `-` handling (if negative ints are not supported by `aterm_int`) or use a signed type with validation.

### [SEV-3] Member variable `character` shadows method parameter `int& character`
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_io_text.h:81` vs lines 47-63
- Why this is a problem: Every parse method takes `int& character` which shadows the member `int character`. In `next_char()` (which takes no such parameter), the member is modified directly. The intent is that `get()` passes the member by reference, so modifications propagate, but this pattern is fragile — any future method that forgets the reference parameter will silently read stale state.
- Evidence type: plausible-issue
- Evidence:
  - Mechanism: A new helper method that accesses `character` without the reference parameter would use the member directly. If `next_char()` updates the member but the caller holds a local copy, state diverges.
  - Current code works because `get()` always passes `this->character` by reference into all methods.
  - Test plan: Static analysis with `-Wshadow` would flag this immediately.
  - Run command: `g++ -Wshadow -c libraries/atermpp/source/aterm_io_text.cpp -I libraries/atermpp/include/`
- Efficiency impact: None
- C++20-specific note: Modern style would use `m_character` for the member to avoid any ambiguity.
- Suggested fix: Rename member to `m_character` (consistent with `m_stream`, `m_line`, `m_column`, `m_history`).

### [SEV-3] Include guard comment mismatch
- Location: `libraries/atermpp/include/mcrl2/atermpp/aterm_io_text.h:83`
- Why this is a problem: The `#endif` comment says `// ATERM_IO_TEXT_H` but the actual guard is `MCRL2_ATERMPP_ATERM_IO_TEXT_H`. Misleading during maintenance.
- Evidence type: static-proof
- Evidence: Line 11: `#ifndef MCRL2_ATERMPP_ATERM_IO_TEXT_H`; last line: `#endif // ATERM_IO_TEXT_H`
- Efficiency impact: None
- C++20-specific note: Could use `#pragma once` (non-standard but universally supported) to eliminate this class of error.
- Suggested fix: `#endif // MCRL2_ATERMPP_ATERM_IO_TEXT_H`

## Proposed Tests
- [ ] `test_parse_long_integer`: Input `"12345678901234567890123456789012"` (32+ digits), expected: either valid parse or `runtime_error`, signal: ASan/UBSan stack-buffer-overflow
- [ ] `test_parse_negative_integer`: Input `"-1"`, verify `aterm_int::value()` round-trips correctly or throws
- [ ] `test_unquoted_with_bracket`: Input `"f[1,2]"`, expected: parse error or correct decomposition, signal: wrong function name containing `[`
- [ ] `test_empty_function_symbol`: Input `""()` (quoted empty string with args), verify round-trip correctness

## Open Questions
- Is negative integer support in `parse_aterm_int` intentional given `aterm_int` is unsigned? If so, what is the expected semantics?
- Should `parse_unquoted_string` truly stop at `[`? The current test `"f([a,f(x),[]],2,[g,g(34566)])"` uses parenthesized arguments for `f`, never bare `f[...]`. Is `f[1,2]` valid syntax in the aterm text format?
- The `next_char()` conversion of trailing `\n` to `EOF` (line 310) — is this documented behavior or a legacy workaround? It means streams ending with `\n` are treated as complete, but it also silently truncates content after the last term if separated only by newlines.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         2m 21s
Total session time:     2m 30s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         152.8k in, 5.0k out, 128.7k cached (Est. 1 Premium request)
 claude-haiku-4.5        148.4k in, 3.0k out, 133.9k cached (Est. 0 Premium requests)

```

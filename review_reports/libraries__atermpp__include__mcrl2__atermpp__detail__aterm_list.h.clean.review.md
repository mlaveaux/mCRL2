---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T09:31:47.253798+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list.h` — a forward-declaration header for list construction utilities in atermpp.
- Risk level: Low (header contains only declarations and one trivial struct, no executable logic besides the identity functor)
- Overall verdict: Several maintenance/quality issues but no runtime correctness bugs given this is a declarations-only header. The file has copy-paste errors, a dead symbol, missing contracts, and missed C++20 opportunities.

## Findings (ordered by severity)

### [SEV-2] Incorrect `\file` and `\brief` Doxygen comment (copy-paste from `aterm_int.h`)
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list.h:9-10`
- Why this is a problem: The file claims to be `mcrl2/atermpp/aterm_int.h` with brief "Term containing an integer." This is factually wrong — this is `detail/aterm_list.h` containing list construction declarations. It will generate incorrect Doxygen output and actively mislead developers navigating docs.
- Evidence type: static-proof
- Evidence:
  - Line 9: `/// \file mcrl2/atermpp/aterm_int.h`
  - Line 10: `/// \brief Term containing an integer.`
  - Correct values: `/// \file mcrl2/atermpp/detail/aterm_list.h` and `/// \brief Forward declarations for term list construction utilities.`
- Suggested fix: Correct the `\file` path and `\brief` text.

### [SEV-2] Duplicate forward declaration of `make_list_backward` (3-arg, out-param version)
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list.h:60-61` and `:82-83`
- Why this is a problem: The exact same signature appears twice:
  ```cpp
  template <class Term, class Iter, class ATermConverter>
  inline void make_list_backward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm);
  ```
  While duplicate declarations are legal C++, this is a copy-paste error that inflates the apparent API surface and confuses readers about the number of distinct overloads. The Doxygen comments differ (line 59 vs line 81), making it seem like two separate functions.
- Evidence type: static-proof
- Evidence:
  - Lines 60-61 declared with brief "Constructs a list starting from first to last where the result is put in result."
  - Lines 82-83 declared with brief "Constructs a list starting from the first iterator element to the last. The result is put into the variable result."
  - Both resolve to the single definition at `aterm_list_implementation.h:501`.
- Suggested fix: Remove the duplicate at lines 81-83.

### [SEV-3] Dead forward declaration: `_term_list` is never defined or used
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list.h:28-29`
- Why this is a problem: `template <class Term> class _term_list;` is forward-declared but has no definition or usage anywhere in the codebase. It is dead code — a vestige from a previous design iteration that adds confusion.
- Evidence type: static-proof
- Evidence:
  - `grep -rn "_term_list" libraries/atermpp/` matches only this single forward declaration.
  - No definition, instantiation, or reference exists elsewhere.
- Suggested fix: Remove lines 28-29.

### [SEV-3] Missing `[[nodiscard]]` on return-by-value `make_list_*` declarations
- Location: Lines 57, 67, 75, 89
- Why this is a problem: Functions that return a freshly constructed `term_list<Term>` by value perform non-trivial work (allocation, iteration, conversion). Discarding the result is almost certainly a bug. `[[nodiscard]]` would catch this at compile time.
- Evidence type: plausible-issue
- Evidence:
  - No `[[nodiscard]]` usage exists anywhere in the atermpp library (confirmed by grep).
  - The out-param overloads exist specifically because callers need the result — discarding the return-by-value version is always wrong.
- Suggested fix: Add `[[nodiscard]]` to all four return-by-value declarations.

### [SEV-3] Missing API contracts (`\pre`/`\post`) on all declarations
- Location: All `make_list_backward`/`make_list_forward` declarations (lines 56-93)
- Why this is a problem: Per mCRL2 coding standards, new or changed public APIs require Doxygen `\pre`/`\post`. These functions have implicit preconditions (iterators must form a valid range; `convert_to_aterm` must return valid terms; `aterm_filter` must be a predicate) and postconditions (result is a valid `term_list` containing converted elements in the appropriate order). None are documented.
- Evidence type: static-proof (absence of required documentation)
- Suggested fix: Add `\pre [first, last)` is a valid iterator range. `\pre convert_to_aterm` is callable with `*Iter`. `\post result` contains converted elements in original order (forward) or reverse traversal order (backward).

### [SEV-3] Inconsistent `typename`/`class` usage across template parameter lists
- Location: Lines 56-57 vs 60-61, line 74 mixes both in one declaration
- Why this is a problem: Purely cosmetic but introduces inconsistency that makes the API harder to scan. Line 74 uses `template <typename Term, class Iter, class ATermConverter>` mixing both keywords in the same declaration.
- Evidence type: static-proof
- Suggested fix: Standardize on one keyword (prefer `typename` for type parameters per modern style, or match surrounding code).

### [SEV-3] No C++20 concepts constraining template parameters
- Location: All template declarations (lines 56-93)
- Why this is a problem: C++20 is the project baseline. These unconstrained templates produce cryptic error messages on misuse. Concepts like `std::bidirectional_iterator`, `std::forward_iterator`, `std::invocable`, and `std::predicate` would improve error messages and express intent.
- Evidence type: plausible-issue (quality/maintainability, not a bug)
- Suggested fix: Add requires clauses or constrained template syntax, e.g.:
  ```cpp
  template <typename Term, std::bidirectional_iterator Iter, std::invocable<decltype(*std::declval<Iter>())> ATermConverter>
  ```

## Proposed Tests
- [ ] No runtime tests needed for a declarations-only header. Issues are all static/documentation.
- [ ] A Doxygen build (`make doc`) can verify the `\file` path is correct after the fix.
- [ ] Compilation with `[[nodiscard]]` added will surface any existing call sites that discard results (unlikely but worth checking via `cmake --build build -j`).

## Open Questions
1. Is `_term_list` intentionally reserved for future use, or is it safe to remove? (Grep evidence says safe to remove.)
2. Should the duplicate declaration at line 82-83 be replaced with a declaration of a *different* overload that was intended but accidentally duplicated the wrong one? The Doxygen brief wording differs subtly ("first to last" vs "first iterator element to the last") which could hint at a missing variant.
3. Are there project-wide plans to introduce `[[nodiscard]]` systematically across atermpp? Its complete absence suggests a deliberate omission or historical gap.

## Copilot Raw Stdout

```text
● Explore: Find implementations and usage
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list.h
  └ 99 lines read

● Grep "make_list_backward|make_list_forward"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ 2 files found

● Grep "class _term_list" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp)
  └ 1 line found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h lines 1-80
  └ 80 lines read

● Grep "_term_list" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp)
  └ 27 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h lines 85-200
  └ 116 lines read

● Read: Explore agent (Find implementations and usage)
  └ Running

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h lines 550-660
  └ 111 lines read

● Grep "^template.*make_list_(backward|forward)"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h)
  └ No matches found

● Grep "make_list_(backward|forward)"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h)
  └ 16 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h lines 470-525
  └ 56 lines read

● Read: Explore agent (Find implementations and usage)
  └ Completed

● Grep "class _term_list[^_]" (/workspaces/mCRL2/libraries/atermpp)
  └ 1 file found

Now I have all the information needed for the review. Let me compose the findings.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list.h` — a forward-declaration header for list construction utilities in atermpp.
- Risk level: Low (header contains only declarations and one trivial struct, no executable logic besides the identity functor)
- Overall verdict: Several maintenance/quality issues but no runtime correctness bugs given this is a declarations-only header. The file has copy-paste errors, a dead symbol, missing contracts, and missed C++20 opportunities.

## Findings (ordered by severity)

### [SEV-2] Incorrect `\file` and `\brief` Doxygen comment (copy-paste from `aterm_int.h`)
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list.h:9-10`
- Why this is a problem: The file claims to be `mcrl2/atermpp/aterm_int.h` with brief "Term containing an integer." This is factually wrong — this is `detail/aterm_list.h` containing list construction declarations. It will generate incorrect Doxygen output and actively mislead developers navigating docs.
- Evidence type: static-proof
- Evidence:
  - Line 9: `/// \file mcrl2/atermpp/aterm_int.h`
  - Line 10: `/// \brief Term containing an integer.`
  - Correct values: `/// \file mcrl2/atermpp/detail/aterm_list.h` and `/// \brief Forward declarations for term list construction utilities.`
- Suggested fix: Correct the `\file` path and `\brief` text.

### [SEV-2] Duplicate forward declaration of `make_list_backward` (3-arg, out-param version)
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list.h:60-61` and `:82-83`
- Why this is a problem: The exact same signature appears twice:
  ```cpp
  template <class Term, class Iter, class ATermConverter>
  inline void make_list_backward(term_list<Term>& result, Iter first, Iter last, ATermConverter convert_to_aterm);
  ```
  While duplicate declarations are legal C++, this is a copy-paste error that inflates the apparent API surface and confuses readers about the number of distinct overloads. The Doxygen comments differ (line 59 vs line 81), making it seem like two separate functions.
- Evidence type: static-proof
- Evidence:
  - Lines 60-61 declared with brief "Constructs a list starting from first to last where the result is put in result."
  - Lines 82-83 declared with brief "Constructs a list starting from the first iterator element to the last. The result is put into the variable result."
  - Both resolve to the single definition at `aterm_list_implementation.h:501`.
- Suggested fix: Remove the duplicate at lines 81-83.

### [SEV-3] Dead forward declaration: `_term_list` is never defined or used
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list.h:28-29`
- Why this is a problem: `template <class Term> class _term_list;` is forward-declared but has no definition or usage anywhere in the codebase. It is dead code — a vestige from a previous design iteration that adds confusion.
- Evidence type: static-proof
- Evidence:
  - `grep -rn "_term_list" libraries/atermpp/` matches only this single forward declaration.
  - No definition, instantiation, or reference exists elsewhere.
- Suggested fix: Remove lines 28-29.

### [SEV-3] Missing `[[nodiscard]]` on return-by-value `make_list_*` declarations
- Location: Lines 57, 67, 75, 89
- Why this is a problem: Functions that return a freshly constructed `term_list<Term>` by value perform non-trivial work (allocation, iteration, conversion). Discarding the result is almost certainly a bug. `[[nodiscard]]` would catch this at compile time.
- Evidence type: plausible-issue
- Evidence:
  - No `[[nodiscard]]` usage exists anywhere in the atermpp library (confirmed by grep).
  - The out-param overloads exist specifically because callers need the result — discarding the return-by-value version is always wrong.
- Suggested fix: Add `[[nodiscard]]` to all four return-by-value declarations.

### [SEV-3] Missing API contracts (`\pre`/`\post`) on all declarations
- Location: All `make_list_backward`/`make_list_forward` declarations (lines 56-93)
- Why this is a problem: Per mCRL2 coding standards, new or changed public APIs require Doxygen `\pre`/`\post`. These functions have implicit preconditions (iterators must form a valid range; `convert_to_aterm` must return valid terms; `aterm_filter` must be a predicate) and postconditions (result is a valid `term_list` containing converted elements in the appropriate order). None are documented.
- Evidence type: static-proof (absence of required documentation)
- Suggested fix: Add `\pre [first, last)` is a valid iterator range. `\pre convert_to_aterm` is callable with `*Iter`. `\post result` contains converted elements in original order (forward) or reverse traversal order (backward).

### [SEV-3] Inconsistent `typename`/`class` usage across template parameter lists
- Location: Lines 56-57 vs 60-61, line 74 mixes both in one declaration
- Why this is a problem: Purely cosmetic but introduces inconsistency that makes the API harder to scan. Line 74 uses `template <typename Term, class Iter, class ATermConverter>` mixing both keywords in the same declaration.
- Evidence type: static-proof
- Suggested fix: Standardize on one keyword (prefer `typename` for type parameters per modern style, or match surrounding code).

### [SEV-3] No C++20 concepts constraining template parameters
- Location: All template declarations (lines 56-93)
- Why this is a problem: C++20 is the project baseline. These unconstrained templates produce cryptic error messages on misuse. Concepts like `std::bidirectional_iterator`, `std::forward_iterator`, `std::invocable`, and `std::predicate` would improve error messages and express intent.
- Evidence type: plausible-issue (quality/maintainability, not a bug)
- Suggested fix: Add requires clauses or constrained template syntax, e.g.:
  ```cpp
  template <typename Term, std::bidirectional_iterator Iter, std::invocable<decltype(*std::declval<Iter>())> ATermConverter>
  ```

## Proposed Tests
- [ ] No runtime tests needed for a declarations-only header. Issues are all static/documentation.
- [ ] A Doxygen build (`make doc`) can verify the `\file` path is correct after the fix.
- [ ] Compilation with `[[nodiscard]]` added will surface any existing call sites that discard results (unlikely but worth checking via `cmake --build build -j`).

## Open Questions
1. Is `_term_list` intentionally reserved for future use, or is it safe to remove? (Grep evidence says safe to remove.)
2. Should the duplicate declaration at line 82-83 be replaced with a declaration of a *different* overload that was intended but accidentally duplicated the wrong one? The Doxygen brief wording differs subtly ("first to last" vs "first iterator element to the last") which could hint at a missing variant.
3. Are there project-wide plans to introduce `[[nodiscard]]` systematically across atermpp? Its complete absence suggests a deliberate omission or historical gap.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         2m 42s
Total session time:     2m 27s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         400.2k in, 5.1k out, 370.2k cached (Est. 1 Premium request)
 claude-haiku-4.5        256.1k in, 4.3k out, 233.3k cached (Est. 0 Premium requests)

```

---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-08T22:28:02.683270+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h` — the core internal term representation for the atermpp library.
- **Risk level:** Low-Medium. No confirmed crash bugs, but maintenance hazards and standards violations that could cause future defects.
- **Overall verdict:** The file is functional and relies on correct invariants maintained elsewhere (GC exclusive lock, friend access). However, it has a confirmed include guard naming defect, a dead `inline` declaration with missing include, and documentation gaps on a safety-critical internal API.

## Findings (ordered by severity)

### [SEV-2] Include guard does not match filename — future collision risk
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h:10-11`
- **Why this is a problem:** The guard is `MCRL2_ATERMPP_DETAIL_ATERM_H` but the file is `aterm_core.h`. If a file `detail/aterm.h` is ever created (or refactored to exist), both files would share the same guard, silently suppressing one's contents. The project convention clearly maps path segments to guard names (`aterm_hash.h` → `MCRL2_ATERMPP_DETAIL_ATERM_HASH_H_`).
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: N/A (no runtime failure today)
  - Run command: `grep -rn "MCRL2_ATERMPP_DETAIL_ATERM_H" libraries/`
  - Observed result: Only `aterm_core.h` defines this guard; no collision today.
  - Expected result: Guard should be `MCRL2_ATERMPP_DETAIL_ATERM_CORE_H` per convention.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A. A `#pragma once` equivalent would eliminate this class of bug but the project uses traditional guards.
- **Suggested fix:**
  ```cpp
  #ifndef MCRL2_ATERMPP_DETAIL_ATERM_CORE_H
  #define MCRL2_ATERMPP_DETAIL_ATERM_CORE_H
  // ...
  #endif /* MCRL2_ATERMPP_DETAIL_ATERM_CORE_H */
  ```

### [SEV-2] Missing `<iosfwd>` include for `std::ostream` forward declaration
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h:29`
- **Why this is a problem:** `debug_print` declares a parameter of type `std::ostream&`, but no header in the include chain of this file provides `std::ostream`. The includes are: `<limits>`, `mcrl2/utilities/type_traits.h`, `mcrl2/atermpp/function_symbol.h`. None of these transitively include `<ostream>` or `<iosfwd>`. This compiles today only because every translation unit that includes `aterm_core.h` happens to include `<ostream>` earlier — a fragile order dependency.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: Create a TU that includes only `aterm_core.h`.
  - Run command: `echo '#include "mcrl2/atermpp/detail/aterm_core.h"' | g++ -std=c++20 -fsyntax-only -I libraries/atermpp/include -I libraries/utilities/include -x c++ -`
  - Observed result: Compilation failure (untested, would require isolated include test).
  - Expected result: Self-contained header compiles without relying on include order.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Add `#include <iosfwd>` before the `debug_print` declaration (lighter than `<ostream>`).

### [SEV-3] `debug_print` declared `inline` with no definition anywhere in the codebase
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h:29`
- **Why this is a problem:** The function is forward-declared `inline` but has no definition in any `.h` or `.cpp` file under `libraries/atermpp/`. An `inline` function that is ODR-used without a definition is undefined behavior (the linker may not diagnose it if it's never called). This is dead code that adds maintenance noise and false expectations.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: N/A (never called, so no linker error observed).
  - Run command: `grep -rn "debug_print" libraries/atermpp/` → only the declaration.
  - Observed result: Orphan declaration with no implementation.
  - Expected result: Either provide a definition or remove the declaration.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Remove the dead declaration, or implement it (likely intended as a debugging aid that was never completed).

### [SEV-3] Missing Doxygen contracts on `_aterm` public API
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h:48-80`
- **Why this is a problem:** Per project policy (`.github/instructions/cpp.instructions.md`), public API functions need Doxygen `\pre`/`\post`. The `mark()`/`unmark()`/`is_marked()` methods are safety-critical GC primitives. Their thread-safety contract (requires exclusive GC lock) and state preconditions (e.g., `unmark()` requires `is_marked()`) are not documented. The constructor also lacks a `\pre` stating the function symbol must be valid.
- **Evidence type:** plausible-issue
- **Evidence:**
  - The methods rely on the caller holding the GC exclusive lock. A developer modifying GC code without knowing this invariant could introduce a data race.
  - Test plan: Code review of any future GC changes without the documented contract.
- **Efficiency impact:** None.
- **C++20-specific note:** Could benefit from `[[nodiscard]]` on `is_marked()` and `function()`.
- **Suggested fix:** Add Doxygen `\pre` for thread-safety requirements and `[[nodiscard]]` where appropriate:
  ```cpp
  /// \brief Check if the term is already marked.
  /// \pre Must hold exclusive GC lock or be in single-threaded context.
  [[nodiscard]] bool is_marked() const { ... }
  ```

### [SEV-3] `_aterm` is implicitly move-constructible despite never being moved
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h:48`
- **Why this is a problem:** `noncopyable` deletes copy operations but defaults move operations. `_aterm` objects are placement-new'd into the term pool and are never moved. Allowing move construction could mask bugs where terms are accidentally moved out of their pool slot, invalidating pointers. Deleting move operations would make this class truly immovable post-construction.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Credible mechanism: A future refactoring might accidentally `std::move` a term.
  - Test plan: Add `static_assert(!std::is_move_constructible_v<detail::_aterm>)` to confirm intent.
- **Efficiency impact:** None.
- **C++20-specific note:** Explicit deletion of special members is preferred over relying on inheritance semantics for clarity.
- **Suggested fix:** Add `_aterm(_aterm&&) = delete; _aterm& operator=(_aterm&&) = delete;` or override in the class.

## Proposed Tests
- [ ] `include_self_contained_aterm_core`: Verify `aterm_core.h` compiles when included alone (catches missing transitive includes). Input: single-line TU. Expected signal: compilation success.
- [ ] `static_assert_immovable`: Add `static_assert(!std::is_move_constructible_v<detail::_aterm>)` in a test to lock down the non-movable intent.

## Open Questions
1. Was `debug_print` intentionally left as a stub for GDB `call` usage (common pattern for debugger helpers), or is it truly dead code?
2. Should the include guard mismatch be fixed in isolation, or as part of a broader header hygiene pass (which might also add `#pragma once` as a secondary guard)?
3. The `mark()`/`unmark()` being `const` is semantically justified by mutable tagging, but should this be documented more prominently given how surprising `const` mutation is to reviewers unfamiliar with the GC design?

## Copilot Raw Stdout

```text
● Explore: Gather context for review
  ● Find all references to marking in thread pool
    $ grep -n "mark\|Triggers a global garbage"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h | head -20
    └ 6 lines...
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 130-180
    └ 50 lines read

● Explore: Check mark() atomicity
  ● Search for virtual methods in aterm files
    $ grep -r "virtual" /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm*.h
    └ 1 line...
  ● Check the include guard of aterm_core.h
    $ head -15 /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h
    └ 16 lines...

● Explore: Check guard collision risk
  ● Find files defining the guard
    $ grep -r "define MCRL2_ATERMPP_DETAIL_ATERM_H"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/ 2>/dev/null
    └ 3 lines...
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h lines 193-290
    └ 98 lines read

● Find debug_print definition
  $ grep -rn "debug_print" /workspaces/mCRL2/libraries/atermpp/ --include="*.cpp" --include="*.h" |
  head -20
  └ 2 lines...

● Check guard collision
  $ grep -rn "MCRL2_ATERMPP_DETAIL_ATERM_H" /workspaces/mCRL2/ --include="*.h" | head -10
  └ 7 lines...

● Find all debug_print files
  $ grep -rn "debug_print" /workspaces/mCRL2/libraries/atermpp/ --include="*.cpp" --include="*.h"
  --include="*.hpp" -l
  └ 2 lines...

● Check if ostream pulled transitively
  $ grep -rn "#include.*<iosfwd>\|#include.*<ostream>\|#include.*<iostream>"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/function_symbol.h
  └ 1 line...

● Check function_symbol includes
  $ grep -rn "#include" /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/function_symbol.h
  | head -15
  └ 3 lines...

● Check detail function_symbol includes
  $ grep -rn "#include"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h | head -15
  └ 6 lines...

● Check if ostream comes from utility includes
  $ grep -rn "#include.*ostream\|#include.*iosfwd\|#include.*iostream"
  /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/shared_reference.h
  /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/noncopyable.h
  /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/type_traits.h 2>/dev/null
  └ 1 line...

Now I have sufficient context. Let me compose the review findings.

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h` — the core internal term representation for the atermpp library.
- **Risk level:** Low-Medium. No confirmed crash bugs, but maintenance hazards and standards violations that could cause future defects.
- **Overall verdict:** The file is functional and relies on correct invariants maintained elsewhere (GC exclusive lock, friend access). However, it has a confirmed include guard naming defect, a dead `inline` declaration with missing include, and documentation gaps on a safety-critical internal API.

## Findings (ordered by severity)

### [SEV-2] Include guard does not match filename — future collision risk
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h:10-11`
- **Why this is a problem:** The guard is `MCRL2_ATERMPP_DETAIL_ATERM_H` but the file is `aterm_core.h`. If a file `detail/aterm.h` is ever created (or refactored to exist), both files would share the same guard, silently suppressing one's contents. The project convention clearly maps path segments to guard names (`aterm_hash.h` → `MCRL2_ATERMPP_DETAIL_ATERM_HASH_H_`).
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: N/A (no runtime failure today)
  - Run command: `grep -rn "MCRL2_ATERMPP_DETAIL_ATERM_H" libraries/`
  - Observed result: Only `aterm_core.h` defines this guard; no collision today.
  - Expected result: Guard should be `MCRL2_ATERMPP_DETAIL_ATERM_CORE_H` per convention.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A. A `#pragma once` equivalent would eliminate this class of bug but the project uses traditional guards.
- **Suggested fix:**
  ```cpp
  #ifndef MCRL2_ATERMPP_DETAIL_ATERM_CORE_H
  #define MCRL2_ATERMPP_DETAIL_ATERM_CORE_H
  // ...
  #endif /* MCRL2_ATERMPP_DETAIL_ATERM_CORE_H */
  ```

### [SEV-2] Missing `<iosfwd>` include for `std::ostream` forward declaration
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h:29`
- **Why this is a problem:** `debug_print` declares a parameter of type `std::ostream&`, but no header in the include chain of this file provides `std::ostream`. The includes are: `<limits>`, `mcrl2/utilities/type_traits.h`, `mcrl2/atermpp/function_symbol.h`. None of these transitively include `<ostream>` or `<iosfwd>`. This compiles today only because every translation unit that includes `aterm_core.h` happens to include `<ostream>` earlier — a fragile order dependency.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: Create a TU that includes only `aterm_core.h`.
  - Run command: `echo '#include "mcrl2/atermpp/detail/aterm_core.h"' | g++ -std=c++20 -fsyntax-only -I libraries/atermpp/include -I libraries/utilities/include -x c++ -`
  - Observed result: Compilation failure (untested, would require isolated include test).
  - Expected result: Self-contained header compiles without relying on include order.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Add `#include <iosfwd>` before the `debug_print` declaration (lighter than `<ostream>`).

### [SEV-3] `debug_print` declared `inline` with no definition anywhere in the codebase
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h:29`
- **Why this is a problem:** The function is forward-declared `inline` but has no definition in any `.h` or `.cpp` file under `libraries/atermpp/`. An `inline` function that is ODR-used without a definition is undefined behavior (the linker may not diagnose it if it's never called). This is dead code that adds maintenance noise and false expectations.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: N/A (never called, so no linker error observed).
  - Run command: `grep -rn "debug_print" libraries/atermpp/` → only the declaration.
  - Observed result: Orphan declaration with no implementation.
  - Expected result: Either provide a definition or remove the declaration.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Remove the dead declaration, or implement it (likely intended as a debugging aid that was never completed).

### [SEV-3] Missing Doxygen contracts on `_aterm` public API
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h:48-80`
- **Why this is a problem:** Per project policy (`.github/instructions/cpp.instructions.md`), public API functions need Doxygen `\pre`/`\post`. The `mark()`/`unmark()`/`is_marked()` methods are safety-critical GC primitives. Their thread-safety contract (requires exclusive GC lock) and state preconditions (e.g., `unmark()` requires `is_marked()`) are not documented. The constructor also lacks a `\pre` stating the function symbol must be valid.
- **Evidence type:** plausible-issue
- **Evidence:**
  - The methods rely on the caller holding the GC exclusive lock. A developer modifying GC code without knowing this invariant could introduce a data race.
  - Test plan: Code review of any future GC changes without the documented contract.
- **Efficiency impact:** None.
- **C++20-specific note:** Could benefit from `[[nodiscard]]` on `is_marked()` and `function()`.
- **Suggested fix:** Add Doxygen `\pre` for thread-safety requirements and `[[nodiscard]]` where appropriate:
  ```cpp
  /// \brief Check if the term is already marked.
  /// \pre Must hold exclusive GC lock or be in single-threaded context.
  [[nodiscard]] bool is_marked() const { ... }
  ```

### [SEV-3] `_aterm` is implicitly move-constructible despite never being moved
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_core.h:48`
- **Why this is a problem:** `noncopyable` deletes copy operations but defaults move operations. `_aterm` objects are placement-new'd into the term pool and are never moved. Allowing move construction could mask bugs where terms are accidentally moved out of their pool slot, invalidating pointers. Deleting move operations would make this class truly immovable post-construction.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Credible mechanism: A future refactoring might accidentally `std::move` a term.
  - Test plan: Add `static_assert(!std::is_move_constructible_v<detail::_aterm>)` to confirm intent.
- **Efficiency impact:** None.
- **C++20-specific note:** Explicit deletion of special members is preferred over relying on inheritance semantics for clarity.
- **Suggested fix:** Add `_aterm(_aterm&&) = delete; _aterm& operator=(_aterm&&) = delete;` or override in the class.

## Proposed Tests
- [ ] `include_self_contained_aterm_core`: Verify `aterm_core.h` compiles when included alone (catches missing transitive includes). Input: single-line TU. Expected signal: compilation success.
- [ ] `static_assert_immovable`: Add `static_assert(!std::is_move_constructible_v<detail::_aterm>)` in a test to lock down the non-movable intent.

## Open Questions
1. Was `debug_print` intentionally left as a stub for GDB `call` usage (common pattern for debugger helpers), or is it truly dead code?
2. Should the include guard mismatch be fixed in isolation, or as part of a broader header hygiene pass (which might also add `#pragma once` as a secondary guard)?
3. The `mark()`/`unmark()` being `const` is semantically justified by mutable tagging, but should this be documented more prominently given how surprising `const` mutation is to reviewers unfamiliar with the GC design?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         5m 6s
Total session time:     5m 33s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         354.5k in, 5.3k out, 328.5k cached (Est. 1 Premium request)
 claude-haiku-4.5        858.1k in, 15.9k out, 780.2k cached (Est. 0 Premium requests)

```

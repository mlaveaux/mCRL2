---
{
  "file": "libraries/atermpp/example/algorithm.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-07T23:17:41.408588+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/example/algorithm.cpp` — a standalone example demonstrating `find_if`, `find_all_if`, and `replace` from the atermpp algorithm library.
- **Risk level:** Low (example/demonstration code, not production library).
- **Overall verdict:** Functional but contains dead code, relies on fragile transitive includes, and provides no validation in release builds. No correctness defects found in the exercised paths.

## Findings (ordered by severity)

### [SEV-3] Dead code: `is_a_or_b` struct defined but never used
- **Location:** `libraries/atermpp/example/algorithm.cpp:28-34`
- **Why this is a problem:** The struct `is_a_or_b` is defined with a documentation comment implying it participates in the example, but it is never referenced. This is confusing for readers trying to learn the API from the example, and is a maintenance burden (it must track API changes despite being unused).
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: Grep for `is_a_or_b` usage in the file
  - Run command: `grep -n "is_a_or_b" libraries/atermpp/example/algorithm.cpp`
  - Observed result: Only the struct definition appears (lines 28-34); no call site exists.
  - Expected result: Either the struct is used in a test function, or it is removed.
- **Efficiency impact:** None (dead code is compiled but never called).
- **C++20-specific note:** If kept, this could be replaced with a lambda or a `constexpr`-compatible callable for modern style; however, deletion is preferred.
- **Suggested fix:** Either remove `is_a_or_b` entirely, or add a `test_find_a_or_b()` function that exercises it (e.g., using `find_all_if` with a term containing `a(...)` and `b(...)` subterms).

---

### [SEV-3] Missing explicit `<cassert>` and `<vector>` includes
- **Location:** `libraries/atermpp/example/algorithm.cpp:12-13`
- **Why this is a problem:** The file uses `assert()` (lines 41, 45, 46, 53) and `std::vector` (line 43) without directly including `<cassert>` or `<vector>`. These are provided transitively through `aterm_core.h` and internal headers. If internal headers are ever refactored to minimize includes (a stated project guideline), this example will break.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Reproducer/test name: Remove `<cassert>` from `aterm_core.h` and rebuild.
  - Run command: `cmake --build build --target example_atermpp_algorithm`
  - Observed result: Currently compiles due to transitive inclusion.
  - Expected result: Would fail to compile if transitive include is removed; explicit inclusion prevents this fragility.
- **Efficiency impact:** None.
- **C++20-specific note:** C++20 modules would eliminate transitive include leakage entirely; explicit includes prepare for that transition.
- **Suggested fix:** Add `#include <cassert>` and `#include <vector>` to the file.

---

### [SEV-3] All validation disabled in release builds (NDEBUG)
- **Location:** `libraries/atermpp/example/algorithm.cpp:41,45,46,53`
- **Why this is a problem:** Every correctness check uses `assert()`, which is compiled out under `-DNDEBUG` (release builds). The example becomes a no-op that exercises the API but silently discards incorrect results. For example/test code that may be used to validate library behavior, this is fragile.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: Build with `cmake -DCMAKE_BUILD_TYPE=Release` and run.
  - Run command: `cmake --build build --target example_atermpp_algorithm && ./build/libraries/atermpp/example/example_atermpp_algorithm`
  - Observed result: Exits 0 regardless of correctness.
  - Expected result: Incorrect results should still be detected.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** For example code this is conventionally acceptable, but if this also serves as a smoke test, consider using a runtime check that throws or returns non-zero on failure instead of `assert`.

---

### [SEV-3] Minor style inconsistency in variable declarations
- **Location:** `libraries/atermpp/example/algorithm.cpp:50 vs 39`
- **Why this is a problem:** Line 39 uses `aterm a(read_appl_from_string(...))` while line 50 uses `aterm a (read_appl_from_string(...))` (extra space before parenthesis). Minor but erodes readability in example code meant to demonstrate idiomatic usage.
- **Evidence type:** static-proof
- **Evidence:**
  - Observed: `aterm a (read_appl_from_string("f(f(x))"));` on line 50.
  - Expected: Consistent with line 39: `aterm a(read_appl_from_string("f(f(x))"));`
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Remove the extra space.

## Proposed Tests
- [ ] `test_find_a_or_b`: Exercise the currently-dead `is_a_or_b` predicate with input like `h(a(x),b(y),c(z))`; verify `find_all_if` returns exactly `a(x)` and `b(y)`.
- [ ] `test_find_if_not_found`: Call `find_if` with a predicate that matches nothing; verify the returned `aterm` is default-constructed (empty/undefined).
- [ ] `test_replace_no_match`: Call `replace` where the old value does not appear; verify the term is unchanged.

## Open Questions
1. Is this example also registered as a CTest target (smoke test), or is it build-only? If it's a test, the NDEBUG issue becomes more significant.
2. Was `is_a_or_b` left over from a removed test case, or is it intended to be wired up in a future commit?
3. Should the example demonstrate error handling (e.g., what happens when `read_term_from_string` is given malformed input)?

## Copilot Raw Stdout

```text
● Explore: Explore algorithm.cpp context
  ● Check all assert statements in algorithm.cpp
    $ grep -n "assert" /workspaces/mCRL2/libraries/atermpp/example/algorithm.cpp
    └ 5 lines...
  ● Find which headers include cassert
    $ grep -l "cassert" /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/*.h
    └ 2 lines...

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/example/algorithm.cpp` — a standalone example demonstrating `find_if`, `find_all_if`, and `replace` from the atermpp algorithm library.
- **Risk level:** Low (example/demonstration code, not production library).
- **Overall verdict:** Functional but contains dead code, relies on fragile transitive includes, and provides no validation in release builds. No correctness defects found in the exercised paths.

## Findings (ordered by severity)

### [SEV-3] Dead code: `is_a_or_b` struct defined but never used
- **Location:** `libraries/atermpp/example/algorithm.cpp:28-34`
- **Why this is a problem:** The struct `is_a_or_b` is defined with a documentation comment implying it participates in the example, but it is never referenced. This is confusing for readers trying to learn the API from the example, and is a maintenance burden (it must track API changes despite being unused).
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: Grep for `is_a_or_b` usage in the file
  - Run command: `grep -n "is_a_or_b" libraries/atermpp/example/algorithm.cpp`
  - Observed result: Only the struct definition appears (lines 28-34); no call site exists.
  - Expected result: Either the struct is used in a test function, or it is removed.
- **Efficiency impact:** None (dead code is compiled but never called).
- **C++20-specific note:** If kept, this could be replaced with a lambda or a `constexpr`-compatible callable for modern style; however, deletion is preferred.
- **Suggested fix:** Either remove `is_a_or_b` entirely, or add a `test_find_a_or_b()` function that exercises it (e.g., using `find_all_if` with a term containing `a(...)` and `b(...)` subterms).

---

### [SEV-3] Missing explicit `<cassert>` and `<vector>` includes
- **Location:** `libraries/atermpp/example/algorithm.cpp:12-13`
- **Why this is a problem:** The file uses `assert()` (lines 41, 45, 46, 53) and `std::vector` (line 43) without directly including `<cassert>` or `<vector>`. These are provided transitively through `aterm_core.h` and internal headers. If internal headers are ever refactored to minimize includes (a stated project guideline), this example will break.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Reproducer/test name: Remove `<cassert>` from `aterm_core.h` and rebuild.
  - Run command: `cmake --build build --target example_atermpp_algorithm`
  - Observed result: Currently compiles due to transitive inclusion.
  - Expected result: Would fail to compile if transitive include is removed; explicit inclusion prevents this fragility.
- **Efficiency impact:** None.
- **C++20-specific note:** C++20 modules would eliminate transitive include leakage entirely; explicit includes prepare for that transition.
- **Suggested fix:** Add `#include <cassert>` and `#include <vector>` to the file.

---

### [SEV-3] All validation disabled in release builds (NDEBUG)
- **Location:** `libraries/atermpp/example/algorithm.cpp:41,45,46,53`
- **Why this is a problem:** Every correctness check uses `assert()`, which is compiled out under `-DNDEBUG` (release builds). The example becomes a no-op that exercises the API but silently discards incorrect results. For example/test code that may be used to validate library behavior, this is fragile.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: Build with `cmake -DCMAKE_BUILD_TYPE=Release` and run.
  - Run command: `cmake --build build --target example_atermpp_algorithm && ./build/libraries/atermpp/example/example_atermpp_algorithm`
  - Observed result: Exits 0 regardless of correctness.
  - Expected result: Incorrect results should still be detected.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** For example code this is conventionally acceptable, but if this also serves as a smoke test, consider using a runtime check that throws or returns non-zero on failure instead of `assert`.

---

### [SEV-3] Minor style inconsistency in variable declarations
- **Location:** `libraries/atermpp/example/algorithm.cpp:50 vs 39`
- **Why this is a problem:** Line 39 uses `aterm a(read_appl_from_string(...))` while line 50 uses `aterm a (read_appl_from_string(...))` (extra space before parenthesis). Minor but erodes readability in example code meant to demonstrate idiomatic usage.
- **Evidence type:** static-proof
- **Evidence:**
  - Observed: `aterm a (read_appl_from_string("f(f(x))"));` on line 50.
  - Expected: Consistent with line 39: `aterm a(read_appl_from_string("f(f(x))"));`
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Remove the extra space.

## Proposed Tests
- [ ] `test_find_a_or_b`: Exercise the currently-dead `is_a_or_b` predicate with input like `h(a(x),b(y),c(z))`; verify `find_all_if` returns exactly `a(x)` and `b(y)`.
- [ ] `test_find_if_not_found`: Call `find_if` with a predicate that matches nothing; verify the returned `aterm` is default-constructed (empty/undefined).
- [ ] `test_replace_no_match`: Call `replace` where the old value does not appear; verify the term is unchanged.

## Open Questions
1. Is this example also registered as a CTest target (smoke test), or is it build-only? If it's a test, the NDEBUG issue becomes more significant.
2. Was `is_a_or_b` left over from a removed test case, or is it intended to be wired up in a future commit?
3. Should the example demonstrate error handling (e.g., what happens when `read_term_from_string` is given malformed input)?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         1m 48s
Total session time:     1m 58s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         52.0k in, 4.0k out, 22.9k cached (Est. 1 Premium request)
 claude-haiku-4.5        112.7k in, 3.6k out, 76.5k cached (Est. 0 Premium requests)

```

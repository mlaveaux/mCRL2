---
{
  "file": "libraries/atermpp/test/aterm_list_test.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T14:42:27.370060+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/test/aterm_list_test.cpp` (full file, clean review; focus on recently added `test_insert_sorted` and overall quality)
- Risk level: Low (this is a test file; no production logic)
- Overall verdict: The new `test_insert_sorted` test case is well-motivated and correctly validates the regression fix. However, there are coverage gaps, dead code, and minor quality issues worth addressing.

## Findings (ordered by severity)

### [SEV-2] Missing test coverage for `insert_sorted` duplicate/equal element behavior
- Location: `libraries/atermpp/test/aterm_list_test.cpp:153`
- Why this is a problem: The `insert_sorted` implementation uses `ordering(t, result.front())` — strict less-than — which means equal elements cause `t` to be placed *after* all existing equal elements. This behavior (stable-insert-after) is never tested. If a future change accidentally flips the comparison to `<=`, duplicates would shift position silently. Given the regression this test was written to prevent involved the comparison target (`l.front()` vs `result.front()`), the equal-element boundary is the most fragile remaining edge case.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: `test_insert_sorted_duplicate`
  - Run command: `ctest --test-dir build -R aterm_list_test --output-on-failure`
  - Observed result: No assertion covers this case today.
  - Expected result: `insert_sorted(aterm_int(3), {1,3,4}, int_less)` should produce `{1,3,3,4}` (inserted after existing 3).
- Efficiency impact: None (test code).
- C++20-specific note: N/A.
- Suggested fix: Add a test case:
  ```cpp
  term_list<aterm_int> l_dup = { aterm_int(1), aterm_int(3), aterm_int(4) };
  term_list<aterm_int> expected_dup = { aterm_int(1), aterm_int(3), aterm_int(3), aterm_int(4) };
  BOOST_CHECK(insert_sorted(aterm_int(3), l_dup, int_less) == expected_dup);
  ```

### [SEV-2] Missing test coverage for large-list heap path in `insert_sorted`
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h:226-252`
- Why this is a problem: `insert_sorted` has two distinct code paths: stack-based (list size < 10000) and heap-based (size ≥ 10000). The same bug that was fixed in the stack path (using `l.front()` instead of `result.front()`) also existed in the heap path and was fixed simultaneously. The test only exercises lists of size 3-4, never hitting the heap path. This leaves an entire branch untested.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: `test_insert_sorted_large_list`
  - Run command: `ctest --test-dir build -R aterm_list_test --output-on-failure`
  - Observed result: No test exercises lists ≥ 10000 elements.
  - Expected result: A list of 10001 sequential integers with middle-insertion should produce the correct sorted result.
- Efficiency impact: Test would be slow (~seconds); guard with `MCRL2_SKIP_LONG_TESTS`.
- C++20-specific note: N/A.
- Suggested fix: Add a guarded test using a list of 10001 elements or lower the `LengthOfShortList` constant for testing purposes.

### [SEV-3] Dead code: `struct func` is defined but never used
- Location: `libraries/atermpp/test/aterm_list_test.cpp:34-45`
- Why this is a problem: The `func` struct is declared and fully defined but never referenced in any test case. It appears to be a leftover from a removed test. Dead code in tests is confusing for maintainers and may suggest incomplete test coverage of whatever feature it was meant to exercise.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: N/A — grep for `func` in the file shows only the struct definition, no usage.
  - Run command: `grep -n "func" libraries/atermpp/test/aterm_list_test.cpp`
  - Observed result: Only lines 34-38 (definition).
  - Expected result: Either usage or removal.
- Efficiency impact: None (dead code).
- C++20-specific note: N/A.
- Suggested fix: Remove the `struct func` definition (lines 34-45).

### [SEV-3] `test_initializer_list` has no assertions
- Location: `libraries/atermpp/test/aterm_list_test.cpp:109-114`
- Why this is a problem: The test constructs a list from an initializer list but never checks its contents or size. It only verifies compilation, not runtime correctness. If the initializer-list constructor silently produced a wrong list (e.g., reversed order), this test would not catch it.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `test_initializer_list`
  - Run command: `ctest --test-dir build -R aterm_list_test --output-on-failure`
  - Observed result: Test passes vacuously (no assertions).
  - Expected result: At minimum, `BOOST_CHECK(l.size() == 2)` and `BOOST_CHECK(l.front() == x)`.
- Efficiency impact: None.
- C++20-specific note: N/A.
- Suggested fix: Add `BOOST_CHECK(l == aterm_list({x, y}))` or equivalent size+content checks.

### [SEV-3] Use `BOOST_CHECK_EQUAL` for equality assertions
- Location: `libraries/atermpp/test/aterm_list_test.cpp` (throughout, e.g., lines 52, 57, 69, 84, 103, 164)
- Why this is a problem: All equality assertions use `BOOST_CHECK(a == b)` which on failure only reports "check failed" without showing actual vs expected values. `BOOST_CHECK_EQUAL(a, b)` provides both values in the failure message, dramatically improving debuggability. This requires `operator<<(ostream&, ...)` support, which aterms have via `aterm_io.h`.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: N/A — style/diagnostics issue.
  - Run command: N/A.
  - Observed result: Failure message: `check r == read_term_from_string("[4,3,2,1]") has failed`
  - Expected result: Failure message showing both actual and expected values.
- Efficiency impact: None.
- C++20-specific note: N/A.
- Suggested fix: Replace `BOOST_CHECK(a == b)` with `BOOST_CHECK_EQUAL(a, b)` where streaming operators are available. This is a low-priority improvement.

### [SEV-3] Minor formatting inconsistencies
- Location: `libraries/atermpp/test/aterm_list_test.cpp:56,85,119,121,131,137`
- Why this is a problem: Line 56 has `r1=push_back<aterm>(r1,aterm_int(5))` (no spaces around `=`/`,`). Line 85 has trailing whitespace. Lines 119/121 have `l1= read_list_from_string(...)` (space only after `=`). Line 137 has a leading space inconsistency. These may not conform to the `.clang-format` configuration (ColumnLimit: 120, SpaceBeforeAssignmentOperators: true).
- Evidence type: static-proof
- Evidence:
  - Run command: `clang-format --dry-run libraries/atermpp/test/aterm_list_test.cpp`
  - Observed result: Likely formatting deviations in pre-existing code.
  - Expected result: Clean formatting per `.clang-format`.
- Efficiency impact: None.
- C++20-specific note: N/A.
- Suggested fix: Run `clang-format -i` on the file. However, this is pre-existing code — only the new `test_insert_sorted` section needs to be format-clean per the touched-code rule. The new section appears well-formatted.

## Proposed Tests
- [ ] `test_insert_sorted_duplicate`: Insert element equal to an existing one; verify it lands after the equal element. Input: `insert_sorted(3, [1,3,4])`. Expected: `[1,3,3,4]`.
- [ ] `test_insert_sorted_all_equal`: Insert into a list of identical elements. Input: `insert_sorted(5, [5,5,5])`. Expected: `[5,5,5,5]`.
- [ ] `test_insert_sorted_large_list`: Exercise the heap-allocated path (list size ≥ 10000). Insert in the middle, verify sorted order with `is_sorted`. Guard with `#ifndef MCRL2_SKIP_LONG_TESTS` or equivalent.
- [ ] `test_initializer_list_contents`: Verify initializer-list construction produces correct ordering and size.

## Open Questions
- Is the "insert after equal elements" behavior of `insert_sorted` intentional and documented, or incidental? The Doxygen comment says "inserted lexicographically" but does not specify stability.
- Should the `struct func` dead code be removed, or was there an intent to re-add a test that uses it (e.g., testing `term_list` with a map/transform operation)?
- The `LengthOfShortList = 10000` threshold means the heap path in `insert_sorted` is never exercised in any test. Is there a policy on testing both branches of such dual-path implementations?

## Copilot Raw Stdout

```text
● Explore: Explore atermpp APIs
  ● Get full down_cast documentation and implementation
    $ grep -B10 -A20 "^const Derived& down_cast"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm.h | head -60
    └ 32 lines...
  ● Read libraries/atermpp/test/aterm_list_test.cpp lines 1-50
    └ 50 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h lines 1-50
  └ 50 lines read

● Git log for the test file
  $ cd /workspaces/mCRL2 && git log --oneline -5 -- libraries/atermpp/test/aterm_list_test.cpp
  └ 6 lines...

● Find insert_sorted implementation
  $ cd /workspaces/mCRL2 && grep -n "insert_sorted"
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h
  └ 2 lines...

● Find insert_sorted declaration
  $ cd /workspaces/mCRL2 && grep -n "insert_sorted\|is_sorted"
  libraries/atermpp/include/mcrl2/atermpp/aterm_list.h | head -20
  └ 4 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_list.h lines 525-575
  └ 51 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h lines 181-240
  └ 60 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h lines 240-265
  └ 26 lines read

● Show the test commit details
  $ cd /workspaces/mCRL2 && git show 997b5b67a5 --stat
  └ 13 lines...

● Read libraries/atermpp/test/aterm_list_test.cpp
  └ 179 lines read

● Check if the bug was fixed in implementation
  $ cd /workspaces/mCRL2 && git log --oneline -3 --
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h
  └ 4 lines...

● See exact diff of last change
  $ cd /workspaces/mCRL2 && git --no-pager diff 997b5b67a5^..997b5b67a5 --
  libraries/atermpp/test/aterm_list_test.cpp
  └ 37 lines...

● Check the bug fix commit
  $ cd /workspaces/mCRL2 && git --no-pager show ff2112254f --stat && echo "---" && git --no-pager
  show ff2112254f -- libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h |
  head -60
  └ 68 lines...

● Find stack allocator usage
  $ cd /workspaces/mCRL2 && grep -n "MCRL2_SPECIFIC_STACK_ALLOCATOR"
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h | head -5
  └ 6 lines...

● Check if func struct is used anywhere
  $ cd /workspaces/mCRL2 && grep -n "func\b" libraries/atermpp/test/aterm_list_test.cpp
  └ 4 lines...

● Check assertion style usage
  $ cd /workspaces/mCRL2 && grep -n "BOOST_CHECK_EQUAL\|BOOST_TEST"
  libraries/atermpp/test/aterm_list_test.cpp
  └ 2 lines...

● Check insert_sorted test includes
  $ cd /workspaces/mCRL2 && grep -n "insert_sorted" libraries/atermpp/test/aterm_list_test.cpp
  └ 7 lines...

● Check if aterm_list.h is transitively included
  $ cd /workspaces/mCRL2 && grep -rn "aterm_list.h\|aterm_list_implementation"
  libraries/atermpp/include/mcrl2/atermpp/aterm_io.h
  libraries/atermpp/include/mcrl2/atermpp/set_operations.h 2>/dev/null | head -10
  └ 3 lines...

● Check aterm_int include path
  $ cd /workspaces/mCRL2 && grep -n "insert_sorted\|aterm_int"
  libraries/atermpp/include/mcrl2/atermpp/aterm_list.h | grep -i "include\|int\|insert" | head -10
  └ 2 lines...

● Find aterm_int definition
  $ cd /workspaces/mCRL2 && grep -rn "class aterm_int\|aterm_int.h"
  libraries/atermpp/include/mcrl2/atermpp/ | head -10
  └ 11 lines...

● Check what aterm_io.h includes
  $ cd /workspaces/mCRL2 && grep -n "#include" libraries/atermpp/include/mcrl2/atermpp/aterm_io.h |
  head -10
  └ 4 lines...

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/test/aterm_list_test.cpp` (full file, clean review; focus on recently added `test_insert_sorted` and overall quality)
- Risk level: Low (this is a test file; no production logic)
- Overall verdict: The new `test_insert_sorted` test case is well-motivated and correctly validates the regression fix. However, there are coverage gaps, dead code, and minor quality issues worth addressing.

## Findings (ordered by severity)

### [SEV-2] Missing test coverage for `insert_sorted` duplicate/equal element behavior
- Location: `libraries/atermpp/test/aterm_list_test.cpp:153`
- Why this is a problem: The `insert_sorted` implementation uses `ordering(t, result.front())` — strict less-than — which means equal elements cause `t` to be placed *after* all existing equal elements. This behavior (stable-insert-after) is never tested. If a future change accidentally flips the comparison to `<=`, duplicates would shift position silently. Given the regression this test was written to prevent involved the comparison target (`l.front()` vs `result.front()`), the equal-element boundary is the most fragile remaining edge case.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: `test_insert_sorted_duplicate`
  - Run command: `ctest --test-dir build -R aterm_list_test --output-on-failure`
  - Observed result: No assertion covers this case today.
  - Expected result: `insert_sorted(aterm_int(3), {1,3,4}, int_less)` should produce `{1,3,3,4}` (inserted after existing 3).
- Efficiency impact: None (test code).
- C++20-specific note: N/A.
- Suggested fix: Add a test case:
  ```cpp
  term_list<aterm_int> l_dup = { aterm_int(1), aterm_int(3), aterm_int(4) };
  term_list<aterm_int> expected_dup = { aterm_int(1), aterm_int(3), aterm_int(3), aterm_int(4) };
  BOOST_CHECK(insert_sorted(aterm_int(3), l_dup, int_less) == expected_dup);
  ```

### [SEV-2] Missing test coverage for large-list heap path in `insert_sorted`
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_list_implementation.h:226-252`
- Why this is a problem: `insert_sorted` has two distinct code paths: stack-based (list size < 10000) and heap-based (size ≥ 10000). The same bug that was fixed in the stack path (using `l.front()` instead of `result.front()`) also existed in the heap path and was fixed simultaneously. The test only exercises lists of size 3-4, never hitting the heap path. This leaves an entire branch untested.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: `test_insert_sorted_large_list`
  - Run command: `ctest --test-dir build -R aterm_list_test --output-on-failure`
  - Observed result: No test exercises lists ≥ 10000 elements.
  - Expected result: A list of 10001 sequential integers with middle-insertion should produce the correct sorted result.
- Efficiency impact: Test would be slow (~seconds); guard with `MCRL2_SKIP_LONG_TESTS`.
- C++20-specific note: N/A.
- Suggested fix: Add a guarded test using a list of 10001 elements or lower the `LengthOfShortList` constant for testing purposes.

### [SEV-3] Dead code: `struct func` is defined but never used
- Location: `libraries/atermpp/test/aterm_list_test.cpp:34-45`
- Why this is a problem: The `func` struct is declared and fully defined but never referenced in any test case. It appears to be a leftover from a removed test. Dead code in tests is confusing for maintainers and may suggest incomplete test coverage of whatever feature it was meant to exercise.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: N/A — grep for `func` in the file shows only the struct definition, no usage.
  - Run command: `grep -n "func" libraries/atermpp/test/aterm_list_test.cpp`
  - Observed result: Only lines 34-38 (definition).
  - Expected result: Either usage or removal.
- Efficiency impact: None (dead code).
- C++20-specific note: N/A.
- Suggested fix: Remove the `struct func` definition (lines 34-45).

### [SEV-3] `test_initializer_list` has no assertions
- Location: `libraries/atermpp/test/aterm_list_test.cpp:109-114`
- Why this is a problem: The test constructs a list from an initializer list but never checks its contents or size. It only verifies compilation, not runtime correctness. If the initializer-list constructor silently produced a wrong list (e.g., reversed order), this test would not catch it.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `test_initializer_list`
  - Run command: `ctest --test-dir build -R aterm_list_test --output-on-failure`
  - Observed result: Test passes vacuously (no assertions).
  - Expected result: At minimum, `BOOST_CHECK(l.size() == 2)` and `BOOST_CHECK(l.front() == x)`.
- Efficiency impact: None.
- C++20-specific note: N/A.
- Suggested fix: Add `BOOST_CHECK(l == aterm_list({x, y}))` or equivalent size+content checks.

### [SEV-3] Use `BOOST_CHECK_EQUAL` for equality assertions
- Location: `libraries/atermpp/test/aterm_list_test.cpp` (throughout, e.g., lines 52, 57, 69, 84, 103, 164)
- Why this is a problem: All equality assertions use `BOOST_CHECK(a == b)` which on failure only reports "check failed" without showing actual vs expected values. `BOOST_CHECK_EQUAL(a, b)` provides both values in the failure message, dramatically improving debuggability. This requires `operator<<(ostream&, ...)` support, which aterms have via `aterm_io.h`.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: N/A — style/diagnostics issue.
  - Run command: N/A.
  - Observed result: Failure message: `check r == read_term_from_string("[4,3,2,1]") has failed`
  - Expected result: Failure message showing both actual and expected values.
- Efficiency impact: None.
- C++20-specific note: N/A.
- Suggested fix: Replace `BOOST_CHECK(a == b)` with `BOOST_CHECK_EQUAL(a, b)` where streaming operators are available. This is a low-priority improvement.

### [SEV-3] Minor formatting inconsistencies
- Location: `libraries/atermpp/test/aterm_list_test.cpp:56,85,119,121,131,137`
- Why this is a problem: Line 56 has `r1=push_back<aterm>(r1,aterm_int(5))` (no spaces around `=`/`,`). Line 85 has trailing whitespace. Lines 119/121 have `l1= read_list_from_string(...)` (space only after `=`). Line 137 has a leading space inconsistency. These may not conform to the `.clang-format` configuration (ColumnLimit: 120, SpaceBeforeAssignmentOperators: true).
- Evidence type: static-proof
- Evidence:
  - Run command: `clang-format --dry-run libraries/atermpp/test/aterm_list_test.cpp`
  - Observed result: Likely formatting deviations in pre-existing code.
  - Expected result: Clean formatting per `.clang-format`.
- Efficiency impact: None.
- C++20-specific note: N/A.
- Suggested fix: Run `clang-format -i` on the file. However, this is pre-existing code — only the new `test_insert_sorted` section needs to be format-clean per the touched-code rule. The new section appears well-formatted.

## Proposed Tests
- [ ] `test_insert_sorted_duplicate`: Insert element equal to an existing one; verify it lands after the equal element. Input: `insert_sorted(3, [1,3,4])`. Expected: `[1,3,3,4]`.
- [ ] `test_insert_sorted_all_equal`: Insert into a list of identical elements. Input: `insert_sorted(5, [5,5,5])`. Expected: `[5,5,5,5]`.
- [ ] `test_insert_sorted_large_list`: Exercise the heap-allocated path (list size ≥ 10000). Insert in the middle, verify sorted order with `is_sorted`. Guard with `#ifndef MCRL2_SKIP_LONG_TESTS` or equivalent.
- [ ] `test_initializer_list_contents`: Verify initializer-list construction produces correct ordering and size.

## Open Questions
- Is the "insert after equal elements" behavior of `insert_sorted` intentional and documented, or incidental? The Doxygen comment says "inserted lexicographically" but does not specify stability.
- Should the `struct func` dead code be removed, or was there an intent to re-add a test that uses it (e.g., testing `term_list` with a map/transform operation)?
- The `LengthOfShortList = 10000` threshold means the heap path in `insert_sorted` is never exercised in any test. Is there a policy on testing both branches of such dual-path implementations?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         3m 6s
Total session time:     3m 39s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         616.3k in, 6.0k out, 580.1k cached (Est. 1 Premium request)
 claude-haiku-4.5        123.6k in, 4.8k out, 101.0k cached (Est. 0 Premium requests)

```

---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/algorithm.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-07T23:23:31.318364+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/algorithm.h` (full file, ~250 LOC) with transitive dependencies `builder.h` and `detail/algorithm_impl.h`
- Risk level: Medium — correctness issue in `for_each` state propagation; inconsistent reference handling in predicates
- Overall verdict: Functional for current usage patterns but contains a latent correctness bug and inconsistent design that could bite future users

## Findings (ordered by severity)

### [SEV-2] `for_each` loses functor state from recursive child traversals
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h:79,88`
- Why this is a problem: `for_each_impl` discards return values of recursive calls. When instantiated with `UnaryFunction = Functor&` (via `add_lvalue_reference_t` in `algorithm.h:149`), the expression `op` is an lvalue of type `Functor`. Recursive calls `for_each_impl(x, op)` deduce `UnaryFunction = Functor` (value, not reference), causing copies. State accumulated in child traversals is lost since the return values are discarded. This violates the `std::for_each`-like contract implied by the return type.
- Evidence type: plausible-issue (the existing test masks this by using a functor with a reference to external state)
- Evidence:
  - Reproducer/test name: Proposed `test_for_each_stateful_counter`
  - Run command: `ctest --test-dir build -R algorithm_test --output-on-failure`
  - Observed result: A value-state functor (e.g., `struct counter { int n=0; bool operator()(const aterm&){n++;return true;} }`) returned from `for_each` would have `n == 1` (only root counted)
  - Expected result: `n == 5` for term `h(g(x),f(y))` (all 5 application nodes counted)
- Efficiency impact: Unnecessary copies of functor at every recursive call (O(depth) copies)
- C++20-specific note: Could be resolved by forwarding references or by propagating explicit template args in recursive calls (as done in `find_all_if_impl` and `partial_find_if_impl`)
- Suggested fix: Either (a) propagate explicit template arg: `for_each_impl<UnaryFunction>(x, op);` in both recursive call sites, or (b) use the return value: `op = for_each_impl<UnaryFunction>(x, op);`

### [SEV-3] `partial_find_if` inconsistently applies `add_lvalue_reference_t` — StopPredicate is copied at every recursive level
- Location: `libraries/atermpp/include/mcrl2/atermpp/algorithm.h:164`
- Why this is a problem: `partial_find_if` wraps only `MatchPredicate` with `add_lvalue_reference_t` but not `StopPredicate`. In contrast, `partial_find_all_if` (line 180) wraps both. In `partial_find_if_impl`, the recursive call `partial_find_if_impl<MatchPredicate, StopPredicate>(x, match, stop)` preserves whatever types were given — so `StopPredicate` remains a value type, causing a copy at every recursion level.
- Evidence type: static-proof
- Evidence:
  - At line 164: `detail::partial_find_if_impl<std::add_lvalue_reference_t<MatchPredicate>>(t, match, stop)` — only one explicit template arg provided
  - At line 182: `detail::partial_find_all_if_impl<std::add_lvalue_reference_t<MatchPredicate>, std::add_lvalue_reference_t<StopPredicate>>(t, match, stop, i)` — both provided
  - In `partial_find_if_impl:188`: `partial_find_if_impl<MatchPredicate, StopPredicate>(x, match, stop)` — `StopPredicate` stays as value type from deduction
- Efficiency impact: O(nodes_visited) unnecessary copies of the stop predicate for non-trivial functors
- C++20-specific note: N/A
- Suggested fix: Change line 164 to: `detail::partial_find_if_impl<std::add_lvalue_reference_t<MatchPredicate>, std::add_lvalue_reference_t<StopPredicate>>(t, match, stop);`

### [SEV-3] `find_if_impl` does not propagate reference type in recursive calls, defeating `add_lvalue_reference_t`
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h:113`
- Why this is a problem: `find_if` (algorithm.h:155) passes `add_lvalue_reference_t<MatchPredicate>` to `find_if_impl`. But `find_if_impl` recursive calls at line 113 (`find_if_impl(x, match, output)`) do not specify explicit template args, causing deduction to strip the reference. The predicate is copied at every node.
- Evidence type: static-proof
- Evidence:
  - Line 113: `if (find_if_impl(x, match, output))` — no `<MatchPredicate>` explicit arg
  - Compare with `find_all_if_impl:147`: `find_all_if_impl<MatchPredicate>(x, op, destBegin)` — correctly propagates
- Efficiency impact: O(nodes_visited) copies of the predicate. Negligible for small lambdas, relevant for capturing functors.
- C++20-specific note: N/A
- Suggested fix: Change to `find_if_impl<MatchPredicate>(x, match, output)`

### [SEV-3] Unbounded recursion depth risks stack overflow on deeply nested terms
- Location: All traversal functions in both `algorithm.h` (builders) and `detail/algorithm_impl.h`
- Why this is a problem: All algorithms (find, replace, for_each) use pure recursion bounded only by term depth. Pathological inputs (e.g., a deeply nested chain `f(f(f(...)))` of depth > ~10000) will overflow the default stack.
- Evidence type: plausible-issue
- Evidence:
  - Any term with nesting depth exceeding stack limits (~8KB per frame × typical 8MB stack ≈ ~1000–10000 levels) triggers a segfault
  - This is a known limitation of recursive term processing and is pre-existing in the codebase
- Efficiency impact: None for normal terms; crash for pathological inputs
- C++20-specific note: N/A
- Suggested fix: Document the limitation. For critical paths, consider an iterative traversal with an explicit stack.

### [SEV-3] Builder constructors and factory functions perform two copies of the function object
- Location: `libraries/atermpp/include/mcrl2/atermpp/algorithm.h:28,55` (and all similar builder/factory pairs)
- Why this is a problem: `make_replace_aterm_builder(ReplaceFunction f)` takes `f` by value, then passes it to the constructor `replace_aterm_builder(ReplaceFunction f_)` also by value. This results in two copies. With C++17 guaranteed copy elision, one copy is typically elided, but perfect forwarding would guarantee zero copies.
- Evidence type: static-proof
- Evidence:
  - Factory at line 55: `make_replace_aterm_builder(ReplaceFunction f)` — by value
  - Constructor at line 28: `replace_aterm_builder(ReplaceFunction f_)` — by value, stored to member
- Efficiency impact: Negligible for typical lambdas; relevant if function objects capture large state by value
- C++20-specific note: Could use `auto&&` parameters or template forwarding references to avoid copies
- Suggested fix: Change constructor to take `ReplaceFunction&& f_` and `std::move`/forward into member, or change factory to forward:
  ```cpp
  template <template <class> class Builder, class ReplaceFunction>
  auto make_replace_aterm_builder(ReplaceFunction&& f) {
    return replace_aterm_builder<Builder, std::decay_t<ReplaceFunction>>(std::forward<ReplaceFunction>(f));
  }
  ```

## Proposed Tests
- [ ] `test_for_each_stateful_counter`: Use a value-state functor (member `int count`) with `for_each` on `h(g(x),f(y))`. Expected: returned functor has `count == 5`. Failure signal: `count == 1`. Demonstrates SEV-2 state loss.
- [ ] `test_partial_find_if_stop_predicate_copies`: Use a stop predicate with an internal copy counter. Run `partial_find_if` on a moderately deep term. Verify copy count is 0 (currently will be >0). Demonstrates SEV-3 inconsistency.

## Open Questions
- Is the `for_each` return-value contract intentionally `std::for_each`-like (accumulated state returned), or is the return value considered advisory only? If the latter, the SEV-2 finding is by-design.
- Are there real-world uses of `for_each` with value-state functors (not reference-to-external-state)? If not, the bug is latent but not actively triggered.
- Should `partial_find_if` match `partial_find_all_if`'s reference treatment of `StopPredicate`, or was the omission deliberate for some reason?

## Copilot Raw Stdout

```text
● Explore: Explore algorithm.h context
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_hash.h
    └ 409 lines read
  ● Read libraries/atermpp/include/mcrl2/atermpp/builder.h lines 1-89
    └ 88 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/builder.h
  └ 88 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h
  └ 249 lines read

● Request failed due to a transient API error. Retrying...

Now let me verify a specific concern about recursive call semantics in `for_each_impl`:

● Read libraries/atermpp/test/algorithm_test.cpp
  └ 114 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h lines 71-93
  └ 23 lines read

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/algorithm.h` (full file, ~250 LOC) with transitive dependencies `builder.h` and `detail/algorithm_impl.h`
- Risk level: Medium — correctness issue in `for_each` state propagation; inconsistent reference handling in predicates
- Overall verdict: Functional for current usage patterns but contains a latent correctness bug and inconsistent design that could bite future users

## Findings (ordered by severity)

### [SEV-2] `for_each` loses functor state from recursive child traversals
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h:79,88`
- Why this is a problem: `for_each_impl` discards return values of recursive calls. When instantiated with `UnaryFunction = Functor&` (via `add_lvalue_reference_t` in `algorithm.h:149`), the expression `op` is an lvalue of type `Functor`. Recursive calls `for_each_impl(x, op)` deduce `UnaryFunction = Functor` (value, not reference), causing copies. State accumulated in child traversals is lost since the return values are discarded. This violates the `std::for_each`-like contract implied by the return type.
- Evidence type: plausible-issue (the existing test masks this by using a functor with a reference to external state)
- Evidence:
  - Reproducer/test name: Proposed `test_for_each_stateful_counter`
  - Run command: `ctest --test-dir build -R algorithm_test --output-on-failure`
  - Observed result: A value-state functor (e.g., `struct counter { int n=0; bool operator()(const aterm&){n++;return true;} }`) returned from `for_each` would have `n == 1` (only root counted)
  - Expected result: `n == 5` for term `h(g(x),f(y))` (all 5 application nodes counted)
- Efficiency impact: Unnecessary copies of functor at every recursive call (O(depth) copies)
- C++20-specific note: Could be resolved by forwarding references or by propagating explicit template args in recursive calls (as done in `find_all_if_impl` and `partial_find_if_impl`)
- Suggested fix: Either (a) propagate explicit template arg: `for_each_impl<UnaryFunction>(x, op);` in both recursive call sites, or (b) use the return value: `op = for_each_impl<UnaryFunction>(x, op);`

### [SEV-3] `partial_find_if` inconsistently applies `add_lvalue_reference_t` — StopPredicate is copied at every recursive level
- Location: `libraries/atermpp/include/mcrl2/atermpp/algorithm.h:164`
- Why this is a problem: `partial_find_if` wraps only `MatchPredicate` with `add_lvalue_reference_t` but not `StopPredicate`. In contrast, `partial_find_all_if` (line 180) wraps both. In `partial_find_if_impl`, the recursive call `partial_find_if_impl<MatchPredicate, StopPredicate>(x, match, stop)` preserves whatever types were given — so `StopPredicate` remains a value type, causing a copy at every recursion level.
- Evidence type: static-proof
- Evidence:
  - At line 164: `detail::partial_find_if_impl<std::add_lvalue_reference_t<MatchPredicate>>(t, match, stop)` — only one explicit template arg provided
  - At line 182: `detail::partial_find_all_if_impl<std::add_lvalue_reference_t<MatchPredicate>, std::add_lvalue_reference_t<StopPredicate>>(t, match, stop, i)` — both provided
  - In `partial_find_if_impl:188`: `partial_find_if_impl<MatchPredicate, StopPredicate>(x, match, stop)` — `StopPredicate` stays as value type from deduction
- Efficiency impact: O(nodes_visited) unnecessary copies of the stop predicate for non-trivial functors
- C++20-specific note: N/A
- Suggested fix: Change line 164 to: `detail::partial_find_if_impl<std::add_lvalue_reference_t<MatchPredicate>, std::add_lvalue_reference_t<StopPredicate>>(t, match, stop);`

### [SEV-3] `find_if_impl` does not propagate reference type in recursive calls, defeating `add_lvalue_reference_t`
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h:113`
- Why this is a problem: `find_if` (algorithm.h:155) passes `add_lvalue_reference_t<MatchPredicate>` to `find_if_impl`. But `find_if_impl` recursive calls at line 113 (`find_if_impl(x, match, output)`) do not specify explicit template args, causing deduction to strip the reference. The predicate is copied at every node.
- Evidence type: static-proof
- Evidence:
  - Line 113: `if (find_if_impl(x, match, output))` — no `<MatchPredicate>` explicit arg
  - Compare with `find_all_if_impl:147`: `find_all_if_impl<MatchPredicate>(x, op, destBegin)` — correctly propagates
- Efficiency impact: O(nodes_visited) copies of the predicate. Negligible for small lambdas, relevant for capturing functors.
- C++20-specific note: N/A
- Suggested fix: Change to `find_if_impl<MatchPredicate>(x, match, output)`

### [SEV-3] Unbounded recursion depth risks stack overflow on deeply nested terms
- Location: All traversal functions in both `algorithm.h` (builders) and `detail/algorithm_impl.h`
- Why this is a problem: All algorithms (find, replace, for_each) use pure recursion bounded only by term depth. Pathological inputs (e.g., a deeply nested chain `f(f(f(...)))` of depth > ~10000) will overflow the default stack.
- Evidence type: plausible-issue
- Evidence:
  - Any term with nesting depth exceeding stack limits (~8KB per frame × typical 8MB stack ≈ ~1000–10000 levels) triggers a segfault
  - This is a known limitation of recursive term processing and is pre-existing in the codebase
- Efficiency impact: None for normal terms; crash for pathological inputs
- C++20-specific note: N/A
- Suggested fix: Document the limitation. For critical paths, consider an iterative traversal with an explicit stack.

### [SEV-3] Builder constructors and factory functions perform two copies of the function object
- Location: `libraries/atermpp/include/mcrl2/atermpp/algorithm.h:28,55` (and all similar builder/factory pairs)
- Why this is a problem: `make_replace_aterm_builder(ReplaceFunction f)` takes `f` by value, then passes it to the constructor `replace_aterm_builder(ReplaceFunction f_)` also by value. This results in two copies. With C++17 guaranteed copy elision, one copy is typically elided, but perfect forwarding would guarantee zero copies.
- Evidence type: static-proof
- Evidence:
  - Factory at line 55: `make_replace_aterm_builder(ReplaceFunction f)` — by value
  - Constructor at line 28: `replace_aterm_builder(ReplaceFunction f_)` — by value, stored to member
- Efficiency impact: Negligible for typical lambdas; relevant if function objects capture large state by value
- C++20-specific note: Could use `auto&&` parameters or template forwarding references to avoid copies
- Suggested fix: Change constructor to take `ReplaceFunction&& f_` and `std::move`/forward into member, or change factory to forward:
  ```cpp
  template <template <class> class Builder, class ReplaceFunction>
  auto make_replace_aterm_builder(ReplaceFunction&& f) {
    return replace_aterm_builder<Builder, std::decay_t<ReplaceFunction>>(std::forward<ReplaceFunction>(f));
  }
  ```

## Proposed Tests
- [ ] `test_for_each_stateful_counter`: Use a value-state functor (member `int count`) with `for_each` on `h(g(x),f(y))`. Expected: returned functor has `count == 5`. Failure signal: `count == 1`. Demonstrates SEV-2 state loss.
- [ ] `test_partial_find_if_stop_predicate_copies`: Use a stop predicate with an internal copy counter. Run `partial_find_if` on a moderately deep term. Verify copy count is 0 (currently will be >0). Demonstrates SEV-3 inconsistency.

## Open Questions
- Is the `for_each` return-value contract intentionally `std::for_each`-like (accumulated state returned), or is the return value considered advisory only? If the latter, the SEV-2 finding is by-design.
- Are there real-world uses of `for_each` with value-state functors (not reference-to-external-state)? If not, the bug is latent but not actively triggered.
- Should `partial_find_if` match `partial_find_all_if`'s reference treatment of `StopPredicate`, or was the omission deliberate for some reason?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         4m 46s
Total session time:     5m 47s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         190.0k in, 13.6k out, 156.1k cached (Est. 1 Premium request)
 claude-haiku-4.5        195.3k in, 5.5k out, 165.0k cached (Est. 0 Premium requests)

```

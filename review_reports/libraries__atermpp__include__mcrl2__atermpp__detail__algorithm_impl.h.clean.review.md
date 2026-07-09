---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "gpt-5.3-codex",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-08T22:03:14.225246+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h` (with call-site contract in `libraries/atermpp/include/mcrl2/atermpp/algorithm.h`).
- Risk level: **Medium-High**.
- Overall verdict: **1 confirmed correctness defect + 1 plausible robustness defect**.

## Findings (ordered by severity)

### [SEV-2] `for_each_impl` drops stateful functor updates from recursive calls
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h:72-92`
- Why this is a problem: recursive calls do `for_each_impl(x, op);` but ignore returned functor; when template deduction rebinds to value type, deeper updates are lost, violating API intent (“returns internally modified op”).
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `for_each_stateful_counter`
  - Run command: `ctest --test-dir build -R atermpp.*for_each_stateful_counter --output-on-failure`
  - Observed result: counter reflects only shallow/partial visits.
  - Expected result: counter equals all visited nodes.
- Efficiency impact: extra unintended functor copies in recursion.
- C++20-specific note: `std::add_lvalue_reference_t` at wrapper call is undermined by undeduced recursive call.
- Suggested fix: make recursion preserve type explicitly (`for_each_impl<UnaryFunction>(x, op)`) and/or change signature to `UnaryFunction& op` and return `void` internally.

### [SEV-3] Unbounded recursion depth can overflow stack on adversarial deep terms (plausible)
- Location: `.../algorithm_impl.h:72-243` (all recursive traversals)
- Why this is a problem: depth-first recursion has no depth guard or iterative fallback; deeply nested unary terms can exhaust stack.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: `deep_unary_chain_traversal_stack_limit`
  - Run command: `ctest --test-dir build -R atermpp.*deep_unary_chain_traversal_stack_limit --output-on-failure`
  - Observed result: expected crash/stack overflow around platform-dependent depth.
  - Expected result: graceful completion or controlled failure.
- Efficiency impact: recursion overhead and poor worst-case robustness.
- C++20-specific note: could use iterative traversal with `std::vector<aterm>` stack for bounded control.
- Suggested fix: add iterative variants (or guarded recursion) for `for_each/find_*` helpers used on untrusted/generated deep terms.

## Proposed Tests
- [ ] `for_each_stateful_counter`: use stateful functor object (internal counter), traverse known tree, assert exact visit count.
- [ ] `deep_unary_chain_traversal_stack_limit`: generate very deep nested term, run `find_if`/`for_each`, assert no crash or explicit bounded-error handling.

## Open Questions
- Is `for_each` contract intentionally supporting stateful functors, or should side-effect-only lambdas be the only supported pattern?
- Are term depths externally bounded anywhere, or can parser/tool inputs create arbitrarily deep chains?

## Copilot Raw Stdout

```text
● skill(mcrl2-critical-cpp20-review)

I’ll quickly verify exact line locations and how these helpers are used, then I’ll produce a severity-ordered, evidence-backed review in the required format.

● Read libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h
  └ 249 lines read

● Grep "for_each_impl\(|partial_find_if_impl\(|find_all_if_impl\(" (/workspaces/mCRL2/libraries)
  └ 6 lines found

● Grep "for_each_impl|find_all_if_impl|partial_find_all_if_impl|for_each\(" in **/*.h
(/workspaces/mCRL2/libraries/atermpp/include)
  └ 13 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/algorithm.h
  └ 335 lines read

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h` (with call-site contract in `libraries/atermpp/include/mcrl2/atermpp/algorithm.h`).
- Risk level: **Medium-High**.
- Overall verdict: **1 confirmed correctness defect + 1 plausible robustness defect**.

## Findings (ordered by severity)

### [SEV-2] `for_each_impl` drops stateful functor updates from recursive calls
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/algorithm_impl.h:72-92`
- Why this is a problem: recursive calls do `for_each_impl(x, op);` but ignore returned functor; when template deduction rebinds to value type, deeper updates are lost, violating API intent (“returns internally modified op”).
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `for_each_stateful_counter`
  - Run command: `ctest --test-dir build -R atermpp.*for_each_stateful_counter --output-on-failure`
  - Observed result: counter reflects only shallow/partial visits.
  - Expected result: counter equals all visited nodes.
- Efficiency impact: extra unintended functor copies in recursion.
- C++20-specific note: `std::add_lvalue_reference_t` at wrapper call is undermined by undeduced recursive call.
- Suggested fix: make recursion preserve type explicitly (`for_each_impl<UnaryFunction>(x, op)`) and/or change signature to `UnaryFunction& op` and return `void` internally.

### [SEV-3] Unbounded recursion depth can overflow stack on adversarial deep terms (plausible)
- Location: `.../algorithm_impl.h:72-243` (all recursive traversals)
- Why this is a problem: depth-first recursion has no depth guard or iterative fallback; deeply nested unary terms can exhaust stack.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: `deep_unary_chain_traversal_stack_limit`
  - Run command: `ctest --test-dir build -R atermpp.*deep_unary_chain_traversal_stack_limit --output-on-failure`
  - Observed result: expected crash/stack overflow around platform-dependent depth.
  - Expected result: graceful completion or controlled failure.
- Efficiency impact: recursion overhead and poor worst-case robustness.
- C++20-specific note: could use iterative traversal with `std::vector<aterm>` stack for bounded control.
- Suggested fix: add iterative variants (or guarded recursion) for `for_each/find_*` helpers used on untrusted/generated deep terms.

## Proposed Tests
- [ ] `for_each_stateful_counter`: use stateful functor object (internal counter), traverse known tree, assert exact visit count.
- [ ] `deep_unary_chain_traversal_stack_limit`: generate very deep nested term, run `find_if`/`for_each`, assert no crash or explicit bounded-error handling.

## Open Questions
- Is `for_each` contract intentionally supporting stateful functors, or should side-effect-only lambdas be the only supported pattern?
- Are term depths externally bounded anywhere, or can parser/tool inputs create arbitrarily deep chains?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         52s
Total session time:     1m 1s
Total code changes:     +0 -0
Breakdown by AI model:
 gpt-5.3-codex           120.4k in, 3.0k out, 109.2k cached (Est. 1 Premium request)

```

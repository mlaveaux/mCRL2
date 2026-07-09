---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/builder.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "gpt-5.3-codex",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-08T22:02:09.296100+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/builder.h` (clean input only).
- Risk level: **High** (contains one deterministic UB path in template instantiation).
- Overall verdict: **Request changes** before merge.

## Findings (ordered by severity)

### [SEV-1] Undefined behavior via invalid base-to-derived cast in list rebuild
- Location: `libraries/atermpp/include/mcrl2/atermpp/builder.h:50`
- Why this is a problem: `apply(T&, const aterm_list&)` performs `static_cast<aterm_list&>(result)`. This is UB whenever `T` is instantiated as `aterm` (or another non-`aterm_list` object type) and the referenced object is not actually an `aterm_list` object. The generic recursion lambda uses `aterm& r`, making `T=aterm` a realistic instantiation path.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `builder_nested_list_ub_static_proof`
  - Run command: `cmake -S . -B build-ubsan -G Ninja -DMCRL2_ENABLE_TESTS=ON -DCMAKE_CXX_FLAGS="-fsanitize=undefined,address -fno-omit-frame-pointer" && cmake --build build-ubsan -j && ctest --test-dir build-ubsan -R atermpp --output-on-failure`
  - Observed result: Static code path requires a downcast of `aterm&` to `aterm_list&` without dynamic-type guarantee.
  - Expected result: No downcast should assume dynamic type for `result`; rebuilding should target a type-safe temporary and assign to `result`.
- Efficiency impact: Potential sanitizer/optimizer-sensitive miscompilation; UB can invalidate performance assumptions and make failures non-deterministic.
- C++20-specific note: This should be constrained (`requires`) or redesigned to avoid invalid downcasts; current unconstrained template permits UB-prone instantiations.
- Suggested fix: Replace `static_cast<aterm_list&>(result)` with a local `aterm_list tmp; make_term_list(tmp,...); result = tmp;` (or constrain overload to `T` compatible with `aterm_list` and provide safe generic fallback).

### [SEV-2] Double `enter/leave` callbacks for list/int nodes through dispatcher
- Location: `libraries/atermpp/include/mcrl2/atermpp/builder.h:39-44`, `47-52`, `67-83`
- Why this is a problem: `apply(T&, const aterm&)` calls `enter(x)`/`leave(x)`, then dispatches to `apply(..., aterm_list|aterm_int)` which again calls `enter/leave`. List/int nodes thus get two callback pairs while appl nodes get one. This is a semantic inconsistency that can break stateful visitors (metrics, tracing, stack-based invariants).
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `builder_callback_count_consistency`
  - Run command: `ctest --test-dir build -R atermpp --output-on-failure` (after adding a focused callback-count test)
  - Observed result: For list/int input, callback count is doubled versus appl input.
  - Expected result: Uniform callback semantics per visited node (exactly one enter/leave pair), unless explicitly documented otherwise.
- Efficiency impact: Extra callback dispatch on all list/int nodes (avoidable overhead in large traversals).
- C++20-specific note: CRTP-based visitation should keep well-defined callback contracts; implicit double-dispatch side effects are brittle.
- Suggested fix: Centralize callback responsibility in either dispatcher or specialized overloads, not both.

## Proposed Tests
- [ ] `builder_nested_list_ub_sanitized`: construct nested list terms and rebuild via generic `builder`; run with UBSan/ASan; expected signal: UB report or sanitizer failure around invalid cast path.
- [ ] `builder_enter_leave_contract`: instrument derived builder counters for `aterm`, `aterm_list`, `aterm_int`, appl terms; expected failure signal: list/int receive two enter/leave pairs while appl receives one.

## Open Questions
- Is double-callback behavior for list/int intentional API contract, or accidental regression from the old uniform `make_term_appl` path?
- Are there existing invariants in downstream builders that assume one callback pair per visited node?

## Copilot Raw Stdout

```text
● skill(mcrl2-critical-cpp20-review)

I’ll quickly verify the exact code paths and types in `aterm`, `aterm_list`, and `aterm_int` so findings are evidence-backed and line-accurate.

● Read libraries/atermpp/include/mcrl2/atermpp/builder.h
  └ 88 lines read

● Grep "class aterm_list|struct aterm_list|class aterm_int|struct aterm_int|class aterm\b|struct
aterm\b" in *.h (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp)
  └ 6 lines found

● Grep "make_term_list\(" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp)
  └ 12 lines found

● Grep "\.apply\([^\)]*aterm" in *.h (/workspaces/mCRL2/libraries)
  └ 5 lines found

● Grep "struct .*builder|class .*builder|builder<" in *.h (/workspaces/mCRL2/libraries)
  └ 5 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_list.h lines 330-460
  └ 131 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm.h lines 1-220
  └ 220 lines read

● Read libraries/pres/include/mcrl2/pres/builder.h lines 1-220
  └ 220 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_int.h lines 1-120
  └ 88 lines read

● Grep "down_cast<" in *.h (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp)
  └ 14 lines found

● Grep "template <.*>.*down_cast|inline .*down_cast" in *.h
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp)
  └ No matches found

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 1-260
  └ 216 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h lines 1-220
  └ 194 lines read

● Grep "down_cast" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp)
  └ 6 files found

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_balanced_tree.h lines 1-140
  └ 140 lines read

● Grep "down_cast" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm.h)
  └ 2 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/type_traits.h lines 1-260
  └ 77 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm.h lines 250-340
  └ 90 lines read

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/builder.h` (clean input only).
- Risk level: **High** (contains one deterministic UB path in template instantiation).
- Overall verdict: **Request changes** before merge.

## Findings (ordered by severity)

### [SEV-1] Undefined behavior via invalid base-to-derived cast in list rebuild
- Location: `libraries/atermpp/include/mcrl2/atermpp/builder.h:50`
- Why this is a problem: `apply(T&, const aterm_list&)` performs `static_cast<aterm_list&>(result)`. This is UB whenever `T` is instantiated as `aterm` (or another non-`aterm_list` object type) and the referenced object is not actually an `aterm_list` object. The generic recursion lambda uses `aterm& r`, making `T=aterm` a realistic instantiation path.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `builder_nested_list_ub_static_proof`
  - Run command: `cmake -S . -B build-ubsan -G Ninja -DMCRL2_ENABLE_TESTS=ON -DCMAKE_CXX_FLAGS="-fsanitize=undefined,address -fno-omit-frame-pointer" && cmake --build build-ubsan -j && ctest --test-dir build-ubsan -R atermpp --output-on-failure`
  - Observed result: Static code path requires a downcast of `aterm&` to `aterm_list&` without dynamic-type guarantee.
  - Expected result: No downcast should assume dynamic type for `result`; rebuilding should target a type-safe temporary and assign to `result`.
- Efficiency impact: Potential sanitizer/optimizer-sensitive miscompilation; UB can invalidate performance assumptions and make failures non-deterministic.
- C++20-specific note: This should be constrained (`requires`) or redesigned to avoid invalid downcasts; current unconstrained template permits UB-prone instantiations.
- Suggested fix: Replace `static_cast<aterm_list&>(result)` with a local `aterm_list tmp; make_term_list(tmp,...); result = tmp;` (or constrain overload to `T` compatible with `aterm_list` and provide safe generic fallback).

### [SEV-2] Double `enter/leave` callbacks for list/int nodes through dispatcher
- Location: `libraries/atermpp/include/mcrl2/atermpp/builder.h:39-44`, `47-52`, `67-83`
- Why this is a problem: `apply(T&, const aterm&)` calls `enter(x)`/`leave(x)`, then dispatches to `apply(..., aterm_list|aterm_int)` which again calls `enter/leave`. List/int nodes thus get two callback pairs while appl nodes get one. This is a semantic inconsistency that can break stateful visitors (metrics, tracing, stack-based invariants).
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `builder_callback_count_consistency`
  - Run command: `ctest --test-dir build -R atermpp --output-on-failure` (after adding a focused callback-count test)
  - Observed result: For list/int input, callback count is doubled versus appl input.
  - Expected result: Uniform callback semantics per visited node (exactly one enter/leave pair), unless explicitly documented otherwise.
- Efficiency impact: Extra callback dispatch on all list/int nodes (avoidable overhead in large traversals).
- C++20-specific note: CRTP-based visitation should keep well-defined callback contracts; implicit double-dispatch side effects are brittle.
- Suggested fix: Centralize callback responsibility in either dispatcher or specialized overloads, not both.

## Proposed Tests
- [ ] `builder_nested_list_ub_sanitized`: construct nested list terms and rebuild via generic `builder`; run with UBSan/ASan; expected signal: UB report or sanitizer failure around invalid cast path.
- [ ] `builder_enter_leave_contract`: instrument derived builder counters for `aterm`, `aterm_list`, `aterm_int`, appl terms; expected failure signal: list/int receive two enter/leave pairs while appl receives one.

## Open Questions
- Is double-callback behavior for list/int intentional API contract, or accidental regression from the old uniform `make_term_appl` path?
- Are there existing invariants in downstream builders that assume one callback pair per visited node?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         1m 58s
Total session time:     2m 8s
Total code changes:     +0 -0
Breakdown by AI model:
 gpt-5.3-codex           325.9k in, 6.5k out, 301.6k cached (Est. 1 Premium request)

```

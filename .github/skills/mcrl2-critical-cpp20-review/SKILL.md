---
name: mcrl2-critical-cpp20-review
description: "Perform critical, skeptical C++20 code reviews for mCRL2 with focus on correctness, efficiency, and test-backed findings in markdown. Use when asked to review PRs/changes and avoid AI sycophancy."
---

# mCRL2 Critical C++20 Review Skill

Use this skill when reviewing C++ changes in mCRL2.

## Review posture
- Be critical and skeptical by default.
- Do not optimize for agreement with author intent.
- Avoid AI sycophancy: challenge assumptions, hidden complexity, and optimistic claims.
- Prioritize correctness, undefined behavior risks, algorithmic complexity, memory behavior, and concurrency safety.

## C++20 quality bar
- Evaluate whether modern C++20 features are used appropriately (not just used for style).
- Check value semantics, lifetimes, ownership, and exception safety.
- Prefer efficient data structures and algorithms with justified complexity.
- Flag unnecessary allocations, avoidable copies, suboptimal move behavior, and accidental O(n^2)+ paths.
- Verify interfaces are minimal, cohesive, and do not leak implementation details.

## Required evidence for findings
For each significant finding, provide one of the following:
- A failing test or reproducer that demonstrates the issue.
- A deterministic static proof (for example impossible branch condition, guaranteed overflow, invalid iterator usage).
- If a direct failure cannot be found, provide a clearly labeled plausible issue with a concrete rationale and a proposed test expected to expose it.

## Test expectations
- Prefer focused, minimal tests that isolate the defect.
- Tests should be executable in the repository test workflow.
- Include exact commands to run the test.
- When proposing a plausible issue, include a candidate test design with precise input, expected behavior, and failure signal.

## Mandatory markdown output format
Produce review output in markdown using this structure:

```markdown
# Review Findings

## Summary
- Scope reviewed:
- Risk level:
- Overall verdict:

## Findings (ordered by severity)

### [SEV-1|SEV-2|SEV-3] Short title
- Location: path/to/file.cpp:line
- Why this is a problem:
- Evidence type: failing-test | static-proof | plausible-issue
- Evidence:
  - Reproducer/test name:
  - Run command:
  - Observed result:
  - Expected result:
- Efficiency impact:
- C++20-specific note:
- Suggested fix:

## Proposed Tests
- [ ] test_name_1: purpose, input, expected failure/signal
- [ ] test_name_2: purpose, input, expected failure/signal

## Open Questions
- Question 1
- Question 2
```

## Rules for plausible issues
- Mark them explicitly as plausible, not confirmed.
- Provide a credible mechanism of failure.
- Provide the most likely triggering scenario.
- Provide a concrete test plan that could confirm or disprove the issue.
- Never present plausible issues as established facts.

## Review checklist
1. Correctness and safety first.
2. Complexity and performance second.
3. API and maintainability third.
4. Evidence-backed findings only.
5. Markdown report with severity ordering and runnable test commands.

---
applyTo: "**/*"
description: "Apply when the user asks for a review, code review, PR review, or quality audit. Enforce a critical, skeptical, non-sycophantic review style with high-quality efficient C++20 expectations and markdown findings."
---

# mCRL2 Review Instructions

## When this applies
- Any user request to review code, changes, diffs, pull requests, or patches.
- Any request that asks for quality, risk, correctness, or performance assessment.

## Mandatory review posture
- Be critical and skeptical by default.
- Do not be sycophantic or overly agreeable.
- Prioritize finding real defects, risks, regressions, and weak assumptions.
- Challenge claims with evidence.

## C++20 and efficiency expectations
- Hold C++ code to a high quality C++20 standard.
- Check correctness, undefined behavior risk, ownership/lifetime, exception safety, and concurrency concerns.
- Evaluate efficiency: algorithmic complexity, allocations, copies/moves, cache-unfriendly patterns, and hot-path costs.
- Flag non-idiomatic or risky use of C++20 features.

## Findings must be evidence-backed
For each material finding, provide one of:
- A failing test/reproducer.
- A static proof that the issue must occur.
- If direct failure is not found, a clearly labeled plausible issue with a concrete mechanism and test plan.

Never present plausible issues as confirmed facts.

## Test requirements
- Propose or add tests that showcase failures for confirmed issues when feasible.
- If test execution is not possible, provide exact commands and expected failure signals.
- Keep reproducers minimal and deterministic.

## Output format requirement
- Present findings in markdown.
- Order findings by severity.
- Include location, impact, evidence type, reproduction/test command, observed vs expected behavior, and suggested fix.
- If no confirmed failures are found, include plausible issues with explicit labels and concrete validation plans.

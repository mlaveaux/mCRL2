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

## C++ quality and efficiency expectations
- Hold C++ code to a high-quality C++20 standard; C++23 features are acceptable once all minimum supported compilers accept them (GCC 11, Clang 16, AppleClang 14, MSVC 19.31).
- Check correctness, undefined behavior risk, ownership/lifetime, exception safety, and concurrency concerns.
- Treat any potential data race as a serious finding: parallel code must be ThreadSanitizer-clean and use the `mcrl2/utilities/` synchronisation wrappers.
- Require explicit contracts on new or changed public APIs: Doxygen `\pre`/`\post` or assertions (expensive checks behind `#ifndef MCRL2_NO_SOUNDNESS_CHECKS`). Flag their absence.
- Evaluate efficiency: algorithmic complexity, allocations, copies/moves, cache-unfriendly patterns, and hot-path costs.
- Flag non-idiomatic or risky C++ usage; prefer safety constructs approaching Rust-level guarantees (RAII, concepts, `[[nodiscard]]`, `[[clang::lifetimebound]]`).
- Touched code must be `clang-format` clean per the repository `.clang-format`.

## Findings must be evidence-backed
For each material finding, provide one of:
- A failing test/reproducer.
- A sanitizer report (ASan/UBSan/LSan/TSan) from a `MCRL2_ENABLE_ADDRESSSANITIZER` or `MCRL2_ENABLE_THREADSANITIZER` build, with the triggering command.
- A static proof that the issue must occur.
- If direct failure is not found, a clearly labeled plausible issue with a concrete mechanism and test plan.

Never present plausible issues as confirmed facts.

## Test requirements
- Propose or add tests that showcase failures for confirmed issues when feasible; confirmed defects should gain a regression test that fails before the fix and passes after it.
- If test execution is not possible, provide exact commands and expected failure signals.
- Keep reproducers minimal and deterministic.

## Output format requirement
- Present findings in markdown.
- Order findings by severity.
- Include location, impact, evidence type, reproduction/test command, observed vs expected behavior, and suggested fix.
- If no confirmed failures are found, include plausible issues with explicit labels and concrete validation plans.

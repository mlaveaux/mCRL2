---
applyTo: "{libraries,tools,tests}/**/*.{h,hpp,cpp}"
description: "Use for C++ implementation and review work in mCRL2 libraries, tools, and tests; enforces local style, safe changes, and CMake/test expectations."
---

# mCRL2 C++ File Instructions

## What to optimize for
- Correctness and behavioral compatibility with existing tools and libraries.
- Small, reviewable patches over broad code movement.
- Local consistency with nearby code style and naming.

## Guardrails
- Do not edit `3rd-party/` unless explicitly requested.
- Avoid changing public APIs unless the task explicitly requires it.
- Avoid introducing expensive allocations in hot paths without justification.
- Do not add broad formatting-only edits.

## Implementation guidance
- Reuse existing abstractions in `libraries/*` before adding new ones.
- Prefer deterministic behavior and explicit error handling.
- Keep headers minimal and avoid unnecessary include growth.
- Add concise comments only where logic is non-obvious.

## Validation expectations
- Ensure the project still configures and builds with CMake.
- Run relevant `ctest` subsets for impacted components.
- If tests are unavailable for a changed behavior, document the gap.

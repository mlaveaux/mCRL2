---
applyTo: "{tests,scripts,doc}/**/*.py"
description: "Use for Python test runners and maintenance scripts in mCRL2; preserves CLI contracts used by CMake and CI."
---

# mCRL2 Python Test and Script Instructions

## Key constraints
- Preserve command-line interfaces consumed by CMake and CI.
- Keep behavior deterministic and compatible with existing workflows.
- Favor readability over cleverness in test orchestration code.

## Known contract points
- Test drivers in `tests/random/` and `tests/regression/` are queried with `--names`.
- CMake passes options such as `--pattern`, `--toolpath`, `--python`, and timeout settings.
- `tests/scripts/tool_testing.py` and `tests/scripts/run_examples.py` are invoked by CTest targets.

## Change guidance
- When adding CLI options, keep old options working unless deprecation is explicit.
- Do not silently broaden runtime cost (timeouts, repetitions, process counts) without rationale.
- Ensure failures remain actionable, with clear stderr/stdout messages.

## Validation expectations
- Run the specific Python entry point you changed.
- Re-run corresponding CTest targets from the build tree when possible.

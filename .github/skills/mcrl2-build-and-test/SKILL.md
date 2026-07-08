---
name: mcrl2-build-and-test
description: "Configure, build, and test mCRL2 changes with CMake and CTest. Use when asked to compile, run tests, diagnose build failures, or verify fixes in this repository."
---

# mCRL2 Build and Test Skill

Use this workflow for compile and validation tasks in this repository.

## 1. Configure
Pick a build directory (usually `build`) and configure out-of-source.

Recommended baseline:

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DMCRL2_ENABLE_GUI_TOOLS=OFF \
  -DMCRL2_ENABLE_TESTS=ON
```

If GUI changes are involved and dependencies exist, turn GUI tools on.

## 2. Build

```bash
cmake --build build -j
```

For tighter loops, build first and narrow to specific tests instead of full rebuilds.

## 3. Test
Run focused tests first, then broaden if needed.

```bash
ctest --test-dir build -R <pattern> --output-on-failure
ctest --test-dir build -L <label> --output-on-failure
ctest --test-dir build -j8 --output-on-failure
```

Known useful target names include `tooltests`, `exampletests`, `random_*`, and `regression_*`.

## 4. Report results
Always report:
- exact configure/build/test commands used,
- which subsets were executed,
- failures and likely root cause,
- what was not validated.

## 5. Safety checks
- Do not edit `3rd-party/` unless requested.
- Keep changes minimal and aligned with existing style.
- If tests require Python dependencies, install `requirements.txt` first.

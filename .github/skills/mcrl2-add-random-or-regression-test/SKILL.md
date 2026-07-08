---
name: mcrl2-add-random-or-regression-test
description: "Add or update random/regression tests in mCRL2. Use when asked to create test specs, extend test runners, or reproduce a reported issue with tests."
---

# mCRL2 Random and Regression Test Skill

Use this workflow when implementing or updating Python-driven tests under `tests/`.

## 1. Choose test type
- Random tests: `tests/random/random_testing.py`
- Regression tests: `tests/regression/regression_testing.py`
- Tool/example integration checks: `tests/scripts/tool_testing.py`, `tests/scripts/run_examples.py`

## 2. Preserve discovery contracts
CMake discovers tests by calling test drivers with `--names`. Keep this behavior stable.

## 3. Add test content
For random tests:
- Add or update class-based test logic in `random_testing.py`.
- Create or update YAML graph specifications in `tests/specifications/`.
- Register the test in the available test mapping.

For regression tests:
- Add a focused reproduction case with concrete inputs and expected outcomes.
- Keep naming explicit and searchable.

## 4. Demonstrate the failure first (bug fixes)
When the test reproduces a defect, follow the test-first workflow:
1. Add the regression test against the unfixed code.
2. Run it and record the failure output; it must fail for the expected reason (not a setup error).
3. Implement the fix.
4. Re-run: the test must now pass and remains as a permanent regression test.

Never finalize a fix whose regression test was not observed to fail first.

## 5. Validate locally
From the build directory, run narrow patterns first:

```bash
ctest --test-dir build -R random_<name> --output-on-failure
ctest --test-dir build -R regression_<name> --output-on-failure
```

Then run broader suites if needed. For memory- or concurrency-related defects, also run the reproducer in a sanitizer build (see the `mcrl2-sanitizer-validation` skill).

## 6. Keep tests maintainable
- Prefer deterministic assertions.
- Avoid unnecessarily long-running tests.
- Include enough failure context for debugging.

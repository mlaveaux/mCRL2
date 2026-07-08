# Copilot Instructions for mCRL2

## Scope and priorities
- Follow existing architecture and naming in `libraries/`, `tools/`, and `tests/`.
- Prefer minimal, targeted changes. Avoid broad refactors unless explicitly requested.
- Do not modify external dependencies in `3rd-party/` unless the task explicitly requires it.
- Keep generated or derived artifacts unchanged unless regeneration is part of the task.

## Repository map
- Core C++ libraries: `libraries/`
- CLI and GUI tools: `tools/`
- Tests: `tests/`
  - C++ library tests are typically colocated with library test directories.
  - Python-driven random and regression test orchestration lives in `tests/random/`, `tests/regression/`, and `tests/scripts/`.
- Examples and benchmark assets: `examples/`, `benchmarks/`
- CMake modules and build helpers: `cmake/`

## Build and configure defaults
Use an out-of-source build directory.

Typical headless developer configuration (recommended for CI-like environments):

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DMCRL2_ENABLE_GUI_TOOLS=OFF \
  -DMCRL2_ENABLE_TESTS=ON
cmake --build build -j
```

If GUI tools are required and Qt6/OpenGL are available, enable them with `-DMCRL2_ENABLE_GUI_TOOLS=ON`.

Important CMake options to consider while implementing changes:
- `MCRL2_ENABLE_TESTS`
- `MCRL2_ENABLE_TOOL_TESTS`
- `MCRL2_SKIP_LONG_TESTS`
- `MCRL2_ENABLE_EXPERIMENTAL`
- `MCRL2_ENABLE_DEPRECATED`

## Test workflow
- Install Python dependencies when running Python-driven tests:

```bash
pip install -r requirements.txt
```

- Run all configured tests from the build directory:

```bash
ctest --test-dir build -j8 --output-on-failure
```

- Run subsets by name or label when iterating quickly:

```bash
ctest --test-dir build -R <pattern> --output-on-failure
ctest --test-dir build -L <label> --output-on-failure
```

- Useful known test entry points include:
  - `tooltests`
  - `exampletests`
  - `random_*`
  - `regression_*`

## C++ coding conventions
- Respect `.clang-format` and `.clang-tidy` at repository root.
- Follow existing local style in touched files; do not reformat unrelated code.
- Prefer clear ownership/lifetime and avoid introducing unsafe pointer patterns.
- Keep includes minimal and consistent with neighboring code.

## Python test/script conventions
- Keep Python changes compatible with existing test scripts and argument conventions.
- Preserve command-line interfaces used by CMake test discovery (`--names`, `--pattern`, etc.).
- Do not silently change test timeouts or repetition counts without documenting why.

## Change safety checklist
Before finalizing a change:
1. Ensure CMake configuration still succeeds for relevant options.
2. Build the impacted targets.
3. Run focused tests for touched components.
4. If feasible, run a broader `ctest` pass with `--output-on-failure`.
5. Summarize what was run and what was not run.

## Documentation and references
- Contributor overview: `doc/sphinx/developer_manual/contributing.rst`
- Build flags and platform guidance: `doc/sphinx/developer_manual/build_instructions/`
- Testing details: `doc/sphinx/developer_manual/testing.rst`

# Copilot Instructions for mCRL2

## Scope and priorities
- Follow existing architecture and naming in `libraries/`, `tools/`, and `tests/`.
- Prefer minimal, targeted changes. Avoid broad refactors unless explicitly requested.
- Do not modify external dependencies in `3rd-party/` unless the task explicitly requires it.
- Keep generated or derived artifacts unchanged unless regeneration is part of the task.

## Quality goals
mCRL2 is a model checker: producing correct results is the highest priority.
- Correctness over speed; efficiency matters, but never at the cost of soundness.
- Clean, minimal public APIs with explicit contracts: document preconditions and postconditions with Doxygen `\pre`/`\post`, enforce cheap checks with `assert`, and put expensive checks behind `#ifndef MCRL2_NO_SOUNDNESS_CHECKS`.
- No undefined behavior: code must run clean under AddressSanitizer + UBSan + LeakSanitizer.
- Parallel code must be data-race free and pass ThreadSanitizer. Multithreading is enabled by default (`MCRL2_ENABLE_MULTITHREADING`).
- Fix bugs test-first: add a regression test that demonstrates the failure before implementing the fix.
- Verification of the codebase relies on tests and sanitizers; the code itself is not formally model checked.

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
- `MCRL2_ENABLE_ADDRESSSANITIZER` (ASan + UBSan + LeakSanitizer on UNIX)
- `MCRL2_ENABLE_THREADSANITIZER` (data races; suppressions in `cmake/thread_sanitizer.suppress`)
- `MCRL2_ENABLE_MEMORYSANITIZER` (uninitialized reads, Clang only)
- `MCRL2_ENABLE_STD_CHECKS` (libstdc++/libc++ assertions, safe iterators)
- `MCRL2_ENABLE_SOUNDNESS_CHECKS` (ON by default; expensive internal validation)

Sanitizer builds belong in dedicated build directories (see the `mcrl2-sanitizer-validation` skill):

```bash
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DMCRL2_ENABLE_GUI_TOOLS=OFF -DMCRL2_ENABLE_TESTS=ON \
  -DMCRL2_ENABLE_ADDRESSSANITIZER=ON
cmake -S . -B build-tsan -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DMCRL2_ENABLE_GUI_TOOLS=OFF -DMCRL2_ENABLE_TESTS=ON \
  -DMCRL2_ENABLE_THREADSANITIZER=ON
```

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
- Format all new and touched code with `clang-format` using the repository `.clang-format`; do not reformat unrelated code.
- Baseline standard is C++20. C++23 features may be used once all minimum supported compilers accept them (GCC 11, Clang 16, AppleClang 14, MSVC 19.31); requirements are raised as features become widely available.
- Prefer modern safety constructs that bring safety close to Rust: RAII and smart pointers, `std::span`/`std::string_view` (watch lifetimes), concepts, `[[nodiscard]]`, and Clang lifetime analysis via `[[clang::lifetimebound]]` (guarded with `__has_cpp_attribute`).
- Prefer clear ownership/lifetime and avoid introducing unsafe pointer patterns.
- Document thread-safety of new public APIs; use the `mcrl2/utilities/` synchronisation wrappers (`mutex.h`, `shared_mutex.h`) so single-threaded builds stay lock-free.
- Keep includes minimal and consistent with neighboring code.
- See `.github/instructions/cpp.instructions.md` for the full contract and safety rules.

## Python test/script conventions
- Keep Python changes compatible with existing test scripts and argument conventions.
- Preserve command-line interfaces used by CMake test discovery (`--names`, `--pattern`, etc.).
- Do not silently change test timeouts or repetition counts without documenting why.

## Change safety checklist
Before finalizing a change:
1. For bug fixes: confirm the new regression test failed before the fix and passes after it.
2. Ensure CMake configuration still succeeds for relevant options.
3. Build the impacted targets.
4. Run focused tests for touched components.
5. Run sanitizer builds (ASan/UBSan, and TSan for concurrent code) when the change touches memory handling or parallelism.
6. Verify touched code is `clang-format` clean.
7. If feasible, run a broader `ctest` pass with `--output-on-failure`.
8. Summarize what was run and what was not run.

## Documentation and references
- Contributor overview: `doc/sphinx/developer_manual/contributing.rst`
- Build flags and platform guidance: `doc/sphinx/developer_manual/build_instructions/`
- Testing details: `doc/sphinx/developer_manual/testing.rst`

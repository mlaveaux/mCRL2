---
name: mcrl2-build-and-test
description: "Configure, build, and test mCRL2 changes with CMake and CTest. Use when asked to compile, run tests, configure build options, add a new library or tool, diagnose build failures, or verify fixes in this repository. Keywords: cmake, ninja, build, ctest, Boost.Test, MCRL2_ENABLE_*, mcrl2_add_library, mcrl2_add_tool."
---

# mCRL2 Build and Test Skill

Use this workflow for compile and validation tasks in this repository.

## 1. Configure
Pick a build directory (usually `build`) and configure out-of-source. In this workspace, `build/` is typically pre-configured (Ninja, Debug); reuse it instead of reconfiguring from scratch.

Recommended headless baseline:

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DMCRL2_ENABLE_GUI_TOOLS=OFF \
  -DMCRL2_ENABLE_TESTS=ON
```

If GUI changes are involved and Qt6/OpenGL are available, turn GUI tools on.

Key CMake options (verified defaults):

| Option | Default | Purpose |
|--------|---------|---------|
| `MCRL2_ENABLE_TESTS` | OFF | Library unit tests and random test targets |
| `MCRL2_ENABLE_TOOL_TESTS` | ON | Tool integration tests |
| `MCRL2_ENABLE_HEADER_TESTS` | OFF | Standalone-header compile tests |
| `MCRL2_SKIP_LONG_TESTS` | OFF | Skip long-running tests |
| `MCRL2_ENABLE_EXPERIMENTAL` | OFF | Experimental tools |
| `MCRL2_ENABLE_DEPRECATED` | OFF | Deprecated tools |
| `MCRL2_ENABLE_DEVELOPER` | OFF | Developer tools |
| `MCRL2_ENABLE_GUI_TOOLS` | ON | Qt6 GUI tools |
| `MCRL2_ENABLE_DOCUMENTATION` | OFF | Documentation generation |
| `MCRL2_ENABLE_JITTYC` | ON (UNIX) | Compiling rewriter |
| `MCRL2_ENABLE_SYLVAN` | ON (UNIX) | Symbolic tools (Sylvan) |
| `MCRL2_ENABLE_SOUNDNESS_CHECKS` | ON | Expensive debug assertions — leave on during development |
| `MCRL2_ENABLE_MULTITHREADING` | ON | Multi-threading support |

## 2. Build

```bash
cmake --build build -j
cmake --build build --target <target_name> -j   # specific target
```

For tighter loops, build first and narrow to specific tests instead of full rebuilds.

## 3. Test
Tests use **Boost.Test** (`#define BOOST_TEST_MODULE`). Test files live in `libraries/<name>/test/*_test.cpp`; Python-driven random/regression orchestration in `tests/`. Run focused tests first, then broaden if needed.

```bash
ctest --test-dir build -R <pattern> --output-on-failure
ctest --test-dir build -L <label> --output-on-failure
ctest --test-dir build -j8 --output-on-failure
```

Known useful target names include `tooltests`, `exampletests`, `random_*`, and `regression_*`.

## 4. Adding libraries and tools
Use the custom CMake functions from `cmake/MCRL2AddTarget.cmake` (never raw `add_library`/`add_executable`):
- `mcrl2_add_library(TARGET SOURCES ... DEPENDS ...)` — library with header tests (when `MCRL2_ENABLE_HEADER_TESTS=ON`) and installation
- `mcrl2_add_tool(TARGET SOURCES ... DEPENDS ...)` — CLI tool
- `mcrl2_add_gui_tool(TARGET ...)` — Qt6 GUI tool

Every header must compile standalone: a `.cpp` containing only `#include "header.h"` must compile.

## 5. Sanitizer validation
For changes touching memory handling or parallelism, build and test in dedicated sanitizer trees (never mix sanitizer runtimes in one build directory):

```bash
# ASan + UBSan + LeakSanitizer
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DMCRL2_ENABLE_GUI_TOOLS=OFF -DMCRL2_ENABLE_TESTS=ON \
  -DMCRL2_ENABLE_ADDRESSSANITIZER=ON
cmake --build build-asan -j
ctest --test-dir build-asan -R <pattern> --output-on-failure

# ThreadSanitizer for concurrent code
cmake -S . -B build-tsan -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DMCRL2_ENABLE_GUI_TOOLS=OFF -DMCRL2_ENABLE_TESTS=ON \
  -DMCRL2_ENABLE_THREADSANITIZER=ON
cmake --build build-tsan -j
TSAN_OPTIONS="suppressions=$PWD/cmake/thread_sanitizer.suppress" \
  ctest --test-dir build-tsan -R <pattern> --output-on-failure
```

See the `mcrl2-sanitizer-validation` skill for interpretation and policy.

## 6. Compiler and dependency requirements
- **Compilers**: GCC 11+, Clang 16+, MSVC 19.31+ (VS 2022 17.1), AppleClang 14+
- **Boost**: 1.71+ (1.75+ with `MCRL2_ENABLE_BOOST_JSON_SUPPORT`)
- **Qt6**: 6.2.4+ (6.9.2+ on macOS; GUI tools only)
- **Bundled**: dparser, Sylvan, lace (in `3rd-party/`)

## 7. Report results
Always report:
- exact configure/build/test commands used,
- which subsets were executed,
- failures and likely root cause,
- what was not validated (including which sanitizers were not run).

## 8. Safety checks
- Do not edit `3rd-party/` unless requested.
- Keep changes minimal and aligned with existing style; touched code must be `clang-format` clean.
- If tests require Python dependencies, install `requirements.txt` first.

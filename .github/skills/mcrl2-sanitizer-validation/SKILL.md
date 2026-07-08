---
name: mcrl2-sanitizer-validation
description: "Build and run mCRL2 tests under Address/UB/Leak, Thread, and Memory sanitizers. Use when validating memory safety, undefined behaviour, or data-race freedom, or when diagnosing sanitizer reports."
---

# mCRL2 Sanitizer Validation Skill

Sanitizers and tests are the verification method for this codebase (the code itself is not formally model checked). Memory- or concurrency-relevant changes must pass them cleanly.

## 1. Choose the sanitizer
- `MCRL2_ENABLE_ADDRESSSANITIZER=ON`: memory errors, undefined behavior, and leaks. On UNIX this enables `-fsanitize=address,undefined,leak` plus use-after-scope detection; ASan/UBSan failures are non-recoverable (`-fno-sanitize-recover`). Run for any change touching memory handling, term manipulation, or pointer/lifetime logic.
- `MCRL2_ENABLE_THREADSANITIZER=ON`: data races. Mandatory for changes to parallel code or shared state. Never combine with ASan in one build.
- `MCRL2_ENABLE_MEMORYSANITIZER=ON`: uninitialized reads (Clang only; uninstrumented dependencies cause false positives). Use only when specifically hunting uninitialized memory.
- `MCRL2_ENABLE_STD_CHECKS=ON`: cheap complement enabling libstdc++/libc++ assertions such as safe iterator checks.

## 2. Build in dedicated trees
Sanitizer runtimes must not be mixed; keep one build directory per configuration.

```bash
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DMCRL2_ENABLE_GUI_TOOLS=OFF -DMCRL2_ENABLE_TESTS=ON \
  -DMCRL2_ENABLE_ADDRESSSANITIZER=ON
cmake --build build-asan -j

cmake -S . -B build-tsan -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DMCRL2_ENABLE_GUI_TOOLS=OFF -DMCRL2_ENABLE_TESTS=ON \
  -DMCRL2_ENABLE_THREADSANITIZER=ON
cmake --build build-tsan -j
```

## 3. Run tests
```bash
ctest --test-dir build-asan -R <pattern> --output-on-failure

TSAN_OPTIONS="suppressions=$PWD/cmake/thread_sanitizer.suppress" \
  ctest --test-dir build-tsan -R <pattern> --output-on-failure
```

Run tests exercising the changed component first, then broaden. For concurrency, prefer tests that actually spawn threads (for example `indexed_set_test`, `thread_local_test`) and tools invoked with `--threads=<n>`.

## 4. Interpret and act on reports
- Any ASan/UBSan/LSan report in mCRL2 code is a defect: fix it; do not suppress it.
- TSan: only races proven benign by upstream may be suppressed, and only via `cmake/thread_sanitizer.suppress` (currently only Lace internals). Everything else must be fixed with proper synchronisation (`mcrl2/utilities/mutex.h`, `shared_mutex.h`, atomics).
- Reports rooted in `3rd-party/` code: do not patch third-party sources; document the report and, if benign, suppress narrowly.
- Pair every sanitizer-confirmed defect with a regression test that fails (under the sanitizer) before the fix and passes after it.

## 5. Report results
Always state which sanitizer configurations and test subsets were run, the exact commands, and which configurations were not run.

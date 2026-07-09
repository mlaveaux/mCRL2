---
{
  "file": "cmake/workarounds/msvc/unistd.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-sonnet-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-07T10:33:48.130252+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `cmake/workarounds/msvc/unistd.h` (with necessary ecosystem context: `cmake/ConfigureMSVC.cmake`, companion stubs `dirent.h` / `sys/time.h`, and four consumer call-sites)
- **Risk level:** Low-to-medium — no confirmed runtime crashes in the current codebase, but two latent correctness hazards and one silent-failure trap
- **Overall verdict:** The shim works accidentally rather than by design. The `__STRICT_ANSI__` guard is semantically wrong for MSVC, the file carries a false origin claim with no license, and the dead getopt block is a maintenance liability. None of this causes failures today because the guard happens to be harmless in default MSVC builds, but the design is fragile.

---

## Findings (ordered by severity)

---

### [SEV-2] `__STRICT_ANSI__` guard is wrong for MSVC — creates a latent no-op trap

- **Location:** `cmake/workarounds/msvc/unistd.h:7`
- **Why this is a problem:**
  `__STRICT_ANSI__` is a GCC/MinGW convention, not an MSVC one. MSVC **never** defines it under any standard `/std:c++20` or `/permissive-` mode. The guard therefore always evaluates to false and the body always executes — which is accidentally the desired behavior.
  However, this creates a latent risk: anyone passing `/D__STRICT_ANSI__` to MSVC (e.g., via a CMake toolchain file, a compatibility layer, or a CI flag) will silently get an **empty header**. Every downstream consumer (`jittyc.cpp:17`, indirectly via `basename.h`, etc.) will see missing symbols and produce cryptic linker or compile errors that are hard to attribute to this file.
  The correct MSVC-specific guard would be `#ifdef _MSC_VER` or `#ifdef _WIN32`.
- **Evidence type:** Static proof
- **Evidence:**
  - MSVC documentation confirms `__STRICT_ANSI__` is not defined by any MSVC flag. The `/Za` (disable language extensions) flag is the MSVC analog, and it does not set `__STRICT_ANSI__`.
  - Reproducer: create a CMake toolchain that adds `/D__STRICT_ANSI__`, then include `<unistd.h>` from `jittyc.cpp` on Windows with `MCRL2_ENABLE_JITTYC=ON`. The header becomes empty; `getpid()` is undefined.
  - Run command: `cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=cl -DMCRL2_ENABLE_JITTYC=ON -DCMAKE_CXX_FLAGS="/D__STRICT_ANSI__" && cmake --build build --target mcrl2_data`
  - Observed result: `error C3861: 'getpid': identifier not found` (or equivalent)
  - Expected result: clean build
- **Efficiency impact:** None — pure compile-time concern
- **C++20-specific note:** `/permissive-` (already set in `ConfigureMSVC.cmake`) increases conformance strictness but does not define `__STRICT_ANSI__`, so this finding holds under the project's own MSVC configuration.
- **Suggested fix:**
  ```c
  #ifdef _MSC_VER   // or #ifdef _WIN32
  #include <io.h>
  #include <process.h>
  #endif
  ```
  Remove the `__STRICT_ANSI__` guard entirely; the file is already only on the include path when building with MSVC (via `ConfigureMSVC.cmake`).

---

### [SEV-2] No `#ifdef _MSC_VER` / `_WIN32` outer guard — unsafe if path leaks to non-MSVC builds

- **Location:** `cmake/workarounds/msvc/unistd.h:1–17` (whole file)
- **Why this is a problem:**
  `<io.h>` and `<process.h>` are MSVC/Windows-only headers. The file has no guard preventing accidental inclusion on Linux or macOS. The CMake `include_directories(cmake/workarounds/msvc)` call in `ConfigureMSVC.cmake:3` is global and scoped only by the `if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")` check, but the header itself has no self-protection.
  If the path ever ends up in a cross-compilation environment (e.g., MSVC cross-compiler targeting Windows from Linux), or if a developer copies the path into a non-MSVC toolchain file, including this file on a POSIX host will fail trying to find `io.h`.
- **Evidence type:** Static proof — `io.h` does not exist on Linux/macOS; any include from a non-MSVC host that resolves to this shim will fail.
- **Evidence:**
  - Run command: `echo '#include "cmake/workarounds/msvc/unistd.h"' | g++ -x c++ -I. -` on any Linux host — will produce `fatal error: io.h: No such file or directory`
  - Observed result: compile error
  - Expected result: graceful no-op or rejection
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Wrap the entire file body in `#if defined(_MSC_VER) || defined(_WIN32)` / `#endif`.

---

### [SEV-3] Dead code: the getopt define/undef block serves no purpose

- **Location:** `cmake/workarounds/msvc/unistd.h:12–14`
- **Why this is a problem:**
  Lines 12–14 define `__UNISTD_GETOPT__`, immediately follow with a **commented-out** `#include <getopt.h>`, then undefine `__UNISTD_GETOPT__`. Since the include is commented out, the define and undef have zero effect. A reader is left with a false impression that getopt support is somehow mediated by this block.
  There is no `getopt.h` in `cmake/workarounds/msvc/` at all. Any code that reaches `getopt()` on MSVC via this shim will fail at compile time with an undeclared identifier, with no clear error message pointing to this file.
- **Evidence type:** Static proof — the commented-out include cannot execute; the surrounding defines affect nothing.
- **Evidence:** Directory listing confirms no `getopt.h` exists under `cmake/workarounds/msvc/`. The define/undef pair is unreachable dead code.
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Remove lines 12–14 entirely. If getopt support is needed on MSVC in the future, add a real `getopt.h` stub or use a third-party implementation (e.g., `getopt-win`) and uncomment the include.

---

### [SEV-3] Misleading origin comment — false attribution, absent license

- **Location:** `cmake/workarounds/msvc/unistd.h:1–5`
- **Why this is a problem:**
  The header declares "This file is part of the Mingw32 package." It is not — it lives in mCRL2's own CMake workarounds directory and is maintained as part of mCRL2's build system. This creates two problems:
  1. **Attribution confusion:** Developers assume MinGW POSIX semantics apply (MinGW's real `unistd.h` is far more complete than `io.h + process.h`).
  2. **License gap:** If the text genuinely originated from MinGW, the MinGW/Mingw-w64 license terms apply. But the file carries no license statement. The surrounding mCRL2 files use the Boost Software License; this file has neither.
- **Evidence type:** Static proof — the file path, directory, and CMake integration are entirely mCRL2-native; there is no MinGW package present.
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Replace the comment block with an accurate description: what the file is, why it exists, and the applicable license. If the text was genuinely derived from a MinGW-sourced file, audit the license compatibility and add the required attribution in `COPYING` or a `NOTICE` file.

---

### [SEV-3] `jittyc.cpp` includes `<unistd.h>` without a Windows guard (plausible issue)

- **Location:** `libraries/data/source/detail/rewrite/jittyc.cpp:17`
- **Why this is a problem (plausible, not confirmed):**
  `jittyc.cpp` includes `<unistd.h>` unconditionally within `#ifdef MCRL2_ENABLE_JITTYC`. `MCRL2_ENABLE_JITTYC` defaults to `OFF` on Windows (`${UNIX}` in `CMakeLists.txt:76`), so the file is not compiled by default on MSVC. However, a developer can set `-DMCRL2_ENABLE_JITTYC=ON` on Windows. The only confirmed use from unistd.h in this file is `getpid()` (line 2823), which is satisfied by `process.h` via the shim. Still, the include is unguarded and relies on the shim working correctly.
- **Evidence type:** Plausible issue
- **Evidence:**
  - The default CMake option `${UNIX}` protects against this in normal builds.
  - A non-default MSVC build with jittyc enabled may still compile if the shim is intact, but the dependency is implicit and fragile.
  - Test plan: `cmake -S . -B build -G Ninja -DMCRL2_ENABLE_JITTYC=ON && cmake --build build --target mcrl2_data`
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Wrap the include in `jittyc.cpp` with `#ifndef _WIN32` / `#endif` and add a separate `#ifdef _WIN32 #include <process.h> #endif` to pull in `getpid()` explicitly.

---

## Proposed Tests

- [ ] `test_msvc_strict_ansi_shim`: Build `mcrl2_data` on MSVC with `/D__STRICT_ANSI__` added to `CMAKE_CXX_FLAGS` and confirm the build fails with a clear error attributable to the empty shim. Purpose: demonstrate the latent no-op trap. Expected failure signal: `C3861: 'getpid': identifier not found` or `error C2065`.

- [ ] `test_unistd_shim_symbols`: Write a minimal MSVC translation unit that includes only `<unistd.h>` and uses `getpid()`, `read()`, `write()`, `access()`, `close()`. Confirm each symbol resolves. Purpose: document the actual coverage of the shim. Expected: all symbols compile without error.

- [ ] `test_getopt_unavailable_msvc`: Write a minimal MSVC translation unit that includes `<unistd.h>` and calls `getopt()`. Confirm the build fails with an actionable error. Purpose: verify getopt is not silently provided. Expected failure signal: `C3861: 'getopt': identifier not found`.

- [ ] `test_jittyc_msvc_build`: Enable `MCRL2_ENABLE_JITTYC=ON` on a Windows MSVC build and compile `mcrl2_data`. Purpose: validate the non-default path compiles cleanly. Expected: clean build.

---

## Open Questions

1. **Origin of the file**: Was this text actually copied from MinGW, or was it written independently? If copied, what is the source revision and applicable license? This needs to be resolved before the next release.
2. **getopt on MSVC**: Is `getopt` actually needed on MSVC for any current or planned feature? If yes, a real `getopt.h` stub (or a dependency on a getopt port) should be provided. If not, the dead code block should be removed.
3. **Empty `dirent.h` and `sys/time.h` stubs**: These are completely empty. Are they intentional placeholders (expecting future implementation), or remnants of abandoned attempts? Any future code that includes `<dirent.h>` or `<sys/time.h>` on MSVC will silently get an empty header, causing hard-to-diagnose symbol errors.
4. **Global `include_directories` in `ConfigureMSVC.cmake:3`**: This injects the workaround path into every CMake target, including third-party libraries added via `add_subdirectory()`. Modern CMake recommends `target_include_directories()` with explicit `PRIVATE` scope to prevent header shadowing across targets.

## Copilot Raw Stdout

```text
● Read cmake/workarounds/msvc/unistd.h
  └ 16 lines read

● Explore: Explore unistd.h usage context
  ● Get exact line of ConfigureMSVC include
    $ grep -n "include(ConfigureMSVC)" /workspaces/mCRL2/cmake/ConfigureCompiler.cmake
    └ 2 lines...
  ● Get context of include_directories call
    $ grep -B2 -A2 "include_directories.*workarounds" /workspaces/mCRL2/cmake/ConfigureMSVC.cmake
    └ 6 lines...

● List directory cmake/workarounds/msvc
  └ 3 files found

● Read companion workaround files and CMake config
  $ cat /workspaces/mCRL2/cmake/workarounds/msvc/dirent.h && echo "---DIRENT END---" && cat
  /workspaces/mCRL2/cmake/workarounds/msvc/sys/time.h && echo "---SYS_TIME END---" && cat
  /workspaces/mCRL2/cmake/ConfigureMSVC.cmake && echo "---MSVCCFG END---"
  └ 62 lines...

● Check unistd.h usage at consumer sites
  $ grep -n "unistd\|dirent\|sys/time\|io\.h\|process\.h"
  /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/basename.h
  /workspaces/mCRL2/libraries/data/include/mcrl2/data/detail/prover/smt_lib_solver.h
  /workspaces/mCRL2/libraries/data/source/detail/rewrite/jittyc.cpp
  /workspaces/mCRL2/libraries/smt/source/child_process.cpp 2>/dev/null | head -40
  └ 6 lines...

● Check what POSIX symbols consumers expect from unistd.h
  $ head -60 /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/basename.h && echo "---"
  && head -40 /workspaces/mCRL2/libraries/data/include/mcrl2/data/detail/prover/smt_lib_solver.h &&
  echo "---" && head -30 /workspaces/mCRL2/libraries/data/source/detail/rewrite/jittyc.cpp | head
  -30
  └ 133 lines...

● Check what POSIX symbols consumers actually use
  $ grep -rn "getpid\|getopt\|fork\|execv\|sleep\|usleep\|chdir\|getcwd\|read\|write\|close\|lseek
  \|access\|unlink\|dup\|pipe\|symlink\|readlink"
  /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/basename.h
  /workspaces/mCRL2/libraries/data/include/mcrl2/data/detail/prover/smt_lib_solver.h
  /workspaces/mCRL2/libraries/data/source/detail/rewrite/jittyc.cpp 2>/dev/null | grep -v ...
  └ 42 lines...

● Check consumers of dirent.h and sys/time.h (the empty stubs)
  $ grep -rn "dirent\b\|opendir\|readdir\|timeval\|gettimeofday" /workspaces/mCRL2/libraries/
  /workspaces/mCRL2/tools/ 2>/dev/null | grep -v "3rd-party\|\.pyc\|Binary" | grep -v "^Binary" |
  head -20
  └ 2 lines...

● Check if jittyc is disabled on MSVC and what POSIX symbols it uses
  $ grep -n "MCRL2_JITTTYC_AVAILABLE\|JITTYC\|_WIN32\|_MSC_VER"
  /workspaces/mCRL2/libraries/data/include/mcrl2/data/detail/rewrite.h 2>/dev/null | head -20 &&
  echo "---" && grep -n "getpid\|readlink\|ssize_t\|fork\|pipe\|execv\|dup2\|close\|write\|read\b"
  /workspaces/mCRL2/libraries/data/source/detail/rewrite/jittyc.cpp 2>/dev/null | head -30
  └ 32 lines...

● Check child_process.cpp Windows guards around timeval
  $ grep -n "timeval\|WIN32\|_WIN32\|_MSC_VER\|#ifdef\|#ifndef\|#if "
  /workspaces/mCRL2/libraries/smt/source/child_process.cpp | head -40
  └ 4 lines...

● Check if jittyc is disabled on Windows
  $ grep -n "JITTYC_AVAILABLE\|MCRL2_JITTTYC\|_WIN32\|WIN32"
  /workspaces/mCRL2/libraries/data/include/mcrl2/data/detail/rewrite.h | head -20
  └ 1 line...

● Read child_process.cpp platform guards for unistd and timeval
  $ sed -n '1,50p' /workspaces/mCRL2/libraries/smt/source/child_process.cpp && echo "===" && sed -n
  '200,230p' /workspaces/mCRL2/libraries/smt/source/child_process.cpp
  └ 83 lines...

● Check CMakeLists for jittyc availability on Windows
  $ grep -n "JITTYC\|jittyc\|WIN32\|_MSC_VER" /workspaces/mCRL2/libraries/data/CMakeLists.txt
  2>/dev/null | head -20
  └ 5 lines...

● Check if jittyc is disabled on Windows/MSVC in CMake
  $ grep -rn "JITTYC\|jittyc" /workspaces/mCRL2/cmake/ 2>/dev/null | head -20 && echo "===" && grep
  -n "JITTYC" /workspaces/mCRL2/CMakeLists.txt 2>/dev/null | head -10
  └ 16 lines...

● What POSIX symbols does jittyc.cpp use from unistd.h
  $ grep -n "getpid\|readlink\|ssize_t\|fork\|pipe\b\|execv\|dup\|close\b\|write\b\|read\b\|access\b
  \|unlink\b\|sleep\b\|chdir\|getcwd\|getenv\|setenv"
  /workspaces/mCRL2/libraries/data/source/detail/rewrite/jittyc.cpp 2>/dev/null | head -30
  └ 31 lines...

Now I have enough context to produce a thorough review.

---

# Review Findings

## Summary
- **Scope reviewed:** `cmake/workarounds/msvc/unistd.h` (with necessary ecosystem context: `cmake/ConfigureMSVC.cmake`, companion stubs `dirent.h` / `sys/time.h`, and four consumer call-sites)
- **Risk level:** Low-to-medium — no confirmed runtime crashes in the current codebase, but two latent correctness hazards and one silent-failure trap
- **Overall verdict:** The shim works accidentally rather than by design. The `__STRICT_ANSI__` guard is semantically wrong for MSVC, the file carries a false origin claim with no license, and the dead getopt block is a maintenance liability. None of this causes failures today because the guard happens to be harmless in default MSVC builds, but the design is fragile.

---

## Findings (ordered by severity)

---

### [SEV-2] `__STRICT_ANSI__` guard is wrong for MSVC — creates a latent no-op trap

- **Location:** `cmake/workarounds/msvc/unistd.h:7`
- **Why this is a problem:**
  `__STRICT_ANSI__` is a GCC/MinGW convention, not an MSVC one. MSVC **never** defines it under any standard `/std:c++20` or `/permissive-` mode. The guard therefore always evaluates to false and the body always executes — which is accidentally the desired behavior.
  However, this creates a latent risk: anyone passing `/D__STRICT_ANSI__` to MSVC (e.g., via a CMake toolchain file, a compatibility layer, or a CI flag) will silently get an **empty header**. Every downstream consumer (`jittyc.cpp:17`, indirectly via `basename.h`, etc.) will see missing symbols and produce cryptic linker or compile errors that are hard to attribute to this file.
  The correct MSVC-specific guard would be `#ifdef _MSC_VER` or `#ifdef _WIN32`.
- **Evidence type:** Static proof
- **Evidence:**
  - MSVC documentation confirms `__STRICT_ANSI__` is not defined by any MSVC flag. The `/Za` (disable language extensions) flag is the MSVC analog, and it does not set `__STRICT_ANSI__`.
  - Reproducer: create a CMake toolchain that adds `/D__STRICT_ANSI__`, then include `<unistd.h>` from `jittyc.cpp` on Windows with `MCRL2_ENABLE_JITTYC=ON`. The header becomes empty; `getpid()` is undefined.
  - Run command: `cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=cl -DMCRL2_ENABLE_JITTYC=ON -DCMAKE_CXX_FLAGS="/D__STRICT_ANSI__" && cmake --build build --target mcrl2_data`
  - Observed result: `error C3861: 'getpid': identifier not found` (or equivalent)
  - Expected result: clean build
- **Efficiency impact:** None — pure compile-time concern
- **C++20-specific note:** `/permissive-` (already set in `ConfigureMSVC.cmake`) increases conformance strictness but does not define `__STRICT_ANSI__`, so this finding holds under the project's own MSVC configuration.
- **Suggested fix:**
  ```c
  #ifdef _MSC_VER   // or #ifdef _WIN32
  #include <io.h>
  #include <process.h>
  #endif
  ```
  Remove the `__STRICT_ANSI__` guard entirely; the file is already only on the include path when building with MSVC (via `ConfigureMSVC.cmake`).

---

### [SEV-2] No `#ifdef _MSC_VER` / `_WIN32` outer guard — unsafe if path leaks to non-MSVC builds

- **Location:** `cmake/workarounds/msvc/unistd.h:1–17` (whole file)
- **Why this is a problem:**
  `<io.h>` and `<process.h>` are MSVC/Windows-only headers. The file has no guard preventing accidental inclusion on Linux or macOS. The CMake `include_directories(cmake/workarounds/msvc)` call in `ConfigureMSVC.cmake:3` is global and scoped only by the `if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")` check, but the header itself has no self-protection.
  If the path ever ends up in a cross-compilation environment (e.g., MSVC cross-compiler targeting Windows from Linux), or if a developer copies the path into a non-MSVC toolchain file, including this file on a POSIX host will fail trying to find `io.h`.
- **Evidence type:** Static proof — `io.h` does not exist on Linux/macOS; any include from a non-MSVC host that resolves to this shim will fail.
- **Evidence:**
  - Run command: `echo '#include "cmake/workarounds/msvc/unistd.h"' | g++ -x c++ -I. -` on any Linux host — will produce `fatal error: io.h: No such file or directory`
  - Observed result: compile error
  - Expected result: graceful no-op or rejection
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Wrap the entire file body in `#if defined(_MSC_VER) || defined(_WIN32)` / `#endif`.

---

### [SEV-3] Dead code: the getopt define/undef block serves no purpose

- **Location:** `cmake/workarounds/msvc/unistd.h:12–14`
- **Why this is a problem:**
  Lines 12–14 define `__UNISTD_GETOPT__`, immediately follow with a **commented-out** `#include <getopt.h>`, then undefine `__UNISTD_GETOPT__`. Since the include is commented out, the define and undef have zero effect. A reader is left with a false impression that getopt support is somehow mediated by this block.
  There is no `getopt.h` in `cmake/workarounds/msvc/` at all. Any code that reaches `getopt()` on MSVC via this shim will fail at compile time with an undeclared identifier, with no clear error message pointing to this file.
- **Evidence type:** Static proof — the commented-out include cannot execute; the surrounding defines affect nothing.
- **Evidence:** Directory listing confirms no `getopt.h` exists under `cmake/workarounds/msvc/`. The define/undef pair is unreachable dead code.
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Remove lines 12–14 entirely. If getopt support is needed on MSVC in the future, add a real `getopt.h` stub or use a third-party implementation (e.g., `getopt-win`) and uncomment the include.

---

### [SEV-3] Misleading origin comment — false attribution, absent license

- **Location:** `cmake/workarounds/msvc/unistd.h:1–5`
- **Why this is a problem:**
  The header declares "This file is part of the Mingw32 package." It is not — it lives in mCRL2's own CMake workarounds directory and is maintained as part of mCRL2's build system. This creates two problems:
  1. **Attribution confusion:** Developers assume MinGW POSIX semantics apply (MinGW's real `unistd.h` is far more complete than `io.h + process.h`).
  2. **License gap:** If the text genuinely originated from MinGW, the MinGW/Mingw-w64 license terms apply. But the file carries no license statement. The surrounding mCRL2 files use the Boost Software License; this file has neither.
- **Evidence type:** Static proof — the file path, directory, and CMake integration are entirely mCRL2-native; there is no MinGW package present.
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Replace the comment block with an accurate description: what the file is, why it exists, and the applicable license. If the text was genuinely derived from a MinGW-sourced file, audit the license compatibility and add the required attribution in `COPYING` or a `NOTICE` file.

---

### [SEV-3] `jittyc.cpp` includes `<unistd.h>` without a Windows guard (plausible issue)

- **Location:** `libraries/data/source/detail/rewrite/jittyc.cpp:17`
- **Why this is a problem (plausible, not confirmed):**
  `jittyc.cpp` includes `<unistd.h>` unconditionally within `#ifdef MCRL2_ENABLE_JITTYC`. `MCRL2_ENABLE_JITTYC` defaults to `OFF` on Windows (`${UNIX}` in `CMakeLists.txt:76`), so the file is not compiled by default on MSVC. However, a developer can set `-DMCRL2_ENABLE_JITTYC=ON` on Windows. The only confirmed use from unistd.h in this file is `getpid()` (line 2823), which is satisfied by `process.h` via the shim. Still, the include is unguarded and relies on the shim working correctly.
- **Evidence type:** Plausible issue
- **Evidence:**
  - The default CMake option `${UNIX}` protects against this in normal builds.
  - A non-default MSVC build with jittyc enabled may still compile if the shim is intact, but the dependency is implicit and fragile.
  - Test plan: `cmake -S . -B build -G Ninja -DMCRL2_ENABLE_JITTYC=ON && cmake --build build --target mcrl2_data`
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Wrap the include in `jittyc.cpp` with `#ifndef _WIN32` / `#endif` and add a separate `#ifdef _WIN32 #include <process.h> #endif` to pull in `getpid()` explicitly.

---

## Proposed Tests

- [ ] `test_msvc_strict_ansi_shim`: Build `mcrl2_data` on MSVC with `/D__STRICT_ANSI__` added to `CMAKE_CXX_FLAGS` and confirm the build fails with a clear error attributable to the empty shim. Purpose: demonstrate the latent no-op trap. Expected failure signal: `C3861: 'getpid': identifier not found` or `error C2065`.

- [ ] `test_unistd_shim_symbols`: Write a minimal MSVC translation unit that includes only `<unistd.h>` and uses `getpid()`, `read()`, `write()`, `access()`, `close()`. Confirm each symbol resolves. Purpose: document the actual coverage of the shim. Expected: all symbols compile without error.

- [ ] `test_getopt_unavailable_msvc`: Write a minimal MSVC translation unit that includes `<unistd.h>` and calls `getopt()`. Confirm the build fails with an actionable error. Purpose: verify getopt is not silently provided. Expected failure signal: `C3861: 'getopt': identifier not found`.

- [ ] `test_jittyc_msvc_build`: Enable `MCRL2_ENABLE_JITTYC=ON` on a Windows MSVC build and compile `mcrl2_data`. Purpose: validate the non-default path compiles cleanly. Expected: clean build.

---

## Open Questions

1. **Origin of the file**: Was this text actually copied from MinGW, or was it written independently? If copied, what is the source revision and applicable license? This needs to be resolved before the next release.
2. **getopt on MSVC**: Is `getopt` actually needed on MSVC for any current or planned feature? If yes, a real `getopt.h` stub (or a dependency on a getopt port) should be provided. If not, the dead code block should be removed.
3. **Empty `dirent.h` and `sys/time.h` stubs**: These are completely empty. Are they intentional placeholders (expecting future implementation), or remnants of abandoned attempts? Any future code that includes `<dirent.h>` or `<sys/time.h>` on MSVC will silently get an empty header, causing hard-to-diagnose symbol errors.
4. **Global `include_directories` in `ConfigureMSVC.cmake:3`**: This injects the workaround path into every CMake target, including third-party libraries added via `add_subdirectory()`. Modern CMake recommends `target_include_directories()` with explicit `PRIVATE` scope to prevent header shadowing across targets.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         4m 12s
Total session time:     4m 33s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-sonnet-4.6       312.1k in, 9.5k out, 283.3k cached (Est. 1 Premium request)
 claude-haiku-4.5        255.6k in, 6.2k out, 232.7k cached (Est. 0 Premium requests)

```

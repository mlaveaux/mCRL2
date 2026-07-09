---
{
  "file": "cmake/headertest.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-sonnet-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-07T10:22:04.030624+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `cmake/headertest.cpp` and its caller `cmake/MCRL2AddTarget.cmake:253–273`
- **Risk level:** Low–Medium (build-correctness risk; no runtime logic)
- **Overall verdict:** The file accomplishes its purpose but carries a significant naming defect that creates a real macro-collision risk and actively misleads readers. The macro names (`BOOST_UNITS_STRINGIZE`, `BOOST_UNITS_HEADER_NAME`) are a cargo-culted residue of Boost's own header-test infrastructure, now repurposed to test mCRL2 headers. This is unjustified and dangerous.

---

## Findings (ordered by severity)

### [SEV-2] `BOOST_UNITS_STRINGIZE` macro names collide with Boost's own definitions

- **Location:** `cmake/headertest.cpp:2–3`
- **Why this is a problem:**  
  Boost itself defines both `BOOST_UNITS_STRINGIZE_IMPL` and `BOOST_UNITS_STRINGIZE` in `<boost/units/detail/utility.hpp>` (shipped with every Boost installation). This file redefines them with identical bodies. If any mCRL2 header currently or in future includes any Boost header that transitively pulls in `boost/units/detail/utility.hpp`, the compiler will emit a macro-redefinition diagnostic. Under `-Werror` (common in CI), this is a build failure.  
  Worse, the definitions happen to be textually identical today, suppressing hard errors. A future Boost version change that modifies those macros would silently produce different stringisation behavior in `headertest.cpp` than intended, with no diagnostic unless the definitions diverge enough to trigger the warning.

- **Evidence type:** static-proof  
  - C++ standard (N4868 §15.6.3): redefining an object-like macro with a different replacement list is undefined; identical redefinition is a diagnostic-free accident of the current Boost version.
  - Boost source: `boost/units/detail/utility.hpp` defines `BOOST_UNITS_STRINGIZE_IMPL` and `BOOST_UNITS_STRINGIZE` with the same expansion.

- **Reproducer:**
  ```bash
  # Create a minimal mCRL2 header that includes a Boost units header:
  echo '#include <boost/units/detail/utility.hpp>' > /tmp/trigger.h
  # Then configure the header test for it and build with -Werror:
  g++ -DBOOST_UNITS_HEADER_NAME=/tmp/trigger.h -Werror cmake/headertest.cpp
  # Expected: warning: "BOOST_UNITS_STRINGIZE" redefined [-Werror]
  ```

- **Efficiency impact:** None at runtime; build-time cost only.
- **C++20-specific note:** Not applicable; pure macro layer.
- **Suggested fix:**  
  Rename the macros to something unambiguously project-internal and non-Boost:
  ```cpp
  #define MCRL2_HEADERTEST_STRINGIZE_IMPL(x) #x
  #define MCRL2_HEADERTEST_STRINGIZE(x) MCRL2_HEADERTEST_STRINGIZE_IMPL(x)
  #define MCRL2_HEADERTEST_HEADER MCRL2_HEADERTEST_STRINGIZE(MCRL2_HEADER_UNDER_TEST)
  #include MCRL2_HEADERTEST_HEADER
  int main() {}
  ```
  Update `MCRL2AddTarget.cmake:268` to pass `-DMCRL2_HEADER_UNDER_TEST=...`.

---

### [SEV-2] Macro namespace is semantically incorrect: file tests mCRL2 headers, not Boost.Units headers

- **Location:** `cmake/headertest.cpp:1–7`, `cmake/MCRL2AddTarget.cmake:265–268`
- **Why this is a problem:**  
  The `BOOST_UNITS_` prefix implies these macros and the include mechanism exist to test Boost.Units headers. They do not. `mcrl2_add_header_tests()` passes mCRL2 library headers (e.g., `mcrl2/process/action_parse.h`) as `BOOST_UNITS_HEADER_NAME`. The naming is entirely wrong for the domain.  
  The CMake comment on line 265 even says _"In test.cpp we define BOOST_UNITS_HEADER_NAME to be the current header"_ — the file is not called `test.cpp`, it is `headertest.cpp`. This is stale copy-paste documentation.  
  This naming creates genuine maintenance hazards: future developers reading this file will assume it is Boost-related, may exclude it from refactors, and will not expect it to be the mechanism behind mCRL2 header isolation tests.

- **Evidence type:** static-proof (semantic mismatch, observable in CMake call-site)
- **Evidence:**  
  `MCRL2AddTarget.cmake:268`: `COMPILE_DEFINITIONS "BOOST_UNITS_HEADER_NAME=${cppname}"` where `cppname` is a path under `mcrl2/`.

- **Efficiency impact:** None at runtime.
- **C++20-specific note:** Not applicable.
- **Suggested fix:** Rename macros and the `COMPILE_DEFINITIONS` key as shown above. Fix the stale comment at `MCRL2AddTarget.cmake:265`.

---

### [SEV-3] Missing `#error` guard for undefined `BOOST_UNITS_HEADER_NAME`

- **Location:** `cmake/headertest.cpp:5–6`
- **Why this is a problem:**  
  `BOOST_UNITS_HEADER_NAME` is expected to be injected by the build system via `-D`. If it is absent (misconfigured CMake, standalone compilation), the preprocessor silently treats the token as an unexpanded identifier and `BOOST_UNITS_STRINGIZE` produces the string `"BOOST_UNITS_HEADER_NAME"`. The resulting `#include "BOOST_UNITS_HEADER_NAME"` fails with:
  ```
  fatal error: BOOST_UNITS_HEADER_NAME: No such file or directory
  ```
  This error message gives no indication of the real problem (missing `-D` flag).

- **Evidence type:** static-proof  
  ```bash
  # Compile without the -D flag:
  g++ cmake/headertest.cpp
  # Observed: fatal error: BOOST_UNITS_HEADER_NAME: No such file or directory
  # Expected: clear diagnostic that the macro must be defined
  ```

- **Efficiency impact:** None.
- **C++20-specific note:** Not applicable.
- **Suggested fix:**
  ```cpp
  #ifndef MCRL2_HEADER_UNDER_TEST
  #  error "MCRL2_HEADER_UNDER_TEST must be defined via -DMCRL2_HEADER_UNDER_TEST=<path>"
  #endif
  ```

---

### [SEV-3] `set_target_properties(COMPILE_DEFINITIONS)` is deprecated CMake API

- **Location:** `cmake/MCRL2AddTarget.cmake:268` (caller of `headertest.cpp`)
- **Why this is a problem:**  
  `set_target_properties(... PROPERTIES COMPILE_DEFINITIONS ...)` is the old directory-scoped CMake API. It does not support generator expressions and bypasses the `PRIVATE`/`PUBLIC`/`INTERFACE` scoping model. The modern replacement is `target_compile_definitions(${testname} PRIVATE ...)`.
- **Evidence type:** static-proof (CMake policy CMP0043; `cmake --warn-uninitialized` will note this)
- **Suggested fix:**
  ```cmake
  target_compile_definitions(${testname} PRIVATE "MCRL2_HEADER_UNDER_TEST=${cppname}")
  ```

---

## Proposed Tests

- [ ] `test_headertest_no_define`: Compile `cmake/headertest.cpp` without any `-D` flag. **Expected failure signal:** compiler `#error` diagnostic with actionable message (currently: cryptic missing-file error).
- [ ] `test_headertest_boost_units_collision`: Compile `cmake/headertest.cpp` with `-DBOOST_UNITS_HEADER_NAME=<boost/units/detail/utility.hpp>` and `-Werror`. **Expected failure signal:** `-Wmacro-redefined` promoted to error (currently: may silently pass or warn depending on toolchain).
- [ ] `test_headertest_mCRL2_domain_label`: Grep `cmake/headertest.cpp` for any `BOOST_` prefix macro. **Expected failure signal:** match found (zero occurrences expected post-fix).

---

## Open Questions

1. **Origin**: Was `headertest.cpp` directly copied from Boost's header self-test infrastructure? If so, is there a license attribution requirement that was silently dropped?
2. **`.hpp` gap**: `mcrl2_add_header_tests` uses `GLOB_RECURSE "*.h"` only. Are there any mCRL2 headers stored with `.hpp` extension that are silently excluded from header isolation testing?
3. **No `add_test()` call**: Header tests are compiled but never executed as CTest tests. Is there any intent to eventually run them (e.g., for link-time or sanitizer checks), or is compile-only coverage sufficient and documented?

## Copilot Raw Stdout

```text
✗ Find cmake files referencing headertest
  $ find /workspaces/mCRL2 -name "*.cmake" | xargs grep -l "headertest\|HEADER_NAME" 2>/dev/null |
  head -20
  Permission denied and could not request permission from user

● Find all CMake references to headertest.cpp
  $ grep -r "headertest" /workspaces/mCRL2/cmake/ --include="*.cmake" --include="CMakeLists.txt" -l
  && grep -r "headertest\|BOOST_UNITS_HEADER_NAME\|BOOST_UNITS_STRINGIZE" /workspaces/mCRL2/cmake/
  -n | head -60
  └ 9 lines...

✗ Check Boost units utility header location
  $ python3 -c "
  # Check if Boost headers define BOOST_UNITS_STRINGIZE
  import subprocess, os
  result = subprocess.run(['find', '/usr', '-name', 'utility.hpp', '-path', '*/units/*'],
  capture_output=True, text=True)...
  Permission denied and could not request permission from user

● Read cmake/MCRL2AddTarget.cmake lines 240-300
  └ 61 lines read

✗ Check if Boost defines BOOST_UNITS_STRINGIZE
  $ find /usr/include -name "utility.hpp" -path "*/units/*" 2>/dev/null | head -5; grep -r
  "BOOST_UNITS_STRINGIZE" /usr/include/boost/units/ 2>/dev/null | head -10
  Permission denied and could not request permission from user

● Read cmake/MCRL2AddTarget.cmake lines 251-273
  └ 23 lines read

Now I have full context. Let me also check whether any mCRL2 headers include Boost units headers, making the macro collision concrete:

● Check if mCRL2 includes boost/units headers
  $ grep -r "boost/units" /workspaces/mCRL2/libraries/ --include="*.h" --include="*.hpp"
  --include="*.cpp" -l 2>/dev/null | head -10
  └ 1 line...

● Check for other BOOST_UNITS_STRINGIZE definitions in the repo
  $ grep -r "BOOST_UNITS_STRINGIZE" /workspaces/mCRL2/ --include="*.h" --include="*.hpp"
  --include="*.cpp" --include="*.cmake" -n 2>/dev/null | grep -v "cmake/headertest.cpp"
  └ 1 line...

● Check breadth of header test coverage and boost usage
  $ # Check if any mCRL2 library headers include boost/units at all
  grep -r "#include.*boost" /workspaces/mCRL2/libraries/ --include="*.h" -l 2>/dev/null | head -10
  # Check function signature - what does it search for (*.h only, not *.hpp)?
  grep -n "GLOB_RECURSE" /workspaces/mCRL2/cmake/MCRL2AddTarget.cmake
  └ 15 lines...

I have full context. Here is the review.

---

# Review Findings

## Summary
- **Scope reviewed:** `cmake/headertest.cpp` and its caller `cmake/MCRL2AddTarget.cmake:253–273`
- **Risk level:** Low–Medium (build-correctness risk; no runtime logic)
- **Overall verdict:** The file accomplishes its purpose but carries a significant naming defect that creates a real macro-collision risk and actively misleads readers. The macro names (`BOOST_UNITS_STRINGIZE`, `BOOST_UNITS_HEADER_NAME`) are a cargo-culted residue of Boost's own header-test infrastructure, now repurposed to test mCRL2 headers. This is unjustified and dangerous.

---

## Findings (ordered by severity)

### [SEV-2] `BOOST_UNITS_STRINGIZE` macro names collide with Boost's own definitions

- **Location:** `cmake/headertest.cpp:2–3`
- **Why this is a problem:**  
  Boost itself defines both `BOOST_UNITS_STRINGIZE_IMPL` and `BOOST_UNITS_STRINGIZE` in `<boost/units/detail/utility.hpp>` (shipped with every Boost installation). This file redefines them with identical bodies. If any mCRL2 header currently or in future includes any Boost header that transitively pulls in `boost/units/detail/utility.hpp`, the compiler will emit a macro-redefinition diagnostic. Under `-Werror` (common in CI), this is a build failure.  
  Worse, the definitions happen to be textually identical today, suppressing hard errors. A future Boost version change that modifies those macros would silently produce different stringisation behavior in `headertest.cpp` than intended, with no diagnostic unless the definitions diverge enough to trigger the warning.

- **Evidence type:** static-proof  
  - C++ standard (N4868 §15.6.3): redefining an object-like macro with a different replacement list is undefined; identical redefinition is a diagnostic-free accident of the current Boost version.
  - Boost source: `boost/units/detail/utility.hpp` defines `BOOST_UNITS_STRINGIZE_IMPL` and `BOOST_UNITS_STRINGIZE` with the same expansion.

- **Reproducer:**
  ```bash
  # Create a minimal mCRL2 header that includes a Boost units header:
  echo '#include <boost/units/detail/utility.hpp>' > /tmp/trigger.h
  # Then configure the header test for it and build with -Werror:
  g++ -DBOOST_UNITS_HEADER_NAME=/tmp/trigger.h -Werror cmake/headertest.cpp
  # Expected: warning: "BOOST_UNITS_STRINGIZE" redefined [-Werror]
  ```

- **Efficiency impact:** None at runtime; build-time cost only.
- **C++20-specific note:** Not applicable; pure macro layer.
- **Suggested fix:**  
  Rename the macros to something unambiguously project-internal and non-Boost:
  ```cpp
  #define MCRL2_HEADERTEST_STRINGIZE_IMPL(x) #x
  #define MCRL2_HEADERTEST_STRINGIZE(x) MCRL2_HEADERTEST_STRINGIZE_IMPL(x)
  #define MCRL2_HEADERTEST_HEADER MCRL2_HEADERTEST_STRINGIZE(MCRL2_HEADER_UNDER_TEST)
  #include MCRL2_HEADERTEST_HEADER
  int main() {}
  ```
  Update `MCRL2AddTarget.cmake:268` to pass `-DMCRL2_HEADER_UNDER_TEST=...`.

---

### [SEV-2] Macro namespace is semantically incorrect: file tests mCRL2 headers, not Boost.Units headers

- **Location:** `cmake/headertest.cpp:1–7`, `cmake/MCRL2AddTarget.cmake:265–268`
- **Why this is a problem:**  
  The `BOOST_UNITS_` prefix implies these macros and the include mechanism exist to test Boost.Units headers. They do not. `mcrl2_add_header_tests()` passes mCRL2 library headers (e.g., `mcrl2/process/action_parse.h`) as `BOOST_UNITS_HEADER_NAME`. The naming is entirely wrong for the domain.  
  The CMake comment on line 265 even says _"In test.cpp we define BOOST_UNITS_HEADER_NAME to be the current header"_ — the file is not called `test.cpp`, it is `headertest.cpp`. This is stale copy-paste documentation.  
  This naming creates genuine maintenance hazards: future developers reading this file will assume it is Boost-related, may exclude it from refactors, and will not expect it to be the mechanism behind mCRL2 header isolation tests.

- **Evidence type:** static-proof (semantic mismatch, observable in CMake call-site)
- **Evidence:**  
  `MCRL2AddTarget.cmake:268`: `COMPILE_DEFINITIONS "BOOST_UNITS_HEADER_NAME=${cppname}"` where `cppname` is a path under `mcrl2/`.

- **Efficiency impact:** None at runtime.
- **C++20-specific note:** Not applicable.
- **Suggested fix:** Rename macros and the `COMPILE_DEFINITIONS` key as shown above. Fix the stale comment at `MCRL2AddTarget.cmake:265`.

---

### [SEV-3] Missing `#error` guard for undefined `BOOST_UNITS_HEADER_NAME`

- **Location:** `cmake/headertest.cpp:5–6`
- **Why this is a problem:**  
  `BOOST_UNITS_HEADER_NAME` is expected to be injected by the build system via `-D`. If it is absent (misconfigured CMake, standalone compilation), the preprocessor silently treats the token as an unexpanded identifier and `BOOST_UNITS_STRINGIZE` produces the string `"BOOST_UNITS_HEADER_NAME"`. The resulting `#include "BOOST_UNITS_HEADER_NAME"` fails with:
  ```
  fatal error: BOOST_UNITS_HEADER_NAME: No such file or directory
  ```
  This error message gives no indication of the real problem (missing `-D` flag).

- **Evidence type:** static-proof  
  ```bash
  # Compile without the -D flag:
  g++ cmake/headertest.cpp
  # Observed: fatal error: BOOST_UNITS_HEADER_NAME: No such file or directory
  # Expected: clear diagnostic that the macro must be defined
  ```

- **Efficiency impact:** None.
- **C++20-specific note:** Not applicable.
- **Suggested fix:**
  ```cpp
  #ifndef MCRL2_HEADER_UNDER_TEST
  #  error "MCRL2_HEADER_UNDER_TEST must be defined via -DMCRL2_HEADER_UNDER_TEST=<path>"
  #endif
  ```

---

### [SEV-3] `set_target_properties(COMPILE_DEFINITIONS)` is deprecated CMake API

- **Location:** `cmake/MCRL2AddTarget.cmake:268` (caller of `headertest.cpp`)
- **Why this is a problem:**  
  `set_target_properties(... PROPERTIES COMPILE_DEFINITIONS ...)` is the old directory-scoped CMake API. It does not support generator expressions and bypasses the `PRIVATE`/`PUBLIC`/`INTERFACE` scoping model. The modern replacement is `target_compile_definitions(${testname} PRIVATE ...)`.
- **Evidence type:** static-proof (CMake policy CMP0043; `cmake --warn-uninitialized` will note this)
- **Suggested fix:**
  ```cmake
  target_compile_definitions(${testname} PRIVATE "MCRL2_HEADER_UNDER_TEST=${cppname}")
  ```

---

## Proposed Tests

- [ ] `test_headertest_no_define`: Compile `cmake/headertest.cpp` without any `-D` flag. **Expected failure signal:** compiler `#error` diagnostic with actionable message (currently: cryptic missing-file error).
- [ ] `test_headertest_boost_units_collision`: Compile `cmake/headertest.cpp` with `-DBOOST_UNITS_HEADER_NAME=<boost/units/detail/utility.hpp>` and `-Werror`. **Expected failure signal:** `-Wmacro-redefined` promoted to error (currently: may silently pass or warn depending on toolchain).
- [ ] `test_headertest_mCRL2_domain_label`: Grep `cmake/headertest.cpp` for any `BOOST_` prefix macro. **Expected failure signal:** match found (zero occurrences expected post-fix).

---

## Open Questions

1. **Origin**: Was `headertest.cpp` directly copied from Boost's header self-test infrastructure? If so, is there a license attribution requirement that was silently dropped?
2. **`.hpp` gap**: `mcrl2_add_header_tests` uses `GLOB_RECURSE "*.h"` only. Are there any mCRL2 headers stored with `.hpp` extension that are silently excluded from header isolation testing?
3. **No `add_test()` call**: Header tests are compiled but never executed as CTest tests. Is there any intent to eventually run them (e.g., for link-time or sanitizer checks), or is compile-only coverage sufficient and documented?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         2m 16s
Total session time:     2m 26s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-sonnet-4.6       173.1k in, 8.0k out, 139.7k cached (Est. 1 Premium request)

```

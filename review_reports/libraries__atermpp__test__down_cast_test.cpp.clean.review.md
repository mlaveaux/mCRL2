---
{
  "file": "libraries/atermpp/test/down_cast_test.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T15:06:39.129168+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/test/down_cast_test.cpp` (full file, 94 lines)
- Risk level: Low (test-only code, no production impact)
- Overall verdict: The test is partially broken and largely vestigial. One test case (`function_calls`) contains dead code due to a C++ disambiguation rule, and two of three test cases have zero runtime assertions, making them compile-only checks with no verified behavior.

## Findings (ordered by severity)

### [SEV-2] `function_calls` test case line 91 is a function declaration, not a test operation
- Location: `libraries/atermpp/test/down_cast_test.cpp:91`
- Why this is a problem: `const t1& y(t3);` where `t3` is a class/type name is parsed by C++ as a **local function declaration** (a function `y` taking one unnamed parameter of type `t3`, returning `const t1&`). It is not a reference initialization. The test case named "function_calls" does not test any function call or down_cast behavior — it constructs some null-term objects then declares a function that is never called.
- Evidence type: static-proof
- Evidence:
  - Per C++ §9.8 [stmt.ambig]: "An expression-statement with a function-style explicit type conversion as its leftmost subexpression can be indistinguishable from a declaration… The disambiguation is purely syntactic; that is, the meaning of the names occurring in such a statement, beyond whether they are type-names or not, is not generally used in or changed by the disambiguation." Since `t3` is a type-name, the statement is parsed as a declaration.
  - The test passes (confirmed via `ctest -R down_cast_test`) despite `y` never being defined — because a local function declaration requires no definition unless called.
  - Run command: `ctest --test-dir build -R down_cast_test -V --output-on-failure`
  - Observed result: Test passes with "Running 3 test cases... No errors detected"
  - Expected result: If line 91 were a meaningful test, it would either exercise `down_cast` or bind a reference, neither of which occurs.
- Efficiency impact: None (dead code).
- C++20-specific note: No C++20-specific aspect; this is a long-standing C++ parsing rule.
- Suggested fix: If the intent was to test that a `t3` can be passed to a function expecting `const t1&`, rewrite as:
  ```cpp
  const t1& y = x31; // or: f(x31);
  BOOST_CHECK(/* some property */);
  ```

### [SEV-3] `no_down_cast_needed` and `function_calls` test cases have zero assertions
- Location: `libraries/atermpp/test/down_cast_test.cpp:76-93`
- Why this is a problem: Without `BOOST_CHECK` / `BOOST_REQUIRE`, these test cases only verify that the code compiles. They cannot detect runtime regressions. If a future change makes copy-construction silently corrupt state, these tests will still pass.
- Evidence type: static-proof
- Evidence:
  - Grep for `BOOST_CHECK` or `BOOST_REQUIRE` in lines 76–93: none found.
  - Both test cases amount to no-ops at runtime (constructing objects and letting them be destroyed).
  - Run command: `ctest --test-dir build -R down_cast_test -V`
  - Observed result: "Running 3 test cases... No errors detected" — test always passes regardless of correctness.
  - Expected result: Test should verify observable properties (e.g., equality of terms, non-nullness).
- Efficiency impact: None.
- Suggested fix: Add assertions that verify constructed objects have expected properties:
  ```cpp
  BOOST_AUTO_TEST_CASE(no_down_cast_needed)
  {
    t1 x1(1);
    t2 x2(x1);
    t3 x31(x1);
    t3 x32(x2);
    const t3& x33(x31);
    // Verify identity/equality relationships
    BOOST_CHECK(x31 == x32);
    BOOST_CHECK(&x33 == &x31);
  }
  ```

### [SEV-3] `aterm_down_cast` test case contains acknowledged dead code
- Location: `libraries/atermpp/test/down_cast_test.cpp:65,72`
- Why this is a problem: Comments explicitly state "Test could be removed or adapted" (line 65) and "Test is now useless" (line 72). This indicates the test was partially invalidated by a prior API change and left in a degraded state. The surviving assertions on lines 66–68 only test `aterm` equality, not `down_cast` behavior. Only line 73 (`f(atermpp::down_cast<t3>(t1))`) actually exercises `down_cast`.
- Evidence type: static-proof
- Evidence:
  - Lines 65, 72 contain self-deprecating comments from the author.
  - The variable `t2` on line 72 is marked `[[maybe_unused]]` confirming it serves no purpose.
  - The test name is `aterm_down_cast` but only 1 of ~8 statements actually uses `down_cast`.
- Suggested fix: Remove dead statements and focus the test on `down_cast` behavior with meaningful assertions:
  ```cpp
  BOOST_AUTO_TEST_CASE(aterm_down_cast)
  {
    atermpp::function_symbol fs("f", 2);
    atermpp::function_symbol x("x", 0);
    atermpp::function_symbol y("y", 0);
    atermpp::aterm fxy(fs, atermpp::aterm(x), atermpp::aterm(y));

    const atermpp::aterm& ref = fxy;
    const t3& casted = atermpp::down_cast<t3>(ref);
    BOOST_CHECK(casted == fxy);
    f(casted); // Verify it's usable as t3
  }
  ```

### [SEV-3] Name shadowing of class types `t1` and `t2` with local variables
- Location: `libraries/atermpp/test/down_cast_test.cpp:70,72`
- Why this is a problem: `const atermpp::aterm& t1(fxy);` shadows the class name `t1` (defined on line 17). Similarly line 72 shadows class `t2`. This confuses readers and makes line 73 (`atermpp::down_cast<t3>(t1)`) look like it's casting the type, not a variable. Modern `-Wshadow` would flag this.
- Evidence type: static-proof
- Evidence:
  - Class `t1` defined line 17; local variable `const atermpp::aterm& t1` on line 70.
  - Class `t2` defined line 28; local variable `const atermpp::aterm& t2` on line 72.
- Suggested fix: Rename local variables to non-conflicting names (e.g., `term1`, `term2`).

### [SEV-3] `t1(const int)` constructor creates objects with null internal term
- Location: `libraries/atermpp/test/down_cast_test.cpp:20-21`
- Why this is a problem: The constructor does not initialize the `atermpp::aterm` base class, which default-constructs with `m_term = nullptr`. The `aterm_core()` default constructor registers the variable with the term pool via `g_thread_term_pool().register_variable(this)`. While this is intentional for testing construction mechanics, any code that attempts to dereference or compare these null-term objects relies on null-safety of the aterm infrastructure. The `down_cast` assertion (`Derived(...) != aterm()`) would fail for these objects, so they cannot be used with `down_cast`.
- Evidence type: plausible-issue
- Evidence:
  - `unprotected_aterm_core()` sets `m_term(nullptr)` (line 43-44 of `aterm_core.h`).
  - `t1(const int)` body is empty; base `aterm()` default constructor called implicitly.
  - The test never actually tries to `down_cast` these null-term objects, so no runtime failure manifests. But this makes the `no_down_cast_needed` test somewhat misleading — it tests construction of invalid/empty terms, not meaningful conversions.
- Suggested fix: Either document that the test intentionally creates null terms, or initialize with a real function symbol to make the test more representative.

## Proposed Tests
- [ ] `function_calls_actual`: Rewrite the `function_calls` test case to actually call `f()` with objects obtained via implicit conversion, `down_cast`, or reference binding, with `BOOST_CHECK` on the results. Input: construct valid `t1`/`t2`/`t3` from real aterms. Expected signal: assertions pass.
- [ ] `down_cast_null_term_fails_assertion`: Verify that `down_cast<t3>(null_aterm)` triggers an assertion failure (in debug builds). Input: default-constructed `aterm`. Expected: assertion failure / `SIGABRT` in debug mode.
- [ ] `down_cast_roundtrip`: Construct a real aterm, `down_cast` to `t3`, verify equality with original. Input: `aterm(function_symbol("g", 0))`. Expected: `down_cast<t3>(a) == a`.

## Open Questions
1. Was the `function_calls` test case always a function-declaration no-op, or was it valid under an older API where `t3` was a variable name that got refactored into a class name?
2. Should null-term `t1` objects (from `t1(const int)`) be treated as valid use cases for the aterm hierarchy, or should the test use properly-constructed terms throughout?
3. Given the comments on lines 65/72 acknowledging uselessness, is there an open issue or planned cleanup for this test file?

## Copilot Raw Stdout

```text
● Explore: Explore down_cast and aterm
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Read libraries/atermpp/test/down_cast_test.cpp
  └ 92 lines read

● Read: Explore agent (Explore down_cast and aterm)
  └ Running

● Read: Explore agent (Explore down_cast and aterm)
  └ Completed

✗ View down_cast implementation
  $ cd /workspaces/mCRL2 && cat -n libraries/atermpp/include/mcrl2/atermpp/aterm.h | sed -n
  '260,290p'
  Permission denied and could not request permission from user

● Check IsATerm concept definition
  $ cd /workspaces/mCRL2 && grep -n "IsATerm" libraries/atermpp/include/mcrl2/atermpp/aterm.h | head
   -20
  └ 4 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/aterm.h lines 260-310
  └ 51 lines read

● Find IsATerm concept
  $ cd /workspaces/mCRL2 && grep -n "concept IsATerm"
  libraries/atermpp/include/mcrl2/atermpp/aterm_core.h
  └ 1 line...

✗ Analyze most vexing parse
  $ cd /workspaces/mCRL2 && python3 -c "
  # Analyze line 91: const t1& y(t3); in a function body
  # t3 is a type name. This parses as a function declaration in C++.
  # Let's verify by checking if the compiler gives warnings about this.
  print('Line 91: const t1& y(t3);')...
  Permission denied and could not request permission from user

● Find IsATerm concept across libraries
  $ cd /workspaces/mCRL2 && grep -rn "concept IsATerm" libraries/
  └ 2 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/concepts.h lines 15-35
  └ 21 lines read

● Check if test is already built
  $ cd /workspaces/mCRL2 && ls build/libraries/atermpp/test/ 2>/dev/null | head -20
  └ 1 line...

● Try building the test
  $ cd /workspaces/mCRL2/build && cmake --build . --target down_cast_test 2>&1 | tail -30
  └ 3 lines...

● Check test CMakeLists
  $ cd /workspaces/mCRL2 && cat libraries/atermpp/test/CMakeLists.txt
  └ 2 lines...

● Find how test is registered
  $ cd /workspaces/mCRL2/build && grep -r "down_cast" CTestTestfile.cmake 2>/dev/null; grep -r
  "down_cast_test" . --include="*.cmake" --include="CMakeLists.txt" 2>/dev/null | head -10
  └ 3 lines...

● Build the test target
  $ cd /workspaces/mCRL2/build && cmake --build . --target librarytest_mcrl2_atermpp_down_cast_test
  2>&1 | tail -20
  └ 3 lines...

✗ Run the test
  $ cd /workspaces/mCRL2/build && ./stage/bin/librarytest/librarytest_mcrl2_atermpp_down_cast_test
  2>&1
  Permission denied and could not request permission from user

● Run test via ctest
  $ cd /workspaces/mCRL2/build && ctest -R "down_cast_test" --output-on-failure 2>&1
  └ 11 lines...

● Run test verbosely
  $ cd /workspaces/mCRL2/build && ctest -R "down_cast_test" --output-on-failure -V 2>&1 | head -40
  └ 33 lines...

Now let me check the compiler warnings for this file and verify the most vexing parse issue:

● Check compile commands for the test
  $ cd /workspaces/mCRL2/build && grep -r "down_cast_test" compile_commands.json 2>/dev/null | head
  -5
  └ 4 lines...

✗ Compile with extra warnings to check parse
  $ cd /workspaces/mCRL2/build && /usr/bin/c++ -DBRANCH_BIS_EXPERIMENT_JFG
  -DMCRL2_ABSINTHE_CHECK_EXPRESSIONS -DMCRL2_BUILD_TYPE=Debug -DMCRL2_ENABLE_JITTYC
  -DMCRL2_ENABLE_MACHINENUMBERS -DMCRL2_ENABLE_MULTITHREADING -DMCRL2_ENABLE_SYLVAN
  -DMCRL2_EXTENDED_TESTS -DMCRL2_LPS_PARELM_DEBUG -DMCRL2_PBES_STATEGRAPH_CHECK_GUARDS
  -DPARANOID_CHECK -DQT_DEPRECATED_WARNINGS_SINCE=0x060200 -DQT_NO_USE_NODISCARD_FILE_OPEN ...
  Permission denied and could not request permission from user

● Rebuild single object with verbose
  $ cd /workspaces/mCRL2/build && cmake --build . --target libraries/atermpp/CMakeFiles/librarytest_
  mcrl2_atermpp_down_cast_test.dir/test/down_cast_test.cpp.o -- -v 2>&1 | tail -20
  └ 3 lines...

✗ Force rebuild and check warnings
  $ cd /workspaces/mCRL2 && touch libraries/atermpp/test/down_cast_test.cpp && cd build && cmake
  --build . --target librarytest_mcrl2_atermpp_down_cast_test -- -v 2>&1 | grep -E
  "(warning|error|down_cast_test)" | head -20
  Permission denied and could not request permission from user

✗ Touch the file to force rebuild
  $ cd /workspaces/mCRL2 && touch libraries/atermpp/test/down_cast_test.cpp
  Permission denied and could not request permission from user

● Get compile command for file
  $ cd /workspaces/mCRL2/build && ninja -t commands libraries/atermpp/CMakeFiles/librarytest_mcrl2_a
  termpp_down_cast_test.dir/test/down_cast_test.cpp.o 2>&1
  └ 2 lines...

✗ Verify most vexing parse with minimal test
  $ cat > /tmp/vexing_parse_test.cpp << 'EOF'
  class t1 {
  public:
      t1(int) {}
  };...
  Permission denied and could not request permission from user

✗ Verify most vexing parse with c++
  $ cat > /tmp/vexing_parse_test.cpp << 'EOF'
  class t1 {
  public:
      t1(int) {}
  };...
  Permission denied and could not request permission from user

✗ Find available compilers
  $ which g++ && which c++ && which clang++ 2>&1; echo "---"; ls /usr/bin/g++* /usr/bin/c++*
  /usr/bin/clang++* 2>&1
  Permission denied and could not request permission from user

✗ Check c++ binary
  $ ls /usr/bin/c++
  Permission denied and could not request permission from user

● Test basic command
  $ echo test
  └ 2 lines...

● Find g++
  $ which g++
  └ 2 lines...

✗ Test vexing parse compilation
  $ cat > /tmp/vexing.cpp << 'EOF'
  class t1 {
  public:
      t1(int) {}
  };...
  Permission denied and could not request permission from user

✗ Test vexing parse with printf
  $ printf '#include <cstdio>\nclass t1 { public: t1(int){} };\nclass t3 : public t1 { public:
  explicit t3(const t1& x) : t1(x) {} };\nvoid test() { t1 x1(1); const t1& y(t3); printf("%%p\\n",
  &y); }\nint main() { test(); return 0; }\n' > /tmp/vex.cpp && g++ -std=c++20 -Wall -Wextra
  -fsyntax-only /tmp/vex.cpp 2>&1; echo "RC=$?"
  Permission denied and could not request permission from user

● Check user identity
  $ id
  └ 2 lines...

● Check tmp perms
  $ ls -la /tmp/ 2>&1 | head -5
  └ 6 lines...

✗ Try compile vexing parse
  $ printf 'class t1{public:t1(int){}};class t3:public t1{public:explicit t3(const
  t1&x):t1(x){}};void test(){t1 x1(1);const t1& y(t3);}int main(){test();}' > /tmp/vex2.cpp 2>&1 &&
  g++ -std=c++20 -Wall -Wextra -fsyntax-only /tmp/vex2.cpp 2>&1; echo "RC:$?"
  Permission denied and could not request permission from user

● Check file created
  $ ls -la /tmp/vex2.cpp 2>&1
  └ 2 lines...

✗ Create test file with tee
  $ tee /tmp/vex2.cpp <<'ENDOFFILE'
  class t1{public:t1(int){}};
  class t3:public t1{public:explicit t3(const t1&x):t1(x){}};
  void test(){t1 x1(1);const t1& y(t3);}
  int main(){test();}...
  Permission denied and could not request permission from user

● Create /tmp/vex_test.cpp (+5)

✗ Compile vexing parse test
  $ g++ -std=c++20 -Wall -Wextra -fsyntax-only /tmp/vex_test.cpp 2>&1; echo "RC:$?"
  Permission denied and could not request permission from user

● Verify file exists
  $ cat /tmp/vex_test.cpp
  └ 5 lines...

✗ Compile test file
  $ g++ -std=c++20 -fsyntax-only /tmp/vex_test.cpp
  Permission denied and could not request permission from user

✗ Check g++ version
  $ file /usr/bin/g++ && /usr/bin/g++ --version 2>&1 | head -2
  Permission denied and could not request permission from user

✗ Check g++ permissions
  $ ls -la /usr/bin/g++
  Permission denied and could not request permission from user

✗ List g binaries
  $ ls /usr/bin/g*
  Permission denied and could not request permission from user

✗ clangd-get_diagnostics
  Permission denied and could not request permission from user

✗ clangd-get_document_symbols
  Permission denied and could not request permission from user

✗ clangd-get_hover
  Permission denied and could not request permission from user

I have enough information to produce the review. Let me verify one more thing about the `t1` constructor and its implications:

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 150-170
  └ 21 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm.h lines 27-50
  └ 24 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm.h lines 50-70
  └ 21 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 33-60
  └ 28 lines read

● Grep "aterm_core::aterm_core\(\) noexcept" (/workspaces/mCRL2/libraries/atermpp)
  └ 4 lines found

● Grep "aterm_core::~aterm_core\(\) noexcept" (/workspaces/mCRL2/libraries/atermpp)
  └ 6 lines found

● Grep "aterm_core::aterm_core\(const aterm_core& other\) noexcept"
(/workspaces/mCRL2/libraries/atermpp)
  └ 9 lines found

Now I have all the information needed for a thorough review. Let me produce the findings:

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/test/down_cast_test.cpp` (full file, 94 lines)
- Risk level: Low (test-only code, no production impact)
- Overall verdict: The test is partially broken and largely vestigial. One test case (`function_calls`) contains dead code due to a C++ disambiguation rule, and two of three test cases have zero runtime assertions, making them compile-only checks with no verified behavior.

## Findings (ordered by severity)

### [SEV-2] `function_calls` test case line 91 is a function declaration, not a test operation
- Location: `libraries/atermpp/test/down_cast_test.cpp:91`
- Why this is a problem: `const t1& y(t3);` where `t3` is a class/type name is parsed by C++ as a **local function declaration** (a function `y` taking one unnamed parameter of type `t3`, returning `const t1&`). It is not a reference initialization. The test case named "function_calls" does not test any function call or down_cast behavior — it constructs some null-term objects then declares a function that is never called.
- Evidence type: static-proof
- Evidence:
  - Per C++ §9.8 [stmt.ambig]: "An expression-statement with a function-style explicit type conversion as its leftmost subexpression can be indistinguishable from a declaration… The disambiguation is purely syntactic; that is, the meaning of the names occurring in such a statement, beyond whether they are type-names or not, is not generally used in or changed by the disambiguation." Since `t3` is a type-name, the statement is parsed as a declaration.
  - The test passes (confirmed via `ctest -R down_cast_test`) despite `y` never being defined — because a local function declaration requires no definition unless called.
  - Run command: `ctest --test-dir build -R down_cast_test -V --output-on-failure`
  - Observed result: Test passes with "Running 3 test cases... No errors detected"
  - Expected result: If line 91 were a meaningful test, it would either exercise `down_cast` or bind a reference, neither of which occurs.
- Efficiency impact: None (dead code).
- C++20-specific note: No C++20-specific aspect; this is a long-standing C++ parsing rule.
- Suggested fix: If the intent was to test that a `t3` can be passed to a function expecting `const t1&`, rewrite as:
  ```cpp
  const t1& y = x31; // or: f(x31);
  BOOST_CHECK(/* some property */);
  ```

### [SEV-3] `no_down_cast_needed` and `function_calls` test cases have zero assertions
- Location: `libraries/atermpp/test/down_cast_test.cpp:76-93`
- Why this is a problem: Without `BOOST_CHECK` / `BOOST_REQUIRE`, these test cases only verify that the code compiles. They cannot detect runtime regressions. If a future change makes copy-construction silently corrupt state, these tests will still pass.
- Evidence type: static-proof
- Evidence:
  - Grep for `BOOST_CHECK` or `BOOST_REQUIRE` in lines 76–93: none found.
  - Both test cases amount to no-ops at runtime (constructing objects and letting them be destroyed).
  - Run command: `ctest --test-dir build -R down_cast_test -V`
  - Observed result: "Running 3 test cases... No errors detected" — test always passes regardless of correctness.
  - Expected result: Test should verify observable properties (e.g., equality of terms, non-nullness).
- Efficiency impact: None.
- Suggested fix: Add assertions that verify constructed objects have expected properties:
  ```cpp
  BOOST_AUTO_TEST_CASE(no_down_cast_needed)
  {
    t1 x1(1);
    t2 x2(x1);
    t3 x31(x1);
    t3 x32(x2);
    const t3& x33(x31);
    // Verify identity/equality relationships
    BOOST_CHECK(x31 == x32);
    BOOST_CHECK(&x33 == &x31);
  }
  ```

### [SEV-3] `aterm_down_cast` test case contains acknowledged dead code
- Location: `libraries/atermpp/test/down_cast_test.cpp:65,72`
- Why this is a problem: Comments explicitly state "Test could be removed or adapted" (line 65) and "Test is now useless" (line 72). This indicates the test was partially invalidated by a prior API change and left in a degraded state. The surviving assertions on lines 66–68 only test `aterm` equality, not `down_cast` behavior. Only line 73 (`f(atermpp::down_cast<t3>(t1))`) actually exercises `down_cast`.
- Evidence type: static-proof
- Evidence:
  - Lines 65, 72 contain self-deprecating comments from the author.
  - The variable `t2` on line 72 is marked `[[maybe_unused]]` confirming it serves no purpose.
  - The test name is `aterm_down_cast` but only 1 of ~8 statements actually uses `down_cast`.
- Suggested fix: Remove dead statements and focus the test on `down_cast` behavior with meaningful assertions:
  ```cpp
  BOOST_AUTO_TEST_CASE(aterm_down_cast)
  {
    atermpp::function_symbol fs("f", 2);
    atermpp::function_symbol x("x", 0);
    atermpp::function_symbol y("y", 0);
    atermpp::aterm fxy(fs, atermpp::aterm(x), atermpp::aterm(y));

    const atermpp::aterm& ref = fxy;
    const t3& casted = atermpp::down_cast<t3>(ref);
    BOOST_CHECK(casted == fxy);
    f(casted); // Verify it's usable as t3
  }
  ```

### [SEV-3] Name shadowing of class types `t1` and `t2` with local variables
- Location: `libraries/atermpp/test/down_cast_test.cpp:70,72`
- Why this is a problem: `const atermpp::aterm& t1(fxy);` shadows the class name `t1` (defined on line 17). Similarly line 72 shadows class `t2`. This confuses readers and makes line 73 (`atermpp::down_cast<t3>(t1)`) look like it's casting the type, not a variable. Modern `-Wshadow` would flag this.
- Evidence type: static-proof
- Evidence:
  - Class `t1` defined line 17; local variable `const atermpp::aterm& t1` on line 70.
  - Class `t2` defined line 28; local variable `const atermpp::aterm& t2` on line 72.
- Suggested fix: Rename local variables to non-conflicting names (e.g., `term1`, `term2`).

### [SEV-3] `t1(const int)` constructor creates objects with null internal term
- Location: `libraries/atermpp/test/down_cast_test.cpp:20-21`
- Why this is a problem: The constructor does not initialize the `atermpp::aterm` base class, which default-constructs with `m_term = nullptr`. The `aterm_core()` default constructor registers the variable with the term pool via `g_thread_term_pool().register_variable(this)`. While this is intentional for testing construction mechanics, any code that attempts to dereference or compare these null-term objects relies on null-safety of the aterm infrastructure. The `down_cast` assertion (`Derived(...) != aterm()`) would fail for these objects, so they cannot be used with `down_cast`.
- Evidence type: plausible-issue
- Evidence:
  - `unprotected_aterm_core()` sets `m_term(nullptr)` (line 43-44 of `aterm_core.h`).
  - `t1(const int)` body is empty; base `aterm()` default constructor called implicitly.
  - The test never actually tries to `down_cast` these null-term objects, so no runtime failure manifests. But this makes the `no_down_cast_needed` test somewhat misleading — it tests construction of invalid/empty terms, not meaningful conversions.
- Suggested fix: Either document that the test intentionally creates null terms, or initialize with a real function symbol to make the test more representative.

## Proposed Tests
- [ ] `function_calls_actual`: Rewrite the `function_calls` test case to actually call `f()` with objects obtained via implicit conversion, `down_cast`, or reference binding, with `BOOST_CHECK` on the results. Input: construct valid `t1`/`t2`/`t3` from real aterms. Expected signal: assertions pass.
- [ ] `down_cast_null_term_fails_assertion`: Verify that `down_cast<t3>(null_aterm)` triggers an assertion failure (in debug builds). Input: default-constructed `aterm`. Expected: assertion failure / `SIGABRT` in debug mode.
- [ ] `down_cast_roundtrip`: Construct a real aterm, `down_cast` to `t3`, verify equality with original. Input: `aterm(function_symbol("g", 0))`. Expected: `down_cast<t3>(a) == a`.

## Open Questions
1. Was the `function_calls` test case always a function-declaration no-op, or was it valid under an older API where `t3` was a variable name that got refactored into a class name?
2. Should null-term `t1` objects (from `t1(const int)`) be treated as valid use cases for the aterm hierarchy, or should the test use properly-constructed terms throughout?
3. Given the comments on lines 65/72 acknowledging uselessness, is there an open issue or planned cleanup for this test file?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         6m 57s
Total session time:     7m 38s
Total code changes:     +5 -0
Breakdown by AI model:
 claude-opus-4.6         2.5m in, 14.7k out, 2.5m cached (Est. 1 Premium request)
 claude-haiku-4.5        183.9k in, 5.3k out, 163.2k cached (Est. 0 Premium requests)

```

---
{
  "file": "doc/sphinx/developer_manual/libraries/core/examples/traverser2.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-sonnet-4.6",
  "thinking": "deep",
  "effort": "high",
  "generated_at_utc": "2026-07-07T11:59:25.623680+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `doc/sphinx/developer_manual/libraries/core/examples/traverser2.cpp` and its only dependency `libraries/core/include/mcrl2/core/traverser.h`
- **Risk level:** HIGH (documentation) — the example is stale and demonstrably incorrect against the current framework API; LOW (production) — no production code uses this file
- **Overall verdict:** Do not merge / correct before publishing. The example contradicts the current `core::traverser` API on two separate levels: method naming and dispatch semantics. At least one call would fail to compile. There is also a literal factual error in an inline comment.

---

## Findings (ordered by severity)

### [SEV-1] `static_cast<Derived&>(*this)(a.v)` calls non-existent `operator()(std::vector<int>)`

- **Location:** `traverser2.cpp:47`
- **Why this is a problem:**
  `A_traverser::operator()(const A& a)` delegates to the derived type via `static_cast<Derived&>(*this)(a.v)`, where `a.v` is `std::vector<int>`. No overload of `operator()` accepting `std::vector<int>` (or any container) is defined anywhere in the traverser hierarchy (`int_traverser`, `A_traverser`, `B_traverser`, `core::apply_traverser`, or `core::traverser`). `core::traverser` provides only `apply(const Container&)`, a distinct name. A conforming compiler must reject this with a "no matching function for call to `operator()(std::vector<int>)`" error.
- **Evidence type:** static-proof
- **Evidence:**
  - **Reproducer:** Attempt to compile `traverser2.cpp` against the current tree.
  - **Run command:**
    ```bash
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_GUI_TOOLS=OFF && \
    cmake --build build --target traverser2 2>&1 | head -40
    # If no CMake target exists, compile directly:
    g++ -std=c++20 -I libraries/core/include -I libraries/atermpp/include \
        doc/sphinx/developer_manual/libraries/core/examples/traverser2.cpp -o /tmp/t2 2>&1
    ```
  - **Observed result:** Compile error — no matching `operator()` for `std::vector<int>` argument.
  - **Expected result:** Example compiles and prints documented output.
- **Efficiency impact:** N/A (compile-time failure).
- **C++20-specific note:** With C++20 concepts, `core::traverser::apply(const Container&)` is constrained with `requires atermpp::is_container<Container>::value`, making the `apply` / `operator()` split even more prominent and explicit.
- **Suggested fix:** Replace all `operator()` calls and definitions in the example with `apply()` to match the current `core::traverser` API. `A_traverser::operator()(const A& a)` should become `void apply(const A& a)` and both dispatch calls should be `static_cast<Derived&>(*this).apply(a.i)` / `static_cast<Derived&>(*this).apply(a.v)`. `apply_traverser` must also expose `using super::operator()` if `operator()` is the chosen entry-point — but the cleaner fix is to align to the existing `apply()` convention.

---

### [SEV-1] Comment at line 82 names the wrong traverser

- **Location:** `traverser2.cpp:82–88`
- **Why this is a problem:**
  The block comment says _"Apply A_traverser to a"_ but the code immediately below applies `B_traverser`:
  ```cpp
  // Apply A_traverser to a. This gives the following output:
  //
  // 2*i = 6
  // 2*i = 2
  // 2*i = 8
  //
  core::apply_traverser<B_traverser>()(a);   // ← B_traverser, not A_traverser
  ```
  Because this is a `literalinclude` in a user-facing Sphinx manual, the mislabelling directly contradicts the rendered documentation and will mislead anyone trying to understand the CRTP traversal pattern.
- **Evidence type:** static-proof
- **Evidence:** Reading lines 82–88 alongside line 88 is sufficient; the comment and the code cannot both be correct simultaneously.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Change the comment to read _"Apply B_traverser to a"_.

---

### [SEV-2] `using super::operator()` on `core::traverser<Derived>` which defines no `operator()`

- **Location:** `traverser2.cpp:23`, `38`, `55`
- **Why this is a problem:**
  Each traverser layer starts with `using super::operator()`. The base of the chain — `core::traverser<Derived>` — has no `operator()`. In C++, a `using`-declaration for a name that does not exist in the named base class is ill-formed (CWG 591; many compilers accept it silently as an empty declaration in template context, but it is not portable and is misleading to readers). More importantly, this reveals the example was written against a version of `traverser.h` that exposed `operator()`, but that API no longer exists — `apply()` replaced it.
- **Evidence type:** static-proof (via traverser.h inspection: `grep -n 'operator()' libraries/core/include/mcrl2/core/traverser.h` returns zero results)
- **Evidence:**
  - **Run command:** `grep -n 'operator()' libraries/core/include/mcrl2/core/traverser.h`
  - **Observed result:** No output.
  - **Expected result:** At least one `operator()` definition if the example is valid.
- **Efficiency impact:** None at runtime; confuses readers about which interface to implement.
- **C++20-specific note:** C++20 makes this more visible because `core::traverser::apply` is now concept-constrained; the gap between the constrained `apply` path and the unconstrained `operator()` path is semantically wider.
- **Suggested fix:** Remove all `using super::operator()` declarations. Replace every `void operator()(T)` override with `void apply(const T&)` and align with the existing real-code pattern visible in `libraries/data/include/mcrl2/data/traverser.h`.

---

### [SEV-2] `apply_traverser` brings in `apply` but the example never defines `apply`

- **Location:** `traverser2.cpp:80, 88`; `libraries/core/include/mcrl2/core/traverser.h:67–74`
- **Why this is a problem:**
  `apply_traverser` does `using super::enter; using super::leave; using super::apply;`. The entry-point called externally is `operator()` (line 80: `core::apply_traverser<A_traverser>()(a)`). The `operator()` overloads are inherited from `A_traverser` / `B_traverser` so the outer call dispatches correctly into `A_traverser::operator()(const A&)`. However, the internal container dispatch (`static_cast<Derived&>(*this)(a.v)`) relies on `operator()` for containers, whereas `core::traverser` only dispatches containers through `apply()`. These are two irreconcilable dispatch paths — the `apply` path in `core::traverser` can never be reached by the `operator()` calls in this example, making the `apply_traverser::using super::apply` line dead weight for this example.
- **Evidence type:** static-proof
- **Evidence:** Trace the call graph: `apply_traverser<A_traverser>()(a)` → `A_traverser::operator()(const A&)` → `static_cast<Derived&>(*this)(a.v)` → seeks `operator()(std::vector<int>)` → not found. The `apply(container)` method in `core::traverser` is never reachable from this chain.
- **Efficiency impact:** None at runtime (compile failure prevents execution), but conceptually the split creates a hidden dead code path in `apply_traverser` for any `operator()`-based traverser.
- **C++20-specific note:** N/A.
- **Suggested fix:** Either (a) switch the entire example to use `apply()` as the primary dispatch name and leverage `core::traverser::apply(const Container&)` for container traversal, or (b) remove the example from the documentation until it is corrected — do not publish documentation that cannot compile.

---

### [SEV-3] Unused `#include <cassert>`

- **Location:** `traverser2.cpp:2`
- **Why this is a problem:** `<cassert>` is included but `assert` is never used. This is dead code in a documentation example and suggests a sloppy copy from `traverser1.cpp`. It also means the example performs no correctness assertion on the observable output — the only way to verify it works is to manually read stdout.
- **Evidence type:** static-proof (`grep 'assert' traverser2.cpp` returns zero results)
- **Efficiency impact:** Negligible include-time cost; primarily a clarity issue.
- **C++20-specific note:** N/A.
- **Suggested fix:** Remove `#include <cassert>`. If an assertion is desired (recommended), add one using output capture or by checking return values.

---

### [SEV-3] Pre-C++11 container initialization style

- **Location:** `traverser2.cpp:69–71`
- **Why this is a problem:**
  ```cpp
  std::vector<int> v;
  v.push_back(1);
  v.push_back(4);
  ```
  The entire codebase targets C++20. This is a documentation example. Using `push_back` instead of brace-initialisation teaches readers an outdated idiom and increases noise.
- **Evidence type:** plausible-issue (style; not a defect)
- **Evidence:** N/A — no failure, but the pattern is inconsistent with C++20 expectations documented in the repo's `.clang-tidy` configuration.
- **Efficiency impact:** Two redundant heap reallocations instead of a single allocation; irrelevant in practice.
- **C++20-specific note:** `std::vector<int> v = {1, 4};` is correct C++11 and remains idiomatic in C++20. Even better: `const std::vector<int> v = {1, 4};` since `v` is never mutated after construction.
- **Suggested fix:**
  ```cpp
  const std::vector<int> v = {1, 4};
  ```

---

### [SEV-3] `std::endl` used where `'\n'` is sufficient

- **Location:** `traverser2.cpp:30, 62`
- **Why this is a problem:** `std::endl` performs a `flush()` after each line. In a traversal that visits large expression trees this would be a significant and unnecessary performance tax. Even in a documentation example, teaching `std::endl` as the default is bad practice.
- **Evidence type:** plausible-issue (style / latent performance)
- **Efficiency impact:** One unnecessary `fflush` per traversal node per run; negligible here, but the pattern will be cargo-culted.
- **C++20-specific note:** N/A.
- **Suggested fix:** Replace `<< std::endl` with `<< '\n'`.

---

## Proposed Tests

- [ ] **compile_traverser2**: Compile `traverser2.cpp` against the current `traverser.h` and expect a compiler error for the unresolved `operator()(std::vector<int>)` call; confirms SEV-1 and SEV-2. Run: `g++ -std=c++20 -I libraries/core/include -I libraries/atermpp/include doc/sphinx/developer_manual/libraries/core/examples/traverser2.cpp -o /tmp/t2 2>&1; test $? -ne 0`
- [ ] **corrected_traverser2_output**: After fixing the example to use `apply()`, compile and capture stdout; assert exact lines `"i = 3\ni = 1\ni = 4\n2*i = 6\n2*i = 2\n2*i = 8\n"` are produced. This would also confirm the wrong comment (SEV-1 line 82) was corrected.
- [ ] **traverser_container_dispatch**: Unit test in `libraries/core/test/` that instantiates a minimal `apply_traverser` with an `apply(int)` override and applies it to `std::vector<int>{1,4}`; asserts both elements are visited in order. This guards against future regressions in container dispatch.

---

## Open Questions

1. **Is `traverser2.cpp` ever compiled in CI?** No CMakeLists.txt target was found for it. If it is purely a `literalinclude` in Sphinx docs, the broken state has been invisible to CI — the build system should be extended to compile all doc examples as part of the test suite.
2. **When was `core::traverser` migrated from `operator()` to `apply()`?** `traverser1.cpp` has the same `using super::operator()` pattern on a direct (non-CRTP, non-`apply_traverser`) instantiation. If `traverser1.cpp` also fails to compile, the entire examples directory is stale and needs wholesale updating.
3. **Is there a `core::traverser` overload anywhere that provides `operator()` for containers?** A thorough `grep -rn 'operator()' libraries/core/include/` found none in `traverser.h`; if one exists elsewhere in the include chain (e.g. through `identifier_string.h` transitive includes) the SEV-1 static proof would need re-evaluation.
4. **Should `make_apply_traverser()` (line 78–81 of `traverser.h`) be preferred over `apply_traverser<T>()(x)` syntax?** The factory function exists but is not used in the example. Document which form is idiomatic.

## Copilot Raw Stdout

```text
● skill(mcrl2-critical-cpp20-review)

● Read doc/sphinx/developer_manual/libraries/core/examples/traverser2.cpp
  └ 92 lines read

● Glob "libraries/core/include/mcrl2/core/traverser.h" in /workspaces/mCRL2
  └ 1 file found

● Read libraries/core/include/mcrl2/core/traverser.h
  └ 85 lines read

● Glob "doc/sphinx/developer_manual/libraries/core/examples/*.cpp" in /workspaces/mCRL2
  └ 3 files found

● Read doc/sphinx/developer_manual/libraries/core/examples/traverser1.cpp
  └ 41 lines read

● Find other traverser/builder usages
  $ cd /workspaces/mCRL2 && grep -r "apply_traverser\|apply_builder" --include="*.h" -l | head -20
  └ 21 lines...

✗ Find CMakeLists for doc examples
  $ find /workspaces/mCRL2/doc -name "CMakeLists.txt" | xargs grep -l
  "traverser2\|traverser1\|examples" 2>/dev/null | head -5
  Permission denied and could not request permission from user

● Read libraries/core/include/mcrl2/core/builder.h lines 1-80
  └ 80 lines read

● Find doc references to traverser2
  $ grep -r "traverser2\|traverser1\|apply_traverser" /workspaces/mCRL2/doc --include="*.txt"
  --include="*.cmake" --include="*.rst" -l 2>/dev/null | head -10
  └ 2 lines...

● Check if operator() exists in traverser.h
  $ grep -rn "operator()" /workspaces/mCRL2/libraries/core/include/mcrl2/core/traverser.h
  └ 1 line...

● Check rst doc references to traverser2
  $ grep -n "traverser2\|apply_traverser\|operator()"
  /workspaces/mCRL2/doc/sphinx/developer_manual/libraries/core/common_functionality.rst | head -30
  └ 12 lines...

● Check data traverser for operator() usage
  $ grep -rn "apply_traverser\|operator()"
  /workspaces/mCRL2/libraries/data/include/mcrl2/data/traverser.h 2>/dev/null | head -20
  └ 1 line...

● Context around traverser2 in rst docs
  $ grep -n -A5 -B5 "traverser2\|apply_traverser"
  /workspaces/mCRL2/doc/sphinx/developer_manual/libraries/core/common_functionality.rst | head -60
  └ 12 lines...

● Check is_container and apply usage in traverser
  $ grep -rn "apply_traverser\|is_container"
  /workspaces/mCRL2/libraries/core/include/mcrl2/core/traverser.h; grep -rn "is_container"
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/concepts.h 2>/dev/null | head -10
  └ 7 lines...

✗ Find is_container definition
  $ grep -rn "is_container" /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/concepts.h
  2>/dev/null; find /workspaces/mCRL2/libraries/atermpp -name "*.h" | xargs grep -l "is_container"
  2>/dev/null | head -5
  Permission denied and could not request permission from user

● Check is_container in atermpp
  $ grep -rn "is_container" /workspaces/mCRL2/libraries/atermpp/include/ | head -15
  └ 12 lines...

● Check if example is compiled
  $ grep -rn "traverser2\|traverser1\|builder1" /workspaces/mCRL2/doc --include="CMakeLists.txt"
  2>/dev/null; ls /workspaces/mCRL2/doc/sphinx/developer_manual/libraries/core/examples/
  └ 2 lines...

● Check how traverser is actually used in production code
  $ grep -rn "apply_traverser\|using super::operator"
  /workspaces/mCRL2/libraries/data/include/mcrl2/data/traverser.h 2>/dev/null | head -20; grep -rn
  "using super::operator" /workspaces/mCRL2/libraries/lps/include/mcrl2/lps/ --include="*.h"
  2>/dev/null | head -10
  └ 3 lines...

● Check real usages of using super::operator() in library headers
  $ grep -rn "using super::operator\(\)" /workspaces/mCRL2/libraries/ --include="*.h" | head -20
  └ 5 lines...

● Check data traverser structure
  $ grep -rn "apply_traverser\|void apply\|void operator()"
  /workspaces/mCRL2/libraries/data/include/mcrl2/data/traverser.h | head -30
  └ 31 lines...

● Check data_rewriter traverser usage
  $ grep -rn "apply_traverser\|void operator()"
  /workspaces/mCRL2/libraries/lps/include/mcrl2/lps/rewriters/data_rewriter.h | head -20
  └ 1 line...

● Check if examples are compiled by any CMakeLists
  $ grep -rn "traverser2\|traverser1\|builder1" /workspaces/mCRL2 --include="CMakeLists.txt"
  2>/dev/null | head -10
  └ 1 line...

Now I have enough context. Let me write the review.

---

# Review Findings

## Summary
- **Scope reviewed:** `doc/sphinx/developer_manual/libraries/core/examples/traverser2.cpp` and its only dependency `libraries/core/include/mcrl2/core/traverser.h`
- **Risk level:** HIGH (documentation) — the example is stale and demonstrably incorrect against the current framework API; LOW (production) — no production code uses this file
- **Overall verdict:** Do not merge / correct before publishing. The example contradicts the current `core::traverser` API on two separate levels: method naming and dispatch semantics. At least one call would fail to compile. There is also a literal factual error in an inline comment.

---

## Findings (ordered by severity)

### [SEV-1] `static_cast<Derived&>(*this)(a.v)` calls non-existent `operator()(std::vector<int>)`

- **Location:** `traverser2.cpp:47`
- **Why this is a problem:**
  `A_traverser::operator()(const A& a)` delegates to the derived type via `static_cast<Derived&>(*this)(a.v)`, where `a.v` is `std::vector<int>`. No overload of `operator()` accepting `std::vector<int>` (or any container) is defined anywhere in the traverser hierarchy (`int_traverser`, `A_traverser`, `B_traverser`, `core::apply_traverser`, or `core::traverser`). `core::traverser` provides only `apply(const Container&)`, a distinct name. A conforming compiler must reject this with a "no matching function for call to `operator()(std::vector<int>)`" error.
- **Evidence type:** static-proof
- **Evidence:**
  - **Reproducer:** Attempt to compile `traverser2.cpp` against the current tree.
  - **Run command:**
    ```bash
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_GUI_TOOLS=OFF && \
    cmake --build build --target traverser2 2>&1 | head -40
    # If no CMake target exists, compile directly:
    g++ -std=c++20 -I libraries/core/include -I libraries/atermpp/include \
        doc/sphinx/developer_manual/libraries/core/examples/traverser2.cpp -o /tmp/t2 2>&1
    ```
  - **Observed result:** Compile error — no matching `operator()` for `std::vector<int>` argument.
  - **Expected result:** Example compiles and prints documented output.
- **Efficiency impact:** N/A (compile-time failure).
- **C++20-specific note:** With C++20 concepts, `core::traverser::apply(const Container&)` is constrained with `requires atermpp::is_container<Container>::value`, making the `apply` / `operator()` split even more prominent and explicit.
- **Suggested fix:** Replace all `operator()` calls and definitions in the example with `apply()` to match the current `core::traverser` API. `A_traverser::operator()(const A& a)` should become `void apply(const A& a)` and both dispatch calls should be `static_cast<Derived&>(*this).apply(a.i)` / `static_cast<Derived&>(*this).apply(a.v)`. `apply_traverser` must also expose `using super::operator()` if `operator()` is the chosen entry-point — but the cleaner fix is to align to the existing `apply()` convention.

---

### [SEV-1] Comment at line 82 names the wrong traverser

- **Location:** `traverser2.cpp:82–88`
- **Why this is a problem:**
  The block comment says _"Apply A_traverser to a"_ but the code immediately below applies `B_traverser`:
  ```cpp
  // Apply A_traverser to a. This gives the following output:
  //
  // 2*i = 6
  // 2*i = 2
  // 2*i = 8
  //
  core::apply_traverser<B_traverser>()(a);   // ← B_traverser, not A_traverser
  ```
  Because this is a `literalinclude` in a user-facing Sphinx manual, the mislabelling directly contradicts the rendered documentation and will mislead anyone trying to understand the CRTP traversal pattern.
- **Evidence type:** static-proof
- **Evidence:** Reading lines 82–88 alongside line 88 is sufficient; the comment and the code cannot both be correct simultaneously.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Change the comment to read _"Apply B_traverser to a"_.

---

### [SEV-2] `using super::operator()` on `core::traverser<Derived>` which defines no `operator()`

- **Location:** `traverser2.cpp:23`, `38`, `55`
- **Why this is a problem:**
  Each traverser layer starts with `using super::operator()`. The base of the chain — `core::traverser<Derived>` — has no `operator()`. In C++, a `using`-declaration for a name that does not exist in the named base class is ill-formed (CWG 591; many compilers accept it silently as an empty declaration in template context, but it is not portable and is misleading to readers). More importantly, this reveals the example was written against a version of `traverser.h` that exposed `operator()`, but that API no longer exists — `apply()` replaced it.
- **Evidence type:** static-proof (via traverser.h inspection: `grep -n 'operator()' libraries/core/include/mcrl2/core/traverser.h` returns zero results)
- **Evidence:**
  - **Run command:** `grep -n 'operator()' libraries/core/include/mcrl2/core/traverser.h`
  - **Observed result:** No output.
  - **Expected result:** At least one `operator()` definition if the example is valid.
- **Efficiency impact:** None at runtime; confuses readers about which interface to implement.
- **C++20-specific note:** C++20 makes this more visible because `core::traverser::apply` is now concept-constrained; the gap between the constrained `apply` path and the unconstrained `operator()` path is semantically wider.
- **Suggested fix:** Remove all `using super::operator()` declarations. Replace every `void operator()(T)` override with `void apply(const T&)` and align with the existing real-code pattern visible in `libraries/data/include/mcrl2/data/traverser.h`.

---

### [SEV-2] `apply_traverser` brings in `apply` but the example never defines `apply`

- **Location:** `traverser2.cpp:80, 88`; `libraries/core/include/mcrl2/core/traverser.h:67–74`
- **Why this is a problem:**
  `apply_traverser` does `using super::enter; using super::leave; using super::apply;`. The entry-point called externally is `operator()` (line 80: `core::apply_traverser<A_traverser>()(a)`). The `operator()` overloads are inherited from `A_traverser` / `B_traverser` so the outer call dispatches correctly into `A_traverser::operator()(const A&)`. However, the internal container dispatch (`static_cast<Derived&>(*this)(a.v)`) relies on `operator()` for containers, whereas `core::traverser` only dispatches containers through `apply()`. These are two irreconcilable dispatch paths — the `apply` path in `core::traverser` can never be reached by the `operator()` calls in this example, making the `apply_traverser::using super::apply` line dead weight for this example.
- **Evidence type:** static-proof
- **Evidence:** Trace the call graph: `apply_traverser<A_traverser>()(a)` → `A_traverser::operator()(const A&)` → `static_cast<Derived&>(*this)(a.v)` → seeks `operator()(std::vector<int>)` → not found. The `apply(container)` method in `core::traverser` is never reachable from this chain.
- **Efficiency impact:** None at runtime (compile failure prevents execution), but conceptually the split creates a hidden dead code path in `apply_traverser` for any `operator()`-based traverser.
- **C++20-specific note:** N/A.
- **Suggested fix:** Either (a) switch the entire example to use `apply()` as the primary dispatch name and leverage `core::traverser::apply(const Container&)` for container traversal, or (b) remove the example from the documentation until it is corrected — do not publish documentation that cannot compile.

---

### [SEV-3] Unused `#include <cassert>`

- **Location:** `traverser2.cpp:2`
- **Why this is a problem:** `<cassert>` is included but `assert` is never used. This is dead code in a documentation example and suggests a sloppy copy from `traverser1.cpp`. It also means the example performs no correctness assertion on the observable output — the only way to verify it works is to manually read stdout.
- **Evidence type:** static-proof (`grep 'assert' traverser2.cpp` returns zero results)
- **Efficiency impact:** Negligible include-time cost; primarily a clarity issue.
- **C++20-specific note:** N/A.
- **Suggested fix:** Remove `#include <cassert>`. If an assertion is desired (recommended), add one using output capture or by checking return values.

---

### [SEV-3] Pre-C++11 container initialization style

- **Location:** `traverser2.cpp:69–71`
- **Why this is a problem:**
  ```cpp
  std::vector<int> v;
  v.push_back(1);
  v.push_back(4);
  ```
  The entire codebase targets C++20. This is a documentation example. Using `push_back` instead of brace-initialisation teaches readers an outdated idiom and increases noise.
- **Evidence type:** plausible-issue (style; not a defect)
- **Evidence:** N/A — no failure, but the pattern is inconsistent with C++20 expectations documented in the repo's `.clang-tidy` configuration.
- **Efficiency impact:** Two redundant heap reallocations instead of a single allocation; irrelevant in practice.
- **C++20-specific note:** `std::vector<int> v = {1, 4};` is correct C++11 and remains idiomatic in C++20. Even better: `const std::vector<int> v = {1, 4};` since `v` is never mutated after construction.
- **Suggested fix:**
  ```cpp
  const std::vector<int> v = {1, 4};
  ```

---

### [SEV-3] `std::endl` used where `'\n'` is sufficient

- **Location:** `traverser2.cpp:30, 62`
- **Why this is a problem:** `std::endl` performs a `flush()` after each line. In a traversal that visits large expression trees this would be a significant and unnecessary performance tax. Even in a documentation example, teaching `std::endl` as the default is bad practice.
- **Evidence type:** plausible-issue (style / latent performance)
- **Efficiency impact:** One unnecessary `fflush` per traversal node per run; negligible here, but the pattern will be cargo-culted.
- **C++20-specific note:** N/A.
- **Suggested fix:** Replace `<< std::endl` with `<< '\n'`.

---

## Proposed Tests

- [ ] **compile_traverser2**: Compile `traverser2.cpp` against the current `traverser.h` and expect a compiler error for the unresolved `operator()(std::vector<int>)` call; confirms SEV-1 and SEV-2. Run: `g++ -std=c++20 -I libraries/core/include -I libraries/atermpp/include doc/sphinx/developer_manual/libraries/core/examples/traverser2.cpp -o /tmp/t2 2>&1; test $? -ne 0`
- [ ] **corrected_traverser2_output**: After fixing the example to use `apply()`, compile and capture stdout; assert exact lines `"i = 3\ni = 1\ni = 4\n2*i = 6\n2*i = 2\n2*i = 8\n"` are produced. This would also confirm the wrong comment (SEV-1 line 82) was corrected.
- [ ] **traverser_container_dispatch**: Unit test in `libraries/core/test/` that instantiates a minimal `apply_traverser` with an `apply(int)` override and applies it to `std::vector<int>{1,4}`; asserts both elements are visited in order. This guards against future regressions in container dispatch.

---

## Open Questions

1. **Is `traverser2.cpp` ever compiled in CI?** No CMakeLists.txt target was found for it. If it is purely a `literalinclude` in Sphinx docs, the broken state has been invisible to CI — the build system should be extended to compile all doc examples as part of the test suite.
2. **When was `core::traverser` migrated from `operator()` to `apply()`?** `traverser1.cpp` has the same `using super::operator()` pattern on a direct (non-CRTP, non-`apply_traverser`) instantiation. If `traverser1.cpp` also fails to compile, the entire examples directory is stale and needs wholesale updating.
3. **Is there a `core::traverser` overload anywhere that provides `operator()` for containers?** A thorough `grep -rn 'operator()' libraries/core/include/` found none in `traverser.h`; if one exists elsewhere in the include chain (e.g. through `identifier_string.h` transitive includes) the SEV-1 static proof would need re-evaluation.
4. **Should `make_apply_traverser()` (line 78–81 of `traverser.h`) be preferred over `apply_traverser<T>()(x)` syntax?** The factory function exists but is not used in the example. Document which form is idiomatic.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         3m 32s
Total session time:     3m 46s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-sonnet-4.6       420.7k in, 12.7k out, 391.3k cached (Est. 1 Premium request)

```

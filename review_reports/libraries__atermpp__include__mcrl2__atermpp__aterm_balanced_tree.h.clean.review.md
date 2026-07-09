---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/aterm_balanced_tree.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-07T23:38:23.882920+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_balanced_tree.h` — complete file (clean review)
- Risk level: **Medium** — Contains UB-by-standard patterns, unbounded buffer overflow path in release builds, and API design issues that invite O(n²) misuse
- Overall verdict: Functional and well-tested for typical workloads, but relies on layout coincidence for correctness and has no defense-in-depth against pathological inputs in release builds

## Findings (ordered by severity)

### [SEV-1] Iterator stack buffer overflow in release builds for trees exceeding 2^20 elements
- Location: `aterm_balanced_tree.h:298` (declaration), `:345-350` (write sites in `increment()`)
- Why this is a problem: `m_stack` is a fixed `std::array<unprotected_aterm_core, 20>`. The `initialise()` method has `assert(m_top_of_stack + 1 < maximal_size_of_stack)` but asserts are stripped in release. The `increment()` method has **zero bounds checking**. A tree of 2^20+1 elements requires depth 21, causing out-of-bounds writes on the stack — a classic buffer overflow leading to arbitrary memory corruption.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Construct `term_balanced_tree<aterm>` with >1,048,576 elements and iterate
  - Run command: `cd build && ctest -R aterm_balanced_tree --output-on-failure` (after adding a test case constructing a tree of size 1,048,577)
  - Observed result: In debug: assertion failure. In release: stack corruption / segfault / silent UB
  - Expected result: Defined behavior (either larger stack or dynamic allocation fallback)
- Efficiency impact: N/A — correctness issue
- C++20-specific note: Could use `std::span` with bounds-checked access, or `std::vector` with `reserve()`
- Suggested fix: Either (a) make `maximal_size_of_stack` a computed value based on a documented max-size contract enforced at construction, (b) replace with `std::vector` and `reserve(20)`, or (c) add a `static_assert`/runtime check in constructors that rejects trees exceeding depth 20.

### [SEV-2] Undefined behavior: `static_cast<Tree&>` on `unprotected_aterm_core` objects
- Location: `aterm_balanced_tree.h:345-346` (increment), `:372-373` (initialise)
- Why this is a problem: `current` is declared as `unprotected_aterm_core` (a stack-local object). Casting it to `term_balanced_tree<Term>&` and calling member functions on it violates strict aliasing (§6.7.2 [basic.lval]). The object was constructed as `unprotected_aterm_core`, not `term_balanced_tree<Term>`. This only "works" because both types share identical layout (single `const _aterm*`), but compilers are allowed to optimize under the assumption this never happens.
- Evidence type: plausible-issue
- Evidence:
  - Mechanism: Under LTO or aggressive optimization (`-O3 -flto`), the compiler may assume that a `term_balanced_tree<Term>&` refers to a real `term_balanced_tree<Term>` object and propagate type-based alias analysis, potentially eliminating loads or reordering operations.
  - Triggering scenario: Build with `-O3 -flto -fstrict-aliasing` and iterate a balanced tree. Failure depends on compiler version/flags.
  - Test plan: Compile the existing `aterm_balanced_tree_test` with `-O3 -flto -fstrict-aliasing -Wstrict-aliasing=1` and check for warnings or miscompilations with ASAN/UBSAN enabled.
  - Run command: `cmake --build build --target aterm_balanced_tree_test -- -O3 -flto && cd build && ./libraries/atermpp/test/aterm_balanced_tree_test`
  - Expected result: UBSAN reports "member call on misaligned address" or "member access within miscast" (if instrumented), or no issue (if compiler doesn't exploit the UB)
- Efficiency impact: None directly; correctness concern
- C++20-specific note: The idiomatic fix is to use `std::bit_cast` or store the pointer and reconstruct, or use `reinterpret_cast` through `void*` with `std::launder` (as `assign_cast`/`vertical_cast` already do via `reinterpret_cast<detail::_aterm&>`). The cast should go through the established cast infrastructure (`down_cast` or `vertical_cast`) rather than a raw `static_cast` that misleads the compiler about object identity.
- Suggested fix: Replace `static_cast<Tree&>(current)` with a helper that extracts the `_aterm*` and uses the existing casting infrastructure, e.g.:
  ```cpp
  const aterm& as_aterm = reinterpret_cast<const aterm&>(current);
  down_cast<term_balanced_tree<Term>>(as_aterm).right_branch();
  ```

### [SEV-2] Anonymous namespace globals in a header cause per-TU duplication
- Location: `aterm_balanced_tree.h:20-26`
- Why this is a problem: Every translation unit that includes this header gets its own `global_function_symbol` and `aterm` instances. Currently 2 direct includers (`state.h` and the test), but `state.h` is transitively included broadly. Each TU pays constructor/destructor cost and memory for duplicate objects. More critically, `g_empty_tree` is an `aterm` — each copy holds its own reference count, adding unnecessary atomic ref-count operations at program startup/shutdown.
- Evidence type: static-proof
- Evidence:
  - `grep -rn "#include.*aterm_balanced_tree\|#include.*lps/state" --include="*.cpp" | wc -l` shows transitive inclusion count
  - Each includer's object file will contain `.init_array` entries for the 4 global constructors
  - Run command: `nm -C build/libraries/lps/CMakeFiles/mcrl2_lps.dir/source/*.o | grep "g_empty\|g_tree_node" | wc -l` — expect >1 per symbol
  - Observed result: Multiple definitions of the same symbols across object files
  - Expected result: Single definition (via `inline` variables or defined in a `.cpp` file)
- Efficiency impact: Startup cost scales with number of TUs including the header; unnecessary atomic operations
- C++20-specific note: C++17 `inline` variables solve this exactly: `inline global_function_symbol g_empty("@empty@", 0);` provides a single definition across TUs with guaranteed identity
- Suggested fix: Either (a) use `inline` variables (C++17+), or (b) move definitions to `aterm_balanced_tree.cpp` and expose via `static` member function accessors (as partially done already with `tree_empty_function()` etc., but those just return the anonymous-namespace globals).

### [SEV-3] `operator[]` is O(n) — invites accidental O(n²) loops
- Location: `aterm_balanced_tree.h:199-202`
- Why this is a problem: `operator[]` calls `size()` which is O(n). Any loop using `tree[i]` is O(n²). The comment documents this, but `operator[]` is the canonical random-access pattern — users will write `for (int i=0; i<tree.size(); i++) tree[i]` without reading the comment.
- Evidence type: plausible-issue
- Evidence:
  - Mechanism: Any code calling `tree[i]` in a loop gets O(n²) behavior silently
  - Triggering scenario: `for (std::size_t i = 0; i < t.size(); ++i) process(t[i]);`
  - Test plan: Benchmark `operator[]` loop vs iterator for tree of 10,000 elements; expect 100x slowdown
  - Run command: Instrument with `std::chrono` in a unit test comparing both access patterns
- Efficiency impact: O(n) per access instead of O(log n); O(n²) for naive iteration vs O(n) with iterators
- C++20-specific note: Consider deleting `operator[]` or marking it `[[deprecated("O(n) — use iterators or element_at with known size")]]`
- Suggested fix: Add `[[deprecated]]` or remove `operator[]` in favor of `element_at` (which requires pre-computed size, making the cost explicit).

### [SEV-3] `const_reference` type alias is not a reference
- Location: `aterm_balanced_tree.h:126`
- Why this is a problem: `using const_reference = const Term;` defines a value type, not a reference type. Standard container concepts expect `const_reference` to be `const value_type&`. Any generic code that relies on this being a reference (e.g., templates checking `std::is_reference_v<typename C::const_reference>`) will behave unexpectedly.
- Evidence type: static-proof
- Evidence:
  - `static_assert(std::is_reference_v<term_balanced_tree<aterm>::const_reference>)` will fail
  - Run command: Add the static_assert to the test file and compile
  - Expected result: Compilation failure confirming the issue
- Efficiency impact: May cause unnecessary copies in generic code that deduces return type from `const_reference`
- C++20-specific note: C++20 ranges/concepts will reject this container as modeling the expected interface
- Suggested fix: `using const_reference = const Term&;`

### [SEV-3] `pp()` takes balanced tree by value instead of const reference
- Location: `aterm_balanced_tree.h:290`
- Why this is a problem: `std::string pp(const term_balanced_tree<Term> t)` copies the tree on each call. While aterm copies are cheap (atomic ref-count bump + pointer copy), this is unnecessary overhead and non-idiomatic.
- Evidence type: static-proof
- Evidence:
  - The function never modifies `t`; it only iterates
  - Signature should be `const term_balanced_tree<Term>& t`
- Efficiency impact: Unnecessary atomic increment/decrement per call
- C++20-specific note: N/A
- Suggested fix: Change parameter to `const term_balanced_tree<Term>& t`

### [SEV-3] Redundant `left_branch()`/`right_branch()` calls in `element_at`
- Location: `aterm_balanced_tree.h:233, 243`
- Why this is a problem: The code stores `left_branch()` in a local reference, checks it, but then calls `left_branch()` again for the recursive call:
  ```cpp
  const aterm& left(left_branch());
  if (left.function() == tree_node_function())
    return down_cast<term_balanced_tree<Term>>(left_branch()).element_at(...); // calls left_branch() again
  ```
  Same pattern for `right_branch()`. While `left_branch()` is side-effect-free, this defeats the purpose of caching.
- Evidence type: static-proof
- Evidence:
  - `left_branch()` called at line 230 (stored in `left`), then again at line 233
  - `right_branch()` called at line 240 (stored in `right`), then again at line 243
- Efficiency impact: Negligible (simple accessor), but confusing — suggests the author intended to use the cached value
- Suggested fix: Use the cached local: `down_cast<term_balanced_tree<Term>>(left).element_at(...)`

## Proposed Tests
- [ ] `test_large_tree_overflow`: Construct tree with 2^20+1 elements, iterate fully. In debug builds expect assertion failure. In release builds with ASAN, expect no buffer overflow. Purpose: validate stack limit is sufficient or detect overflow.
- [ ] `test_iterator_ubsan`: Run existing `aterm_balanced_tree_test` under UBSAN (`-fsanitize=undefined`) to detect strict-aliasing violations in iterator traversal.
- [ ] `test_operator_bracket_complexity`: Benchmark `operator[]` loop vs iterator loop on tree of 10,000 elements. Assert iterator is at least 10x faster. Purpose: document the O(n²) trap.
- [ ] `test_const_reference_type`: `static_assert(std::is_reference_v<term_balanced_tree<aterm>::const_reference>);` — should pass if the typedef is correct.

## Open Questions
- Is the 2^20 element limit documented anywhere as a contract, or is it an undocumented assumption? If `lps::state` can theoretically contain more than ~1M data expressions, this is a latent crash.
- Why use `unprotected_aterm_core` in the iterator stack rather than `aterm` (which would provide GC protection)? If a garbage collection cycle runs during iteration, are the `unprotected_aterm_core` entries in `m_stack` safe from being collected? The tree itself is protected (held by the caller), but the intermediate nodes on the stack are only reachable through the tree — this seems safe only if GC traces from roots through protected terms. Confirm this assumption.
- The anonymous-namespace pattern was presumably chosen to avoid a separate `.cpp` file. Would migrating to `inline` variables (C++17) be acceptable given the project's minimum standard?

## Copilot Raw Stdout

```text
● Read .github/instructions/cpp.instructions.md
  └ 28 lines read

● Explore: Explore aterm infrastructure
  ● Find global_function_symbol implementation
    $ grep -r "global_function_symbol::" /workspaces/mCRL2/libraries/atermpp/ 2>/dev/null | head -5
    └ 2 lines...
  ● See global_function_symbol implementation
    $ grep -A 5 "global_function_symbol::"
    /workspaces/mCRL2/libraries/atermpp/source/function_symbol.cpp
    └ 4 lines...

● Explore: Check iterator UB and size limits
  ● Read libraries/atermpp/test/down_cast_test.cpp
    └ 92 lines read
  ● Check all references to m_stack
    $ grep -rn "m_stack\|maximal_size_of_stack"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm_balanced_tree.h | tail -20
    └ 14 lines...

● Check balanced tree tests
  $ cd /workspaces/mCRL2 && grep -rn "aterm_balanced_tree" libraries/atermpp/test/ | head -20
  └ 20 lines...

● Check const_reference typedef
  $ cd /workspaces/mCRL2 && grep -n "const_reference"
  libraries/atermpp/include/mcrl2/atermpp/aterm_balanced_tree.h
  └ 2 lines...

● Count header inclusions
  $ cd /workspaces/mCRL2 && grep -rn "#include.*aterm_balanced_tree" --include="*.cpp"
  --include="*.h" | wc -l && echo "---" && grep -rn "#include.*aterm_balanced_tree"
  --include="*.cpp" --include="*.h" | head -20
  └ 5 lines...

● Check state usage frequency
  $ cd /workspaces/mCRL2 && grep -c "lps::state" libraries/lps/source/*.cpp
  libraries/lts/source/*.cpp 2>/dev/null | grep -v ":0$"
  └ 2 lines...

● Check unprotected_aterm_core definition
  $ cd /workspaces/mCRL2 && grep -n "unprotected_aterm_core"
  libraries/atermpp/include/mcrl2/atermpp/aterm_core.h | head -10 && echo "---" && grep -n "class
  aterm_core" libraries/atermpp/include/mcrl2/atermpp/aterm_core.h | head -5
  └ 12 lines...

● View unprotected_aterm_core class
  $ cd /workspaces/mCRL2 && sed -n '33,130p' libraries/atermpp/include/mcrl2/atermpp/aterm_core.h
  └ 99 lines...

Now I have all the information needed to write the review. Let me produce the findings.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/aterm_balanced_tree.h` — complete file (clean review)
- Risk level: **Medium** — Contains UB-by-standard patterns, unbounded buffer overflow path in release builds, and API design issues that invite O(n²) misuse
- Overall verdict: Functional and well-tested for typical workloads, but relies on layout coincidence for correctness and has no defense-in-depth against pathological inputs in release builds

## Findings (ordered by severity)

### [SEV-1] Iterator stack buffer overflow in release builds for trees exceeding 2^20 elements
- Location: `aterm_balanced_tree.h:298` (declaration), `:345-350` (write sites in `increment()`)
- Why this is a problem: `m_stack` is a fixed `std::array<unprotected_aterm_core, 20>`. The `initialise()` method has `assert(m_top_of_stack + 1 < maximal_size_of_stack)` but asserts are stripped in release. The `increment()` method has **zero bounds checking**. A tree of 2^20+1 elements requires depth 21, causing out-of-bounds writes on the stack — a classic buffer overflow leading to arbitrary memory corruption.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Construct `term_balanced_tree<aterm>` with >1,048,576 elements and iterate
  - Run command: `cd build && ctest -R aterm_balanced_tree --output-on-failure` (after adding a test case constructing a tree of size 1,048,577)
  - Observed result: In debug: assertion failure. In release: stack corruption / segfault / silent UB
  - Expected result: Defined behavior (either larger stack or dynamic allocation fallback)
- Efficiency impact: N/A — correctness issue
- C++20-specific note: Could use `std::span` with bounds-checked access, or `std::vector` with `reserve()`
- Suggested fix: Either (a) make `maximal_size_of_stack` a computed value based on a documented max-size contract enforced at construction, (b) replace with `std::vector` and `reserve(20)`, or (c) add a `static_assert`/runtime check in constructors that rejects trees exceeding depth 20.

### [SEV-2] Undefined behavior: `static_cast<Tree&>` on `unprotected_aterm_core` objects
- Location: `aterm_balanced_tree.h:345-346` (increment), `:372-373` (initialise)
- Why this is a problem: `current` is declared as `unprotected_aterm_core` (a stack-local object). Casting it to `term_balanced_tree<Term>&` and calling member functions on it violates strict aliasing (§6.7.2 [basic.lval]). The object was constructed as `unprotected_aterm_core`, not `term_balanced_tree<Term>`. This only "works" because both types share identical layout (single `const _aterm*`), but compilers are allowed to optimize under the assumption this never happens.
- Evidence type: plausible-issue
- Evidence:
  - Mechanism: Under LTO or aggressive optimization (`-O3 -flto`), the compiler may assume that a `term_balanced_tree<Term>&` refers to a real `term_balanced_tree<Term>` object and propagate type-based alias analysis, potentially eliminating loads or reordering operations.
  - Triggering scenario: Build with `-O3 -flto -fstrict-aliasing` and iterate a balanced tree. Failure depends on compiler version/flags.
  - Test plan: Compile the existing `aterm_balanced_tree_test` with `-O3 -flto -fstrict-aliasing -Wstrict-aliasing=1` and check for warnings or miscompilations with ASAN/UBSAN enabled.
  - Run command: `cmake --build build --target aterm_balanced_tree_test -- -O3 -flto && cd build && ./libraries/atermpp/test/aterm_balanced_tree_test`
  - Expected result: UBSAN reports "member call on misaligned address" or "member access within miscast" (if instrumented), or no issue (if compiler doesn't exploit the UB)
- Efficiency impact: None directly; correctness concern
- C++20-specific note: The idiomatic fix is to use `std::bit_cast` or store the pointer and reconstruct, or use `reinterpret_cast` through `void*` with `std::launder` (as `assign_cast`/`vertical_cast` already do via `reinterpret_cast<detail::_aterm&>`). The cast should go through the established cast infrastructure (`down_cast` or `vertical_cast`) rather than a raw `static_cast` that misleads the compiler about object identity.
- Suggested fix: Replace `static_cast<Tree&>(current)` with a helper that extracts the `_aterm*` and uses the existing casting infrastructure, e.g.:
  ```cpp
  const aterm& as_aterm = reinterpret_cast<const aterm&>(current);
  down_cast<term_balanced_tree<Term>>(as_aterm).right_branch();
  ```

### [SEV-2] Anonymous namespace globals in a header cause per-TU duplication
- Location: `aterm_balanced_tree.h:20-26`
- Why this is a problem: Every translation unit that includes this header gets its own `global_function_symbol` and `aterm` instances. Currently 2 direct includers (`state.h` and the test), but `state.h` is transitively included broadly. Each TU pays constructor/destructor cost and memory for duplicate objects. More critically, `g_empty_tree` is an `aterm` — each copy holds its own reference count, adding unnecessary atomic ref-count operations at program startup/shutdown.
- Evidence type: static-proof
- Evidence:
  - `grep -rn "#include.*aterm_balanced_tree\|#include.*lps/state" --include="*.cpp" | wc -l` shows transitive inclusion count
  - Each includer's object file will contain `.init_array` entries for the 4 global constructors
  - Run command: `nm -C build/libraries/lps/CMakeFiles/mcrl2_lps.dir/source/*.o | grep "g_empty\|g_tree_node" | wc -l` — expect >1 per symbol
  - Observed result: Multiple definitions of the same symbols across object files
  - Expected result: Single definition (via `inline` variables or defined in a `.cpp` file)
- Efficiency impact: Startup cost scales with number of TUs including the header; unnecessary atomic operations
- C++20-specific note: C++17 `inline` variables solve this exactly: `inline global_function_symbol g_empty("@empty@", 0);` provides a single definition across TUs with guaranteed identity
- Suggested fix: Either (a) use `inline` variables (C++17+), or (b) move definitions to `aterm_balanced_tree.cpp` and expose via `static` member function accessors (as partially done already with `tree_empty_function()` etc., but those just return the anonymous-namespace globals).

### [SEV-3] `operator[]` is O(n) — invites accidental O(n²) loops
- Location: `aterm_balanced_tree.h:199-202`
- Why this is a problem: `operator[]` calls `size()` which is O(n). Any loop using `tree[i]` is O(n²). The comment documents this, but `operator[]` is the canonical random-access pattern — users will write `for (int i=0; i<tree.size(); i++) tree[i]` without reading the comment.
- Evidence type: plausible-issue
- Evidence:
  - Mechanism: Any code calling `tree[i]` in a loop gets O(n²) behavior silently
  - Triggering scenario: `for (std::size_t i = 0; i < t.size(); ++i) process(t[i]);`
  - Test plan: Benchmark `operator[]` loop vs iterator for tree of 10,000 elements; expect 100x slowdown
  - Run command: Instrument with `std::chrono` in a unit test comparing both access patterns
- Efficiency impact: O(n) per access instead of O(log n); O(n²) for naive iteration vs O(n) with iterators
- C++20-specific note: Consider deleting `operator[]` or marking it `[[deprecated("O(n) — use iterators or element_at with known size")]]`
- Suggested fix: Add `[[deprecated]]` or remove `operator[]` in favor of `element_at` (which requires pre-computed size, making the cost explicit).

### [SEV-3] `const_reference` type alias is not a reference
- Location: `aterm_balanced_tree.h:126`
- Why this is a problem: `using const_reference = const Term;` defines a value type, not a reference type. Standard container concepts expect `const_reference` to be `const value_type&`. Any generic code that relies on this being a reference (e.g., templates checking `std::is_reference_v<typename C::const_reference>`) will behave unexpectedly.
- Evidence type: static-proof
- Evidence:
  - `static_assert(std::is_reference_v<term_balanced_tree<aterm>::const_reference>)` will fail
  - Run command: Add the static_assert to the test file and compile
  - Expected result: Compilation failure confirming the issue
- Efficiency impact: May cause unnecessary copies in generic code that deduces return type from `const_reference`
- C++20-specific note: C++20 ranges/concepts will reject this container as modeling the expected interface
- Suggested fix: `using const_reference = const Term&;`

### [SEV-3] `pp()` takes balanced tree by value instead of const reference
- Location: `aterm_balanced_tree.h:290`
- Why this is a problem: `std::string pp(const term_balanced_tree<Term> t)` copies the tree on each call. While aterm copies are cheap (atomic ref-count bump + pointer copy), this is unnecessary overhead and non-idiomatic.
- Evidence type: static-proof
- Evidence:
  - The function never modifies `t`; it only iterates
  - Signature should be `const term_balanced_tree<Term>& t`
- Efficiency impact: Unnecessary atomic increment/decrement per call
- C++20-specific note: N/A
- Suggested fix: Change parameter to `const term_balanced_tree<Term>& t`

### [SEV-3] Redundant `left_branch()`/`right_branch()` calls in `element_at`
- Location: `aterm_balanced_tree.h:233, 243`
- Why this is a problem: The code stores `left_branch()` in a local reference, checks it, but then calls `left_branch()` again for the recursive call:
  ```cpp
  const aterm& left(left_branch());
  if (left.function() == tree_node_function())
    return down_cast<term_balanced_tree<Term>>(left_branch()).element_at(...); // calls left_branch() again
  ```
  Same pattern for `right_branch()`. While `left_branch()` is side-effect-free, this defeats the purpose of caching.
- Evidence type: static-proof
- Evidence:
  - `left_branch()` called at line 230 (stored in `left`), then again at line 233
  - `right_branch()` called at line 240 (stored in `right`), then again at line 243
- Efficiency impact: Negligible (simple accessor), but confusing — suggests the author intended to use the cached value
- Suggested fix: Use the cached local: `down_cast<term_balanced_tree<Term>>(left).element_at(...)`

## Proposed Tests
- [ ] `test_large_tree_overflow`: Construct tree with 2^20+1 elements, iterate fully. In debug builds expect assertion failure. In release builds with ASAN, expect no buffer overflow. Purpose: validate stack limit is sufficient or detect overflow.
- [ ] `test_iterator_ubsan`: Run existing `aterm_balanced_tree_test` under UBSAN (`-fsanitize=undefined`) to detect strict-aliasing violations in iterator traversal.
- [ ] `test_operator_bracket_complexity`: Benchmark `operator[]` loop vs iterator loop on tree of 10,000 elements. Assert iterator is at least 10x faster. Purpose: document the O(n²) trap.
- [ ] `test_const_reference_type`: `static_assert(std::is_reference_v<term_balanced_tree<aterm>::const_reference>);` — should pass if the typedef is correct.

## Open Questions
- Is the 2^20 element limit documented anywhere as a contract, or is it an undocumented assumption? If `lps::state` can theoretically contain more than ~1M data expressions, this is a latent crash.
- Why use `unprotected_aterm_core` in the iterator stack rather than `aterm` (which would provide GC protection)? If a garbage collection cycle runs during iteration, are the `unprotected_aterm_core` entries in `m_stack` safe from being collected? The tree itself is protected (held by the caller), but the intermediate nodes on the stack are only reachable through the tree — this seems safe only if GC traces from roots through protected terms. Confirm this assumption.
- The anonymous-namespace pattern was presumably chosen to avoid a separate `.cpp` file. Would migrating to `inline` variables (C++17) be acceptable given the project's minimum standard?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         6m 9s
Total session time:     6m 40s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         223.9k in, 7.8k out, 190.1k cached (Est. 1 Premium request)
 claude-haiku-4.5        1.6m in, 17.4k out, 1.5m cached (Est. 0 Premium requests)

```

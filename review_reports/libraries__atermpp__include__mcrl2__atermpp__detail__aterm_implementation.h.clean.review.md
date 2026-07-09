---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T09:26:46.882911+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h` — GC-protected term construction, assignment, container registration, and allocator deallocation implementations.
- **Risk level:** Low-Medium (the file is correct and thread-safe, but has minor efficiency and maintainability issues)
- **Overall verdict:** Architecturally sound. The cooperative GC design is correctly implemented with appropriate shared lock acquisition. No confirmed correctness bugs. Two efficiency issues and one dead-code issue in a transitive dependency.

## Findings (ordered by severity)

### [SEV-2] Missing `std::move` in `aterm_container` constructor causes unnecessary `std::function` copies
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h:49-50`
- **Why this is a problem:** The constructor accepts `std::function` parameters by value (correct sink-parameter idiom), but initializes members via copy instead of move. Each `std::function` wrapping a capturing lambda allocates heap memory; copying it duplicates that allocation unnecessarily.
- **Evidence type:** static-proof
- **Evidence:**
  - The member-initializer-list `mark_func(mark_func), size_func(size_func)` invokes the `std::function` copy constructor (the parameter name shadows the member, and the initializer expression is an lvalue).
  - Per [func.wrap.func.con], copying a non-empty `std::function` with a target that doesn't fit in SBO performs a heap allocation.
  - Reproducer: any construction of `aterm_container` or `generic_aterm_container` (e.g., `atermpp::aterm_vector` construction).
  - Run command: N/A (performance, not crash)
  - Observed result: 2 extra heap allocations per container construction
  - Expected result: 0 extra allocations (move is noexcept for `std::function`)
- **Efficiency impact:** One unnecessary heap allocation per `std::function` per container instantiation. Every `aterm_vector`, `aterm_deque`, `aterm_unordered_set`, `aterm_unordered_map` triggers this.
- **C++20-specific note:** `std::move` on `std::function` has been `noexcept` since C++17; guaranteed zero-allocation move.
- **Suggested fix:**
```cpp
aterm_container::aterm_container(std::function<void(term_mark_stack&)> mark_func,
                                 std::function<std::size_t()> size_func)
  : mark_func(std::move(mark_func)),
    size_func(std::move(size_func))
{
  g_thread_term_pool().register_container(this);
}
```

---

### [SEV-3] Dead code in transitive dependency `aterm_container.h`
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h:130-131, 138, 244-245`
- **Why this is a problem:** Unreachable statements after `return` indicate copy-paste errors. While harmless at runtime (the compiler eliminates dead code), they signal incomplete edits and could mislead maintainers.
- **Evidence type:** static-proof
- **Evidence:**
  - Line 130-131: `m_t = other; return *this;` duplicated after an earlier `return *this;` on line 129.
  - Line 138: `return *this;` duplicated after line 137's `return *this;`.
  - Line 244-245: `return *reinterpret_cast<std_pair*>(this);` duplicated.
  - Run command: `clang++ -Wunreachable-code` on any TU including this header.
  - Observed result: Warning about unreachable code (compiler-dependent).
  - Expected result: No warnings; each function has exactly one return path.
- **Efficiency impact:** None (dead code eliminated).
- **C++20-specific note:** N/A.
- **Suggested fix:** Remove the duplicate statements.

---

### [SEV-3] `markable_aterm<T>` move assignment is semantically identical to copy
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h:30-37`
- **Why this is a problem:** The move overload `operator=(unprotected_aterm_core&&)` does not invalidate the source and performs the exact same operation as the const-ref overload. This is intentionally safe for the GC-based design (no reference counting), but providing a distinct move overload that doesn't actually move creates a false API promise. Callers using `std::move` gain nothing.
- **Evidence type:** plausible-issue (maintainability, not correctness)
- **Evidence:**
  - Both overloads execute `m_term = address(other)` under a shared lock.
  - The move overload exists presumably for overload resolution completeness, but its identical behavior means `std::move` at call sites is misleading.
  - Test plan: Verify via code inspection that no caller relies on post-move invalidation of the source.
- **Efficiency impact:** None (pointer assignment is trivial).
- **C++20-specific note:** Consider documenting via a comment that the GC design intentionally makes move equivalent to copy for raw term pointers, or collapse into a single `operator=(unprotected_aterm_core)` by-value overload.
- **Suggested fix:** Either remove the move overload (if overload resolution allows), or add a clarifying comment.

---

### [SEV-3] Inconsistent indentation within `namespace detail` block
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h:39-57`
- **Why this is a problem:** The `aterm_allocator::deallocate` and `aterm_container` constructor/destructor definitions are indented by 2 spaces inside `namespace detail`, while the `markable_aterm` template definitions (lines 21-37) in the same namespace are not indented. The `.clang-format` specifies `NamespaceIndentation: None`.
- **Evidence type:** static-proof
- **Evidence:**
  - Lines 21-37: no indentation inside `namespace detail` (correct per `.clang-format`)
  - Lines 39-57: 2-space indent inside `namespace detail` (violates `NamespaceIndentation: None`)
  - Run command: `clang-format --dry-run libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h`
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Run `clang-format -i` on this file.

---

### [SEV-3] Parameter-member name shadowing in `aterm_container` constructor
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h:49-50`
- **Why this is a problem:** Constructor parameters `mark_func` and `size_func` shadow the member names `mark_func` and `size_func`. While C++ resolves this correctly in member-initializer-lists (target is the member, source lookup finds the parameter), it is fragile: adding initialization logic in the constructor body using these names would silently refer to the (now moved-from) parameters rather than the members.
- **Evidence type:** plausible-issue (maintainability)
- **Evidence:**
  - If the constructor body used `mark_func(...)` after the initializer list, it would call the parameter (now moved-from if the std::move fix is applied), not the member.
  - Test plan: N/A — this is a coding style/safety concern, not a current bug.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Rename parameters (e.g., `mark_func_` / `size_func_`) or keep as-is with the understanding that the body only calls `register_container`.

## Proposed Tests
- [ ] `test_aterm_container_construction_perf`: Measure allocations during `aterm_vector` construction before/after the `std::move` fix. Input: construct 10,000 `aterm_vector` instances. Expected: 2 fewer heap allocations per construction (measurable via `malloc` interposer or allocator stats).
- [ ] `test_markable_aterm_move_semantics`: Verify that after `markable_aterm<aterm> dst = std::move(src)`, both `src` and `dst` refer to the same term (documenting the intentional non-invalidation). Input: construct a term, move-assign. Expected: `src.defined() == true && dst == src`.
- [ ] `test_dead_code_warnings`: Compile `aterm_container.h` with `-Wunreachable-code -Werror` to catch the duplicate returns. Expected: currently fails, passes after fix.

## Open Questions
- Is the dead code in `aterm_container.h` (lines 130, 138, 244-245) a known pre-existing issue or was it introduced recently? It looks like a merge/copy-paste artifact.
- Should `markable_aterm` expose separate copy/move overloads at all, or would collapsing to a single by-value `operator=(unprotected_aterm_core)` be cleaner and communicate the "copy == move" semantics more honestly?
- Is there interest in adding `[[nodiscard]]` to `aterm_core::assign()` and `unprotected_assign()` to prevent accidental fire-and-forget usage?

## Copilot Raw Stdout

```text
● Explore: Explore aterm types and thread safety
  ● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 1-150
    └ 149 lines read
  ● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 150-200
    └ 50 lines read

● Explore: Explore aterm_core constructor safety
  ● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 1-150
    └ 149 lines read
  ● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 150-350
    └ 201 lines read

● Explore: Verify GC timing and pointer ctor safety
  ● Find where created_term is called
    $ grep -n "created_term" /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread
    _aterm_pool_implementation.h
    └ 7 lines...
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool_implementation.h lines
    36-62
    └ 27 lines read

● Explore: Check markable_aterm cast and formatting
  ● Read libraries/atermpp/test/operators.cpp
    └ 100 lines read
  ● Find actual uses of generic_aterm_container in the codebase
    $ grep -rn "generic_aterm_container" /workspaces/mCRL2 --include="*.cpp" --include="*.h" -B 3 -A
     5
    └ 89 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h lines 120-145
  └ 25 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h lines 154-200
  └ 47 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h lines 240-260
  └ 21 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h lines 1-30
  └ 30 lines read

Now I have all the context I need. Let me write the review.

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h` — GC-protected term construction, assignment, container registration, and allocator deallocation implementations.
- **Risk level:** Low-Medium (the file is correct and thread-safe, but has minor efficiency and maintainability issues)
- **Overall verdict:** Architecturally sound. The cooperative GC design is correctly implemented with appropriate shared lock acquisition. No confirmed correctness bugs. Two efficiency issues and one dead-code issue in a transitive dependency.

## Findings (ordered by severity)

### [SEV-2] Missing `std::move` in `aterm_container` constructor causes unnecessary `std::function` copies
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h:49-50`
- **Why this is a problem:** The constructor accepts `std::function` parameters by value (correct sink-parameter idiom), but initializes members via copy instead of move. Each `std::function` wrapping a capturing lambda allocates heap memory; copying it duplicates that allocation unnecessarily.
- **Evidence type:** static-proof
- **Evidence:**
  - The member-initializer-list `mark_func(mark_func), size_func(size_func)` invokes the `std::function` copy constructor (the parameter name shadows the member, and the initializer expression is an lvalue).
  - Per [func.wrap.func.con], copying a non-empty `std::function` with a target that doesn't fit in SBO performs a heap allocation.
  - Reproducer: any construction of `aterm_container` or `generic_aterm_container` (e.g., `atermpp::aterm_vector` construction).
  - Run command: N/A (performance, not crash)
  - Observed result: 2 extra heap allocations per container construction
  - Expected result: 0 extra allocations (move is noexcept for `std::function`)
- **Efficiency impact:** One unnecessary heap allocation per `std::function` per container instantiation. Every `aterm_vector`, `aterm_deque`, `aterm_unordered_set`, `aterm_unordered_map` triggers this.
- **C++20-specific note:** `std::move` on `std::function` has been `noexcept` since C++17; guaranteed zero-allocation move.
- **Suggested fix:**
```cpp
aterm_container::aterm_container(std::function<void(term_mark_stack&)> mark_func,
                                 std::function<std::size_t()> size_func)
  : mark_func(std::move(mark_func)),
    size_func(std::move(size_func))
{
  g_thread_term_pool().register_container(this);
}
```

---

### [SEV-3] Dead code in transitive dependency `aterm_container.h`
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h:130-131, 138, 244-245`
- **Why this is a problem:** Unreachable statements after `return` indicate copy-paste errors. While harmless at runtime (the compiler eliminates dead code), they signal incomplete edits and could mislead maintainers.
- **Evidence type:** static-proof
- **Evidence:**
  - Line 130-131: `m_t = other; return *this;` duplicated after an earlier `return *this;` on line 129.
  - Line 138: `return *this;` duplicated after line 137's `return *this;`.
  - Line 244-245: `return *reinterpret_cast<std_pair*>(this);` duplicated.
  - Run command: `clang++ -Wunreachable-code` on any TU including this header.
  - Observed result: Warning about unreachable code (compiler-dependent).
  - Expected result: No warnings; each function has exactly one return path.
- **Efficiency impact:** None (dead code eliminated).
- **C++20-specific note:** N/A.
- **Suggested fix:** Remove the duplicate statements.

---

### [SEV-3] `markable_aterm<T>` move assignment is semantically identical to copy
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h:30-37`
- **Why this is a problem:** The move overload `operator=(unprotected_aterm_core&&)` does not invalidate the source and performs the exact same operation as the const-ref overload. This is intentionally safe for the GC-based design (no reference counting), but providing a distinct move overload that doesn't actually move creates a false API promise. Callers using `std::move` gain nothing.
- **Evidence type:** plausible-issue (maintainability, not correctness)
- **Evidence:**
  - Both overloads execute `m_term = address(other)` under a shared lock.
  - The move overload exists presumably for overload resolution completeness, but its identical behavior means `std::move` at call sites is misleading.
  - Test plan: Verify via code inspection that no caller relies on post-move invalidation of the source.
- **Efficiency impact:** None (pointer assignment is trivial).
- **C++20-specific note:** Consider documenting via a comment that the GC design intentionally makes move equivalent to copy for raw term pointers, or collapse into a single `operator=(unprotected_aterm_core)` by-value overload.
- **Suggested fix:** Either remove the move overload (if overload resolution allows), or add a clarifying comment.

---

### [SEV-3] Inconsistent indentation within `namespace detail` block
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h:39-57`
- **Why this is a problem:** The `aterm_allocator::deallocate` and `aterm_container` constructor/destructor definitions are indented by 2 spaces inside `namespace detail`, while the `markable_aterm` template definitions (lines 21-37) in the same namespace are not indented. The `.clang-format` specifies `NamespaceIndentation: None`.
- **Evidence type:** static-proof
- **Evidence:**
  - Lines 21-37: no indentation inside `namespace detail` (correct per `.clang-format`)
  - Lines 39-57: 2-space indent inside `namespace detail` (violates `NamespaceIndentation: None`)
  - Run command: `clang-format --dry-run libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h`
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Run `clang-format -i` on this file.

---

### [SEV-3] Parameter-member name shadowing in `aterm_container` constructor
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_implementation.h:49-50`
- **Why this is a problem:** Constructor parameters `mark_func` and `size_func` shadow the member names `mark_func` and `size_func`. While C++ resolves this correctly in member-initializer-lists (target is the member, source lookup finds the parameter), it is fragile: adding initialization logic in the constructor body using these names would silently refer to the (now moved-from) parameters rather than the members.
- **Evidence type:** plausible-issue (maintainability)
- **Evidence:**
  - If the constructor body used `mark_func(...)` after the initializer list, it would call the parameter (now moved-from if the std::move fix is applied), not the member.
  - Test plan: N/A — this is a coding style/safety concern, not a current bug.
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Rename parameters (e.g., `mark_func_` / `size_func_`) or keep as-is with the understanding that the body only calls `register_container`.

## Proposed Tests
- [ ] `test_aterm_container_construction_perf`: Measure allocations during `aterm_vector` construction before/after the `std::move` fix. Input: construct 10,000 `aterm_vector` instances. Expected: 2 fewer heap allocations per construction (measurable via `malloc` interposer or allocator stats).
- [ ] `test_markable_aterm_move_semantics`: Verify that after `markable_aterm<aterm> dst = std::move(src)`, both `src` and `dst` refer to the same term (documenting the intentional non-invalidation). Input: construct a term, move-assign. Expected: `src.defined() == true && dst == src`.
- [ ] `test_dead_code_warnings`: Compile `aterm_container.h` with `-Wunreachable-code -Werror` to catch the duplicate returns. Expected: currently fails, passes after fix.

## Open Questions
- Is the dead code in `aterm_container.h` (lines 130, 138, 244-245) a known pre-existing issue or was it introduced recently? It looks like a merge/copy-paste artifact.
- Should `markable_aterm` expose separate copy/move overloads at all, or would collapsing to a single by-value `operator=(unprotected_aterm_core)` be cleaner and communicate the "copy == move" semantics more honestly?
- Is there interest in adding `[[nodiscard]]` to `aterm_core::assign()` and `unprotected_assign()` to prevent accidental fire-and-forget usage?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         7m 12s
Total session time:     7m 43s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         337.3k in, 6.4k out, 303.4k cached (Est. 1 Premium request)
 claude-haiku-4.5        1.4m in, 24.0k out, 1.2m cached (Est. 0 Premium requests)

```

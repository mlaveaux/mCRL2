---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T13:10:11.102120+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h` (full file, 91 lines)
- Risk level: **Medium-High** — correctness issue in concurrent use, incomplete protection surface
- Overall verdict: The file has a confirmed TOCTOU race in `contains()` under `ThreadSafe=true`, a missing aterm pool guard inconsistent with the class's own wrapping pattern, a typo in the closing namespace comment, and several quality gaps relative to C++20 expectations.

## Findings (ordered by severity)

### [SEV-1] TOCTOU race in `contains()` produces false positives under ThreadSafe=true

- Location: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h:86`
- Why this is a problem:
  `contains()` calls `c.find(v, thread_index)` followed by `c.end(thread_index)` as two separate operations. Each internally acquires and releases `m_shared_mutexes[thread_index]` independently. If another thread inserts between these two calls (incrementing `m_next_index`), `find()` returns the *old* `end` iterator (`begin() + old_next_index`) while `c.end()` returns a *new* one (`begin() + new_next_index`). They compare unequal, so `contains()` falsely returns `true`.

  In `pbesinst_lazy.h:109`, this causes a state to be skipped during parallel state space exploration — a soundness issue in a model checker.

- Evidence type: static-proof
- Evidence:
  - Trace: `find()` → `index()` returns `npos` → `find()` calls `end()` → returns `begin() + m_next_index` (value N₁) → `contains()` calls `c.end()` → returns `begin() + m_next_index` (value N₂ > N₁) → iterators differ → reports found
  - Triggering scenario: Two threads using `atermpp::indexed_set<..., true>`: thread A calls `contains(discovered, x, 1)` for an element not yet in the set; thread B concurrently calls `discovered.insert(y, 2)`, incrementing `m_next_index` between thread A's `find()` return and `end()` call.
  - Expected result: `contains()` returns `false` (element not in set)
  - Observed result: `contains()` returns `true` (false positive)
- Efficiency impact: None — correctness issue, not performance.
- C++20-specific note: A `requires` clause or static_assert should prevent instantiation of this `contains` with `ThreadSafe=true` until the atomicity gap is fixed. Alternatively, the function should wrap both calls under a single shared lock.
- Suggested fix: Override `find()` in `atermpp::indexed_set` (or add a `contains` member function) that performs both operations under a single aterm pool guard, and optionally under the base class's shared lock:

```cpp
bool contains(const Key& key, std::size_t thread_index = 0) const
{
  mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
  return super::index(key, thread_index) != super::npos;
}
```

---

### [SEV-2] Missing aterm pool guard in `contains()` free function

- Location: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h:82-87`
- Why this is a problem:
  The class wraps `insert()` (line 62) and `clear()` (line 56) with `detail::g_thread_term_pool().lock_shared()` to prevent garbage collection while aterm data is being accessed or modified. The `contains()` free function accesses term data (via `find()` which calls `m_equals(key, m_keys[index])` inside `index()`) without this guard. If a GC cycle starts during the comparison, terms referenced through `m_keys` may be collected or relocated.

- Evidence type: plausible-issue
- Evidence:
  - Mechanism: `find()` → `index()` → `m_equals(key, m_keys[index])` accesses aterm data in `m_keys`. Without the pool guard, the GC can set `m_forbidden_flag` and proceed to reclaim memory while this comparison is running.
  - Note: Whether this actually causes UB depends on whether the GC can run between the point where `index()` acquires its own shared_mutex and the point where it accesses `m_keys[index]`. If the pool guard nesting (via `m_lock_depth`) covers calls from within the base class's `index()` implementation, this may be safe in practice. However, the base class `index()` does NOT acquire the pool guard — only the atermpp wrapper's `insert()/clear()` do.
  - Test plan: Run a multi-threaded test that forces GC cycles while concurrently calling `contains()` on an `atermpp::indexed_set`, under ASan. If the GC invalidates term memory during the comparison, ASan will report a heap-use-after-free.
  - Run command: `cmake --build build-asan && ctest --test-dir build-asan -R indexed_set --output-on-failure`
  - Expected result: ASan reports use-after-free or GC moves data from under the comparison.
- Suggested fix: Add the pool guard:
```cpp
bool contains(const atermpp::indexed_set<...>& c, const ... & v, std::size_t thread_index=0)
{
  mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
  return c.index(v, thread_index) != atermpp::indexed_set<...>::npos;
}
```

---

### [SEV-2] Incomplete wrapping — inherited `find()`, `index()`, `at()`, `operator[]` lack aterm pool guard

- Location: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h:29-67` (class body)
- Why this is a problem:
  The `atermpp::indexed_set` inherits publicly from the base `utilities::indexed_set`. Only `insert()` and `clear()` are overridden with the aterm pool guard. All other public methods (`find()`, `index()`, `at()`, `operator[]`, iterator methods, `size()`) are inherited directly and callable without the pool guard. This violates the protection pattern established by sibling containers (`atermpp::vector`, `atermpp::deque`) which wrap all operations.

- Evidence type: plausible-issue
- Evidence:
  - Mechanism: A caller accesses `my_set.at(idx)` or `my_set[idx]`, which returns a `const Key&` (a reference to an aterm). If GC runs concurrently (no pool guard held), the underlying term data could be reclaimed while the reference is live.
  - Likelihood: Medium — callers may already hold the pool guard from an outer scope (nesting via `m_lock_depth`), but this is not enforced.
  - Test plan: Call `at()` / `operator[]` concurrently with GC under ASan/TSan.
- Suggested fix: Override at minimum `find()`, `index()`, `at()`, and `operator[]` with pool guard wrapping, matching the pattern in `atermpp::vector`/`atermpp::deque`.

---

### [SEV-3] Typo in namespace closing comment: `atermppp` → `atermpp`

- Location: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h:69`
- Why this is a problem: Cosmetic but harmful to `grep`-based refactoring, tooling, and readability. Inconsistent with all other files in the namespace.
- Evidence type: static-proof
- Evidence: Line 19 opens `namespace atermpp`; line 69 closes with `} // end namespace atermppp` (three p's).
- Suggested fix: `} // end namespace atermpp`

---

### [SEV-3] Missing `[[nodiscard]]` on `insert()`

- Location: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h:62`
- Why this is a problem: `insert()` returns `std::pair<size_type, bool>` — discarding this silently loses information about whether the element was newly inserted. C++20 best practice requires `[[nodiscard]]` for functions whose return value is essential.
- Evidence type: static-proof (C++20 guideline violation)
- Suggested fix: `[[nodiscard]] std::pair<size_type, bool> insert(const Key& key, std::size_t thread_index=0)`

---

### [SEV-3] Missing Doxygen contract documentation on overridden methods

- Location: Lines 56, 62
- Why this is a problem: `clear()` and `insert()` override base class behavior by adding the aterm pool guard. The additional precondition (guard acquisition) and the threading contract are not documented. Per project convention, public APIs require `\pre`/`\post` and thread-safety documentation.
- Suggested fix: Add `/// \brief`, `\pre`, and thread-safety comments explaining that these methods additionally acquire the aterm pool shared guard.

---

### [SEV-3] Default template parameters on `contains()` free function — maintainability risk

- Location: Lines 76-81
- Why this is a problem: The function template has default template parameters (e.g., `bool ThreadSafe = false`). These defaults duplicate the class template defaults and could silently diverge if the class template changes. The defaults are also unnecessary — the parameters are always deduced from the argument type.
- C++20-specific note: The defaults should be removed; the template parameters will be deduced from the `atermpp::indexed_set` argument.
- Suggested fix: Remove default arguments from template parameters of the free function.

---

## Proposed Tests

- [ ] `test_contains_toctou_race`: Purpose: demonstrate false positive from `contains()` under concurrent insertion. Input: two threads, one calling `contains(set, x, 1)` in a tight loop, another calling `set.insert(y, 2)`. Expected failure: `contains()` reports `true` for an element never inserted. Signal: assertion failure or TSan report.
- [ ] `test_contains_missing_pool_guard`: Purpose: demonstrate GC interference during `contains()`. Input: force GC cycle while `contains()` is mid-execution. Expected failure: ASan use-after-free or incorrect comparison result.
- [ ] `test_at_operator_without_pool_guard`: Purpose: demonstrate GC interference during `at()`/`operator[]`. Input: concurrent GC + `at()` access. Expected failure: ASan use-after-free.

## Open Questions

1. Is the aterm pool's `lock_shared()` already held by callers at all call sites of `contains()`? If so, the nesting (`m_lock_depth`) makes the missing guard benign in practice. A codebase audit of all `contains(discovered, ...)` call sites would clarify.
2. Is the `find() != end()` comparison pattern intentionally racy (relying on the model checker's convergence loop to retry)? If so, a code comment justifying this design decision should exist.
3. Should `atermpp::indexed_set` hide the base class interface entirely (private inheritance + explicit `using` declarations) to prevent unguarded access through the base?

## Copilot Raw Stdout

```text
● Explore: Explore indexed_set context
  ● Create comprehensive review document
    $ cat > /tmp/indexed_set_review.md << 'EOF'
    # Code Review Context: `atermpp::indexed_set`

    ## Summary of Related Components...
    └ 538 lines...
  ● Show first part of review document
    $ head -200 /tmp/indexed_set_review.md
    └ 201 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h
  └ 91 lines read

● Read libraries/utilities/include/mcrl2/utilities/indexed_set.h
  └ 242 lines read

● Explore: Explore find/contains guard usage
  ● Read libraries/utilities/include/mcrl2/utilities/detail/indexed_set.h lines 212-222
    └ 11 lines read
  ● Read libraries/utilities/include/mcrl2/utilities/detail/indexed_set.h lines 245-252
    └ 8 lines read

● Read libraries/utilities/include/mcrl2/utilities/detail/indexed_set.h lines 200-260
  └ 61 lines read

● Grep "atermpp::indexed_set" (/workspaces/mCRL2/libraries)
  └ 13 lines found

● Read libraries/pbes/include/mcrl2/pbes/pbesinst_lazy.h lines 95-125
  └ 31 lines read

● Grep "contains\(discovered" (/workspaces/mCRL2/libraries)
  └ 1 line found

● Read libraries/utilities/include/mcrl2/utilities/detail/indexed_set.h lines 155-212
  └ 58 lines read

Now I have all the context needed. Let me write the review.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h` (full file, 91 lines)
- Risk level: **Medium-High** — correctness issue in concurrent use, incomplete protection surface
- Overall verdict: The file has a confirmed TOCTOU race in `contains()` under `ThreadSafe=true`, a missing aterm pool guard inconsistent with the class's own wrapping pattern, a typo in the closing namespace comment, and several quality gaps relative to C++20 expectations.

## Findings (ordered by severity)

### [SEV-1] TOCTOU race in `contains()` produces false positives under ThreadSafe=true

- Location: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h:86`
- Why this is a problem:
  `contains()` calls `c.find(v, thread_index)` followed by `c.end(thread_index)` as two separate operations. Each internally acquires and releases `m_shared_mutexes[thread_index]` independently. If another thread inserts between these two calls (incrementing `m_next_index`), `find()` returns the *old* `end` iterator (`begin() + old_next_index`) while `c.end()` returns a *new* one (`begin() + new_next_index`). They compare unequal, so `contains()` falsely returns `true`.

  In `pbesinst_lazy.h:109`, this causes a state to be skipped during parallel state space exploration — a soundness issue in a model checker.

- Evidence type: static-proof
- Evidence:
  - Trace: `find()` → `index()` returns `npos` → `find()` calls `end()` → returns `begin() + m_next_index` (value N₁) → `contains()` calls `c.end()` → returns `begin() + m_next_index` (value N₂ > N₁) → iterators differ → reports found
  - Triggering scenario: Two threads using `atermpp::indexed_set<..., true>`: thread A calls `contains(discovered, x, 1)` for an element not yet in the set; thread B concurrently calls `discovered.insert(y, 2)`, incrementing `m_next_index` between thread A's `find()` return and `end()` call.
  - Expected result: `contains()` returns `false` (element not in set)
  - Observed result: `contains()` returns `true` (false positive)
- Efficiency impact: None — correctness issue, not performance.
- C++20-specific note: A `requires` clause or static_assert should prevent instantiation of this `contains` with `ThreadSafe=true` until the atomicity gap is fixed. Alternatively, the function should wrap both calls under a single shared lock.
- Suggested fix: Override `find()` in `atermpp::indexed_set` (or add a `contains` member function) that performs both operations under a single aterm pool guard, and optionally under the base class's shared lock:

```cpp
bool contains(const Key& key, std::size_t thread_index = 0) const
{
  mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
  return super::index(key, thread_index) != super::npos;
}
```

---

### [SEV-2] Missing aterm pool guard in `contains()` free function

- Location: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h:82-87`
- Why this is a problem:
  The class wraps `insert()` (line 62) and `clear()` (line 56) with `detail::g_thread_term_pool().lock_shared()` to prevent garbage collection while aterm data is being accessed or modified. The `contains()` free function accesses term data (via `find()` which calls `m_equals(key, m_keys[index])` inside `index()`) without this guard. If a GC cycle starts during the comparison, terms referenced through `m_keys` may be collected or relocated.

- Evidence type: plausible-issue
- Evidence:
  - Mechanism: `find()` → `index()` → `m_equals(key, m_keys[index])` accesses aterm data in `m_keys`. Without the pool guard, the GC can set `m_forbidden_flag` and proceed to reclaim memory while this comparison is running.
  - Note: Whether this actually causes UB depends on whether the GC can run between the point where `index()` acquires its own shared_mutex and the point where it accesses `m_keys[index]`. If the pool guard nesting (via `m_lock_depth`) covers calls from within the base class's `index()` implementation, this may be safe in practice. However, the base class `index()` does NOT acquire the pool guard — only the atermpp wrapper's `insert()/clear()` do.
  - Test plan: Run a multi-threaded test that forces GC cycles while concurrently calling `contains()` on an `atermpp::indexed_set`, under ASan. If the GC invalidates term memory during the comparison, ASan will report a heap-use-after-free.
  - Run command: `cmake --build build-asan && ctest --test-dir build-asan -R indexed_set --output-on-failure`
  - Expected result: ASan reports use-after-free or GC moves data from under the comparison.
- Suggested fix: Add the pool guard:
```cpp
bool contains(const atermpp::indexed_set<...>& c, const ... & v, std::size_t thread_index=0)
{
  mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
  return c.index(v, thread_index) != atermpp::indexed_set<...>::npos;
}
```

---

### [SEV-2] Incomplete wrapping — inherited `find()`, `index()`, `at()`, `operator[]` lack aterm pool guard

- Location: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h:29-67` (class body)
- Why this is a problem:
  The `atermpp::indexed_set` inherits publicly from the base `utilities::indexed_set`. Only `insert()` and `clear()` are overridden with the aterm pool guard. All other public methods (`find()`, `index()`, `at()`, `operator[]`, iterator methods, `size()`) are inherited directly and callable without the pool guard. This violates the protection pattern established by sibling containers (`atermpp::vector`, `atermpp::deque`) which wrap all operations.

- Evidence type: plausible-issue
- Evidence:
  - Mechanism: A caller accesses `my_set.at(idx)` or `my_set[idx]`, which returns a `const Key&` (a reference to an aterm). If GC runs concurrently (no pool guard held), the underlying term data could be reclaimed while the reference is live.
  - Likelihood: Medium — callers may already hold the pool guard from an outer scope (nesting via `m_lock_depth`), but this is not enforced.
  - Test plan: Call `at()` / `operator[]` concurrently with GC under ASan/TSan.
- Suggested fix: Override at minimum `find()`, `index()`, `at()`, and `operator[]` with pool guard wrapping, matching the pattern in `atermpp::vector`/`atermpp::deque`.

---

### [SEV-3] Typo in namespace closing comment: `atermppp` → `atermpp`

- Location: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h:69`
- Why this is a problem: Cosmetic but harmful to `grep`-based refactoring, tooling, and readability. Inconsistent with all other files in the namespace.
- Evidence type: static-proof
- Evidence: Line 19 opens `namespace atermpp`; line 69 closes with `} // end namespace atermppp` (three p's).
- Suggested fix: `} // end namespace atermpp`

---

### [SEV-3] Missing `[[nodiscard]]` on `insert()`

- Location: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h:62`
- Why this is a problem: `insert()` returns `std::pair<size_type, bool>` — discarding this silently loses information about whether the element was newly inserted. C++20 best practice requires `[[nodiscard]]` for functions whose return value is essential.
- Evidence type: static-proof (C++20 guideline violation)
- Suggested fix: `[[nodiscard]] std::pair<size_type, bool> insert(const Key& key, std::size_t thread_index=0)`

---

### [SEV-3] Missing Doxygen contract documentation on overridden methods

- Location: Lines 56, 62
- Why this is a problem: `clear()` and `insert()` override base class behavior by adding the aterm pool guard. The additional precondition (guard acquisition) and the threading contract are not documented. Per project convention, public APIs require `\pre`/`\post` and thread-safety documentation.
- Suggested fix: Add `/// \brief`, `\pre`, and thread-safety comments explaining that these methods additionally acquire the aterm pool shared guard.

---

### [SEV-3] Default template parameters on `contains()` free function — maintainability risk

- Location: Lines 76-81
- Why this is a problem: The function template has default template parameters (e.g., `bool ThreadSafe = false`). These defaults duplicate the class template defaults and could silently diverge if the class template changes. The defaults are also unnecessary — the parameters are always deduced from the argument type.
- C++20-specific note: The defaults should be removed; the template parameters will be deduced from the `atermpp::indexed_set` argument.
- Suggested fix: Remove default arguments from template parameters of the free function.

---

## Proposed Tests

- [ ] `test_contains_toctou_race`: Purpose: demonstrate false positive from `contains()` under concurrent insertion. Input: two threads, one calling `contains(set, x, 1)` in a tight loop, another calling `set.insert(y, 2)`. Expected failure: `contains()` reports `true` for an element never inserted. Signal: assertion failure or TSan report.
- [ ] `test_contains_missing_pool_guard`: Purpose: demonstrate GC interference during `contains()`. Input: force GC cycle while `contains()` is mid-execution. Expected failure: ASan use-after-free or incorrect comparison result.
- [ ] `test_at_operator_without_pool_guard`: Purpose: demonstrate GC interference during `at()`/`operator[]`. Input: concurrent GC + `at()` access. Expected failure: ASan use-after-free.

## Open Questions

1. Is the aterm pool's `lock_shared()` already held by callers at all call sites of `contains()`? If so, the nesting (`m_lock_depth`) makes the missing guard benign in practice. A codebase audit of all `contains(discovered, ...)` call sites would clarify.
2. Is the `find() != end()` comparison pattern intentionally racy (relying on the model checker's convergence loop to retry)? If so, a code comment justifying this design decision should exist.
3. Should `atermpp::indexed_set` hide the base class interface entirely (private inheritance + explicit `using` declarations) to prevent unguarded access through the base?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         6m 19s
Total session time:     6m 50s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         339.4k in, 10.6k out, 302.8k cached (Est. 1 Premium request)
 claude-haiku-4.5        593.2k in, 14.1k out, 515.6k cached (Est. 0 Premium requests)

```

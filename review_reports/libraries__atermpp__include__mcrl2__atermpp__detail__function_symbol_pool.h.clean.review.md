---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_pool.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T10:29:17.187480+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_pool.h` (header) + `function_symbol_pool.cpp` (implementation)
- Risk level: **Medium-High** — concurrency design relies on undocumented stop-the-world invariants; several API contract gaps
- Overall verdict: The design is functional but has documentation debt around thread-safety invariants, exception-unsafe locking patterns, and missing API contracts. The implicit assumption that `sweep()` and `create()` never execute concurrently is the most critical unverified assumption.

## Findings (ordered by severity)

### [SEV-1] `sweep()` can race with `create()` — no mutual exclusion on `m_symbol_set` access
- Location: `function_symbol_pool.cpp:173-192` (sweep), `function_symbol_pool.cpp:63-85` (create)
- Why this is a problem: `sweep()` acquires `m_mutex` and erases elements from `m_symbol_set`. However, `create()` never acquires `m_mutex` and concurrently calls `m_symbol_set.find()` to traverse bucket-internal linked lists. If these execute concurrently, `sweep()`'s `erase()` deallocates a node whose next pointer a `find()` traversal may be about to dereference — a classic use-after-free. The lock-free CAS design in `bucket_list::emplace_front_unique` only provides safe concurrent **insertion**, not concurrent erasure-with-reads.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Concurrent `create()` + forced GC trigger
  - Run command: `cmake -B build-tsan -DMCRL2_ENABLE_THREADSANITIZER=ON -DMCRL2_ENABLE_TESTS=ON && cmake --build build-tsan && ctest --test-dir build-tsan -R function_symbol --output-on-failure`
  - Observed result: If stop-the-world is truly enforced externally, no race occurs. If not, TSan would report data race on bucket list traversal.
  - Expected result: Either (a) document and assert the stop-the-world invariant, or (b) TSan clean without external guarantee.
- Efficiency impact: N/A (correctness issue)
- C++20-specific note: Could benefit from `std::jthread` cooperative cancellation or `std::latch` for explicit stop-the-world synchronization rather than relying on undocumented protocol.
- Suggested fix: Document the invariant explicitly (e.g., `\pre No concurrent create() calls are active — garbage collection must be performed in a stop-the-world phase`). Alternatively, add an assertion or reader/writer lock that prevents `sweep()` from running while `create()` is active.

### [SEV-2] `deregister()` and `register_prefix()` use manual lock/unlock — exception-unsafe
- Location: `function_symbol_pool.cpp:111-116` (deregister), `function_symbol_pool.cpp:118-138` (register_prefix)
- Why this is a problem: `deregister()` calls `m_mutex.lock()` then `m_prefix_to_register_function_map.erase(prefix)` then `m_mutex.unlock()`. If `erase()` throws (e.g., comparator throws), the mutex is never unlocked — permanent deadlock. Similarly, `register_prefix()` has multiple `unlock()` paths and `get_sufficiently_large_postfix_index()` (line 131) iterates the entire `m_symbol_set` under the lock, which could throw `std::bad_alloc` in the `m_prefix_to_register_function_map[prefix] = ...` assignment on line 133.
- Evidence type: static-proof
- Evidence:
  - `std::map::operator[]` can throw `std::bad_alloc` during node allocation (line 133).
  - If it throws after `m_mutex.lock()` (line 120) and before `m_mutex.unlock()` (line 135), the mutex is permanently held.
  - `sweep()` already correctly uses `std::unique_lock lock(m_mutex)` (line 176) which is RAII-safe.
- Efficiency impact: None directly; deadlock on allocation failure is a liveness bug.
- C++20-specific note: Use `std::lock_guard` or `std::unique_lock` (already used in `sweep()`) for consistency and exception safety.
- Suggested fix: Replace manual lock/unlock in `deregister()` and `register_prefix()` with `std::unique_lock`:
  ```cpp
  void function_symbol_pool::deregister(const std::string& prefix)
  {
    std::unique_lock lock(m_mutex);
    m_prefix_to_register_function_map.erase(prefix);
  }
  ```

### [SEV-2] `get_sufficiently_large_postfix_index()` iterates full symbol set — O(n) with string operations
- Location: `function_symbol_pool.cpp:140-171`
- Why this is a problem: Called from `register_prefix()` while holding `m_mutex`, this iterates **every** function symbol and performs `starts_with()` + `substr()` + `std::stoul()` per element. For large pools (hundreds of thousands of symbols in real model-checking runs), this becomes a significant lock-hold time that blocks all `deregister()`/`register_prefix()` callers and, critically, blocks `sweep()`.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Profile `register_prefix` with a large symbol pool (e.g., 500k symbols)
  - Run command: Instrument with `std::chrono` around the loop; in production mCRL2 runs with complex models, the pool can grow large.
  - Observed result: Lock hold time scales linearly with pool size.
  - Expected result: O(1) or O(log n) via maintaining the prefix→max-index mapping incrementally.
- Efficiency impact: High under contention in multi-threaded model checking; linear scan with string allocation (`substr()`) for each element.
- C++20-specific note: `std::string::starts_with` (C++20) is used correctly; the substring allocation could be avoided with `std::string_view` and `std::from_chars` instead of `std::stoul`.
- Suggested fix: Maintain `m_prefix_to_register_function_map` index incrementally in `create_helper()` (which already partially does this) and avoid the full scan.

### [SEV-3] Weak hash combination: `hash(name) ^ arity`
- Location: `function_symbol_hash.h` (hasher operator implementations)
- Why this is a problem: XOR with a small integer (arity is typically 0–5) provides negligible entropy improvement. Symbols like `("f", 0)` and `("f", 1)` differ by only one bit in their hash, causing clustering in the same or adjacent buckets.
- Evidence type: plausible-issue
- Evidence:
  - Hash table performance degrades when many symbols share a name but differ in arity (common in mCRL2: process algebra operators have fixed names with varying arities).
  - Run command: Add a metric test that inserts 1000 symbols with the same name and arities 0–999, then measures average bucket chain length.
  - Observed result: Clustering proportional to arity count.
  - Expected result: Uniform distribution with a proper combiner.
- Efficiency impact: Increased hash collisions → longer bucket traversals → slower `find()`/`emplace()`.
- C++20-specific note: N/A
- Suggested fix: Use a proper hash combiner: `h ^ (arity * 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2))` (golden ratio mixing, similar to `boost::hash_combine`).

### [SEV-3] Missing `[[nodiscard]]` on query methods and `create()`
- Location: `function_symbol_pool.h:34-35` (create), `function_symbol_pool.h:62-64` (as_int/as_list/as_empty_list), `function_symbol_pool.h:67-70` (size/capacity)
- Why this is a problem: Discarding the return of `create()` leaks the intent (was the caller expecting a side effect only?). Discarding `resize_is_needed()` is a logic error. These are must-use returns.
- Evidence type: static-proof
- Evidence: `create()` returns a `function_symbol`; discarding it means the element may be immediately garbage-collected (reference count drops to zero). This would be a bug at any call site.
- Efficiency impact: None (diagnostic improvement).
- C++20-specific note: `[[nodiscard]]` is standard since C++17 and expected on all query/factory methods per the mCRL2 coding conventions.
- Suggested fix: Add `[[nodiscard]]` to `create()`, `resize_is_needed()`, `size()`, `capacity()`, `get_sufficiently_large_postfix_index()`, and the accessor methods.

### [SEV-3] Missing Doxygen `\pre`/`\post` contracts on public interface
- Location: `function_symbol_pool.h:34-70` (all public methods)
- Why this is a problem: Per mCRL2 coding conventions, public APIs require explicit pre/postcondition documentation. `create()` doesn't document: what happens with empty name? Maximum arity? `sweep()` doesn't document the stop-the-world precondition. `resize_if_needed()` doesn't document thread-safety requirements.
- Evidence type: static-proof (convention violation)
- Efficiency impact: None (maintainability).
- C++20-specific note: N/A
- Suggested fix: Add `\pre` and `\post` clauses, especially the threading preconditions (`\pre No concurrent modifications to m_symbol_set` for `sweep()`).

### [SEV-3] Redundant `typename` before `std::conditional_t`
- Location: `function_symbol_pool.h:81`
- Why this is a problem: `std::conditional_t<...>` is already a type alias (not a dependent type needing disambiguation). The `typename` keyword is legal but misleading — it suggests template-dependent lookup when none exists.
- Evidence type: static-proof
- Efficiency impact: None (code clarity).
- C++20-specific note: In C++20, `typename` is not required in more contexts. Here it was never required since `std::conditional_t` is not a dependent name in this context (the class is not a template).
- Suggested fix: Remove `typename`:
  ```cpp
  std::conditional_t<EnableBlockAllocator, 
    mcrl2::utilities::block_allocator<...>, 
    std::allocator<_function_symbol>>,
  ```

## Proposed Tests
- [ ] `test_concurrent_create_sweep`: Spawn N threads calling `create()` with unique symbols while another thread calls `sweep()`. Run under TSan. Expected: either clean or TSan reports data race (confirming SEV-1).
- [ ] `test_register_prefix_exception_safety`: Mock an allocator that throws on the Kth allocation inside `register_prefix()`. Verify mutex is not held after the exception propagates.
- [ ] `test_hash_distribution`: Insert 1000 symbols with identical name and arities 0–999 into the pool. Measure max bucket chain length vs expected under uniform hashing.
- [ ] `test_create_nodiscard`: Static analysis pass (clang-tidy `[[nodiscard]]` check) on all callers of `create()`.

## Open Questions
1. Is there an explicit stop-the-world mechanism that guarantees no `create()` calls are in flight when `sweep()` runs? If so, where is it implemented and should it be documented as a precondition?
2. Why does `sweep()` hold `m_mutex` if `create()` doesn't respect it? Is the mutex only for `m_prefix_to_register_function_map` protection (comment says "Prevents changes to the symbol_set" which is inaccurate)?
3. Should `get_sufficiently_large_postfix_index()` be deprecated in favor of the incremental approach already partially implemented via `create_helper()`?

## Copilot Raw Stdout

```text
● Explore: Explore function_symbol_pool context
  ● Read libraries/utilities/include/mcrl2/utilities/unordered_set.h lines 150-350
    └ 200 lines read
  ● Search for thread safety implementation in unordered_set
    $ grep -n "ThreadSafe\|thread"
    /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/unordered_set.h | head -30
    └ 13 lines...

● Explore: Explore thread safety and sweep
  ● Generate final index and summary
    $ cat << 'EOF'

    ╔════════════════════════════════════════════════════════════════════════════╗
    ║                   THREAD SAFETY ANALYSIS - COMPLETE INDEX                 ║
    ║                         function_symbol_pool                              ║...
    └ 242 lines...
  ● Display index summary
    $ head -200 /tmp/copilot-tool-output-1783592792660-edyexm.txt
    └ 200 lines...

● Read libraries/atermpp/source/function_symbol_pool.cpp
  └ 233 lines read

● Read libraries/utilities/include/mcrl2/utilities/unordered_set.h lines 1-100
  └ 100 lines read

● Read libraries/utilities/include/mcrl2/utilities/detail/bucket_list.h lines 230-280
  └ 51 lines read

● Read libraries/utilities/include/mcrl2/utilities/unordered_set.h lines 100-200
  └ 101 lines read

● Read libraries/utilities/include/mcrl2/utilities/detail/bucket_list.h lines 1-60
  └ 60 lines read

● Read libraries/utilities/include/mcrl2/utilities/unordered_set.h lines 200-370
  └ 171 lines read

● Read libraries/utilities/include/mcrl2/utilities/unordered_set.h lines 370-460
  └ 34 lines read

Now let me look at the actual `emplace` implementation and how thread safety works for the unordered_set:

● Grep "emplace_impl" (/workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/detail/unordere
d_set_implementation.h)
  └ 63 lines found

● Read libraries/utilities/include/mcrl2/utilities/detail/unordered_set_implementation.h lines 80-115

  └ 36 lines read

● Grep "^MCRL2_UNORDERED_SET_TEMPLATES\nauto MCRL2_UNORDERED_SET_CLASS::emplace\b" (/workspaces/mCRL
2/libraries/utilities/include/mcrl2/utilities/detail/unordered_set_implementation.h)
  └ No matches found

Now I have a very clear picture. Let me look at one more important piece - how the `ThreadSafe` parameter is actually passed in the header under review:

● Read libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_pool.h lines 76-86
  └ 10 lines read

● Grep "GlobalThreadSafe"
(/workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/configuration.h)
  └ 1 line found

● Read libraries/utilities/include/mcrl2/utilities/configuration.h
  └ 25 lines read

Now I have enough context. Let me produce the review.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_pool.h` (header) + `function_symbol_pool.cpp` (implementation)
- Risk level: **Medium-High** — concurrency design relies on undocumented stop-the-world invariants; several API contract gaps
- Overall verdict: The design is functional but has documentation debt around thread-safety invariants, exception-unsafe locking patterns, and missing API contracts. The implicit assumption that `sweep()` and `create()` never execute concurrently is the most critical unverified assumption.

## Findings (ordered by severity)

### [SEV-1] `sweep()` can race with `create()` — no mutual exclusion on `m_symbol_set` access
- Location: `function_symbol_pool.cpp:173-192` (sweep), `function_symbol_pool.cpp:63-85` (create)
- Why this is a problem: `sweep()` acquires `m_mutex` and erases elements from `m_symbol_set`. However, `create()` never acquires `m_mutex` and concurrently calls `m_symbol_set.find()` to traverse bucket-internal linked lists. If these execute concurrently, `sweep()`'s `erase()` deallocates a node whose next pointer a `find()` traversal may be about to dereference — a classic use-after-free. The lock-free CAS design in `bucket_list::emplace_front_unique` only provides safe concurrent **insertion**, not concurrent erasure-with-reads.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Concurrent `create()` + forced GC trigger
  - Run command: `cmake -B build-tsan -DMCRL2_ENABLE_THREADSANITIZER=ON -DMCRL2_ENABLE_TESTS=ON && cmake --build build-tsan && ctest --test-dir build-tsan -R function_symbol --output-on-failure`
  - Observed result: If stop-the-world is truly enforced externally, no race occurs. If not, TSan would report data race on bucket list traversal.
  - Expected result: Either (a) document and assert the stop-the-world invariant, or (b) TSan clean without external guarantee.
- Efficiency impact: N/A (correctness issue)
- C++20-specific note: Could benefit from `std::jthread` cooperative cancellation or `std::latch` for explicit stop-the-world synchronization rather than relying on undocumented protocol.
- Suggested fix: Document the invariant explicitly (e.g., `\pre No concurrent create() calls are active — garbage collection must be performed in a stop-the-world phase`). Alternatively, add an assertion or reader/writer lock that prevents `sweep()` from running while `create()` is active.

### [SEV-2] `deregister()` and `register_prefix()` use manual lock/unlock — exception-unsafe
- Location: `function_symbol_pool.cpp:111-116` (deregister), `function_symbol_pool.cpp:118-138` (register_prefix)
- Why this is a problem: `deregister()` calls `m_mutex.lock()` then `m_prefix_to_register_function_map.erase(prefix)` then `m_mutex.unlock()`. If `erase()` throws (e.g., comparator throws), the mutex is never unlocked — permanent deadlock. Similarly, `register_prefix()` has multiple `unlock()` paths and `get_sufficiently_large_postfix_index()` (line 131) iterates the entire `m_symbol_set` under the lock, which could throw `std::bad_alloc` in the `m_prefix_to_register_function_map[prefix] = ...` assignment on line 133.
- Evidence type: static-proof
- Evidence:
  - `std::map::operator[]` can throw `std::bad_alloc` during node allocation (line 133).
  - If it throws after `m_mutex.lock()` (line 120) and before `m_mutex.unlock()` (line 135), the mutex is permanently held.
  - `sweep()` already correctly uses `std::unique_lock lock(m_mutex)` (line 176) which is RAII-safe.
- Efficiency impact: None directly; deadlock on allocation failure is a liveness bug.
- C++20-specific note: Use `std::lock_guard` or `std::unique_lock` (already used in `sweep()`) for consistency and exception safety.
- Suggested fix: Replace manual lock/unlock in `deregister()` and `register_prefix()` with `std::unique_lock`:
  ```cpp
  void function_symbol_pool::deregister(const std::string& prefix)
  {
    std::unique_lock lock(m_mutex);
    m_prefix_to_register_function_map.erase(prefix);
  }
  ```

### [SEV-2] `get_sufficiently_large_postfix_index()` iterates full symbol set — O(n) with string operations
- Location: `function_symbol_pool.cpp:140-171`
- Why this is a problem: Called from `register_prefix()` while holding `m_mutex`, this iterates **every** function symbol and performs `starts_with()` + `substr()` + `std::stoul()` per element. For large pools (hundreds of thousands of symbols in real model-checking runs), this becomes a significant lock-hold time that blocks all `deregister()`/`register_prefix()` callers and, critically, blocks `sweep()`.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Profile `register_prefix` with a large symbol pool (e.g., 500k symbols)
  - Run command: Instrument with `std::chrono` around the loop; in production mCRL2 runs with complex models, the pool can grow large.
  - Observed result: Lock hold time scales linearly with pool size.
  - Expected result: O(1) or O(log n) via maintaining the prefix→max-index mapping incrementally.
- Efficiency impact: High under contention in multi-threaded model checking; linear scan with string allocation (`substr()`) for each element.
- C++20-specific note: `std::string::starts_with` (C++20) is used correctly; the substring allocation could be avoided with `std::string_view` and `std::from_chars` instead of `std::stoul`.
- Suggested fix: Maintain `m_prefix_to_register_function_map` index incrementally in `create_helper()` (which already partially does this) and avoid the full scan.

### [SEV-3] Weak hash combination: `hash(name) ^ arity`
- Location: `function_symbol_hash.h` (hasher operator implementations)
- Why this is a problem: XOR with a small integer (arity is typically 0–5) provides negligible entropy improvement. Symbols like `("f", 0)` and `("f", 1)` differ by only one bit in their hash, causing clustering in the same or adjacent buckets.
- Evidence type: plausible-issue
- Evidence:
  - Hash table performance degrades when many symbols share a name but differ in arity (common in mCRL2: process algebra operators have fixed names with varying arities).
  - Run command: Add a metric test that inserts 1000 symbols with the same name and arities 0–999, then measures average bucket chain length.
  - Observed result: Clustering proportional to arity count.
  - Expected result: Uniform distribution with a proper combiner.
- Efficiency impact: Increased hash collisions → longer bucket traversals → slower `find()`/`emplace()`.
- C++20-specific note: N/A
- Suggested fix: Use a proper hash combiner: `h ^ (arity * 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2))` (golden ratio mixing, similar to `boost::hash_combine`).

### [SEV-3] Missing `[[nodiscard]]` on query methods and `create()`
- Location: `function_symbol_pool.h:34-35` (create), `function_symbol_pool.h:62-64` (as_int/as_list/as_empty_list), `function_symbol_pool.h:67-70` (size/capacity)
- Why this is a problem: Discarding the return of `create()` leaks the intent (was the caller expecting a side effect only?). Discarding `resize_is_needed()` is a logic error. These are must-use returns.
- Evidence type: static-proof
- Evidence: `create()` returns a `function_symbol`; discarding it means the element may be immediately garbage-collected (reference count drops to zero). This would be a bug at any call site.
- Efficiency impact: None (diagnostic improvement).
- C++20-specific note: `[[nodiscard]]` is standard since C++17 and expected on all query/factory methods per the mCRL2 coding conventions.
- Suggested fix: Add `[[nodiscard]]` to `create()`, `resize_is_needed()`, `size()`, `capacity()`, `get_sufficiently_large_postfix_index()`, and the accessor methods.

### [SEV-3] Missing Doxygen `\pre`/`\post` contracts on public interface
- Location: `function_symbol_pool.h:34-70` (all public methods)
- Why this is a problem: Per mCRL2 coding conventions, public APIs require explicit pre/postcondition documentation. `create()` doesn't document: what happens with empty name? Maximum arity? `sweep()` doesn't document the stop-the-world precondition. `resize_if_needed()` doesn't document thread-safety requirements.
- Evidence type: static-proof (convention violation)
- Efficiency impact: None (maintainability).
- C++20-specific note: N/A
- Suggested fix: Add `\pre` and `\post` clauses, especially the threading preconditions (`\pre No concurrent modifications to m_symbol_set` for `sweep()`).

### [SEV-3] Redundant `typename` before `std::conditional_t`
- Location: `function_symbol_pool.h:81`
- Why this is a problem: `std::conditional_t<...>` is already a type alias (not a dependent type needing disambiguation). The `typename` keyword is legal but misleading — it suggests template-dependent lookup when none exists.
- Evidence type: static-proof
- Efficiency impact: None (code clarity).
- C++20-specific note: In C++20, `typename` is not required in more contexts. Here it was never required since `std::conditional_t` is not a dependent name in this context (the class is not a template).
- Suggested fix: Remove `typename`:
  ```cpp
  std::conditional_t<EnableBlockAllocator, 
    mcrl2::utilities::block_allocator<...>, 
    std::allocator<_function_symbol>>,
  ```

## Proposed Tests
- [ ] `test_concurrent_create_sweep`: Spawn N threads calling `create()` with unique symbols while another thread calls `sweep()`. Run under TSan. Expected: either clean or TSan reports data race (confirming SEV-1).
- [ ] `test_register_prefix_exception_safety`: Mock an allocator that throws on the Kth allocation inside `register_prefix()`. Verify mutex is not held after the exception propagates.
- [ ] `test_hash_distribution`: Insert 1000 symbols with identical name and arities 0–999 into the pool. Measure max bucket chain length vs expected under uniform hashing.
- [ ] `test_create_nodiscard`: Static analysis pass (clang-tidy `[[nodiscard]]` check) on all callers of `create()`.

## Open Questions
1. Is there an explicit stop-the-world mechanism that guarantees no `create()` calls are in flight when `sweep()` runs? If so, where is it implemented and should it be documented as a precondition?
2. Why does `sweep()` hold `m_mutex` if `create()` doesn't respect it? Is the mutex only for `m_prefix_to_register_function_map` protection (comment says "Prevents changes to the symbol_set" which is inaccurate)?
3. Should `get_sufficiently_large_postfix_index()` be deprecated in favor of the incremental approach already partially implemented via `create_helper()`?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         9m 51s
Total session time:     11m 35s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         612.8k in, 7.9k out, 531.6k cached (Est. 1 Premium request)
 claude-haiku-4.5        1.0m in, 38.2k out, 926.0k cached (Est. 0 Premium requests)

```

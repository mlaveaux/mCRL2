---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_set_implementation.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T13:03:17.728018+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_set_implementation.h` — GC-guard wrappers for `atermpp::unordered_set` (wrapping `std::unordered_set`) and `atermpp::utilities::unordered_set` (wrapping `mcrl2::utilities::unordered_set`)
- Risk level: **High** — concurrency correctness issues in the lock-free data structure interaction
- Overall verdict: The shared-guard pattern correctly prevents GC during individual operations, but the `atermpp::utilities` variant has a latent use-after-free race between `erase` and lock-free `emplace`, and several operations have gaps where the thread is unprotected.

## Findings (ordered by severity)

### [SEV-1] Use-after-free: `erase_after` destroys node while `emplace_front_unique` iterates it

- Location: `libraries/utilities/include/mcrl2/utilities/detail/bucket_list.h:284-292` (exposed via this file's `erase()` wrappers at lines wrapping `super::erase()`)
- Why this is a problem: When `ThreadSafe=true`, `emplace_front_unique` performs a lock-free CAS loop that iterates bucket nodes (line 258: `for (auto it = iterator(old_head); it != old_end; ++it)`). Concurrently, `erase_after` unlinks a node with `set_next(next_node)` using `memory_order_relaxed` (line 288), then immediately destroys and deallocates the node (lines 291-292). A thread in `emplace_front_unique` that loaded `old_head` pointing to the to-be-erased node (via `m_head.next()` with acquire at line 244) can still be iterating that node when it is freed. The relaxed store in erase does not establish a happens-before relationship with the acquire load that preceded the iteration.
  
  In the atermpp::utilities wrapper, both `erase()` and `emplace()` acquire only a **shared** guard (lines in the `atermpp::utilities` namespace block), which provides GC protection but NOT mutual exclusion between threads. Two threads can concurrently execute erase+emplace on the same container.

- Evidence type: static-proof
- Evidence:
  - Interleaving:
    1. Thread A: `emplace_front_unique` → `old_head = m_head.next()` loads node X (acquire)
    2. Thread A: begins iterating, reads `X->key()` 
    3. Thread B: `erase_after` → `current_node->set_next(X->next())` (relaxed store)
    4. Thread B: `destroy(X)`, `deallocate(X)` → memory freed
    5. Thread A: dereferences X → **use-after-free**
  - Run command: `cmake -S . -B build-tsan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_THREADSANITIZER=ON -DMCRL2_ENABLE_TESTS=ON && cmake --build build-tsan -j$(nproc) && ctest --test-dir build-tsan -R pbesinst --output-on-failure`
  - Expected result: TSan reports a data race on the bucket node pointer / use-after-free from ASan
- Efficiency impact: N/A — this is a correctness/safety issue
- C++20-specific note: The `memory_order_relaxed` store in `set_next` is insufficient to publish the unlink before deallocation. A `release` store + corresponding `acquire` fence in iterating code, or an epoch-based reclamation / hazard pointer scheme, is needed for safe lock-free erase.
- Suggested fix: Either (a) protect `erase` with an exclusive lock (stop-the-world) rather than a shared guard, or (b) use deferred reclamation (e.g., epoch-based or RCU) so freed nodes are not deallocated while readers may still reference them, or (c) acquire the per-bucket mutex in `erase_after` when `ThreadSafe=true`.

### [SEV-2] `rehash_if_needed()` leaves thread unprotected in gap between unlock and exclusive lock

- Location: `unordered_set_implementation.h` (utilities namespace), `rehash_if_needed()` — the `guard.unlock_shared(); rehash(count * 2);` sequence
- Why this is a problem: Between `unlock_shared()` (clears the busy flag) and `rehash()`'s `lock()` (sets forbidden flags and waits), this thread is not in any shared section. The garbage collector can trigger during this window. While the container's `generic_aterm_container` registration ensures marking, there is a secondary issue: if GC triggers between unlock and lock, and the marking function traverses the container while `rehash()` is about to restructure it, there's a potential for the marker to see inconsistent state — though the exclusive lock in `rehash` should prevent this since GC also needs to stop all threads.
  
  The more concrete problem is a **TOCTOU**: the load factor is checked under shared lock, then rehash targets `count * 2`. Between the check and the rehash, other threads can insert elements, making the target bucket count insufficient.

- Evidence type: plausible-issue
- Evidence:
  - Mechanism: Thread A checks load_factor ≥ max_load_factor under shared lock, unlocks. Threads B, C insert N elements. Thread A rehashes to `count*2`, but load_factor might now exceed max_load_factor again immediately, causing another rehash on the next emplace. Under heavy contention this creates a rehash storm.
  - Run command: A concurrent stress test with many threads calling `emplace` simultaneously on a `ThreadSafe=true` instance
  - Expected result: Performance degradation from redundant rehashes; potentially observable as excessive time in exclusive lock acquisition
- Efficiency impact: Redundant rehashes under contention
- C++20-specific note: N/A
- Suggested fix: Re-check load factor inside `rehash()` after acquiring the exclusive lock; skip rehash if no longer needed.

### [SEV-2] `clear()` under shared guard is not safe for concurrent access

- Location: `unordered_set_implementation.h`, utilities namespace `clear()` method
- Why this is a problem: `clear()` iterates all buckets calling `erase` on each element, destroying and deallocating nodes. It acquires only a shared guard, which does not exclude other threads. If another thread is concurrently calling `find()` or `emplace()` on the same container, the nodes being iterated by the other thread may be freed. This has the same use-after-free mechanism as SEV-1.
- Evidence type: static-proof (same race mechanism as SEV-1, amplified to all buckets)
- Evidence:
  - Reproducer: Thread A calls `clear()`, Thread B calls `find()` — B iterates a bucket whose nodes A is destroying
  - Run command: TSan build + concurrent test with clear+find
  - Expected result: TSan data race report
- Efficiency impact: N/A — correctness issue
- C++20-specific note: N/A
- Suggested fix: `clear()` should use an exclusive lock (`detail::g_thread_term_pool().lock()`) when `ThreadSafe=true`, since it is a destructive bulk operation.

### [SEV-3] Iterators escape lock scope — lifetime validity unguarded

- Location: All `begin()`, `end()`, `cbegin()`, `cend()`, and `find()` methods in both namespaces
- Why this is a problem: These methods acquire a shared guard, obtain an iterator from the base class, then release the guard (on scope exit) before returning the iterator to the caller. The caller then operates on the iterator without any GC protection. For the `atermpp::unordered_set` (std-based, single-threaded), this is acceptable because the container is GC-registered and elements won't be collected. For `atermpp::utilities::unordered_set` with `ThreadSafe=true`, a concurrent `erase` can invalidate the returned iterator before the caller dereferences it.
- Evidence type: plausible-issue
- Evidence:
  - Mechanism: Thread A calls `find()`, gets iterator to element E (lock released). Thread B calls `erase(E)`. Thread A dereferences iterator → dangling pointer.
  - Validation plan: Stress test with concurrent find+erase; ASan should report heap-use-after-free.
- C++20-specific note: This is a case where `[[clang::lifetimebound]]` annotations cannot help because the lifetime issue is cross-thread, not syntactic.
- Suggested fix: Document that callers must hold their own shared guard for the duration of iterator use, or provide a `find_and_act` callback-style API that holds the lock.

### [SEV-3] `size()` inconsistently bypasses the lock pattern

- Location: `unordered_set.h` lines 172-175 and 277-280 (both `size()` definitions)
- Why this is a problem: Every other method in this file acquires a shared guard before delegating to `super`. The `size()` method does not. While for the utilities variant (where `m_number_of_elements` is `std::atomic` when `ThreadSafe=true`) this is technically safe as an atomic read, the inconsistency is confusing and for the `atermpp::unordered_set` wrapping `std::unordered_set` there is no atomicity guarantee if size is called while another thread mutates (though single-threaded use is assumed).
- Evidence type: plausible-issue (inconsistency; actual race depends on usage)
- Suggested fix: Either add the shared guard for consistency, or add a code comment documenting why `size()` intentionally skips it.

### [SEV-3] No API contracts or thread-safety documentation

- Location: All methods in this file
- Why this is a problem: Per mCRL2 coding standards (`.github/instructions/cpp.instructions.md`), public APIs require Doxygen `\pre`/`\post` documentation. None of these methods document their thread-safety guarantees, the relationship between the shared guard and GC protection, or the caller's responsibilities for iterator lifetime.
- Evidence type: static-proof (absence of documentation)
- C++20-specific note: This is a missed opportunity to use `[[nodiscard]]` on iterator-returning methods and to document the GC-protection contract.
- Suggested fix: Add `/// \brief` + `\pre`/`\post` for each method, documenting that the shared guard prevents GC (not concurrent access), and that the container must not be concurrently modified unless `ThreadSafe=true` in the underlying type.

## Proposed Tests
- [ ] `test_concurrent_erase_emplace`: Two threads, one continuously emplacing keys, one continuously erasing. Container is `atermpp::utilities::unordered_set<aterm, ..., true>`. Run under TSan. Expected: TSan reports data race on bucket node pointers.
- [ ] `test_concurrent_clear_find`: One thread calls `clear()` in a loop, another calls `find()`. Run under ASan. Expected: heap-use-after-free.
- [ ] `test_rehash_under_contention`: N threads emplacing simultaneously triggering rehash. Measure rehash count. Expected: no correctness failure, but validates the TOCTOU concern.
- [ ] `test_iterator_outlives_erase`: Thread A gets iterator from `find()`, Thread B erases that key, Thread A dereferences. Run under ASan. Expected: use-after-free if bug is real.

## Open Questions
1. Is `atermpp::utilities::unordered_set` ever instantiated with `ThreadSafe=true` in production code? The only visible `ThreadSafe=true` usage is through `atermpp::indexed_set` in `pbesinst_lazy.h`, which takes a different code path (own hash table in `indexed_set`). If `ThreadSafe=true` is never used for this specific unordered_set wrapper, SEV-1/SEV-2/SEV-3 iterator issues are latent only.
2. Is `erase` on `mcrl2::utilities::unordered_set<..., ThreadSafe=true>` intended to be called only under an exclusive (stop-the-world) lock at a higher level? If so, the lock-free CAS in emplace and the non-atomic erase are architecturally consistent, but this invariant is completely undocumented and unenforced.
3. The `set_next` uses `memory_order_relaxed` — was this a deliberate performance optimization assuming erase is always serialized, or an oversight?

## Copilot Raw Stdout

```text
● Explore: Explore locking and container context
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h lines 1-80
    └ 80 lines read
  ● Get file sizes
    $ wc -l /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/shared_mutex.h
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_set.h
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h
    └ 5 lines...

● Read explore agent output
  $ cat /tmp/1783601686503-copilot-tool-output-sxgic8.txt
  └ 8 lines...

● Read first part of explore output
  $ head -200 /tmp/copilot-tool-output-1783601690914-49kl5u.txt
  └ 201 lines...

● Read remaining explore output
  $ sed -n '200,400p' /tmp/copilot-tool-output-1783601690914-49kl5u.txt
  └ 201 lines...

● Read more explore output
  $ sed -n '400,500p' /tmp/copilot-tool-output-1783601690914-49kl5u.txt
  └ 101 lines...

● Read final part of explore output
  $ sed -n '500,600p' /tmp/copilot-tool-output-1783601690914-49kl5u.txt
  └ 102 lines...

● Read last part of explore output
  $ sed -n '600,700p' /tmp/copilot-tool-output-1783601690914-49kl5u.txt
  └ 102 lines...

Now let me check the underlying `mcrl2::utilities::unordered_set` to understand its thread-safety model, and also check `size()` which appears unprotected:

● Explore: Check utilities unordered_set base
  ● Grep "ThreadSafe|fetch_add|lock_guard|emplace_front_unique" (/workspaces/mCRL2/libraries/utiliti
  es/include/mcrl2/utilities/detail/unordered_set_implementation.h)
    └ 7 lines found
  ● Read
    libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_set_implementation.h
     lines 1-80
    └ 80 lines read

● Check m_number_of_elements usage
  $ grep -n "m_number_of_elements"
  /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/unordered_set.h | head -20
  └ 3 lines...

● Check erase declarations
  $ grep -n "erase" /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/unordered_set.h |
  head -20
  └ 4 lines...

● Read libraries/utilities/include/mcrl2/utilities/unordered_set.h lines 285-300
  └ 15 lines read

● Read libraries/utilities/include/mcrl2/utilities/unordered_set.h lines 356-395
  └ 40 lines read

✗ Find erase implementation file
  $ find /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities -name "*.h" | xargs grep -l
  "erase_impl\|erase.*const_iterator" | head -5
  Permission denied and could not request permission from user

● Glob "libraries/utilities/include/mcrl2/utilities/detail/unordered_set*"
  └ 1 file found

● Check erase impl and counter usage
  $ grep -n "erase\|m_number_of_elements" /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilit
  ies/detail/unordered_set_implementation.h | head -30
  └ 16 lines...

● Read libraries/utilities/include/mcrl2/utilities/detail/unordered_set_implementation.h lines
  118-165
  └ 48 lines read

● Read libraries/utilities/include/mcrl2/utilities/detail/unordered_set_implementation.h lines
  315-345
  └ 29 lines read

● Read libraries/utilities/include/mcrl2/utilities/detail/unordered_set_implementation.h lines
  250-320
  └ 71 lines read

Now let me check the bucket type's `erase_after` and `emplace_front_unique` for thread safety:

● Find bucket implementations
  $ grep -rn "emplace_front_unique\|erase_after"
  /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/ --include="*.h" | grep -v
  "unordered_set_implementation" | head -20
  └ 7 lines...

● Read libraries/utilities/include/mcrl2/utilities/detail/bucket_list.h lines 237-310
  └ 73 lines read

● Check atomic ops in bucket_list
  $ grep -n "exchange\|atomic\|compare_exchange\|CAS\|memory_order"
  /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/detail/bucket_list.h | head -20
  └ 8 lines...

● Read libraries/utilities/include/mcrl2/utilities/detail/unordered_set_implementation.h lines 65-95
  └ 31 lines read

● Read libraries/utilities/include/mcrl2/utilities/detail/bucket_list.h lines 50-75
  └ 25 lines read

● Read libraries/utilities/include/mcrl2/utilities/detail/bucket_list.h lines 200-220
  └ 21 lines read

● Find usages of utilities::unordered_set
  $ grep -rn "atermpp::utilities::unordered_set\|utilities::unordered_set<"
  /workspaces/mCRL2/libraries/ --include="*.h" --include="*.cpp" | grep -v "test/" | grep -v
  "/detail/" | head -15
  └ 6 lines...

● Find ThreadSafe=true usages
  $ grep -rn "atermpp::utilities::unordered_set" /workspaces/mCRL2/libraries/
  /workspaces/mCRL2/tools/ --include="*.h" --include="*.cpp" | grep -v "test/" | grep "true" | head
  -10
  └ 1 line...

● Find all usages of utilities unordered_set
  $ grep -rn "utilities::unordered_set<" /workspaces/mCRL2/libraries/ /workspaces/mCRL2/tools/
  --include="*.h" --include="*.cpp" | grep -v "test/" | grep -v "detail/" | grep -v
  "standard_containers" | head -15
  └ 3 lines...

● Find ThreadSafe=true instantiations
  $ grep -rn "unordered_set.*ThreadSafe\|ThreadSafe.*true\|unordered_set.*true"
  /workspaces/mCRL2/libraries/atermpp/ --include="*.h" --include="*.cpp" | grep -v "detail/unordered
  _set_implementation\|standard_containers/unordered_set.h\|standard_containers/detail" | head -10
  └ 1 line...

● Search for ThreadSafe=true instantiation
  $ grep -rn "unordered_set.*true" /workspaces/mCRL2/libraries/ /workspaces/mCRL2/tools/
  --include="*.h" --include="*.cpp" | grep -i "thread\|safe" | head -10
  └ 1 line...

● Search broader for atermpp utilities usage
  $ grep -rn "atermpp.*utilities.*unordered_set\|atermpp::utilities::unordered_set"
  /workspaces/mCRL2/ --include="*.h" --include="*.cpp" | grep -v "standard_containers" | grep -v
  "detail/unordered" | head -10
  └ 1 line...

● Check if indexed_set uses utilities unordered_set
  $ grep -rn "indexed_set" /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/
  --include="*.h" | head -10
  └ 11 lines...

● Find ThreadSafe indexed_set usage
  $ grep -rn "indexed_set.*true\|indexed_set<.*,.*true" /workspaces/mCRL2/libraries/
  /workspaces/mCRL2/tools/ --include="*.h" --include="*.cpp" | grep -v "test/" | head -10
  └ 2 lines...

● Check if indexed_set uses unordered_set
  $ grep -n "unordered_set"
  /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/indexed_set.h | head -5
  └ 1 line...

● Check indexed_set header
  $ head -40 /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/indexed_set.h
  └ 41 lines...

Now I have all the context needed. Let me verify one more thing about the `set_next` memory ordering issue:

● Read libraries/utilities/include/mcrl2/utilities/detail/bucket_list.h lines 278-296
  └ 18 lines read

● Read libraries/utilities/include/mcrl2/utilities/detail/bucket_list.h lines 237-275
  └ 39 lines read

I now have complete context for the review. Here is my analysis:

---

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_set_implementation.h` — GC-guard wrappers for `atermpp::unordered_set` (wrapping `std::unordered_set`) and `atermpp::utilities::unordered_set` (wrapping `mcrl2::utilities::unordered_set`)
- Risk level: **High** — concurrency correctness issues in the lock-free data structure interaction
- Overall verdict: The shared-guard pattern correctly prevents GC during individual operations, but the `atermpp::utilities` variant has a latent use-after-free race between `erase` and lock-free `emplace`, and several operations have gaps where the thread is unprotected.

## Findings (ordered by severity)

### [SEV-1] Use-after-free: `erase_after` destroys node while `emplace_front_unique` iterates it

- Location: `libraries/utilities/include/mcrl2/utilities/detail/bucket_list.h:284-292` (exposed via this file's `erase()` wrappers at lines wrapping `super::erase()`)
- Why this is a problem: When `ThreadSafe=true`, `emplace_front_unique` performs a lock-free CAS loop that iterates bucket nodes (line 258: `for (auto it = iterator(old_head); it != old_end; ++it)`). Concurrently, `erase_after` unlinks a node with `set_next(next_node)` using `memory_order_relaxed` (line 288), then immediately destroys and deallocates the node (lines 291-292). A thread in `emplace_front_unique` that loaded `old_head` pointing to the to-be-erased node (via `m_head.next()` with acquire at line 244) can still be iterating that node when it is freed. The relaxed store in erase does not establish a happens-before relationship with the acquire load that preceded the iteration.
  
  In the atermpp::utilities wrapper, both `erase()` and `emplace()` acquire only a **shared** guard (lines in the `atermpp::utilities` namespace block), which provides GC protection but NOT mutual exclusion between threads. Two threads can concurrently execute erase+emplace on the same container.

- Evidence type: static-proof
- Evidence:
  - Interleaving:
    1. Thread A: `emplace_front_unique` → `old_head = m_head.next()` loads node X (acquire)
    2. Thread A: begins iterating, reads `X->key()` 
    3. Thread B: `erase_after` → `current_node->set_next(X->next())` (relaxed store)
    4. Thread B: `destroy(X)`, `deallocate(X)` → memory freed
    5. Thread A: dereferences X → **use-after-free**
  - Run command: `cmake -S . -B build-tsan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_THREADSANITIZER=ON -DMCRL2_ENABLE_TESTS=ON && cmake --build build-tsan -j$(nproc) && ctest --test-dir build-tsan -R pbesinst --output-on-failure`
  - Expected result: TSan reports a data race on the bucket node pointer / use-after-free from ASan
- Efficiency impact: N/A — this is a correctness/safety issue
- C++20-specific note: The `memory_order_relaxed` store in `set_next` is insufficient to publish the unlink before deallocation. A `release` store + corresponding `acquire` fence in iterating code, or an epoch-based reclamation / hazard pointer scheme, is needed for safe lock-free erase.
- Suggested fix: Either (a) protect `erase` with an exclusive lock (stop-the-world) rather than a shared guard, or (b) use deferred reclamation (e.g., epoch-based or RCU) so freed nodes are not deallocated while readers may still reference them, or (c) acquire the per-bucket mutex in `erase_after` when `ThreadSafe=true`.

### [SEV-2] `rehash_if_needed()` leaves thread unprotected in gap between unlock and exclusive lock

- Location: `unordered_set_implementation.h` (utilities namespace), `rehash_if_needed()` — the `guard.unlock_shared(); rehash(count * 2);` sequence
- Why this is a problem: Between `unlock_shared()` (clears the busy flag) and `rehash()`'s `lock()` (sets forbidden flags and waits), this thread is not in any shared section. The garbage collector can trigger during this window. While the container's `generic_aterm_container` registration ensures marking, there is a secondary issue: if GC triggers between unlock and lock, and the marking function traverses the container while `rehash()` is about to restructure it, there's a potential for the marker to see inconsistent state — though the exclusive lock in `rehash` should prevent this since GC also needs to stop all threads.
  
  The more concrete problem is a **TOCTOU**: the load factor is checked under shared lock, then rehash targets `count * 2`. Between the check and the rehash, other threads can insert elements, making the target bucket count insufficient.

- Evidence type: plausible-issue
- Evidence:
  - Mechanism: Thread A checks load_factor ≥ max_load_factor under shared lock, unlocks. Threads B, C insert N elements. Thread A rehashes to `count*2`, but load_factor might now exceed max_load_factor again immediately, causing another rehash on the next emplace. Under heavy contention this creates a rehash storm.
  - Run command: A concurrent stress test with many threads calling `emplace` simultaneously on a `ThreadSafe=true` instance
  - Expected result: Performance degradation from redundant rehashes; potentially observable as excessive time in exclusive lock acquisition
- Efficiency impact: Redundant rehashes under contention
- C++20-specific note: N/A
- Suggested fix: Re-check load factor inside `rehash()` after acquiring the exclusive lock; skip rehash if no longer needed.

### [SEV-2] `clear()` under shared guard is not safe for concurrent access

- Location: `unordered_set_implementation.h`, utilities namespace `clear()` method
- Why this is a problem: `clear()` iterates all buckets calling `erase` on each element, destroying and deallocating nodes. It acquires only a shared guard, which does not exclude other threads. If another thread is concurrently calling `find()` or `emplace()` on the same container, the nodes being iterated by the other thread may be freed. This has the same use-after-free mechanism as SEV-1.
- Evidence type: static-proof (same race mechanism as SEV-1, amplified to all buckets)
- Evidence:
  - Reproducer: Thread A calls `clear()`, Thread B calls `find()` — B iterates a bucket whose nodes A is destroying
  - Run command: TSan build + concurrent test with clear+find
  - Expected result: TSan data race report
- Efficiency impact: N/A — correctness issue
- C++20-specific note: N/A
- Suggested fix: `clear()` should use an exclusive lock (`detail::g_thread_term_pool().lock()`) when `ThreadSafe=true`, since it is a destructive bulk operation.

### [SEV-3] Iterators escape lock scope — lifetime validity unguarded

- Location: All `begin()`, `end()`, `cbegin()`, `cend()`, and `find()` methods in both namespaces
- Why this is a problem: These methods acquire a shared guard, obtain an iterator from the base class, then release the guard (on scope exit) before returning the iterator to the caller. The caller then operates on the iterator without any GC protection. For the `atermpp::unordered_set` (std-based, single-threaded), this is acceptable because the container is GC-registered and elements won't be collected. For `atermpp::utilities::unordered_set` with `ThreadSafe=true`, a concurrent `erase` can invalidate the returned iterator before the caller dereferences it.
- Evidence type: plausible-issue
- Evidence:
  - Mechanism: Thread A calls `find()`, gets iterator to element E (lock released). Thread B calls `erase(E)`. Thread A dereferences iterator → dangling pointer.
  - Validation plan: Stress test with concurrent find+erase; ASan should report heap-use-after-free.
- C++20-specific note: This is a case where `[[clang::lifetimebound]]` annotations cannot help because the lifetime issue is cross-thread, not syntactic.
- Suggested fix: Document that callers must hold their own shared guard for the duration of iterator use, or provide a `find_and_act` callback-style API that holds the lock.

### [SEV-3] `size()` inconsistently bypasses the lock pattern

- Location: `unordered_set.h` lines 172-175 and 277-280 (both `size()` definitions)
- Why this is a problem: Every other method in this file acquires a shared guard before delegating to `super`. The `size()` method does not. While for the utilities variant (where `m_number_of_elements` is `std::atomic` when `ThreadSafe=true`) this is technically safe as an atomic read, the inconsistency is confusing and for the `atermpp::unordered_set` wrapping `std::unordered_set` there is no atomicity guarantee if size is called while another thread mutates (though single-threaded use is assumed).
- Evidence type: plausible-issue (inconsistency; actual race depends on usage)
- Suggested fix: Either add the shared guard for consistency, or add a code comment documenting why `size()` intentionally skips it.

### [SEV-3] No API contracts or thread-safety documentation

- Location: All methods in this file
- Why this is a problem: Per mCRL2 coding standards (`.github/instructions/cpp.instructions.md`), public APIs require Doxygen `\pre`/`\post` documentation. None of these methods document their thread-safety guarantees, the relationship between the shared guard and GC protection, or the caller's responsibilities for iterator lifetime.
- Evidence type: static-proof (absence of documentation)
- C++20-specific note: This is a missed opportunity to use `[[nodiscard]]` on iterator-returning methods and to document the GC-protection contract.
- Suggested fix: Add `/// \brief` + `\pre`/`\post` for each method, documenting that the shared guard prevents GC (not concurrent access), and that the container must not be concurrently modified unless `ThreadSafe=true` in the underlying type.

## Proposed Tests
- [ ] `test_concurrent_erase_emplace`: Two threads, one continuously emplacing keys, one continuously erasing. Container is `atermpp::utilities::unordered_set<aterm, ..., true>`. Run under TSan. Expected: TSan reports data race on bucket node pointers.
- [ ] `test_concurrent_clear_find`: One thread calls `clear()` in a loop, another calls `find()`. Run under ASan. Expected: heap-use-after-free.
- [ ] `test_rehash_under_contention`: N threads emplacing simultaneously triggering rehash. Measure rehash count. Expected: no correctness failure, but validates the TOCTOU concern.
- [ ] `test_iterator_outlives_erase`: Thread A gets iterator from `find()`, Thread B erases that key, Thread A dereferences. Run under ASan. Expected: use-after-free if bug is real.

## Open Questions
1. Is `atermpp::utilities::unordered_set` ever instantiated with `ThreadSafe=true` in production code? The only visible `ThreadSafe=true` usage is through `atermpp::indexed_set` in `pbesinst_lazy.h`, which takes a different code path (own hash table in `indexed_set`). If `ThreadSafe=true` is never used for this specific unordered_set wrapper, SEV-1/SEV-2/SEV-3 iterator issues are latent only.
2. Is `erase` on `mcrl2::utilities::unordered_set<..., ThreadSafe=true>` intended to be called only under an exclusive (stop-the-world) lock at a higher level? If so, the lock-free CAS in emplace and the non-atomic erase are architecturally consistent, but this invariant is completely undocumented and unenforced.
3. The `set_next` uses `memory_order_relaxed` — was this a deliberate performance optimization assuming erase is always serialized, or an oversight?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         10m 11s
Total session time:     10m 47s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         2.2m in, 21.3k out, 2.1m cached (Est. 1 Premium request)
 claude-haiku-4.5        367.9k in, 12.7k out, 301.6k cached (Est. 0 Premium requests)

```

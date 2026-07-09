---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T12:52:26.287373+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h` — full file, clean review
- **Risk level:** Medium-High
- **Overall verdict:** The file implements GC-coordination wrappers for two `unordered_map` variants. It has a TOCTOU race in `rehash_if_needed()`, significant inconsistencies in which operations acquire the GC protection lock (several `std::unordered_map` methods are inherited without protection), and missing API contracts. The code is functional for its intended single-owner usage pattern but the protection boundaries are incomplete and undocumented.

## Findings (ordered by severity)

### [SEV-1] `operator[]`, `at()`, `find()`, `count()` bypass GC protection in `atermpp::unordered_map`
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h` (class declaration, lines 31–229)
- **Why this is a problem:** The class publicly inherits from `std::unordered_map` and overrides mutating methods (`insert`, `erase`, `clear`) and iterators (`begin`/`end`) with `lock_shared()` GC protection. However, `operator[]`, `at()`, `find()`, `count()`, `contains()`, `equal_range()`, and `rehash()`/`reserve()` are NOT overridden and are directly callable on the derived class without any lock. `operator[]` is particularly dangerous because it can insert a default-constructed element (mutation) entirely bypassing the GC protection that every other mutating method acquires. If the GC runs concurrently (via another thread acquiring exclusive), this thread would be modifying the container while the GC is scanning it.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: Grep of class declaration confirms no `operator[]`/`at()`/`find()` override.
  - Run command: `grep -n 'operator\[\]\|at(\|find(' libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h`
  - Observed result: Only `find()` for the `utilities` variant (line 334); none for `atermpp::unordered_map`.
  - Expected result: All methods that access or modify aterm-containing elements should hold `lock_shared()`.
- **Efficiency impact:** None (this is a correctness issue).
- **C++20-specific note:** Public inheritance from standard library containers is generally discouraged precisely because you cannot intercept all calls without virtual dispatch. A composition-with-forwarding approach or a `private` inheritance with explicit `using` declarations for safe methods would prevent this class of bugs.
- **Suggested fix:** Either (a) override `operator[]`, `at()`, `find()`, `count()`, `contains()` with lock-protected versions, or (b) switch to private/protected inheritance and explicitly expose only the wrapped API, or (c) document prominently that these inherited methods MUST NOT be called directly.

---

### [SEV-2] TOCTOU race in `rehash_if_needed()` — unprotected window between unlock and exclusive lock
- **Location:** `unordered_map_implementation.h`, lines ~255–264 (in `atermpp::utilities` section)
- **Why this is a problem:** `rehash_if_needed()` acquires a shared lock, reads `bucket_count()`, calls `guard.unlock_shared()`, then calls `rehash(count * 2)`. Between `unlock_shared()` and `rehash()` acquiring the exclusive lock, this thread holds NO GC protection. A concurrent GC cycle can execute during this window. While `container_wrapper` keeps map contents reachable, the thread may be executing code that reads internal state (the `count` variable is already read, but the map's internal pointers are accessible without protection). Additionally, multiple threads can simultaneously detect `rehash_is_needed()`, each read the same stale `count`, and queue redundant exclusive-lock acquisitions.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Thread A reads `bucket_count()=64`, unlocks. Thread B also reads 64, unlocks. Both call `rehash(128)`. First succeeds; second either no-ops (if idempotent) or does redundant work. The unprotected window between unlock and exclusive acquisition is the real concern for GC safety.
  - Test plan: Run a TSan build with N threads concurrently inserting into a shared `atermpp::utilities::unordered_map<aterm, aterm, ..., true>` near load-factor threshold. Check for TSan reports.
  - Run command: `cmake --build build-tsan && ctest --test-dir build-tsan -R unordered_map`
  - Expected failure: TSan data-race report on bucket_count read vs. rehash modification, or redundant rehash observed under contention.
- **Efficiency impact:** Redundant rehash operations under contention (O(n) per redundant rehash).
- **Suggested fix:** Use a compare-and-swap or check-after-exclusive-acquire pattern:
  ```cpp
  void rehash_if_needed() {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    if (rehash_is_needed()) {
      std::size_t count = super::bucket_count();
      guard.unlock_shared();
      mcrl2::utilities::lock_guard ex = detail::g_thread_term_pool().lock();
      // Re-check under exclusive lock
      if (super::load_factor() >= super::max_load_factor()) {
        super::rehash(super::bucket_count() * 2);
      }
    }
  }
  ```

---

### [SEV-2] `size()` has no lock while `empty()` and `max_size()` do — inconsistent protection boundary
- **Location:** `unordered_map.h:204–207` (inline `size()`) vs `unordered_map_implementation.h` (`empty()`, `max_size()` with locks)
- **Why this is a problem:** `size()` is defined inline in the header as `return super::size();` with NO lock. `empty()` (which is semantically `size() == 0`) IS wrapped with `lock_shared()`. Either `size()` needs the lock (GC could run during the call, though the risk for a simple read is minimal) or `empty()`/`max_size()` don't need it. The inconsistency indicates confusion about the protection invariant and makes reasoning about safety impossible for maintainers.
- **Evidence type:** static-proof
- **Evidence:**
  - `unordered_map.h:204`: `std::size_t size() const { return super::size(); }` — no lock.
  - `unordered_map_implementation.h`: `bool empty() const noexcept { lock_shared(); return super::empty(); }` — locked.
  - These are semantically equivalent operations with different protection.
- **Efficiency impact:** Negligible.
- **Suggested fix:** Make `size()` consistent with `empty()` — either both lock or neither locks. Given the overall design philosophy of "always lock during container access," add the lock to `size()`.

---

### [SEV-2] `swap()` only protects `*this`, not `other` — potential GC inconsistency
- **Location:** `unordered_map_implementation.h`, both `atermpp::unordered_map::swap` and `atermpp::utilities::unordered_map::swap`
- **Why this is a problem:** `swap(other)` acquires `lock_shared()` on `*this` thread's term pool, then calls `super::swap(other)`. If `other` is owned by a different thread (rare but possible), `other`'s terms are being moved into `*this` without any lock held on `other`'s thread. More practically, `swap` modifies BOTH containers but only signals GC protection for one thread. If GC triggers on `other`'s owning thread mid-swap, the container contents are in an inconsistent state.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Thread A calls `map1.swap(map2)` while Thread B is iterating `map2` (both hold shared locks on their own thread pools). The swap modifies `map2` without Thread B's knowledge.
  - Test plan: Two threads sharing map references, one swapping while the other iterates, under TSan.
  - Expected failure: TSan report on concurrent modification of the swapped container.
- **Efficiency impact:** None (correctness issue).
- **Suggested fix:** Document that `swap` requires exclusive external synchronization on both maps, or acquire locks on both thread pools.

---

### [SEV-2] `find()` in `utilities` variant conditionally locks but all other read-only operations unconditionally lock
- **Location:** `unordered_map_implementation.h`, `find()` vs `begin()`/`end()`/`empty()`
- **Why this is a problem:** When `ThreadSafe=false` in a multithreaded build (`MCRL2_ENABLE_MULTITHREADING=ON`), `find()` skips the shared lock entirely. But `begin()`, `end()`, `empty()`, `max_size()` still acquire it unconditionally. Either `find()` is unsafe without the lock (GC can run during find, corrupting term access) or the other read-only methods are paying unnecessary lock overhead. Given that `find()` dereferences aterms for hash/equality comparison, it has the same GC-interaction requirements as iteration.
- **Evidence type:** static-proof
- **Evidence:**
  - `find()`: `if constexpr (ThreadSafe) { lock; find; } return super::find(args...);` — no lock when `ThreadSafe=false`
  - `begin()`: always locks regardless of `ThreadSafe`
  - In a multithreaded build with `ThreadSafe=false` instantiation, find() exposes the thread to GC interference.
- **C++20-specific note:** The `if constexpr` pattern is correct C++17/20 usage, but the semantic inconsistency between methods is the issue.
- **Suggested fix:** Either gate ALL methods on `ThreadSafe` (consistent performance optimization) or gate NONE (uniform safety). The current mix is the worst of both worlds.

---

### [SEV-3] No Doxygen `\pre`/`\post` on any public method
- **Location:** Entire file — only one `\brief` comment (line 32)
- **Why this is a problem:** Per mCRL2 coding standards, all public API methods require Doxygen documentation with `\pre`/`\post` clauses. The thread-safety requirements (caller must ensure no concurrent access to the same map instance from multiple threads) are completely undocumented. The relationship between the GC lock and the container's thread-safety model is nowhere explained.
- **Evidence type:** static-proof
- **Suggested fix:** Add at minimum:
  - `\pre` The calling thread must not already hold an exclusive lock on its term pool.
  - `\post` The GC has been prevented from running during the operation.
  - Class-level documentation explaining the locking model and that external synchronization is required for multi-threaded access to the same map.

---

### [SEV-3] Trailing whitespace and tab characters in several lines
- **Location:** Multiple lines (e.g., line with `erase(const_iterator pos)` has a tab character after the closing brace)
- **Why this is a problem:** Does not comply with `.clang-format` clean requirement for touched code.
- **Evidence type:** static-proof
- **Suggested fix:** Run `clang-format` on the file.

---

## Proposed Tests
- [ ] `test_unordered_map_operator_bracket_gc_safety`: Multi-threaded test calling `operator[]` on `atermpp::unordered_map` while forcing GC. Expected: TSan clean (or reproduces a race). Command: build with `-DMCRL2_ENABLE_THREADSANITIZER=ON`, run test.
- [ ] `test_rehash_if_needed_toctou`: N threads inserting into `atermpp::utilities::unordered_map<aterm,aterm,...,true>` at load-factor boundary simultaneously. Verify no double-rehash or TSan report.
- [ ] `test_swap_cross_thread`: Two threads, each owning a map, swap maps between threads under TSan. Verify no data race.
- [ ] `test_find_no_lock_gc_interference`: Instantiate `atermpp::utilities::unordered_map<..., false>` in a multithreaded build, call `find()` while another thread triggers GC. Verify no crash or TSan report.

## Open Questions
1. Is the omission of `operator[]`/`at()`/`find()` overrides in `atermpp::unordered_map` intentional? If so, what invariant makes them safe without the lock?
2. Is `size()` deliberately unlocked for performance while `empty()` is locked? What is the reasoning?
3. For `atermpp::utilities::unordered_map` with `ThreadSafe=false` in a multithreaded build — is this combination expected to be used, and if so, what external guarantees does the caller provide?
4. Should `rehash_if_needed()` re-verify the condition after acquiring exclusive access to avoid redundant work under contention?

## Copilot Raw Stdout

```text
● Explore: Explore unordered_map declarations
  ● Get line count of implementation file
    $ wc -l /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/uno
    rdered_map_implementation.h
    └ 2 lines...
  ● Read
    libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h
     lines 1-50
    └ 50 lines read

● Explore: Explore utilities unordered_map base
  ● Check find method signature
    $ grep -A 10 "inline.*iterator find"
    /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/unordered_set.h
    └ 1 line...
  ● Create comprehensive analysis document
    $ cat > /tmp/analysis.md << 'EOF'
    # Analysis of mcrl2::utilities::unordered_map and shared_mutex

    ## 1. Full Class Declaration of mcrl2::utilities::unordered_map...
    └ 449 lines...

Now I have full context. Let me verify one more thing about the `swap` semantics and iterator patterns:

● Grep "operator\[\]|T&\s+at\("
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h)
  └ No matches found

● Grep "generic_aterm_container"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h)
  └ 9 lines found

● Grep "operator\[|at\(|count\(|contains\(|size\(\)"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h)
  └ 6 lines found

● Grep "find"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h)
  └ 4 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h lines 195-230
  └ 35 lines read

Now I have complete context. Here is the review:

---

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h` — full file, clean review
- **Risk level:** Medium-High
- **Overall verdict:** The file implements GC-coordination wrappers for two `unordered_map` variants. It has a TOCTOU race in `rehash_if_needed()`, significant inconsistencies in which operations acquire the GC protection lock (several `std::unordered_map` methods are inherited without protection), and missing API contracts. The code is functional for its intended single-owner usage pattern but the protection boundaries are incomplete and undocumented.

## Findings (ordered by severity)

### [SEV-1] `operator[]`, `at()`, `find()`, `count()` bypass GC protection in `atermpp::unordered_map`
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h` (class declaration, lines 31–229)
- **Why this is a problem:** The class publicly inherits from `std::unordered_map` and overrides mutating methods (`insert`, `erase`, `clear`) and iterators (`begin`/`end`) with `lock_shared()` GC protection. However, `operator[]`, `at()`, `find()`, `count()`, `contains()`, `equal_range()`, and `rehash()`/`reserve()` are NOT overridden and are directly callable on the derived class without any lock. `operator[]` is particularly dangerous because it can insert a default-constructed element (mutation) entirely bypassing the GC protection that every other mutating method acquires. If the GC runs concurrently (via another thread acquiring exclusive), this thread would be modifying the container while the GC is scanning it.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer/test name: Grep of class declaration confirms no `operator[]`/`at()`/`find()` override.
  - Run command: `grep -n 'operator\[\]\|at(\|find(' libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h`
  - Observed result: Only `find()` for the `utilities` variant (line 334); none for `atermpp::unordered_map`.
  - Expected result: All methods that access or modify aterm-containing elements should hold `lock_shared()`.
- **Efficiency impact:** None (this is a correctness issue).
- **C++20-specific note:** Public inheritance from standard library containers is generally discouraged precisely because you cannot intercept all calls without virtual dispatch. A composition-with-forwarding approach or a `private` inheritance with explicit `using` declarations for safe methods would prevent this class of bugs.
- **Suggested fix:** Either (a) override `operator[]`, `at()`, `find()`, `count()`, `contains()` with lock-protected versions, or (b) switch to private/protected inheritance and explicitly expose only the wrapped API, or (c) document prominently that these inherited methods MUST NOT be called directly.

---

### [SEV-2] TOCTOU race in `rehash_if_needed()` — unprotected window between unlock and exclusive lock
- **Location:** `unordered_map_implementation.h`, lines ~255–264 (in `atermpp::utilities` section)
- **Why this is a problem:** `rehash_if_needed()` acquires a shared lock, reads `bucket_count()`, calls `guard.unlock_shared()`, then calls `rehash(count * 2)`. Between `unlock_shared()` and `rehash()` acquiring the exclusive lock, this thread holds NO GC protection. A concurrent GC cycle can execute during this window. While `container_wrapper` keeps map contents reachable, the thread may be executing code that reads internal state (the `count` variable is already read, but the map's internal pointers are accessible without protection). Additionally, multiple threads can simultaneously detect `rehash_is_needed()`, each read the same stale `count`, and queue redundant exclusive-lock acquisitions.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Thread A reads `bucket_count()=64`, unlocks. Thread B also reads 64, unlocks. Both call `rehash(128)`. First succeeds; second either no-ops (if idempotent) or does redundant work. The unprotected window between unlock and exclusive acquisition is the real concern for GC safety.
  - Test plan: Run a TSan build with N threads concurrently inserting into a shared `atermpp::utilities::unordered_map<aterm, aterm, ..., true>` near load-factor threshold. Check for TSan reports.
  - Run command: `cmake --build build-tsan && ctest --test-dir build-tsan -R unordered_map`
  - Expected failure: TSan data-race report on bucket_count read vs. rehash modification, or redundant rehash observed under contention.
- **Efficiency impact:** Redundant rehash operations under contention (O(n) per redundant rehash).
- **Suggested fix:** Use a compare-and-swap or check-after-exclusive-acquire pattern:
  ```cpp
  void rehash_if_needed() {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    if (rehash_is_needed()) {
      std::size_t count = super::bucket_count();
      guard.unlock_shared();
      mcrl2::utilities::lock_guard ex = detail::g_thread_term_pool().lock();
      // Re-check under exclusive lock
      if (super::load_factor() >= super::max_load_factor()) {
        super::rehash(super::bucket_count() * 2);
      }
    }
  }
  ```

---

### [SEV-2] `size()` has no lock while `empty()` and `max_size()` do — inconsistent protection boundary
- **Location:** `unordered_map.h:204–207` (inline `size()`) vs `unordered_map_implementation.h` (`empty()`, `max_size()` with locks)
- **Why this is a problem:** `size()` is defined inline in the header as `return super::size();` with NO lock. `empty()` (which is semantically `size() == 0`) IS wrapped with `lock_shared()`. Either `size()` needs the lock (GC could run during the call, though the risk for a simple read is minimal) or `empty()`/`max_size()` don't need it. The inconsistency indicates confusion about the protection invariant and makes reasoning about safety impossible for maintainers.
- **Evidence type:** static-proof
- **Evidence:**
  - `unordered_map.h:204`: `std::size_t size() const { return super::size(); }` — no lock.
  - `unordered_map_implementation.h`: `bool empty() const noexcept { lock_shared(); return super::empty(); }` — locked.
  - These are semantically equivalent operations with different protection.
- **Efficiency impact:** Negligible.
- **Suggested fix:** Make `size()` consistent with `empty()` — either both lock or neither locks. Given the overall design philosophy of "always lock during container access," add the lock to `size()`.

---

### [SEV-2] `swap()` only protects `*this`, not `other` — potential GC inconsistency
- **Location:** `unordered_map_implementation.h`, both `atermpp::unordered_map::swap` and `atermpp::utilities::unordered_map::swap`
- **Why this is a problem:** `swap(other)` acquires `lock_shared()` on `*this` thread's term pool, then calls `super::swap(other)`. If `other` is owned by a different thread (rare but possible), `other`'s terms are being moved into `*this` without any lock held on `other`'s thread. More practically, `swap` modifies BOTH containers but only signals GC protection for one thread. If GC triggers on `other`'s owning thread mid-swap, the container contents are in an inconsistent state.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Thread A calls `map1.swap(map2)` while Thread B is iterating `map2` (both hold shared locks on their own thread pools). The swap modifies `map2` without Thread B's knowledge.
  - Test plan: Two threads sharing map references, one swapping while the other iterates, under TSan.
  - Expected failure: TSan report on concurrent modification of the swapped container.
- **Efficiency impact:** None (correctness issue).
- **Suggested fix:** Document that `swap` requires exclusive external synchronization on both maps, or acquire locks on both thread pools.

---

### [SEV-2] `find()` in `utilities` variant conditionally locks but all other read-only operations unconditionally lock
- **Location:** `unordered_map_implementation.h`, `find()` vs `begin()`/`end()`/`empty()`
- **Why this is a problem:** When `ThreadSafe=false` in a multithreaded build (`MCRL2_ENABLE_MULTITHREADING=ON`), `find()` skips the shared lock entirely. But `begin()`, `end()`, `empty()`, `max_size()` still acquire it unconditionally. Either `find()` is unsafe without the lock (GC can run during find, corrupting term access) or the other read-only methods are paying unnecessary lock overhead. Given that `find()` dereferences aterms for hash/equality comparison, it has the same GC-interaction requirements as iteration.
- **Evidence type:** static-proof
- **Evidence:**
  - `find()`: `if constexpr (ThreadSafe) { lock; find; } return super::find(args...);` — no lock when `ThreadSafe=false`
  - `begin()`: always locks regardless of `ThreadSafe`
  - In a multithreaded build with `ThreadSafe=false` instantiation, find() exposes the thread to GC interference.
- **C++20-specific note:** The `if constexpr` pattern is correct C++17/20 usage, but the semantic inconsistency between methods is the issue.
- **Suggested fix:** Either gate ALL methods on `ThreadSafe` (consistent performance optimization) or gate NONE (uniform safety). The current mix is the worst of both worlds.

---

### [SEV-3] No Doxygen `\pre`/`\post` on any public method
- **Location:** Entire file — only one `\brief` comment (line 32)
- **Why this is a problem:** Per mCRL2 coding standards, all public API methods require Doxygen documentation with `\pre`/`\post` clauses. The thread-safety requirements (caller must ensure no concurrent access to the same map instance from multiple threads) are completely undocumented. The relationship between the GC lock and the container's thread-safety model is nowhere explained.
- **Evidence type:** static-proof
- **Suggested fix:** Add at minimum:
  - `\pre` The calling thread must not already hold an exclusive lock on its term pool.
  - `\post` The GC has been prevented from running during the operation.
  - Class-level documentation explaining the locking model and that external synchronization is required for multi-threaded access to the same map.

---

### [SEV-3] Trailing whitespace and tab characters in several lines
- **Location:** Multiple lines (e.g., line with `erase(const_iterator pos)` has a tab character after the closing brace)
- **Why this is a problem:** Does not comply with `.clang-format` clean requirement for touched code.
- **Evidence type:** static-proof
- **Suggested fix:** Run `clang-format` on the file.

---

## Proposed Tests
- [ ] `test_unordered_map_operator_bracket_gc_safety`: Multi-threaded test calling `operator[]` on `atermpp::unordered_map` while forcing GC. Expected: TSan clean (or reproduces a race). Command: build with `-DMCRL2_ENABLE_THREADSANITIZER=ON`, run test.
- [ ] `test_rehash_if_needed_toctou`: N threads inserting into `atermpp::utilities::unordered_map<aterm,aterm,...,true>` at load-factor boundary simultaneously. Verify no double-rehash or TSan report.
- [ ] `test_swap_cross_thread`: Two threads, each owning a map, swap maps between threads under TSan. Verify no data race.
- [ ] `test_find_no_lock_gc_interference`: Instantiate `atermpp::utilities::unordered_map<..., false>` in a multithreaded build, call `find()` while another thread triggers GC. Verify no crash or TSan report.

## Open Questions
1. Is the omission of `operator[]`/`at()`/`find()` overrides in `atermpp::unordered_map` intentional? If so, what invariant makes them safe without the lock?
2. Is `size()` deliberately unlocked for performance while `empty()` is locked? What is the reasoning?
3. For `atermpp::utilities::unordered_map` with `ThreadSafe=false` in a multithreaded build — is this combination expected to be used, and if so, what external guarantees does the caller provide?
4. Should `rehash_if_needed()` re-verify the condition after acquiring exclusive access to avoid redundant work under contention?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         7m 53s
Total session time:     8m 21s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         421.2k in, 15.0k out, 367.2k cached (Est. 1 Premium request)
 claude-haiku-4.5        422.3k in, 13.8k out, 363.4k cached (Est. 0 Premium requests)

```

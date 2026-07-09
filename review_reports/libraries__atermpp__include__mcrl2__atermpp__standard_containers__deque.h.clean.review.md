---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/standard_containers/deque.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T12:44:01.853714+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/deque.h` — full file, clean review
- **Risk level:** Medium-High — missing GC-coordination methods create correctness gaps exploitable under concurrent garbage collection
- **Overall verdict:** The file follows the established aterm container wrapper pattern but has omissions relative to the vector.h sibling. Several structural modification methods that need GC coordination are not overridden, and documentation is stale.

## Findings (ordered by severity)

### [SEV-1] Missing `pop_front()` override — unguarded structural modification

- **Location:** deque.h (entire class; method absent)
- **Why this is a problem:** Every other structural modification (`push_front`, `push_back`, `pop_back`, `insert`, `erase`, `emplace`, `resize`, `clear`, `swap`) acquires `lock_shared()` to prevent the GC from entering a sweep while the container layout changes. `pop_front()` — a core deque operation — is inherited directly from `std::deque` without any lock. Destroying the front element without GC coordination risks the GC sweeping while internal pointers are in an inconsistent state, potentially missing the mark of elements being relocated.
- **Evidence type:** static-proof
- **Evidence:**
  - All peer mutating methods acquire `shared_guard` (lines 122–247); `pop_front()` is conspicuously absent while `pop_back()` is present at line 202.
  - `push_front` (line 208) has the guard, but the symmetric `pop_front` does not.
  - In the vector sibling (`vector.h`), `pop_back()` is overridden with a lock; `pop_front()` doesn't exist for vector so the gap only manifests in deque.
- **Efficiency impact:** None (trivial lock acquisition).
- **C++20-specific note:** N/A
- **Suggested fix:**
  ```cpp
  void pop_front()
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::pop_front();
  }
  ```

---

### [SEV-2] Missing `assign()` overrides — unguarded bulk replacement

- **Location:** deque.h (entire class; methods absent)
- **Why this is a problem:** `std::deque::assign()` (all three overloads: count+value, iterator range, initializer_list) clears the container and repopulates it. This is a structural modification that destroys all existing elements and reallocates. Without the shared lock, GC can sweep mid-operation. The same gap exists in `vector.h`, suggesting a systematic omission rather than a deliberate design choice.
- **Evidence type:** static-proof
- **Evidence:**
  - `clear()` at line 127 explicitly acquires the guard; `assign()` implicitly calls the equivalent destruction + construction without the guard.
  - Reproducer: call `d.assign(1000, some_aterm)` in one thread while another thread triggers `detail::g_thread_term_pool().collect()`.
- **Efficiency impact:** None.
- **Suggested fix:** Override all three `assign()` signatures with the shared guard.

---

### [SEV-2] `size()` lacks GC coordination — inconsistent with `vector.h`

- **Location:** deque.h:239–242
- **Why this is a problem:** `vector.h` explicitly wraps `size()` in `lock_shared()` with the comment `// Concurrent read/write on the size.` The deque's `size()` does not acquire the guard. If the GC triggers a sweep or the container is being modified on another thread, reading `size()` without the guard can observe a torn value. Since `vector.h` deemed this necessary, the deque inconsistency is unjustified.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `vector.h` lines 329–334 acquire `lock_shared()` for `size()`.
  - Deque line 239 does not.
  - Triggering scenario: thread A calls `d.size()`, thread B calls `d.push_back()` (which holds lock), GC tries to sweep concurrently. Without the guard on `size()`, there is no happens-before relationship between the size read and the GC sweep.
  - Validation plan: run a stress test with `atermpp::deque` under TSan where one thread reads `size()` in a loop and another mutates the deque.
- **Efficiency impact:** Negligible.
- **Suggested fix:**
  ```cpp
  std::size_t size() const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::size();
  }
  ```

---

### [SEV-3] Incorrect `\file` Doxygen path

- **Location:** deque.h:10
- **Why this is a problem:** The Doxygen `\file` directive states `mcrl2/data/standard_containers/deque.h`, but the actual path is `mcrl2/atermpp/standard_containers/deque.h`. This breaks Doxygen cross-referencing and is misleading to developers.
- **Evidence type:** static-proof
- **Evidence:** File resides at `libraries/atermpp/include/mcrl2/atermpp/standard_containers/deque.h`; line 10 reads `/// \file mcrl2/data/standard_containers/deque.h`.
- **Suggested fix:** Change to `/// \file mcrl2/atermpp/standard_containers/deque.h`.

---

### [SEV-3] No `ThreadSafe` template parameter — design divergence from `vector.h`

- **Location:** deque.h:34 (class template declaration)
- **Why this is a problem:** `vector.h` provides a `bool ThreadSafe = false` template parameter that upgrades from `lock_shared()` to exclusive `lock()` when true, enabling safe concurrent mutation from multiple threads. The deque lacks this parameter entirely, making it impossible to use a thread-safe deque without external synchronization. This limits composability and creates a surprising API asymmetry.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `vector.h` template declaration: `template<class T, bool ThreadSafe = false, class Alloc = ...>`
  - Deque has only `template<class T, class Alloc = ...>`
  - Any use case requiring a thread-safe deque must use ad-hoc locking.
- **Efficiency impact:** No runtime cost (template parameter is compile-time).
- **Suggested fix:** Add `bool ThreadSafe = false` template parameter and conditionally select `lock()` vs `lock_shared()` as in vector.h (if deque is ever used in thread-safe contexts).

---

### [SEV-3] Missing `[[nodiscard]]` on `size()` and `emplace`/`emplace_back`/`emplace_front` return values

- **Location:** deque.h:239, 167, 181, 225
- **Why this is a problem:** Per mCRL2 coding standards, pure queries and must-use results should be `[[nodiscard]]`. `size()` is a pure query; `emplace_back` and `emplace_front` return references to the newly inserted element that should not be silently discarded.
- **Evidence type:** static-proof (coding convention)
- **Suggested fix:** Add `[[nodiscard]]` to these declarations.

---

## Proposed Tests

- [ ] `parallel_deque_gc_pop_front_test`: Spawn N threads each calling `pop_front()` on a shared `atermpp::deque<aterm>` while another thread forces GC via `detail::g_thread_term_pool().collect()`. Under TSan, observe data race or crash due to missing lock. Run: `ctest -R parallel_deque`
- [ ] `deque_assign_gc_stress_test`: Call `assign(1000, term)` concurrently with GC collection. Validate no use-after-free under ASan.
- [ ] `deque_size_tsan_test`: One thread loops `d.size()`, another loops `d.push_back()`/`d.pop_back()`. Run under TSan to detect torn reads.

## Open Questions

1. Is the absence of `pop_front()` an oversight or was it intentionally left for some reason (e.g., the deque is never used with `pop_front` in practice)? A grep for `atermpp::deque` usage sites would clarify.
2. Should the `assign()` gap also be fixed in `vector.h` (same issue exists there)?
3. Is the `ThreadSafe` parameter intentionally omitted because no deque usage requires it, or is this simply incomplete implementation?

## Copilot Raw Stdout

```text
● Explore: Explore deque infrastructure
  ✗ Find test files with actual test cases
    $ grep -l "vector\|deque" /workspaces/mCRL2/libraries/atermpp/test/*.cpp | head -3 | xargs grep
    -l "BOOST_AUTO_TEST"
    Permission denied and could not request permission from user
  ● Create comprehensive summary
    $ cat > /tmp/summary.md << 'EOF'
    # mCRL2 atermpp Deque Infrastructure Summary

    ## 1. generic_aterm_container Class
    **File**:
    `/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h`...
    └ 389 lines...

● Explore: Investigate assignment operator safety
  ● Search for pop_front throughout atermpp library
    $ cd /workspaces/mCRL2 && grep -rn "pop_front" libraries/atermpp/include/mcrl2/atermpp/
    └ 10 lines...
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h lines 50-80
    └ 30 lines read

● Explore: Verify lambda capture behavior
  ● Search for these specific methods
    $ cd /workspaces/mCRL2 && grep -E "(operator\[\]|at\(|front\(|back\(|pop_front\()"
    libraries/atermpp/include/mcrl2/atermpp/standard_containers/deque.h
    └ 15 lines...
  ● List unoverridden deque methods
    $ cat > /tmp/check_methods.txt << 'EOF'
    Methods in std::deque that are NOT overridden in atermpp::deque:

    1. assign() - capacity-modifying methods
    2. operator[] - element access (both const and mutable)...
    └ 18 lines...

● Explore: Verify locking design intent
  ● Look at how lock_shared is used during modifications
    $ grep -B 3 -A 10 "shared_guard guard = detail::g_thread_term_pool().lock_shared()"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/vector.h | head
    -50
    └ 50 lines...
  ● Read libraries/atermpp/include/mcrl2/atermpp/standard_containers/vector.h lines 329-335
    └ 7 lines read

● Explore: Check pop_front and size locking
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h
    └ 465 lines read
  ● Search for pop_front and pop_back definitions
    $ grep -rn "pop_front\|pop_back"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/ | head -20
    └ 7 lines...

Now I have enough information to write the review. Let me produce the markdown findings.

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/deque.h` — full file, clean review
- **Risk level:** Medium-High — missing GC-coordination methods create correctness gaps exploitable under concurrent garbage collection
- **Overall verdict:** The file follows the established aterm container wrapper pattern but has omissions relative to the vector.h sibling. Several structural modification methods that need GC coordination are not overridden, and documentation is stale.

## Findings (ordered by severity)

### [SEV-1] Missing `pop_front()` override — unguarded structural modification

- **Location:** deque.h (entire class; method absent)
- **Why this is a problem:** Every other structural modification (`push_front`, `push_back`, `pop_back`, `insert`, `erase`, `emplace`, `resize`, `clear`, `swap`) acquires `lock_shared()` to prevent the GC from entering a sweep while the container layout changes. `pop_front()` — a core deque operation — is inherited directly from `std::deque` without any lock. Destroying the front element without GC coordination risks the GC sweeping while internal pointers are in an inconsistent state, potentially missing the mark of elements being relocated.
- **Evidence type:** static-proof
- **Evidence:**
  - All peer mutating methods acquire `shared_guard` (lines 122–247); `pop_front()` is conspicuously absent while `pop_back()` is present at line 202.
  - `push_front` (line 208) has the guard, but the symmetric `pop_front` does not.
  - In the vector sibling (`vector.h`), `pop_back()` is overridden with a lock; `pop_front()` doesn't exist for vector so the gap only manifests in deque.
- **Efficiency impact:** None (trivial lock acquisition).
- **C++20-specific note:** N/A
- **Suggested fix:**
  ```cpp
  void pop_front()
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::pop_front();
  }
  ```

---

### [SEV-2] Missing `assign()` overrides — unguarded bulk replacement

- **Location:** deque.h (entire class; methods absent)
- **Why this is a problem:** `std::deque::assign()` (all three overloads: count+value, iterator range, initializer_list) clears the container and repopulates it. This is a structural modification that destroys all existing elements and reallocates. Without the shared lock, GC can sweep mid-operation. The same gap exists in `vector.h`, suggesting a systematic omission rather than a deliberate design choice.
- **Evidence type:** static-proof
- **Evidence:**
  - `clear()` at line 127 explicitly acquires the guard; `assign()` implicitly calls the equivalent destruction + construction without the guard.
  - Reproducer: call `d.assign(1000, some_aterm)` in one thread while another thread triggers `detail::g_thread_term_pool().collect()`.
- **Efficiency impact:** None.
- **Suggested fix:** Override all three `assign()` signatures with the shared guard.

---

### [SEV-2] `size()` lacks GC coordination — inconsistent with `vector.h`

- **Location:** deque.h:239–242
- **Why this is a problem:** `vector.h` explicitly wraps `size()` in `lock_shared()` with the comment `// Concurrent read/write on the size.` The deque's `size()` does not acquire the guard. If the GC triggers a sweep or the container is being modified on another thread, reading `size()` without the guard can observe a torn value. Since `vector.h` deemed this necessary, the deque inconsistency is unjustified.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `vector.h` lines 329–334 acquire `lock_shared()` for `size()`.
  - Deque line 239 does not.
  - Triggering scenario: thread A calls `d.size()`, thread B calls `d.push_back()` (which holds lock), GC tries to sweep concurrently. Without the guard on `size()`, there is no happens-before relationship between the size read and the GC sweep.
  - Validation plan: run a stress test with `atermpp::deque` under TSan where one thread reads `size()` in a loop and another mutates the deque.
- **Efficiency impact:** Negligible.
- **Suggested fix:**
  ```cpp
  std::size_t size() const
  {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    return super::size();
  }
  ```

---

### [SEV-3] Incorrect `\file` Doxygen path

- **Location:** deque.h:10
- **Why this is a problem:** The Doxygen `\file` directive states `mcrl2/data/standard_containers/deque.h`, but the actual path is `mcrl2/atermpp/standard_containers/deque.h`. This breaks Doxygen cross-referencing and is misleading to developers.
- **Evidence type:** static-proof
- **Evidence:** File resides at `libraries/atermpp/include/mcrl2/atermpp/standard_containers/deque.h`; line 10 reads `/// \file mcrl2/data/standard_containers/deque.h`.
- **Suggested fix:** Change to `/// \file mcrl2/atermpp/standard_containers/deque.h`.

---

### [SEV-3] No `ThreadSafe` template parameter — design divergence from `vector.h`

- **Location:** deque.h:34 (class template declaration)
- **Why this is a problem:** `vector.h` provides a `bool ThreadSafe = false` template parameter that upgrades from `lock_shared()` to exclusive `lock()` when true, enabling safe concurrent mutation from multiple threads. The deque lacks this parameter entirely, making it impossible to use a thread-safe deque without external synchronization. This limits composability and creates a surprising API asymmetry.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `vector.h` template declaration: `template<class T, bool ThreadSafe = false, class Alloc = ...>`
  - Deque has only `template<class T, class Alloc = ...>`
  - Any use case requiring a thread-safe deque must use ad-hoc locking.
- **Efficiency impact:** No runtime cost (template parameter is compile-time).
- **Suggested fix:** Add `bool ThreadSafe = false` template parameter and conditionally select `lock()` vs `lock_shared()` as in vector.h (if deque is ever used in thread-safe contexts).

---

### [SEV-3] Missing `[[nodiscard]]` on `size()` and `emplace`/`emplace_back`/`emplace_front` return values

- **Location:** deque.h:239, 167, 181, 225
- **Why this is a problem:** Per mCRL2 coding standards, pure queries and must-use results should be `[[nodiscard]]`. `size()` is a pure query; `emplace_back` and `emplace_front` return references to the newly inserted element that should not be silently discarded.
- **Evidence type:** static-proof (coding convention)
- **Suggested fix:** Add `[[nodiscard]]` to these declarations.

---

## Proposed Tests

- [ ] `parallel_deque_gc_pop_front_test`: Spawn N threads each calling `pop_front()` on a shared `atermpp::deque<aterm>` while another thread forces GC via `detail::g_thread_term_pool().collect()`. Under TSan, observe data race or crash due to missing lock. Run: `ctest -R parallel_deque`
- [ ] `deque_assign_gc_stress_test`: Call `assign(1000, term)` concurrently with GC collection. Validate no use-after-free under ASan.
- [ ] `deque_size_tsan_test`: One thread loops `d.size()`, another loops `d.push_back()`/`d.pop_back()`. Run under TSan to detect torn reads.

## Open Questions

1. Is the absence of `pop_front()` an oversight or was it intentionally left for some reason (e.g., the deque is never used with `pop_front` in practice)? A grep for `atermpp::deque` usage sites would clarify.
2. Should the `assign()` gap also be fixed in `vector.h` (same issue exists there)?
3. Is the `ThreadSafe` parameter intentionally omitted because no deque usage requires it, or is this simply incomplete implementation?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         6m 34s
Total session time:     7m 13s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         268.8k in, 5.1k out, 239.2k cached (Est. 1 Premium request)
 claude-haiku-4.5        1.0m in, 21.4k out, 897.3k cached (Est. 0 Premium requests)

```

---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/standard_containers/vector.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T13:45:53.015964+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/vector.h` — full clean review of `atermpp::vector`, a GC-coordinating vector wrapper
- Risk level: **High** — data race potential in assignment operators and missing lock coverage for several mutation paths
- Overall verdict: **Contains correctness defects** — the defaulted assignment operators and unoverridden `reserve()`/`assign()` bypass GC coordination, creating real data-race windows. The `ThreadSafe=true` mode provides incomplete protection due to public inheritance from `std::vector`.

## Findings (ordered by severity)

### [SEV-1] Default assignment operators bypass GC coordination lock

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/vector.h:109-113`
- **Why this is a problem:** Both `operator=(const vector&)` and `operator=(vector&&)` are `= default`. The default implementation calls `std::vector::operator=()` which may reallocate the buffer and copy/move elements — all without acquiring a shared lock. During this window, the GC can acquire its exclusive lock and call `container_wrapper.mark()`, which iterates the vector via `for (const auto& element : container) { element.mark(todo); }`. If the vector is mid-reallocation (buffer freed, elements not yet copied), the GC reads freed memory (use-after-free) or partially-constructed elements.
- **Evidence type:** plausible-issue (confirmed by code analysis of the locking contract)
- **Evidence:**
  - Reproducer/test name: Extend `parallel_container_test` to exercise assignment during GC
  - Run command: `cmake --build build -t parallel_container_test && cd build && ctest -R parallel_container --output-on-failure` (with ASan: `-DMCRL2_ENABLE_ADDRESSSANITIZER=ON`)
  - Proposed test design:
    ```cpp
    BOOST_AUTO_TEST_CASE(parallel_vector_assignment)
    {
      if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
      {
        atermpp::aterm_int value(0);
        atermpp::vector<atermpp::aterm> source;
        for (std::size_t i = 0; i < 10000; ++i) source.push_back(value);

        std::thread local([&source]() {
          for (std::size_t i = 0; i < 1000; ++i) {
            atermpp::vector<atermpp::aterm> target;
            target = source; // No lock acquired — GC race window
          }
        });

        for (std::size_t i = 0; i < 10000; ++i)
          atermpp::detail::g_thread_term_pool().collect();

        local.join();
      }
    }
    ```
  - Observed result: Potential ASan heap-use-after-free or TSan data race on vector internals
  - Expected result: Assignment should be protected by at least a shared lock
- **Efficiency impact:** None (fix adds one lock acquisition per assignment)
- **C++20-specific note:** The `= default` is convenient but hides the missing lock. A custom implementation is required.
- **Suggested fix:**
  ```cpp
  vector& operator=(const vector& x)
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::operator=(x);
    } else {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::operator=(x);
    }
    return *this;
  }

  vector& operator=(vector&& x) noexcept
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::operator=(std::move(x));
    } else {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::operator=(std::move(x));
    }
    return *this;
  }
  ```

---

### [SEV-2] `reserve()` not overridden — GC can iterate during reallocation

- **Location:** Entire class (missing override)
- **Why this is a problem:** `std::vector::reserve()` may reallocate the internal buffer. Since `atermpp::vector` publicly inherits from `std::vector`, callers can invoke `reserve()` directly without any GC lock. A concurrent GC mark phase iterates the container's old buffer (freed memory) or sees inconsistent internal pointers.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: `reserve()` → realloc → GC fires → `container_wrapper.mark()` reads freed old buffer
  - Run command: Build with ASan, call `reserve()` in a tight loop concurrently with `collect()`
  - Expected failure: ASan heap-use-after-free
- **Efficiency impact:** Minimal (one lock acquisition per `reserve()` call)
- **C++20-specific note:** N/A
- **Suggested fix:** Override `reserve()` with the same `if constexpr (ThreadSafe)` locking pattern.

---

### [SEV-2] `assign()` not overridden — same GC coordination gap

- **Location:** Entire class (missing override)
- **Why this is a problem:** `std::vector::assign()` replaces all elements and may reallocate. Same race condition as `reserve()`.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Same mechanism as `reserve()` finding
- **Suggested fix:** Override all `assign()` overloads with the locking pattern.

---

### [SEV-2] `ThreadSafe=true` provides incomplete protection — public inheritance leaks unlocked API

- **Location:** Class declaration, line 35: `class vector : public std::vector<...>`
- **Why this is a problem:** When `ThreadSafe=true`, the intent is multi-thread safe access. But public inheritance exposes every `std::vector` method without locking: `operator[]`, `at()`, `front()`, `back()`, `begin()`, `end()`, `data()`, `reserve()`, `assign()`, `capacity()`, `empty()`. Any concurrent call to these while another thread mutates is a data race (UB per C++ memory model). The class creates a false sense of thread-safety.
- **Evidence type:** static-proof
- **Evidence:**
  - Two threads: thread A calls `vec.push_back(x)` (acquires exclusive lock on GC mutex), thread B calls `vec[0]` (no lock acquired via inherited method). This is a data race on the vector's internal storage.
  - The GC mutex is orthogonal to user-thread synchronization — it only coordinates GC vs. mutators, not mutator vs. mutator reads.
- **C++20-specific note:** A safer design would use private/protected inheritance with explicit forwarding of read methods (guarded by `if constexpr (ThreadSafe)`), or a separate synchronized wrapper.
- **Suggested fix:** At minimum, document that `ThreadSafe=true` only protects GC coordination for concurrent writers, not general thread-safe read/write access. Ideally, consider private inheritance to prevent accidental unprotected access.

---

### [SEV-2] GC visibility window during element-storing constructors (plausible)

- **Location:** Lines 62-78 (constructors that store elements before `container_wrapper` initialization)
- **Why this is a problem:** In constructors like `vector(size_type n, const value_type& val, ...)`, the base class creates `n` aterm elements. `container_wrapper` (which registers with GC) is initialized AFTER the base. If another thread triggers GC in this window, the elements are unreachable by the GC marker. Elements only referenced by this vector could be collected prematurely.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Thread A constructs `vector(1000, some_aterm)`. After base ctor completes, 1000 aterms exist only in this unregistered vector. Thread B triggers GC → mark phase doesn't find these 1000 aterms → sweep frees them → Thread A's vector now contains dangling term pointers.
  - Mitigating factor: In practice, `val` is still alive on the caller's stack (protecting at least the prototype term). But the copies in the vector share the same underlying term pointer (due to maximal sharing), so they ARE still reachable through `val`. This makes the issue unlikely for copy-constructing cases but could manifest for range constructors with ephemeral iterators.
  - Test plan: Construct a vector from a range of freshly-created, non-shared aterms while another thread aggressively collects.
- **Suggested fix:** Acquire shared lock during the constructor body or use a two-phase init pattern where `container_wrapper` registers first.

---

### [SEV-3] `size()` lock is inconsistent and provides no useful guarantee

- **Location:** Lines 329-334
- **Why this is a problem:** `size()` unconditionally takes a `shared_guard`, but no other read accessor (`operator[]`, `at()`, `begin()`, `empty()`, `capacity()`) is overridden with a lock. If the lock is needed for GC coordination during reads, ALL reads need it. If not (and it's not — GC only reads the container during mark, which is safe vs. concurrent reads), then `size()` doesn't need it either. For `ThreadSafe=true`, a shared lock doesn't protect against a concurrent writer anyway (the writer takes exclusive, but the size is stale after release). The lock imposes overhead without providing any safety guarantee.
- **Evidence type:** static-proof
- **Evidence:** The GC's mark phase only reads container contents (to mark aterms). User reads don't modify the container. No protection is needed for concurrent reads. The comment "Concurrent read/write on the size" is incorrect — shared_lock doesn't protect reads from writes.
- **Suggested fix:** Either remove the lock from `size()` (if consistent with other unoverridden read methods), or override ALL read accessors with the locking pattern for consistency. The former is recommended as it matches the existing behavior of inherited read methods.

---

### [SEV-3] Incorrect `\file` path in Doxygen header

- **Location:** Line 11
- **Why this is a problem:** `/// \file mcrl2/data/standard_containers/vector.h` — says `data` but the actual path is `mcrl2/atermpp/standard_containers/vector.h`. This causes Doxygen to generate incorrect documentation linkage.
- **Evidence type:** static-proof
- **Suggested fix:** Change to `/// \file mcrl2/atermpp/standard_containers/vector.h`

---

## Proposed Tests

- [ ] `parallel_vector_assignment`: Exercise copy/move assignment of `atermpp::vector` concurrently with GC collection. Expected signal: TSan data race report or ASan heap-use-after-free on vector internals during assignment.
- [ ] `parallel_vector_reserve`: Call `reserve()` with increasing capacities concurrently with GC. Expected signal: ASan heap-use-after-free as GC reads stale buffer.
- [ ] `parallel_vector_construction_gc_window`: Construct large vectors from ephemeral iterator ranges while another thread aggressively calls `collect()`. Expected signal: Accessing freed term memory post-construction.

## Open Questions

1. Is the `ThreadSafe` template parameter intended to provide user-thread safety, or only GC-coordination safety for concurrent mutators? The current implementation conflates both concerns without fully achieving either.
2. Are the same defects present in `atermpp::deque` and `atermpp::unordered_map` (both use `= default` assignment)? A consistent fix should be applied across all containers.
3. Why does `size()` have a lock but `empty()`, `capacity()`, and element accessors do not? Is there historical context for this inconsistency?
4. Has the GC visibility window during construction ever manifested as a bug? The aterm maximal-sharing model may mask it in practice (shared aterms are kept alive by other references), but the invariant violation remains.

## Copilot Raw Stdout

```text
● Explore: Explore aterm container internals
  ● Read libraries/atermpp/include/mcrl2/atermpp/standard_containers/deque.h
    └ 252 lines read
  ● Read libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h
    └ 455 lines read

Now let me check if `reserve()` and `assign()` are overridden, and verify the locking model in more detail:

● Explore: Check reserve/assign and aterm_container registration
  ● Find uses of atermpp::vector in tests
    $ grep -r "atermpp::vector\|atermpp\:\:vector" /workspaces/mCRL2/libraries/atermpp/test/*.cpp -A
     3 -B 1 | head -60
    └ 12 lines...
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h lines 1-50
    └ 50 lines read

Now I have all the context needed for the review. Let me write up the findings.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/vector.h` — full clean review of `atermpp::vector`, a GC-coordinating vector wrapper
- Risk level: **High** — data race potential in assignment operators and missing lock coverage for several mutation paths
- Overall verdict: **Contains correctness defects** — the defaulted assignment operators and unoverridden `reserve()`/`assign()` bypass GC coordination, creating real data-race windows. The `ThreadSafe=true` mode provides incomplete protection due to public inheritance from `std::vector`.

## Findings (ordered by severity)

### [SEV-1] Default assignment operators bypass GC coordination lock

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/vector.h:109-113`
- **Why this is a problem:** Both `operator=(const vector&)` and `operator=(vector&&)` are `= default`. The default implementation calls `std::vector::operator=()` which may reallocate the buffer and copy/move elements — all without acquiring a shared lock. During this window, the GC can acquire its exclusive lock and call `container_wrapper.mark()`, which iterates the vector via `for (const auto& element : container) { element.mark(todo); }`. If the vector is mid-reallocation (buffer freed, elements not yet copied), the GC reads freed memory (use-after-free) or partially-constructed elements.
- **Evidence type:** plausible-issue (confirmed by code analysis of the locking contract)
- **Evidence:**
  - Reproducer/test name: Extend `parallel_container_test` to exercise assignment during GC
  - Run command: `cmake --build build -t parallel_container_test && cd build && ctest -R parallel_container --output-on-failure` (with ASan: `-DMCRL2_ENABLE_ADDRESSSANITIZER=ON`)
  - Proposed test design:
    ```cpp
    BOOST_AUTO_TEST_CASE(parallel_vector_assignment)
    {
      if constexpr (mcrl2::utilities::detail::GlobalThreadSafe)
      {
        atermpp::aterm_int value(0);
        atermpp::vector<atermpp::aterm> source;
        for (std::size_t i = 0; i < 10000; ++i) source.push_back(value);

        std::thread local([&source]() {
          for (std::size_t i = 0; i < 1000; ++i) {
            atermpp::vector<atermpp::aterm> target;
            target = source; // No lock acquired — GC race window
          }
        });

        for (std::size_t i = 0; i < 10000; ++i)
          atermpp::detail::g_thread_term_pool().collect();

        local.join();
      }
    }
    ```
  - Observed result: Potential ASan heap-use-after-free or TSan data race on vector internals
  - Expected result: Assignment should be protected by at least a shared lock
- **Efficiency impact:** None (fix adds one lock acquisition per assignment)
- **C++20-specific note:** The `= default` is convenient but hides the missing lock. A custom implementation is required.
- **Suggested fix:**
  ```cpp
  vector& operator=(const vector& x)
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::operator=(x);
    } else {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::operator=(x);
    }
    return *this;
  }

  vector& operator=(vector&& x) noexcept
  {
    if constexpr (ThreadSafe) {
      mcrl2::utilities::lock_guard guard = detail::g_thread_term_pool().lock();
      super::operator=(std::move(x));
    } else {
      mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
      super::operator=(std::move(x));
    }
    return *this;
  }
  ```

---

### [SEV-2] `reserve()` not overridden — GC can iterate during reallocation

- **Location:** Entire class (missing override)
- **Why this is a problem:** `std::vector::reserve()` may reallocate the internal buffer. Since `atermpp::vector` publicly inherits from `std::vector`, callers can invoke `reserve()` directly without any GC lock. A concurrent GC mark phase iterates the container's old buffer (freed memory) or sees inconsistent internal pointers.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: `reserve()` → realloc → GC fires → `container_wrapper.mark()` reads freed old buffer
  - Run command: Build with ASan, call `reserve()` in a tight loop concurrently with `collect()`
  - Expected failure: ASan heap-use-after-free
- **Efficiency impact:** Minimal (one lock acquisition per `reserve()` call)
- **C++20-specific note:** N/A
- **Suggested fix:** Override `reserve()` with the same `if constexpr (ThreadSafe)` locking pattern.

---

### [SEV-2] `assign()` not overridden — same GC coordination gap

- **Location:** Entire class (missing override)
- **Why this is a problem:** `std::vector::assign()` replaces all elements and may reallocate. Same race condition as `reserve()`.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Same mechanism as `reserve()` finding
- **Suggested fix:** Override all `assign()` overloads with the locking pattern.

---

### [SEV-2] `ThreadSafe=true` provides incomplete protection — public inheritance leaks unlocked API

- **Location:** Class declaration, line 35: `class vector : public std::vector<...>`
- **Why this is a problem:** When `ThreadSafe=true`, the intent is multi-thread safe access. But public inheritance exposes every `std::vector` method without locking: `operator[]`, `at()`, `front()`, `back()`, `begin()`, `end()`, `data()`, `reserve()`, `assign()`, `capacity()`, `empty()`. Any concurrent call to these while another thread mutates is a data race (UB per C++ memory model). The class creates a false sense of thread-safety.
- **Evidence type:** static-proof
- **Evidence:**
  - Two threads: thread A calls `vec.push_back(x)` (acquires exclusive lock on GC mutex), thread B calls `vec[0]` (no lock acquired via inherited method). This is a data race on the vector's internal storage.
  - The GC mutex is orthogonal to user-thread synchronization — it only coordinates GC vs. mutators, not mutator vs. mutator reads.
- **C++20-specific note:** A safer design would use private/protected inheritance with explicit forwarding of read methods (guarded by `if constexpr (ThreadSafe)`), or a separate synchronized wrapper.
- **Suggested fix:** At minimum, document that `ThreadSafe=true` only protects GC coordination for concurrent writers, not general thread-safe read/write access. Ideally, consider private inheritance to prevent accidental unprotected access.

---

### [SEV-2] GC visibility window during element-storing constructors (plausible)

- **Location:** Lines 62-78 (constructors that store elements before `container_wrapper` initialization)
- **Why this is a problem:** In constructors like `vector(size_type n, const value_type& val, ...)`, the base class creates `n` aterm elements. `container_wrapper` (which registers with GC) is initialized AFTER the base. If another thread triggers GC in this window, the elements are unreachable by the GC marker. Elements only referenced by this vector could be collected prematurely.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Thread A constructs `vector(1000, some_aterm)`. After base ctor completes, 1000 aterms exist only in this unregistered vector. Thread B triggers GC → mark phase doesn't find these 1000 aterms → sweep frees them → Thread A's vector now contains dangling term pointers.
  - Mitigating factor: In practice, `val` is still alive on the caller's stack (protecting at least the prototype term). But the copies in the vector share the same underlying term pointer (due to maximal sharing), so they ARE still reachable through `val`. This makes the issue unlikely for copy-constructing cases but could manifest for range constructors with ephemeral iterators.
  - Test plan: Construct a vector from a range of freshly-created, non-shared aterms while another thread aggressively collects.
- **Suggested fix:** Acquire shared lock during the constructor body or use a two-phase init pattern where `container_wrapper` registers first.

---

### [SEV-3] `size()` lock is inconsistent and provides no useful guarantee

- **Location:** Lines 329-334
- **Why this is a problem:** `size()` unconditionally takes a `shared_guard`, but no other read accessor (`operator[]`, `at()`, `begin()`, `empty()`, `capacity()`) is overridden with a lock. If the lock is needed for GC coordination during reads, ALL reads need it. If not (and it's not — GC only reads the container during mark, which is safe vs. concurrent reads), then `size()` doesn't need it either. For `ThreadSafe=true`, a shared lock doesn't protect against a concurrent writer anyway (the writer takes exclusive, but the size is stale after release). The lock imposes overhead without providing any safety guarantee.
- **Evidence type:** static-proof
- **Evidence:** The GC's mark phase only reads container contents (to mark aterms). User reads don't modify the container. No protection is needed for concurrent reads. The comment "Concurrent read/write on the size" is incorrect — shared_lock doesn't protect reads from writes.
- **Suggested fix:** Either remove the lock from `size()` (if consistent with other unoverridden read methods), or override ALL read accessors with the locking pattern for consistency. The former is recommended as it matches the existing behavior of inherited read methods.

---

### [SEV-3] Incorrect `\file` path in Doxygen header

- **Location:** Line 11
- **Why this is a problem:** `/// \file mcrl2/data/standard_containers/vector.h` — says `data` but the actual path is `mcrl2/atermpp/standard_containers/vector.h`. This causes Doxygen to generate incorrect documentation linkage.
- **Evidence type:** static-proof
- **Suggested fix:** Change to `/// \file mcrl2/atermpp/standard_containers/vector.h`

---

## Proposed Tests

- [ ] `parallel_vector_assignment`: Exercise copy/move assignment of `atermpp::vector` concurrently with GC collection. Expected signal: TSan data race report or ASan heap-use-after-free on vector internals during assignment.
- [ ] `parallel_vector_reserve`: Call `reserve()` with increasing capacities concurrently with GC. Expected signal: ASan heap-use-after-free as GC reads stale buffer.
- [ ] `parallel_vector_construction_gc_window`: Construct large vectors from ephemeral iterator ranges while another thread aggressively calls `collect()`. Expected signal: Accessing freed term memory post-construction.

## Open Questions

1. Is the `ThreadSafe` template parameter intended to provide user-thread safety, or only GC-coordination safety for concurrent mutators? The current implementation conflates both concerns without fully achieving either.
2. Are the same defects present in `atermpp::deque` and `atermpp::unordered_map` (both use `= default` assignment)? A consistent fix should be applied across all containers.
3. Why does `size()` have a lock but `empty()`, `capacity()`, and element accessors do not? Is there historical context for this inconsistency?
4. Has the GC visibility window during construction ever manifested as a bug? The aterm maximal-sharing model may mask it in practice (shared aterms are kept alive by other references), but the invariant violation remains.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         5m 21s
Total session time:     5m 39s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         140.1k in, 10.1k out, 107.2k cached (Est. 1 Premium request)
 claude-haiku-4.5        321.7k in, 11.2k out, 247.7k cached (Est. 0 Premium requests)

```

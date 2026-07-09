---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T13:22:03.236215+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h` — both the `atermpp::unordered_map` (wrapping `std::unordered_map`) and `atermpp::utilities::unordered_map` (wrapping `mcrl2::utilities::unordered_map`) classes; plus implementation in `detail/unordered_map_implementation.h` and the `generic_aterm_container` / `markable_aterm` infrastructure in `detail/aterm_container.h`.
- **Risk level:** Medium-High
- **Overall verdict:** The GC-protection design (bulk marking via `container_wrapper`) is sound in principle — captured references to the base sub-object remain valid across assignment and swap. However, the public interface leaks inherited mutating operations that bypass the GC shared lock, and one constructor is structurally broken. Several missing API overrides create thread-safety holes.

## Findings (ordered by severity)

### [SEV-1] Defaulted copy/move assignment operators bypass GC lock — data race with GC marking

- **Location:** `unordered_map.h:118–121` (first class), `unordered_map.h:289–292` (utilities class)
- **Why this is a problem:** The defaulted `operator=` directly assigns the base `std::unordered_map` (or `mcrl2::utilities::unordered_map`) without acquiring the GC shared lock (`lock_shared()`). Every other mutating operation in the implementation file acquires `shared_guard guard = detail::g_thread_term_pool().lock_shared()` before modifying the container. During a defaulted assignment, the thread's "busy flag" is not set, so the GC is free to run concurrently. The GC mark lambda (`for (const auto& element : container) { element.mark(todo); }`) iterates the container via iterators — if the container's internal hash table is being reassigned (copy) or pointer-swapped (move) simultaneously, iterators are invalidated → **undefined behavior (data race, use-after-free)**.
- **Evidence type:** static-proof
- **Evidence:**
  - The protocol requires `lock_shared()` for any container structural modification to prevent concurrent GC marking. Defaulted `operator=` calls `std::unordered_map::operator=` which is a structural modification.
  - Compare: `clear()`, `insert()`, `erase()`, `swap()` all acquire `shared_guard` in `detail/unordered_map_implementation.h`.
  - Reproducer design: two threads — thread A performs `map_a = map_b` in a loop; thread B forces GC collection via `g_thread_term_pool().collect()`. Under TSan, this should report a data race on the hash table internals.
  - Run command: `cmake -B build-tsan -DMCRL2_ENABLE_THREADSANITIZER=ON && cmake --build build-tsan && ctest --test-dir build-tsan -R parallel_container`
  - Expected result: TSan reports data race between assignment and GC mark iteration.
- **Efficiency impact:** None — fix adds one lock acquisition per assignment.
- **C++20-specific note:** Consider deleting the defaulted operators and providing explicit ones. With C++20, you could use `= default` inside the function body after the lock guard, but that changes semantics. Custom implementations are clearer.
- **Suggested fix:** Replace `= default` with explicit implementations that acquire `shared_guard` before delegating to `super::operator=`:
  ```cpp
  unordered_map& operator=(const unordered_map& other) {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::operator=(other);
    return *this;
  }
  ```

---

### [SEV-2] Inherited `operator[]` and `at()` bypass GC lock — mutating operations exposed without synchronization

- **Location:** `unordered_map.h:40` (public inheritance from `std::unordered_map`)
- **Why this is a problem:** The first class (`atermpp::unordered_map`) publicly inherits from `std::unordered_map`. While mutating operations like `insert`, `emplace`, `erase`, `clear` are overridden to acquire the GC lock, `operator[]` and non-const `at()` are NOT overridden. `operator[]` is particularly dangerous: it may insert a default-constructed element and trigger a rehash — all without the shared lock. This creates the same race condition as SEV-1 (GC mark lambda iterating concurrently with a structural modification).
- **Evidence type:** plausible-issue
- **Evidence:**
  - `grep -n "operator\[" unordered_map.h` returns no matches — `operator[]` is not overridden.
  - `std::unordered_map::operator[]` inserts a default element when key is absent (per [unord.map.elem]).
  - Any user calling `my_atermpp_map[key] = value;` bypasses the GC lock.
  - Candidate test: thread A calls `map[aterm_key]` repeatedly with new keys; thread B forces GC. TSan should report race.
- **Efficiency impact:** None.
- **C++20-specific note:** Consider making the inheritance `private` or `protected` and exposing only the safe API explicitly, or at minimum override `operator[]` and `at()`.
- **Suggested fix:** Override `operator[]` with lock acquisition, or (better) switch to private inheritance with explicit forwarding of the safe subset.

---

### [SEV-2] Constructor `(size_type, value_type, allocator)` is structurally broken — will fail to compile when instantiated

- **Location:** `unordered_map.h:74–77` (first class), `unordered_map.h:277–280` (utilities class)
- **Why this is a problem:** This constructor:
  ```cpp
  unordered_map(size_type n, const value_type& val, const allocator_type& alloc = allocator_type())
   : super::unordered_map(n, detail::markable_aterm<std::pair<const Key, T>>(val), alloc),
  ```
  attempts to pass a `markable_aterm<std::pair<...>>` as the second argument to `std::unordered_map`'s constructor. `std::unordered_map` has no fill-constructor; its relevant constructor is `(size_type bucket_count, const hasher& hf, const allocator_type& a)`. The second argument would need to be a `Hash` object. `markable_aterm<std::pair<...>>` is not convertible to `std::hash<markable_aterm<Key>>`. This is dead code that will produce a hard compile error if ever instantiated. It appears to be a copy-paste from a `std::vector`-style API.
- **Evidence type:** static-proof
- **Evidence:**
  - `std::unordered_map` (C++11–C++23) has no constructor with signature `(size_type, value_type, allocator)`.
  - The 3-argument form is `(size_type n, const hasher& hf, const allocator_type& a)` per [unord.map.cnstr].
  - This code survives only because it's a template — SFINAE-by-non-instantiation.
  - Reproducer: `atermpp::unordered_map<aterm, aterm> m(10, std::make_pair(aterm(), aterm()));` → compile error.
- **Efficiency impact:** None (dead code).
- **C++20-specific note:** N/A.
- **Suggested fix:** Delete this constructor entirely. `std::unordered_map` has no fill constructor — this concept doesn't apply to associative containers.

---

### [SEV-3] Incorrect Doxygen `\file` path

- **Location:** `unordered_map.h:11`
- **Why this is a problem:** The `\file` directive says `mcrl2/data/standard_containers/unordered_map.h` but the actual path is `mcrl2/atermpp/standard_containers/unordered_map.h`. This causes Doxygen to create a phantom documentation page for a non-existent file.
- **Evidence type:** static-proof
- **Evidence:** Compare line 11 (`mcrl2/data/...`) with the actual filesystem path.
- **Suggested fix:** Change to `\file mcrl2/atermpp/standard_containers/unordered_map.h`.

---

### [SEV-3] Dead code in `markable_aterm` fundamental-type specialization (transitive dependency)

- **Location:** `detail/aterm_container.h:128–131` and `detail/aterm_container.h:137–138`
- **Why this is a problem:** Duplicated unreachable statements after `return`:
  ```cpp
  markable_aterm& operator=(const T& other) noexcept {
    m_t = other;
    return *this;
    m_t = other;    // ← dead code
    return *this;   // ← dead code
  }
  ```
  Same pattern in the move-assignment operator and in `operator std_pair&()` (line 245). While harmless at runtime, it suggests a faulty merge or accidental duplication that should be cleaned up.
- **Evidence type:** static-proof
- **Suggested fix:** Remove the duplicated lines after the first `return *this;`.

---

### [SEV-3] Missing `[[nodiscard]]` on query functions

- **Location:** `unordered_map.h:196` (`size()`), lines 219–233 (`begin`, `end`, `cbegin`, `cend`, `empty`, `max_size`)
- **Why this is a problem:** Per the project's C++ instructions, pure queries and must-use results should be marked `[[nodiscard]]`. Discarding the return value of `empty()` or `size()` is almost certainly a bug.
- **Evidence type:** plausible-issue (no runtime failure, but violates project convention for new/touched code).
- **Suggested fix:** Add `[[nodiscard]]` to `size()`, `empty()`, `max_size()`, and the iterator accessors.

---

### [SEV-3] Redundant `#pragma once` with traditional include guard

- **Location:** `unordered_map.h:20–21`
- **Why this is a problem:** Both `#ifndef MCRL2_ATERMPP_...` guard and `#pragma once` are present. Harmless but redundant. The project appears to use include guards as the standard mechanism.
- **Evidence type:** style observation.
- **Suggested fix:** Remove `#pragma once` for consistency, or adopt it project-wide.

---

## Proposed Tests

- [ ] **`test_assignment_race`**: Two threads — thread A performs copy/move assignment of `atermpp::unordered_map` in a loop (source with 1000 elements); thread B calls `g_thread_term_pool().collect()` in a loop. Build with TSan. Expected: TSan reports data race (currently) → after fix, no report.
- [ ] **`test_operator_bracket_race`**: Thread A calls `map[new_key]` (inserting) in a loop; thread B forces GC. Build with TSan. Expected: TSan reports race on hash table internals.
- [ ] **`test_bogus_fill_constructor`**: Compile-time test that explicitly instantiates the `(size_type, value_type, allocator)` constructor. Expected: compilation failure confirming it's broken dead code.

## Open Questions

1. Is the public inheritance from `std::unordered_map` intentional for performance (avoiding virtual dispatch) or a legacy design choice? Private inheritance with explicit forwarding would be safer and prevent the lock-bypass problem for any future std::unordered_map additions (e.g., C++23 `insert_range`).
2. Are there any callers relying on the `(size_type, value_type, allocator)` constructor? A codebase search found no instantiations — safe to remove?
3. The `utilities::unordered_map` variant conditionally acquires the lock only in `find()` when `ThreadSafe=true`, but always acquires it in `insert()`. Is the intention that non-ThreadSafe instances are never accessed concurrently with GC? If so, this relies on the single-threaded code never triggering GC, which may be fragile.

## Copilot Raw Stdout

```text
● Explore: Explore atermpp container context
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Explore: Check assignment semantics
  └ Agent started in background with agent_id: agent-1. You can use read_agent to...

● Explore: Check thread safety
  └ Agent started in background with agent_id: agent-2. You can use read_agent to...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h
  └ 465 lines read

● Glob "libraries/atermpp/source/*unordered_map*"
  └ No matches found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h lines 126-138
  └ 13 lines read

● Grep "unordered_map" (libraries/atermpp/source)
  └ No matches found

● Glob "libraries/atermpp/source/*container*"
  └ No matches found

● Grep "atermpp::unordered_map" in *.cpp (libraries/atermpp)
  └ No matches found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h lines 55-75
  └ 21 lines read

● Grep "void.*unordered_map.*clear|unordered_map.*::insert|unordered_map.*::emplace" in *.cpp
(libraries/atermpp)
  └ No matches found

● Grep "aterm_container::" in *.cpp (libraries/atermpp)
  └ No matches found

● Read libraries/utilities/include/mcrl2/utilities/noncopyable.h
  └ 34 lines read

● Grep "unordered_map.*::"" in *.h (libraries/atermpp/include)
  └ No matches found

● Grep "aterm_container"" in *.cpp (libraries/atermpp)
  └ No matches found

● Glob "libraries/atermpp/**/*unordered_map*"
  └ 2 files found

● Read: Explore agent (Explore atermpp container context)
  └ Completed

● Read: Explore agent (Check assignment semantics)
  └ Completed

● Read
  libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h
  └ 1 line read

● Read: Explore agent (Check thread safety)
  └ Running

● Read
  libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h
  lines 1-100
  └ 100 lines read

● Grep "void.*swap|noexcept.*swap" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standa
rd_containers/detail/unordered_map_implementation.h)
  └ 23 lines found

● Read
  libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h
  lines 250-330
  └ 81 lines read

● Grep "operator\[\]|at\(" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_conta
iners/detail/unordered_map_implementation.h)
  └ No matches found

● Read
  libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h
  lines 100-200
  └ 101 lines read

● Grep "class shared_guard|shared_guard\(|lock_shared"
(/workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/shared_mutex.h)
  └ 86 lines found

● Read: Explore agent (Check thread safety)
  └ Completed

● Grep "operator\[|find\(|at\(|count\(|contains\("
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h)
  └ 2 lines found

● Grep "atermpp::(unordered_map|vector|unordered_set)" in *.{h,cpp} (/workspaces/mCRL2/libraries)
  └ 10 files found

● Grep "find\(" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/detai
l/unordered_map_implementation.h)
  └ 6 lines found

● Read
  libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_map_implementation.h
  lines 350-380
  └ 31 lines read

● Grep "rehash_if_needed" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_contai
ners/detail/unordered_map_implementation.h)
  └ 5 lines found

Now I have all the information needed to produce the review. Let me verify one more detail about whether `operator[]` is accessible through the public inheritance.

● Check operator[] usage on atermpp maps
  $ cd /workspaces/mCRL2 && grep -rn "atermpp::unordered_map.*\[" --include="*.h" --include="*.cpp"
  libraries/ | head -10
  └ 1 line...

● Check the bogus constructor argument
  $ cd /workspaces/mCRL2 && grep -rn "detail::markable_aterm<std::pair" --include="*.h"
  libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h
  └ 3 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h lines 70-80
  └ 11 lines read

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h` — both the `atermpp::unordered_map` (wrapping `std::unordered_map`) and `atermpp::utilities::unordered_map` (wrapping `mcrl2::utilities::unordered_map`) classes; plus implementation in `detail/unordered_map_implementation.h` and the `generic_aterm_container` / `markable_aterm` infrastructure in `detail/aterm_container.h`.
- **Risk level:** Medium-High
- **Overall verdict:** The GC-protection design (bulk marking via `container_wrapper`) is sound in principle — captured references to the base sub-object remain valid across assignment and swap. However, the public interface leaks inherited mutating operations that bypass the GC shared lock, and one constructor is structurally broken. Several missing API overrides create thread-safety holes.

## Findings (ordered by severity)

### [SEV-1] Defaulted copy/move assignment operators bypass GC lock — data race with GC marking

- **Location:** `unordered_map.h:118–121` (first class), `unordered_map.h:289–292` (utilities class)
- **Why this is a problem:** The defaulted `operator=` directly assigns the base `std::unordered_map` (or `mcrl2::utilities::unordered_map`) without acquiring the GC shared lock (`lock_shared()`). Every other mutating operation in the implementation file acquires `shared_guard guard = detail::g_thread_term_pool().lock_shared()` before modifying the container. During a defaulted assignment, the thread's "busy flag" is not set, so the GC is free to run concurrently. The GC mark lambda (`for (const auto& element : container) { element.mark(todo); }`) iterates the container via iterators — if the container's internal hash table is being reassigned (copy) or pointer-swapped (move) simultaneously, iterators are invalidated → **undefined behavior (data race, use-after-free)**.
- **Evidence type:** static-proof
- **Evidence:**
  - The protocol requires `lock_shared()` for any container structural modification to prevent concurrent GC marking. Defaulted `operator=` calls `std::unordered_map::operator=` which is a structural modification.
  - Compare: `clear()`, `insert()`, `erase()`, `swap()` all acquire `shared_guard` in `detail/unordered_map_implementation.h`.
  - Reproducer design: two threads — thread A performs `map_a = map_b` in a loop; thread B forces GC collection via `g_thread_term_pool().collect()`. Under TSan, this should report a data race on the hash table internals.
  - Run command: `cmake -B build-tsan -DMCRL2_ENABLE_THREADSANITIZER=ON && cmake --build build-tsan && ctest --test-dir build-tsan -R parallel_container`
  - Expected result: TSan reports data race between assignment and GC mark iteration.
- **Efficiency impact:** None — fix adds one lock acquisition per assignment.
- **C++20-specific note:** Consider deleting the defaulted operators and providing explicit ones. With C++20, you could use `= default` inside the function body after the lock guard, but that changes semantics. Custom implementations are clearer.
- **Suggested fix:** Replace `= default` with explicit implementations that acquire `shared_guard` before delegating to `super::operator=`:
  ```cpp
  unordered_map& operator=(const unordered_map& other) {
    mcrl2::utilities::shared_guard guard = detail::g_thread_term_pool().lock_shared();
    super::operator=(other);
    return *this;
  }
  ```

---

### [SEV-2] Inherited `operator[]` and `at()` bypass GC lock — mutating operations exposed without synchronization

- **Location:** `unordered_map.h:40` (public inheritance from `std::unordered_map`)
- **Why this is a problem:** The first class (`atermpp::unordered_map`) publicly inherits from `std::unordered_map`. While mutating operations like `insert`, `emplace`, `erase`, `clear` are overridden to acquire the GC lock, `operator[]` and non-const `at()` are NOT overridden. `operator[]` is particularly dangerous: it may insert a default-constructed element and trigger a rehash — all without the shared lock. This creates the same race condition as SEV-1 (GC mark lambda iterating concurrently with a structural modification).
- **Evidence type:** plausible-issue
- **Evidence:**
  - `grep -n "operator\[" unordered_map.h` returns no matches — `operator[]` is not overridden.
  - `std::unordered_map::operator[]` inserts a default element when key is absent (per [unord.map.elem]).
  - Any user calling `my_atermpp_map[key] = value;` bypasses the GC lock.
  - Candidate test: thread A calls `map[aterm_key]` repeatedly with new keys; thread B forces GC. TSan should report race.
- **Efficiency impact:** None.
- **C++20-specific note:** Consider making the inheritance `private` or `protected` and exposing only the safe API explicitly, or at minimum override `operator[]` and `at()`.
- **Suggested fix:** Override `operator[]` with lock acquisition, or (better) switch to private inheritance with explicit forwarding of the safe subset.

---

### [SEV-2] Constructor `(size_type, value_type, allocator)` is structurally broken — will fail to compile when instantiated

- **Location:** `unordered_map.h:74–77` (first class), `unordered_map.h:277–280` (utilities class)
- **Why this is a problem:** This constructor:
  ```cpp
  unordered_map(size_type n, const value_type& val, const allocator_type& alloc = allocator_type())
   : super::unordered_map(n, detail::markable_aterm<std::pair<const Key, T>>(val), alloc),
  ```
  attempts to pass a `markable_aterm<std::pair<...>>` as the second argument to `std::unordered_map`'s constructor. `std::unordered_map` has no fill-constructor; its relevant constructor is `(size_type bucket_count, const hasher& hf, const allocator_type& a)`. The second argument would need to be a `Hash` object. `markable_aterm<std::pair<...>>` is not convertible to `std::hash<markable_aterm<Key>>`. This is dead code that will produce a hard compile error if ever instantiated. It appears to be a copy-paste from a `std::vector`-style API.
- **Evidence type:** static-proof
- **Evidence:**
  - `std::unordered_map` (C++11–C++23) has no constructor with signature `(size_type, value_type, allocator)`.
  - The 3-argument form is `(size_type n, const hasher& hf, const allocator_type& a)` per [unord.map.cnstr].
  - This code survives only because it's a template — SFINAE-by-non-instantiation.
  - Reproducer: `atermpp::unordered_map<aterm, aterm> m(10, std::make_pair(aterm(), aterm()));` → compile error.
- **Efficiency impact:** None (dead code).
- **C++20-specific note:** N/A.
- **Suggested fix:** Delete this constructor entirely. `std::unordered_map` has no fill constructor — this concept doesn't apply to associative containers.

---

### [SEV-3] Incorrect Doxygen `\file` path

- **Location:** `unordered_map.h:11`
- **Why this is a problem:** The `\file` directive says `mcrl2/data/standard_containers/unordered_map.h` but the actual path is `mcrl2/atermpp/standard_containers/unordered_map.h`. This causes Doxygen to create a phantom documentation page for a non-existent file.
- **Evidence type:** static-proof
- **Evidence:** Compare line 11 (`mcrl2/data/...`) with the actual filesystem path.
- **Suggested fix:** Change to `\file mcrl2/atermpp/standard_containers/unordered_map.h`.

---

### [SEV-3] Dead code in `markable_aterm` fundamental-type specialization (transitive dependency)

- **Location:** `detail/aterm_container.h:128–131` and `detail/aterm_container.h:137–138`
- **Why this is a problem:** Duplicated unreachable statements after `return`:
  ```cpp
  markable_aterm& operator=(const T& other) noexcept {
    m_t = other;
    return *this;
    m_t = other;    // ← dead code
    return *this;   // ← dead code
  }
  ```
  Same pattern in the move-assignment operator and in `operator std_pair&()` (line 245). While harmless at runtime, it suggests a faulty merge or accidental duplication that should be cleaned up.
- **Evidence type:** static-proof
- **Suggested fix:** Remove the duplicated lines after the first `return *this;`.

---

### [SEV-3] Missing `[[nodiscard]]` on query functions

- **Location:** `unordered_map.h:196` (`size()`), lines 219–233 (`begin`, `end`, `cbegin`, `cend`, `empty`, `max_size`)
- **Why this is a problem:** Per the project's C++ instructions, pure queries and must-use results should be marked `[[nodiscard]]`. Discarding the return value of `empty()` or `size()` is almost certainly a bug.
- **Evidence type:** plausible-issue (no runtime failure, but violates project convention for new/touched code).
- **Suggested fix:** Add `[[nodiscard]]` to `size()`, `empty()`, `max_size()`, and the iterator accessors.

---

### [SEV-3] Redundant `#pragma once` with traditional include guard

- **Location:** `unordered_map.h:20–21`
- **Why this is a problem:** Both `#ifndef MCRL2_ATERMPP_...` guard and `#pragma once` are present. Harmless but redundant. The project appears to use include guards as the standard mechanism.
- **Evidence type:** style observation.
- **Suggested fix:** Remove `#pragma once` for consistency, or adopt it project-wide.

---

## Proposed Tests

- [ ] **`test_assignment_race`**: Two threads — thread A performs copy/move assignment of `atermpp::unordered_map` in a loop (source with 1000 elements); thread B calls `g_thread_term_pool().collect()` in a loop. Build with TSan. Expected: TSan reports data race (currently) → after fix, no report.
- [ ] **`test_operator_bracket_race`**: Thread A calls `map[new_key]` (inserting) in a loop; thread B forces GC. Build with TSan. Expected: TSan reports race on hash table internals.
- [ ] **`test_bogus_fill_constructor`**: Compile-time test that explicitly instantiates the `(size_type, value_type, allocator)` constructor. Expected: compilation failure confirming it's broken dead code.

## Open Questions

1. Is the public inheritance from `std::unordered_map` intentional for performance (avoiding virtual dispatch) or a legacy design choice? Private inheritance with explicit forwarding would be safer and prevent the lock-bypass problem for any future std::unordered_map additions (e.g., C++23 `insert_range`).
2. Are there any callers relying on the `(size_type, value_type, allocator)` constructor? A codebase search found no instantiations — safe to remove?
3. The `utilities::unordered_map` variant conditionally acquires the lock only in `find()` when `ThreadSafe=true`, but always acquires it in `insert()`. Is the intention that non-ThreadSafe instances are never accessed concurrently with GC? If so, this relies on the single-threaded code never triggering GC, which may be fragile.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         10m 9s
Total session time:     7m 20s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         1.1m in, 16.5k out, 1.0m cached (Est. 1 Premium request)
 claude-haiku-4.5        1.7m in, 19.7k out, 1.6m cached (Est. 0 Premium requests)

```

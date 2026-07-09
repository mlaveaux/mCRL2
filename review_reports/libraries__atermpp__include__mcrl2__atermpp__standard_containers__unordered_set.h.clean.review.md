---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_set.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T13:40:09.292360+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_set.h` and its implementation file `detail/unordered_set_implementation.h`
- Risk level: **Medium** — one potential GC-safety gap under specific concurrency conditions; several style/contract violations
- Overall verdict: Functionally sound for the primary `atermpp::unordered_set` (std-based) class. The `atermpp::utilities::unordered_set` variant has a plausible GC protection gap during rehash and a declared-but-undefined member function.

## Findings (ordered by severity)

### [SEV-2] GC protection gap during `rehash()` in non-ThreadSafe `atermpp::utilities::unordered_set::emplace()`
- Location: `detail/unordered_set_implementation.h:284-294` (rehash_if_needed) and `detail/unordered_set_implementation.h:329-337` (emplace)
- Why this is a problem: `rehash_if_needed()` acquires the shared GC guard, checks load factor, then **releases the guard** (`guard.unlock_shared()`) before calling `rehash()`. When `ThreadSafe=false`, `rehash()` executes without any GC protection (line 278-280). During rehash, elements are spliced between bucket lists via the base class implementation. If another thread triggers GC in this window, the GC mark phase iterates the container (via `container_wrapper`'s lambda) while it is structurally inconsistent — elements temporarily reside in a stack-local `old_keys` list unreachable from the container. This could cause aterms to be missed during marking and prematurely collected.
- Evidence type: **plausible-issue**
- Evidence:
  - Triggering scenario: Multi-threaded program using `atermpp::utilities::unordered_set<..., ThreadSafe=false>` where thread A calls `emplace()` triggering rehash, and thread B triggers GC during the rehash window.
  - Expected failure: Dangling aterm reference after GC incorrectly frees a term that was in the set but unreachable during the mark due to ongoing rehash.
  - Validation plan: Build with `MCRL2_ENABLE_THREADSANITIZER=ON`. Create a test with a shared `atermpp::utilities::unordered_set<aterm, ..., false>`, one thread calling `emplace()` in a tight loop (to trigger frequent rehash), another thread triggering `detail::g_thread_term_pool().collect()`. TSan should report a data race on the bucket list internals.
  - Run command: `cmake -DMCRL2_ENABLE_THREADSANITIZER=ON ... && ctest -R <test_name>`
- Efficiency impact: None (design issue, not performance)
- C++20-specific note: N/A
- Suggested fix: For the non-ThreadSafe case, `rehash()` should still acquire at minimum the shared GC guard (or exclusive) to prevent GC from marking during structural modification. Alternatively, avoid releasing the shared guard before rehash by restructuring `rehash_if_needed()` to keep protection throughout.

### [SEV-3] `rehash_is_needed()` declared but never defined — shadows base class
- Location: `unordered_set.h:307`
- Why this is a problem: `atermpp::utilities::unordered_set` declares `bool rehash_is_needed() const;` in its `protected` section. No definition exists in the implementation file. The base class `mcrl2::utilities::unordered_set` has a public `rehash_is_needed()` with a working definition, but the derived class's declaration hides it. Any call through the derived type would produce a linker error (`undefined reference`).
- Evidence type: **static-proof**
- Evidence:
  - `grep -rn "rehash_is_needed" libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_set_implementation.h` — returns no hits
  - The declaration at line 307 has no corresponding out-of-line definition anywhere
- Efficiency impact: None (dead code)
- C++20-specific note: N/A
- Suggested fix: Remove the `rehash_is_needed()` declaration from the derived class to un-hide the base class's working version, or provide a definition that delegates to `super::rehash_is_needed()`.

### [SEV-3] Unconstrained universal-reference `insert(P&&)` templates
- Location: `unordered_set.h:122-123` and `unordered_set.h:127-128`
- Why this is a problem: `template<class P> std::pair<iterator, bool> insert(P&& value);` accepts any type without a concept or `requires` clause. Per mCRL2 coding guidelines, SFINAE and unconstrained templates should be replaced by concepts. This produces poor error messages on misuse and violates the project's C++20 style requirements.
- Evidence type: **static-proof** (style rule violation)
- Evidence: No `requires` clause, no concept, no `static_assert` in the implementation (lines 39-46 of impl file simply forward to `super::insert`).
- Efficiency impact: None
- C++20-specific note: Should use `requires std::constructible_from<value_type, P>` or similar constraint.
- Suggested fix:
  ```cpp
  template<class P> requires std::constructible_from<value_type, P>
  std::pair<iterator, bool> insert(P&& value);
  ```

### [SEV-3] Inconsistent GC protection: `size()` inline without lock vs `empty()`/`max_size()` with lock
- Location: `unordered_set.h:172-175` (inline `size()`) vs `detail/unordered_set_implementation.h:251-255` (`empty()` with lock)
- Why this is a problem: `size()` calls `super::size()` without acquiring the shared GC guard. `empty()` and `max_size()` — which are semantically equivalent read-only queries — acquire the lock. While `size()` is likely safe (reads a stored counter, doesn't dereference aterm pointers), the inconsistency makes the protection contract unclear and error-prone for maintainers. Applies to both the `atermpp::unordered_set` and `atermpp::utilities::unordered_set` variants.
- Evidence type: **plausible-issue**
- Evidence: If `std::unordered_set::size()` implementation accesses internal data that GC's exclusive lock protects, this would be a data race. In practice, libstdc++/libc++ `size()` reads a single `size_t` member, so it's safe. But this is an implementation assumption, not a contract guarantee.
- Efficiency impact: Negligible (one fewer lock acquire per `size()` call)
- C++20-specific note: N/A
- Suggested fix: Either add the lock to `size()` for consistency, or document with a comment why `size()` is intentionally lock-free.

### [SEV-3] Missing `[[nodiscard]]` on pure query methods
- Location: `unordered_set.h:172` (`size()`), `unordered_set.h:190` (`empty()`), `unordered_set.h:193` (`max_size()`)
- Why this is a problem: Per mCRL2 C++ guidelines, pure queries and must-use results should be marked `[[nodiscard]]`. Discarding the result of `size()` or `empty()` is almost certainly a bug.
- Evidence type: **static-proof** (style rule)
- Efficiency impact: None
- C++20-specific note: `[[nodiscard]]` has been available since C++17.
- Suggested fix: Add `[[nodiscard]]` to `size()`, `empty()`, `max_size()`, and `find()`.

### [SEV-3] Missing Doxygen `\pre`/`\post` contracts on public API
- Location: Throughout both classes
- Why this is a problem: Per mCRL2 requirements, new or changed public APIs need Doxygen `\pre`/`\post` documentation. Most methods have only `/// \brief Constructor.` or no documentation. Key contracts that should be documented: what types `Key` must satisfy, thread-safety guarantees, iterator invalidation rules, and the GC protection semantics.
- Evidence type: **static-proof** (documentation requirement)
- Suggested fix: Add at minimum `\pre` for `insert`/`emplace` (element must be a valid aterm), `\post` for modifiers (element is GC-protected), and thread-safety annotations.

## Proposed Tests
- [ ] `test_rehash_gc_race`: Purpose: expose data race during `atermpp::utilities::unordered_set` rehash with concurrent GC. Input: ThreadSafe=false set, one thread emplacing elements until rehash triggers, another thread forcing GC. Expected failure: TSan data race report on bucket internals during rehash.
- [ ] `test_rehash_is_needed_link`: Purpose: verify `rehash_is_needed()` is callable through derived type. Input: `atermpp::utilities::unordered_set<aterm> s; s.rehash_is_needed();` — should either compile+link or be removed. Expected failure: linker error (confirms the bug).

## Open Questions
1. Is the `rehash_if_needed()` → unlock → `rehash()` pattern intentional (to avoid holding the shared guard during potentially expensive rehash), and if so, is the GC-safety argument for the non-ThreadSafe case formally justified somewhere?
2. Why does `find()` in `atermpp::utilities::unordered_set` skip the GC lock for `ThreadSafe=false`? If the `find` constructs a temporary `Key(args...)`, that temporary aterm is unprotected during the search. Is this safe only because single-threaded contexts won't have concurrent GC?
3. The `erase(const_iterator pos)` overload is commented out with a TODO about iterator/const_iterator distinction. Is there a tracking issue for this, and does the absence of this overload break standard container compatibility?

## Copilot Raw Stdout

```text
● Explore: Explore atermpp container internals
  ● Read libraries/utilities/include/mcrl2/utilities/unordered_set.h
    └ 403 lines read
  ● Read libraries/utilities/include/mcrl2/utilities/detail/unordered_set_implementation.h
    └ 400 lines read

● Explore: Check defaulted assignment safety
  ● Generate final detailed verification document
    $ cat << 'EOF'

    ╔═══════════════════════════════════════════════════════════════════════════════╗
    ║         DETAILED VERIFICATION: ANSWER TO YOUR THREE QUESTIONS                ║
    ╚═══════════════════════════════════════════════════════════════════════════════╝...
    └ 277 lines...
  ● Final concise summary
    $ cat << 'EOF'

    ╔════════════════════════════════════════════════════════════════════════════════╗
    ║           FINAL ANSWER: Your Analysis is CORRECT & SAFE ✓✓✓                  ║
    ╚════════════════════════════════════════════════════════════════════════════════╝...
    └ 255 lines...

● Explore: Check locking and thread safety
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 1-100
    └ 99 lines read
  ● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 100-250
    └ 99 lines read

● Explore: Check inconsistency and SFINAE issues
  ● Get rehash_is_needed implementation
    $ cat /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/detail/unordered_set_impleme
    ntation.h | grep -A 6 "bool MCRL2_UNORDERED_SET_CLASS::rehash_is_needed"
    └ 8 lines...
  ● Find insert template lines
    $ grep -n "template<class P>"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_set.h |
    head -5
    └ 3 lines...

● Read
  libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_set_implementation.h
  └ 427 lines read

● Grep "unlock_shared"
(/workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/shared_mutex.h)
  └ 4 lines found

● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 1-50
  └ 50 lines read

● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 215-310
  └ 93 lines read

● Grep "lock_shared"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h)
  └ 7 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 145-165
  └ 20 lines read

● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 155-215
  └ 61 lines read

● Check exact content of lock_shared line
  $ sed -n '150,156p'
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h | cat -A
  └ 8 lines...

● Find all lock method signatures
  $ grep -n "lock_shared\|lock()"
  /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/shared_mutex.h | head -20
  └ 19 lines...

● Check if build exists
  $ cd /workspaces/mCRL2 && ls build/compile_commands.json 2>/dev/null && head -5
  build/compile_commands.json 2>/dev/null || echo "No compile commands"
  └ 7 lines...

● Find compilation of thread_aterm_pool
  $ cd /workspaces/mCRL2 && grep -l "thread_aterm_pool.h" build/compile_commands.json | head -1 &&
  grep "thread_aterm_pool" build/compile_commands.json | head -3
  └ 1 line...

✗ Find atermpp cpp files in build
  $ cd /workspaces/mCRL2 && python3 -c "
  import json
  with open('build/compile_commands.json') as f:
      cmds = json.load(f)
  # Find any .cpp that includes atermpp...
  Permission denied and could not request permission from user

● Find atermpp source files
  $ find /workspaces/mCRL2/libraries/atermpp/source -name "*.cpp" | head -5
  └ 6 lines...

● Check lock_shared usage in source
  $ cd /workspaces/mCRL2 && grep -rn "lock_shared\(\)" libraries/atermpp/source/ | head -10
  └ 3 lines...

● Read libraries/atermpp/source/aterm_implementation.cpp lines 55-65
  └ 11 lines read

✗ Verify exact bytes of line 152
  $ cd /workspaces/mCRL2 && cat -A
  libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h | sed -n '152p'
  Permission denied and could not request permission from user

✗ Check hex content of line 152
  $ sed -n '152p'
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h | xxd | head
  -5
  Permission denied and could not request permission from user

✗ Get exact line 152
  $ awk 'NR==152'
  /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h
  Permission denied and could not request permission from user

● Read libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h lines 151-156
  └ 5 lines read

● Check recent changes to shared_mutex files
  $ cd /workspaces/mCRL2 && git --no-pager log --oneline -5
  libraries/utilities/include/mcrl2/utilities/shared_mutex.h && echo "---" && git --no-pager log
  --oneline -5 libraries/atermpp/include/mcrl2/atermpp/detail/thread_aterm_pool.h
  └ 12 lines...

● Check what changed in shared_mutex
  $ cd /workspaces/mCRL2 && git --no-pager show 16ae8d0c5b --stat --
  libraries/utilities/include/mcrl2/utilities/shared_mutex.h | tail -5
  └ 5 lines...

● Check structural improvement commit
  $ cd /workspaces/mCRL2 && git --no-pager show a69ae98eb5 --
  libraries/utilities/include/mcrl2/utilities/shared_mutex.h | head -80
  └ 81 lines...

● Full history of shared_mutex.h
  $ cd /workspaces/mCRL2 && git --no-pager log --oneline -10 --all --
  libraries/utilities/include/mcrl2/utilities/shared_mutex.h
  └ 11 lines...

✗ Check dates of shared_mutex changes
  $ cd /workspaces/mCRL2 && git --no-pager log --oneline --format="%h %ai %s" -3 --
  libraries/utilities/include/mcrl2/utilities/shared_mutex.h
  Permission denied and could not request permission from user

● Check the "belongs to previous commit" change
  $ cd /workspaces/mCRL2 && git --no-pager show --stat 16ae8d0c5b | head -10
  └ 9 lines...

● See what changed between the two commits
  $ cd /workspaces/mCRL2 && git --no-pager diff a69ae98eb5..16ae8d0c5b --
  libraries/utilities/include/mcrl2/utilities/shared_mutex.h
  └ 15 lines...

● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 25-100
  └ 76 lines read

● Check file length
  $ wc -l /workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/shared_mutex.h
  └ 2 lines...

● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 100-160
  └ 60 lines read

● Read libraries/utilities/include/mcrl2/utilities/shared_mutex.h lines 160-307
  └ 148 lines read

● Check if shared_mutex.h has local changes
  $ cd /workspaces/mCRL2 && git --no-pager status --
  libraries/utilities/include/mcrl2/utilities/shared_mutex.h && echo "---" && git --no-pager diff
  HEAD -- libraries/utilities/include/mcrl2/utilities/shared_mutex.h | head -50
  └ 59 lines...

● Check git status and branch
  $ cd /workspaces/mCRL2 && git --no-pager status --short && echo "---BRANCH---" && git --no-pager
  branch --show-current
  └ 12 lines...

● View committed version of shared_mutex.h
  $ cd /workspaces/mCRL2 && git --no-pager show
  HEAD:libraries/utilities/include/mcrl2/utilities/shared_mutex.h | head -100
  └ 101 lines...

● Find lock_shared return type in committed version
  $ cd /workspaces/mCRL2 && git --no-pager show
  HEAD:libraries/utilities/include/mcrl2/utilities/shared_mutex.h | grep -n
  "lock_shared\|shared_guard" | head -20
  └ 21 lines...

✗ Verify shared_mutex::lock_shared return type
  $ cd /workspaces/mCRL2 && git --no-pager show
  HEAD:libraries/utilities/include/mcrl2/utilities/shared_mutex.h | sed -n '200,215p'
  Permission denied and could not request permission from user

● Verify lock_shared return type in committed version
  $ cd /workspaces/mCRL2 && git --no-pager show
  HEAD:libraries/utilities/include/mcrl2/utilities/shared_mutex.h | grep -n "shared_guard
  lock_shared\|lock_shared_impl"
  └ 5 lines...

● Find rehash_is_needed definition
  $ cd /workspaces/mCRL2 && grep -rn "rehash_is_needed" libraries/ | head -10
  └ 10 lines...

● Check base class rehash_is_needed visibility
  $ cd /workspaces/mCRL2 && grep -B5 -A2 "rehash_is_needed"
  libraries/utilities/include/mcrl2/utilities/unordered_set.h
  └ 9 lines...

● Check access specifier for base rehash_is_needed
  $ cd /workspaces/mCRL2 && grep -B2 "rehash_is_needed"
  libraries/utilities/include/mcrl2/utilities/unordered_set.h | grep "public\|private\|protected" |
  tail -1 && echo "---" && grep -n "public\|private\|protected"
  libraries/utilities/include/mcrl2/utilities/unordered_set.h | head -20
  └ 8 lines...

● Find access level of rehash_is_needed in base
  $ cd /workspaces/mCRL2 && grep -n "rehash_is_needed\|public\|private\|protected"
  libraries/utilities/include/mcrl2/utilities/unordered_set.h | grep -B1 "rehash_is_needed"
  └ 3 lines...

● Read libraries/utilities/include/mcrl2/utilities/unordered_set.h lines 340-360
  └ 21 lines read

● Check context around line 151
  $ cd /workspaces/mCRL2 && sed -n '147,155p'
  libraries/utilities/include/mcrl2/utilities/unordered_set.h
  └ 10 lines...

● Check access level context near rehash_is_needed
  $ cd /workspaces/mCRL2 && sed -n '320,355p'
  libraries/utilities/include/mcrl2/utilities/unordered_set.h
  └ 37 lines...

● Check type of m_term_set
  $ cd /workspaces/mCRL2 && grep -B5 "m_term_set"
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h | grep -A5
  "m_term_set.rehash_is_needed" | head -3 && echo "---" && grep "m_term_set"
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage.h | head -5
  └ 9 lines...

● Check unordered_set typedef in aterm_pool_storage
  $ cd /workspaces/mCRL2 && grep -n "using unordered_set\|typedef.*unordered_set\|unordered_set
  m_term_set" libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage.h | head -5
  └ 3 lines...

● Check rehash function definitions in impl
  $ cd /workspaces/mCRL2 && grep -n "rehash_is_needed\|rehash_if_needed"
  libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_set_implementation.h
  └ 3 lines...

● Check if redundant guards are common
  $ cd /workspaces/mCRL2 && head -22
  libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_map.h | tail -5
  └ 6 lines...

Now I have enough information to produce the review. Here's my analysis:

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/standard_containers/unordered_set.h` and its implementation file `detail/unordered_set_implementation.h`
- Risk level: **Medium** — one potential GC-safety gap under specific concurrency conditions; several style/contract violations
- Overall verdict: Functionally sound for the primary `atermpp::unordered_set` (std-based) class. The `atermpp::utilities::unordered_set` variant has a plausible GC protection gap during rehash and a declared-but-undefined member function.

## Findings (ordered by severity)

### [SEV-2] GC protection gap during `rehash()` in non-ThreadSafe `atermpp::utilities::unordered_set::emplace()`
- Location: `detail/unordered_set_implementation.h:284-294` (rehash_if_needed) and `detail/unordered_set_implementation.h:329-337` (emplace)
- Why this is a problem: `rehash_if_needed()` acquires the shared GC guard, checks load factor, then **releases the guard** (`guard.unlock_shared()`) before calling `rehash()`. When `ThreadSafe=false`, `rehash()` executes without any GC protection (line 278-280). During rehash, elements are spliced between bucket lists via the base class implementation. If another thread triggers GC in this window, the GC mark phase iterates the container (via `container_wrapper`'s lambda) while it is structurally inconsistent — elements temporarily reside in a stack-local `old_keys` list unreachable from the container. This could cause aterms to be missed during marking and prematurely collected.
- Evidence type: **plausible-issue**
- Evidence:
  - Triggering scenario: Multi-threaded program using `atermpp::utilities::unordered_set<..., ThreadSafe=false>` where thread A calls `emplace()` triggering rehash, and thread B triggers GC during the rehash window.
  - Expected failure: Dangling aterm reference after GC incorrectly frees a term that was in the set but unreachable during the mark due to ongoing rehash.
  - Validation plan: Build with `MCRL2_ENABLE_THREADSANITIZER=ON`. Create a test with a shared `atermpp::utilities::unordered_set<aterm, ..., false>`, one thread calling `emplace()` in a tight loop (to trigger frequent rehash), another thread triggering `detail::g_thread_term_pool().collect()`. TSan should report a data race on the bucket list internals.
  - Run command: `cmake -DMCRL2_ENABLE_THREADSANITIZER=ON ... && ctest -R <test_name>`
- Efficiency impact: None (design issue, not performance)
- C++20-specific note: N/A
- Suggested fix: For the non-ThreadSafe case, `rehash()` should still acquire at minimum the shared GC guard (or exclusive) to prevent GC from marking during structural modification. Alternatively, avoid releasing the shared guard before rehash by restructuring `rehash_if_needed()` to keep protection throughout.

### [SEV-3] `rehash_is_needed()` declared but never defined — shadows base class
- Location: `unordered_set.h:307`
- Why this is a problem: `atermpp::utilities::unordered_set` declares `bool rehash_is_needed() const;` in its `protected` section. No definition exists in the implementation file. The base class `mcrl2::utilities::unordered_set` has a public `rehash_is_needed()` with a working definition, but the derived class's declaration hides it. Any call through the derived type would produce a linker error (`undefined reference`).
- Evidence type: **static-proof**
- Evidence:
  - `grep -rn "rehash_is_needed" libraries/atermpp/include/mcrl2/atermpp/standard_containers/detail/unordered_set_implementation.h` — returns no hits
  - The declaration at line 307 has no corresponding out-of-line definition anywhere
- Efficiency impact: None (dead code)
- C++20-specific note: N/A
- Suggested fix: Remove the `rehash_is_needed()` declaration from the derived class to un-hide the base class's working version, or provide a definition that delegates to `super::rehash_is_needed()`.

### [SEV-3] Unconstrained universal-reference `insert(P&&)` templates
- Location: `unordered_set.h:122-123` and `unordered_set.h:127-128`
- Why this is a problem: `template<class P> std::pair<iterator, bool> insert(P&& value);` accepts any type without a concept or `requires` clause. Per mCRL2 coding guidelines, SFINAE and unconstrained templates should be replaced by concepts. This produces poor error messages on misuse and violates the project's C++20 style requirements.
- Evidence type: **static-proof** (style rule violation)
- Evidence: No `requires` clause, no concept, no `static_assert` in the implementation (lines 39-46 of impl file simply forward to `super::insert`).
- Efficiency impact: None
- C++20-specific note: Should use `requires std::constructible_from<value_type, P>` or similar constraint.
- Suggested fix:
  ```cpp
  template<class P> requires std::constructible_from<value_type, P>
  std::pair<iterator, bool> insert(P&& value);
  ```

### [SEV-3] Inconsistent GC protection: `size()` inline without lock vs `empty()`/`max_size()` with lock
- Location: `unordered_set.h:172-175` (inline `size()`) vs `detail/unordered_set_implementation.h:251-255` (`empty()` with lock)
- Why this is a problem: `size()` calls `super::size()` without acquiring the shared GC guard. `empty()` and `max_size()` — which are semantically equivalent read-only queries — acquire the lock. While `size()` is likely safe (reads a stored counter, doesn't dereference aterm pointers), the inconsistency makes the protection contract unclear and error-prone for maintainers. Applies to both the `atermpp::unordered_set` and `atermpp::utilities::unordered_set` variants.
- Evidence type: **plausible-issue**
- Evidence: If `std::unordered_set::size()` implementation accesses internal data that GC's exclusive lock protects, this would be a data race. In practice, libstdc++/libc++ `size()` reads a single `size_t` member, so it's safe. But this is an implementation assumption, not a contract guarantee.
- Efficiency impact: Negligible (one fewer lock acquire per `size()` call)
- C++20-specific note: N/A
- Suggested fix: Either add the lock to `size()` for consistency, or document with a comment why `size()` is intentionally lock-free.

### [SEV-3] Missing `[[nodiscard]]` on pure query methods
- Location: `unordered_set.h:172` (`size()`), `unordered_set.h:190` (`empty()`), `unordered_set.h:193` (`max_size()`)
- Why this is a problem: Per mCRL2 C++ guidelines, pure queries and must-use results should be marked `[[nodiscard]]`. Discarding the result of `size()` or `empty()` is almost certainly a bug.
- Evidence type: **static-proof** (style rule)
- Efficiency impact: None
- C++20-specific note: `[[nodiscard]]` has been available since C++17.
- Suggested fix: Add `[[nodiscard]]` to `size()`, `empty()`, `max_size()`, and `find()`.

### [SEV-3] Missing Doxygen `\pre`/`\post` contracts on public API
- Location: Throughout both classes
- Why this is a problem: Per mCRL2 requirements, new or changed public APIs need Doxygen `\pre`/`\post` documentation. Most methods have only `/// \brief Constructor.` or no documentation. Key contracts that should be documented: what types `Key` must satisfy, thread-safety guarantees, iterator invalidation rules, and the GC protection semantics.
- Evidence type: **static-proof** (documentation requirement)
- Suggested fix: Add at minimum `\pre` for `insert`/`emplace` (element must be a valid aterm), `\post` for modifiers (element is GC-protected), and thread-safety annotations.

## Proposed Tests
- [ ] `test_rehash_gc_race`: Purpose: expose data race during `atermpp::utilities::unordered_set` rehash with concurrent GC. Input: ThreadSafe=false set, one thread emplacing elements until rehash triggers, another thread forcing GC. Expected failure: TSan data race report on bucket internals during rehash.
- [ ] `test_rehash_is_needed_link`: Purpose: verify `rehash_is_needed()` is callable through derived type. Input: `atermpp::utilities::unordered_set<aterm> s; s.rehash_is_needed();` — should either compile+link or be removed. Expected failure: linker error (confirms the bug).

## Open Questions
1. Is the `rehash_if_needed()` → unlock → `rehash()` pattern intentional (to avoid holding the shared guard during potentially expensive rehash), and if so, is the GC-safety argument for the non-ThreadSafe case formally justified somewhere?
2. Why does `find()` in `atermpp::utilities::unordered_set` skip the GC lock for `ThreadSafe=false`? If the `find` constructs a temporary `Key(args...)`, that temporary aterm is unprotected during the search. Is this safe only because single-threaded contexts won't have concurrent GC?
3. The `erase(const_iterator pos)` overload is commented out with a TODO about iterator/const_iterator distinction. Is there a tracking issue for this, and does the absence of this overload break standard container compatibility?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         16m 28s
Total session time:     18m 2s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         3.8m in, 28.0k out, 3.8m cached (Est. 1 Premium request)
 claude-haiku-4.5        1.3m in, 30.3k out, 1.2m cached (Est. 0 Premium requests)

```

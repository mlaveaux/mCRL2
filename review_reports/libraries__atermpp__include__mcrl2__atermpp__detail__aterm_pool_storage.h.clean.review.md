---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T10:03:25.833674+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage.h` — header declaring the template class that manages unique term storage and garbage collection.
- **Risk level:** Medium — contains dead code with latent UB, a missing overload definition, and several correctness/maintainability issues.
- **Overall verdict:** The header is functional in its current usage but contains technical debt and latent undefined behavior that would trigger if dormant code paths were ever instantiated.

## Findings (ordered by severity)

### [SEV-2] Fake copy constructor performs silent copy instead of move on const reference
- **Location:** `aterm_pool_storage.h:158-161`
- **Why this is a problem:** `std::move(other.m_term_set)` where `other` is `const aterm_pool_storage&` produces `const unordered_set&&`, which binds to the copy constructor (not the move constructor). The intent is a move but the result is a copy. Additionally, this constructor bypasses `aterm_core`-level reference management and contradicts the `noncopyable` base class.
- **Evidence type:** static-proof
- **Evidence:**
  - `std::move(const T&)` → `const T&&` → binds to `T(const T&)` (copy ctor), not `T(T&&)` (move ctor). See [dcl.init.ref]/5 and `mcrl2::utilities::unordered_set` at `unordered_set.h:246` which declares a copy constructor.
  - The class inherits `private mcrl2::utilities::noncopyable` which deletes copy construction, yet this user-defined constructor compiles because (a) it default-initializes the base (not copy-initializes), and (b) as a member of a class template it is only instantiated on use.
  - The comment references "GCC 4 and 5" which are far below the minimum supported GCC 11.
- **Efficiency impact:** If ever instantiated, silently copies the entire hash table instead of moving it — O(n) when O(1) was intended.
- **C++20-specific note:** Dead code in a template is never checked until instantiated. A `static_assert(false, ...)` or deletion would be safer than leaving latent UB.
- **Suggested fix:** Delete this constructor entirely. If it was needed for container emplacement in ancient GCC, it is no longer needed with move semantics and C++17 guaranteed copy elision. If a move-like constructor is needed, declare `aterm_pool_storage(aterm_pool_storage&&) = default;` (already permitted by the `noncopyable` base).

### [SEV-2] Unconstrained `create_appl_dynamic` overload (lines 123-128) has no definition
- **Location:** `aterm_pool_storage.h:121-128`
- **Why this is a problem:** The header declares a `create_appl_dynamic` overload with template parameters `<InputIterator, TermConverter>` and NO `requires` clause. The implementation file (`aterm_pool_storage_implementation.h`) only defines the two constrained overloads (lines 236-258 and 260-286). If overload resolution ever selects this unconstrained overload (i.e., when the converter matches neither constraint), the result is an undefined symbol at link time.
- **Evidence type:** static-proof
- **Evidence:**
  - The implementation file has exactly 3 definitions of `create_appl_dynamic`: one without converter (line 228), one with `requires std::is_convertible_v<...>` (line 239), and one with `requires std::is_same_v<..., void>` (line 263). No unconstrained with-converter definition exists.
  - Run command to confirm: `grep -n "create_appl_dynamic" libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h`
- **Efficiency impact:** None at runtime (linker catches it), but adds confusion and maintenance burden.
- **C++20-specific note:** In C++20 partial ordering, constrained overloads are preferred over unconstrained ones. So this overload is only selected as a last resort — when the converter doesn't satisfy either constraint. This makes it function as an accidental "negative constraint" trap that produces unhelpful linker errors instead of clear SFINAE/concept failures.
- **Suggested fix:** Remove the unconstrained overload entirely. If a catch-all is needed, add a `static_assert(false, "TermConverter must either return aterm or accept (aterm&, value_type)")` inside a constrained negation overload.

### [SEV-2] `reinterpret_cast` in `call_deletion_hook` acknowledged as questionable via TODO
- **Location:** `aterm_pool_storage_implementation.h:358`
- **Why this is a problem:** `callback(reinterpret_cast<const aterm&>(term))` casts `unprotected_aterm_core` to `const aterm&`. Since `aterm` inherits from `aterm_core` → `unprotected_aterm_core`, this is a downcast through `reinterpret_cast` which is UB unless the object was actually constructed as an `aterm`. The code has a TODO: `"Check whether this cast is OK. It was a static cast."` The original `static_cast` would have been equally undefined but at least required an inheritance relationship check. With `reinterpret_cast`, even that compile-time check is gone.
- **Evidence type:** plausible-issue
- **Evidence:**
  - The actual objects in the set ARE `Element` types (e.g., `_aterm_appl<N>`) which inherit from `_aterm`. The pointer `&term` stored in `unprotected_aterm_core` points to an `_aterm`-derived object. The cast reinterprets the `unprotected_aterm_core` wrapper (containing the pointer) as an `aterm` wrapper. This works because both are single-pointer wrappers with identical layout, but the C++ object model does not guarantee this is defined.
  - Triggering scenario: A callback function that uses dynamic_cast or RTTI on the passed `const aterm&` would observe UB.
  - Test plan: Call `typeid()` on the reference inside a deletion hook callback; compare against expected `aterm` typeinfo.
- **C++20-specific note:** C++20's `std::bit_cast` cannot help here (different sizes). `std::start_lifetime_as` (C++23) would be the correct solution. Current workaround: construct a proper `aterm` from the pointer.
- **Suggested fix:** Replace with explicit construction: `aterm temp(detail::address(term)); callback(temp);` — this goes through the proper `aterm` constructor path.

### [SEV-3] `\threadsafe` annotation on `call_deletion_hook` is misleading
- **Location:** `aterm_pool_storage.h:173`
- **Why this is a problem:** The docstring claims `\threadsafe` but `call_deletion_hook` iterates `m_deletion_hooks` without any synchronization. While `add_deletion_hook` (line 122 of impl) appends to the same vector. If hooks are only added before threads start and `call_deletion_hook` only runs during stop-the-world GC, this is safe in practice — but the `\threadsafe` annotation implies the function can be called concurrently at any time.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `add_deletion_hook` does not acquire any lock (impl line 129: `m_deletion_hooks.emplace_back(...)`).
  - `call_deletion_hook` iterates without locks (impl line 353).
  - If a future change adds deletion hooks at runtime (not inconceivable), this becomes a data race.
  - Test plan: Register a deletion hook from a worker thread while GC is in progress; run under TSan (`cmake -DMCRL2_ENABLE_THREADSANITIZER=ON`).
- **Suggested fix:** Either (a) add a `static_assert` or runtime check that hooks are only added during single-threaded init, or (b) remove the `\threadsafe` annotation and document the actual contract ("must be called only during stop-the-world GC"), or (c) protect with `mcrl2::utilities::shared_mutex`.

### [SEV-3] Placement new in `emplace` bypasses `aterm_core` destructor
- **Location:** `aterm_pool_storage_implementation.h:414`
- **Why this is a problem:** `new (&term) atermpp::unprotected_aterm_core(&*it)` overwrites a live `aterm_core` object (which has a non-trivial destructor per `aterm_core.h:158`) without calling its destructor first. Per [basic.life]/1, this ends the old object's lifetime without invoking its destructor — technically UB if the destructor has side effects (which `~aterm_core()` does: reference counting / protection management).
- **Evidence type:** plausible-issue
- **Evidence:**
  - `aterm_core.h:169` explicitly states: "This class has a non-trivial destructor."
  - This pattern is pervasive in atermpp and appears intentional (term creation with busy flag held bypasses normal protection). The protection system likely accounts for this.
  - However, relying on this pattern makes the code fragile to future refactoring.
  - Test plan: Add ASan annotations to `~aterm_core` to track missed destructions; run existing tests under `MCRL2_ENABLE_ADDRESSSANITIZER=ON`.
- **Suggested fix:** Consider using a designated "raw assign" method (`unprotected_assign`) which is already available in the API (line 196-197 of `aterm_core.h`) instead of placement new.

### [SEV-3] Missing `[[nodiscard]]` on pure query functions
- **Location:** `aterm_pool_storage.h:70,149,155,164,167`
- **Why this is a problem:** `capacity()`, `resize_is_needed()`, `size()`, `verify_mark()`, `verify_sweep()` are pure queries whose return values should never be discarded. Discarding them is always a bug.
- **Evidence type:** static-proof
- **Suggested fix:** Add `[[nodiscard]]` to all query functions.

### [SEV-3] Redundant `typename` before `std::conditional_t` in type alias
- **Location:** `aterm_pool_storage.h:50,52`
- **Why this is a problem:** In C++20 (P0634R3), `typename` is optional in template argument contexts where a type is syntactically required. `typename std::conditional_t<...>` is redundant since `conditional_t` is already a type alias. Not harmful but adds noise.
- **C++20-specific note:** This is a C++17 leftover; C++20 relaxed the requirement.
- **Suggested fix:** Remove the `typename` keywords.

## Proposed Tests
- [ ] `test_deletion_hook_thread_race`: Register deletion hook from worker thread during GC cycle, run under TSan. Expected: TSan report if race exists.
- [ ] `test_reinterpret_cast_deletion_hook`: Install deletion hook that uses `typeid()` on the `const aterm&` argument. Expected: correct type identity OR expose UB.
- [ ] `test_unconstrained_overload_linker_error`: Construct a converter type that satisfies neither constraint; call `create_appl_dynamic`. Expected: compile error (currently would be linker error — demonstrating the overload's existence is problematic).

## Open Questions
1. Is the fake copy constructor (lines 157-161) still instantiated anywhere, or is it truly dead code? If dead, can it be removed immediately?
2. What is the actual threading contract for `add_deletion_hook` — is it guaranteed to be called only during single-threaded initialization? If so, the `\threadsafe` annotation should be removed and a precondition documented.
3. The placement-new pattern in `emplace` appears intentional but technically undefined — is there an appetite to migrate to the existing `unprotected_assign` API, or is the current pattern considered acceptable given the protection model?

## Copilot Raw Stdout

```text
● Explore: Explore aterm_pool_storage context
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Explore: Explore dependencies and tests
  └ Agent started in background with agent_id: agent-1. You can use read_agent to...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage.h
  └ 214 lines read

● Read: Explore agent (Explore dependencies and tests)
  └ Completed

● Read: Explore agent (Explore aterm_pool_storage context)
  └ Running

● Read: Explore agent (Explore aterm_pool_storage context)
  └ Completed

● Grep "create_appl_dynamic" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ 20 lines found

● Grep "std::move\(other\." (/workspaces/mCRL2/libraries/atermpp)
  └ 5 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
  225-290
  └ 66 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h lines 315-380
  └ 66 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
  345-365
  └ 21 lines read

● Grep "class unordered_set|noncopyable|unordered_set\(const|unordered_set\(unordered_set&&"
(/workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/unordered_set.h)
  └ 5 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
  115-135
  └ 21 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_implementation.h lines 40-90
  └ 51 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
  395-430
  └ 36 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines 1-30
  └ 30 lines read

● Grep "class (aterm_core|unprotected_aterm_core)\b"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp)
  └ 4 lines found

● Grep "~aterm_core|destructor|trivial"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm_core.h)
  └ 3 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 33-80
  └ 48 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 150-200
  └ 51 lines read

Now I have all the context needed. Let me produce the critical review.

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage.h` — header declaring the template class that manages unique term storage and garbage collection.
- **Risk level:** Medium — contains dead code with latent UB, a missing overload definition, and several correctness/maintainability issues.
- **Overall verdict:** The header is functional in its current usage but contains technical debt and latent undefined behavior that would trigger if dormant code paths were ever instantiated.

## Findings (ordered by severity)

### [SEV-2] Fake copy constructor performs silent copy instead of move on const reference
- **Location:** `aterm_pool_storage.h:158-161`
- **Why this is a problem:** `std::move(other.m_term_set)` where `other` is `const aterm_pool_storage&` produces `const unordered_set&&`, which binds to the copy constructor (not the move constructor). The intent is a move but the result is a copy. Additionally, this constructor bypasses `aterm_core`-level reference management and contradicts the `noncopyable` base class.
- **Evidence type:** static-proof
- **Evidence:**
  - `std::move(const T&)` → `const T&&` → binds to `T(const T&)` (copy ctor), not `T(T&&)` (move ctor). See [dcl.init.ref]/5 and `mcrl2::utilities::unordered_set` at `unordered_set.h:246` which declares a copy constructor.
  - The class inherits `private mcrl2::utilities::noncopyable` which deletes copy construction, yet this user-defined constructor compiles because (a) it default-initializes the base (not copy-initializes), and (b) as a member of a class template it is only instantiated on use.
  - The comment references "GCC 4 and 5" which are far below the minimum supported GCC 11.
- **Efficiency impact:** If ever instantiated, silently copies the entire hash table instead of moving it — O(n) when O(1) was intended.
- **C++20-specific note:** Dead code in a template is never checked until instantiated. A `static_assert(false, ...)` or deletion would be safer than leaving latent UB.
- **Suggested fix:** Delete this constructor entirely. If it was needed for container emplacement in ancient GCC, it is no longer needed with move semantics and C++17 guaranteed copy elision. If a move-like constructor is needed, declare `aterm_pool_storage(aterm_pool_storage&&) = default;` (already permitted by the `noncopyable` base).

### [SEV-2] Unconstrained `create_appl_dynamic` overload (lines 123-128) has no definition
- **Location:** `aterm_pool_storage.h:121-128`
- **Why this is a problem:** The header declares a `create_appl_dynamic` overload with template parameters `<InputIterator, TermConverter>` and NO `requires` clause. The implementation file (`aterm_pool_storage_implementation.h`) only defines the two constrained overloads (lines 236-258 and 260-286). If overload resolution ever selects this unconstrained overload (i.e., when the converter matches neither constraint), the result is an undefined symbol at link time.
- **Evidence type:** static-proof
- **Evidence:**
  - The implementation file has exactly 3 definitions of `create_appl_dynamic`: one without converter (line 228), one with `requires std::is_convertible_v<...>` (line 239), and one with `requires std::is_same_v<..., void>` (line 263). No unconstrained with-converter definition exists.
  - Run command to confirm: `grep -n "create_appl_dynamic" libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h`
- **Efficiency impact:** None at runtime (linker catches it), but adds confusion and maintenance burden.
- **C++20-specific note:** In C++20 partial ordering, constrained overloads are preferred over unconstrained ones. So this overload is only selected as a last resort — when the converter doesn't satisfy either constraint. This makes it function as an accidental "negative constraint" trap that produces unhelpful linker errors instead of clear SFINAE/concept failures.
- **Suggested fix:** Remove the unconstrained overload entirely. If a catch-all is needed, add a `static_assert(false, "TermConverter must either return aterm or accept (aterm&, value_type)")` inside a constrained negation overload.

### [SEV-2] `reinterpret_cast` in `call_deletion_hook` acknowledged as questionable via TODO
- **Location:** `aterm_pool_storage_implementation.h:358`
- **Why this is a problem:** `callback(reinterpret_cast<const aterm&>(term))` casts `unprotected_aterm_core` to `const aterm&`. Since `aterm` inherits from `aterm_core` → `unprotected_aterm_core`, this is a downcast through `reinterpret_cast` which is UB unless the object was actually constructed as an `aterm`. The code has a TODO: `"Check whether this cast is OK. It was a static cast."` The original `static_cast` would have been equally undefined but at least required an inheritance relationship check. With `reinterpret_cast`, even that compile-time check is gone.
- **Evidence type:** plausible-issue
- **Evidence:**
  - The actual objects in the set ARE `Element` types (e.g., `_aterm_appl<N>`) which inherit from `_aterm`. The pointer `&term` stored in `unprotected_aterm_core` points to an `_aterm`-derived object. The cast reinterprets the `unprotected_aterm_core` wrapper (containing the pointer) as an `aterm` wrapper. This works because both are single-pointer wrappers with identical layout, but the C++ object model does not guarantee this is defined.
  - Triggering scenario: A callback function that uses dynamic_cast or RTTI on the passed `const aterm&` would observe UB.
  - Test plan: Call `typeid()` on the reference inside a deletion hook callback; compare against expected `aterm` typeinfo.
- **C++20-specific note:** C++20's `std::bit_cast` cannot help here (different sizes). `std::start_lifetime_as` (C++23) would be the correct solution. Current workaround: construct a proper `aterm` from the pointer.
- **Suggested fix:** Replace with explicit construction: `aterm temp(detail::address(term)); callback(temp);` — this goes through the proper `aterm` constructor path.

### [SEV-3] `\threadsafe` annotation on `call_deletion_hook` is misleading
- **Location:** `aterm_pool_storage.h:173`
- **Why this is a problem:** The docstring claims `\threadsafe` but `call_deletion_hook` iterates `m_deletion_hooks` without any synchronization. While `add_deletion_hook` (line 122 of impl) appends to the same vector. If hooks are only added before threads start and `call_deletion_hook` only runs during stop-the-world GC, this is safe in practice — but the `\threadsafe` annotation implies the function can be called concurrently at any time.
- **Evidence type:** plausible-issue
- **Evidence:**
  - `add_deletion_hook` does not acquire any lock (impl line 129: `m_deletion_hooks.emplace_back(...)`).
  - `call_deletion_hook` iterates without locks (impl line 353).
  - If a future change adds deletion hooks at runtime (not inconceivable), this becomes a data race.
  - Test plan: Register a deletion hook from a worker thread while GC is in progress; run under TSan (`cmake -DMCRL2_ENABLE_THREADSANITIZER=ON`).
- **Suggested fix:** Either (a) add a `static_assert` or runtime check that hooks are only added during single-threaded init, or (b) remove the `\threadsafe` annotation and document the actual contract ("must be called only during stop-the-world GC"), or (c) protect with `mcrl2::utilities::shared_mutex`.

### [SEV-3] Placement new in `emplace` bypasses `aterm_core` destructor
- **Location:** `aterm_pool_storage_implementation.h:414`
- **Why this is a problem:** `new (&term) atermpp::unprotected_aterm_core(&*it)` overwrites a live `aterm_core` object (which has a non-trivial destructor per `aterm_core.h:158`) without calling its destructor first. Per [basic.life]/1, this ends the old object's lifetime without invoking its destructor — technically UB if the destructor has side effects (which `~aterm_core()` does: reference counting / protection management).
- **Evidence type:** plausible-issue
- **Evidence:**
  - `aterm_core.h:169` explicitly states: "This class has a non-trivial destructor."
  - This pattern is pervasive in atermpp and appears intentional (term creation with busy flag held bypasses normal protection). The protection system likely accounts for this.
  - However, relying on this pattern makes the code fragile to future refactoring.
  - Test plan: Add ASan annotations to `~aterm_core` to track missed destructions; run existing tests under `MCRL2_ENABLE_ADDRESSSANITIZER=ON`.
- **Suggested fix:** Consider using a designated "raw assign" method (`unprotected_assign`) which is already available in the API (line 196-197 of `aterm_core.h`) instead of placement new.

### [SEV-3] Missing `[[nodiscard]]` on pure query functions
- **Location:** `aterm_pool_storage.h:70,149,155,164,167`
- **Why this is a problem:** `capacity()`, `resize_is_needed()`, `size()`, `verify_mark()`, `verify_sweep()` are pure queries whose return values should never be discarded. Discarding them is always a bug.
- **Evidence type:** static-proof
- **Suggested fix:** Add `[[nodiscard]]` to all query functions.

### [SEV-3] Redundant `typename` before `std::conditional_t` in type alias
- **Location:** `aterm_pool_storage.h:50,52`
- **Why this is a problem:** In C++20 (P0634R3), `typename` is optional in template argument contexts where a type is syntactically required. `typename std::conditional_t<...>` is redundant since `conditional_t` is already a type alias. Not harmful but adds noise.
- **C++20-specific note:** This is a C++17 leftover; C++20 relaxed the requirement.
- **Suggested fix:** Remove the `typename` keywords.

## Proposed Tests
- [ ] `test_deletion_hook_thread_race`: Register deletion hook from worker thread during GC cycle, run under TSan. Expected: TSan report if race exists.
- [ ] `test_reinterpret_cast_deletion_hook`: Install deletion hook that uses `typeid()` on the `const aterm&` argument. Expected: correct type identity OR expose UB.
- [ ] `test_unconstrained_overload_linker_error`: Construct a converter type that satisfies neither constraint; call `create_appl_dynamic`. Expected: compile error (currently would be linker error — demonstrating the overload's existence is problematic).

## Open Questions
1. Is the fake copy constructor (lines 157-161) still instantiated anywhere, or is it truly dead code? If dead, can it be removed immediately?
2. What is the actual threading contract for `add_deletion_hook` — is it guaranteed to be called only during single-threaded initialization? If so, the `\threadsafe` annotation should be removed and a precondition documented.
3. The placement-new pattern in `emplace` appears intentional but technically undefined — is there an appetite to migrate to the existing `unprotected_assign` API, or is the current pattern considered acceptable given the protection model?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         6m 7s
Total session time:     5m 30s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         614.6k in, 12.0k out, 571.2k cached (Est. 1 Premium request)
 claude-haiku-4.5        422.5k in, 11.2k out, 378.6k cached (Est. 0 Premium requests)

```

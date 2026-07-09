---
{
  "file": "libraries/atermpp/source/function_symbol.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T14:16:17.065064+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/source/function_symbol.cpp` (full file, 31 lines) with transitive context of global_aterm_pool.h, function_symbol.h, shared_reference.h, function_symbol_pool.cpp
- Risk level: **Medium** — foundational infrastructure with a subtle lifetime violation and an unguarded lazy-init path
- Overall verdict: Correct in practice due to known initialization ordering (static init is single-threaded, symbols are permanent), but contains **technically undefined behavior** from operating on objects before their lifetime begins, and a **reference count leak** per built-in symbol.

## Findings (ordered by severity)

### [SEV-2] Undefined behavior: `operator=` on `g_as_int` before its lifetime begins
- Location: `libraries/atermpp/source/function_symbol.cpp:15–17` (definitions) + `libraries/atermpp/source/function_symbol_pool.cpp:27–29` (assignments)
- Why this is a problem: The dynamic initializer of `g_as_int` evaluates `g_term_pool<true>().as_int()`. During that expression, `aterm_pool` is placement-constructed, which invokes `function_symbol_pool::function_symbol_pool()`, which calls `g_as_int = m_as_int;`. At that point, `g_as_int` has only been zero-initialized — its constructor has **not** completed. Calling `operator=` on a non-trivially-destructible object whose lifetime has not yet begun is UB per [basic.life]§6.8/6. After the pool constructor returns, the copy-constructor of `g_as_int` fires on the same storage, overwriting `m_reference` without first decrementing the value stored by the prior `operator=` — leaking one reference count per symbol.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Static analysis trace of initialization order
  - Run command: N/A (compile-time reasoning)
  - Observed result: `shared_reference::operator=` increments the underlying `_function_symbol`'s ref count at line `function_symbol_pool.cpp:27`. Then `shared_reference` copy constructor (from the direct-initialization of `g_as_int`) sets `m_reference` to the same pointer and increments again, without decrementing the stale write.
  - Expected result: Each constructed `function_symbol` holding a reference to a `_function_symbol` should contribute exactly +1/−1 to its ref count over the object's lifetime.
- Efficiency impact: Permanent +1 ref-count leak per built-in symbol prevents GC from reclaiming them even if no live references remain (benign for these specific symbols since they are conceptually immortal, but violates invariants for sanitizer/leak detection).
- C++20-specific note: Could be addressed by making `g_as_int` etc. `inline constinit` placeholders that are only ever assigned-to (never double-constructed), or by removing the redundant definition-with-initializer and relying solely on the pool-constructor assignment after a placement-new of the `function_symbol` objects.
- Suggested fix: Define the globals without an initializer expression (rely on zero-initialization + the assignment inside `function_symbol_pool` constructor), OR construct them via `std::construct_at`/placement-new inside the pool constructor rather than using `operator=` on not-yet-alive objects.

### [SEV-2] Non-atomic lazy initialization guard in `g_term_pool<true>()`
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:45–51`
- Why this is a problem: The `static bool initialized` variable plus the unguarded `if (!initialized)` check constitutes a data race if `g_term_pool<true>()` is ever called concurrently. Two threads could both read `false`, enter the body, and perform placement-new on the same storage simultaneously — corrupting the pool. While today all call sites are during static initialization (guaranteed single-threaded), no mechanism enforces this: any future `global_function_symbol` constructed from a thread after `main()` begins (before the pool is initialized by a prior static init in another TU) would trigger the race.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Hypothetical concurrent `global_function_symbol` construction
  - Run command: Build with `MCRL2_ENABLE_THREADSANITIZER=ON`, add a test that constructs a `global_function_symbol` from two threads simultaneously before the pool is otherwise initialized.
  - Observed result: TSan would report a data race on `initialized`.
  - Expected result: Initialization should be atomic (e.g., `std::call_once`, `std::atomic_flag`, or relying on C++11 function-local static with non-trivial init).
- Efficiency impact: Minimal — initialization occurs once.
- C++20-specific note: `std::call_once` or wrapping the pool reference in a function-local static with a non-trivial initializer would leverage the C++11 thread-safe static initialization guarantee at zero ongoing cost.
- Suggested fix: Replace the bare `bool` with `std::call_once` / `std::once_flag`, or restructure so that `g_term_pool<true>()` uses a Meyers singleton (function-local static reference to an object with the lazy construction inside).

### [SEV-3] Missing `[[nodiscard]]` and Doxygen contracts on private constructors
- Location: `libraries/atermpp/source/function_symbol.cpp:20–26`
- Why this is a problem: The three-argument private constructors are the main creation path for runtime function symbols but lack `\pre`/`\post` documentation specifying: (1) `name` must not be empty, (2) the returned symbol is always `defined()`, (3) the `check_for_registered_functions` semantics are undocumented in any header-visible Doxygen. The `global_function_symbol` constructor similarly lacks a `\pre` stating it must only be called during static initialization.
- Evidence type: code-inspection
- Evidence:
  - Reproducer/test name: N/A
  - Run command: N/A
  - Observed result: No documentation on preconditions/thread-safety of these constructors.
  - Expected result: Per project policy, public/semi-public creation paths should have Doxygen `\pre`/`\post` and assertions.
- Efficiency impact: None.
- C++20-specific note: Concepts or `static_assert` could enforce compile-time invariants for `arity` bounds if desired.
- Suggested fix: Add Doxygen comments and `assert(!name.empty())` (or document that empty names are valid).

### [SEV-3] `const` by-value parameters in definitions
- Location: `libraries/atermpp/source/function_symbol.cpp:20, 25`
- Why this is a problem: `const std::size_t arity` and `const bool check_for_registered_functions` are by-value. The top-level `const` is ignored in declarations and is a meaningless implementation detail to callers. It's inconsistent with modern C++ style but consistent with existing mCRL2 code, so this is informational only.
- Evidence type: code-inspection
- Suggested fix: No action needed (consistency with codebase style prevails).

## Proposed Tests
- [ ] `test_global_symbol_refcount`: After pool initialization, verify that `detail::g_as_int`'s underlying `_function_symbol` has the expected reference count (exactly the number of live `function_symbol` objects pointing to it). Failure signal: ref count mismatch indicating the leak.
- [ ] `test_concurrent_global_function_symbol_init` (TSan): Spawn two threads that simultaneously construct `global_function_symbol` instances; run under TSan. Failure signal: TSan data-race report on the `initialized` variable in `g_term_pool<true>()`.

## Open Questions
1. Is the +1 ref-count leak on built-in symbols intentional (to guarantee they are never GC'd)? If so, it should be documented. If not, the double-initialization should be eliminated.
2. Should `g_term_pool<true>()` be hardened for concurrent use, or should a `static_assert` / runtime check enforce it is only called before threads are spawned?
3. The `global_function_symbol` constructor calls `create_function_symbol` on the global pool directly (bypassing `g_thread_term_pool()`). Is there a lock protecting the pool's symbol table during concurrent `function_symbol` creation from other threads that have already registered via `g_thread_term_pool()`?

## Copilot Raw Stdout

```text
● Explore: Explore function_symbol context
  ● Create comprehensive analysis
    $ cat > /tmp/analysis.md << 'EOF'
    # Static Initialization Order Analysis: function_symbol.cpp

    ## File Locations
    - **Main Source**: `/workspaces/mCRL2/libraries/atermpp/source/function_symbol.cpp`...
    └ 549 lines...
  ● Create detailed code snippet document
    $ cat > /tmp/code_snippets.md << 'EOF'
    # Key Code Snippets and Context

    ## 1. FILE: function_symbol.cpp (THE PROBLEM)
    **Location**: `/workspaces/mCRL2/libraries/atermpp/source/function_symbol.cpp`...
    └ 631 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h
  └ 65 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/function_symbol.h
  └ 166 lines read

● Explore: Explore thread pool and init order
  ● Locate shared_reference.h
    $ find /workspaces/mCRL2 -name "shared_reference.h" -type f
    └ 2 lines...
  ● Read libraries/utilities/include/mcrl2/utilities/shared_reference.h lines 1-150
    └ 150 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h lines 40-55
  └ 16 lines read

● Grep "g_term_pool<true>" (/workspaces/mCRL2/libraries/atermpp/source)
  └ 4 lines found

● Grep "function_symbol_pool|g_as_int|g_as_list|g_as_empty"
(/workspaces/mCRL2/libraries/atermpp/source/function_symbol_pool.cpp)
  └ 18 lines found

● Grep "EnableGCStressThread" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ 2 lines found

● Read libraries/atermpp/source/function_symbol_pool.cpp lines 19-30
  └ 12 lines read

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_configuration.h lines 30-38
  └ 9 lines read

Now I have enough context to produce the review. Let me verify one more detail about the initialization ordering:

● Read libraries/utilities/include/mcrl2/utilities/shared_reference.h lines 99-170
  └ 72 lines read

I have all the information needed. Here is the review:

---

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/source/function_symbol.cpp` (full file, 31 lines) with transitive context of global_aterm_pool.h, function_symbol.h, shared_reference.h, function_symbol_pool.cpp
- Risk level: **Medium** — foundational infrastructure with a subtle lifetime violation and an unguarded lazy-init path
- Overall verdict: Correct in practice due to known initialization ordering (static init is single-threaded, symbols are permanent), but contains **technically undefined behavior** from operating on objects before their lifetime begins, and a **reference count leak** per built-in symbol.

## Findings (ordered by severity)

### [SEV-2] Undefined behavior: `operator=` on `g_as_int` before its lifetime begins
- Location: `libraries/atermpp/source/function_symbol.cpp:15–17` (definitions) + `libraries/atermpp/source/function_symbol_pool.cpp:27–29` (assignments)
- Why this is a problem: The dynamic initializer of `g_as_int` evaluates `g_term_pool<true>().as_int()`. During that expression, `aterm_pool` is placement-constructed, which invokes `function_symbol_pool::function_symbol_pool()`, which calls `g_as_int = m_as_int;`. At that point, `g_as_int` has only been zero-initialized — its constructor has **not** completed. Calling `operator=` on a non-trivially-destructible object whose lifetime has not yet begun is UB per [basic.life]§6.8/6. After the pool constructor returns, the copy-constructor of `g_as_int` fires on the same storage, overwriting `m_reference` without first decrementing the value stored by the prior `operator=` — leaking one reference count per symbol.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Static analysis trace of initialization order
  - Run command: N/A (compile-time reasoning)
  - Observed result: `shared_reference::operator=` increments the underlying `_function_symbol`'s ref count at line `function_symbol_pool.cpp:27`. Then `shared_reference` copy constructor (from the direct-initialization of `g_as_int`) sets `m_reference` to the same pointer and increments again, without decrementing the stale write.
  - Expected result: Each constructed `function_symbol` holding a reference to a `_function_symbol` should contribute exactly +1/−1 to its ref count over the object's lifetime.
- Efficiency impact: Permanent +1 ref-count leak per built-in symbol prevents GC from reclaiming them even if no live references remain (benign for these specific symbols since they are conceptually immortal, but violates invariants for sanitizer/leak detection).
- C++20-specific note: Could be addressed by making `g_as_int` etc. `inline constinit` placeholders that are only ever assigned-to (never double-constructed), or by removing the redundant definition-with-initializer and relying solely on the pool-constructor assignment after a placement-new of the `function_symbol` objects.
- Suggested fix: Define the globals without an initializer expression (rely on zero-initialization + the assignment inside `function_symbol_pool` constructor), OR construct them via `std::construct_at`/placement-new inside the pool constructor rather than using `operator=` on not-yet-alive objects.

### [SEV-2] Non-atomic lazy initialization guard in `g_term_pool<true>()`
- Location: `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:45–51`
- Why this is a problem: The `static bool initialized` variable plus the unguarded `if (!initialized)` check constitutes a data race if `g_term_pool<true>()` is ever called concurrently. Two threads could both read `false`, enter the body, and perform placement-new on the same storage simultaneously — corrupting the pool. While today all call sites are during static initialization (guaranteed single-threaded), no mechanism enforces this: any future `global_function_symbol` constructed from a thread after `main()` begins (before the pool is initialized by a prior static init in another TU) would trigger the race.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Hypothetical concurrent `global_function_symbol` construction
  - Run command: Build with `MCRL2_ENABLE_THREADSANITIZER=ON`, add a test that constructs a `global_function_symbol` from two threads simultaneously before the pool is otherwise initialized.
  - Observed result: TSan would report a data race on `initialized`.
  - Expected result: Initialization should be atomic (e.g., `std::call_once`, `std::atomic_flag`, or relying on C++11 function-local static with non-trivial init).
- Efficiency impact: Minimal — initialization occurs once.
- C++20-specific note: `std::call_once` or wrapping the pool reference in a function-local static with a non-trivial initializer would leverage the C++11 thread-safe static initialization guarantee at zero ongoing cost.
- Suggested fix: Replace the bare `bool` with `std::call_once` / `std::once_flag`, or restructure so that `g_term_pool<true>()` uses a Meyers singleton (function-local static reference to an object with the lazy construction inside).

### [SEV-3] Missing `[[nodiscard]]` and Doxygen contracts on private constructors
- Location: `libraries/atermpp/source/function_symbol.cpp:20–26`
- Why this is a problem: The three-argument private constructors are the main creation path for runtime function symbols but lack `\pre`/`\post` documentation specifying: (1) `name` must not be empty, (2) the returned symbol is always `defined()`, (3) the `check_for_registered_functions` semantics are undocumented in any header-visible Doxygen. The `global_function_symbol` constructor similarly lacks a `\pre` stating it must only be called during static initialization.
- Evidence type: code-inspection
- Evidence:
  - Reproducer/test name: N/A
  - Run command: N/A
  - Observed result: No documentation on preconditions/thread-safety of these constructors.
  - Expected result: Per project policy, public/semi-public creation paths should have Doxygen `\pre`/`\post` and assertions.
- Efficiency impact: None.
- C++20-specific note: Concepts or `static_assert` could enforce compile-time invariants for `arity` bounds if desired.
- Suggested fix: Add Doxygen comments and `assert(!name.empty())` (or document that empty names are valid).

### [SEV-3] `const` by-value parameters in definitions
- Location: `libraries/atermpp/source/function_symbol.cpp:20, 25`
- Why this is a problem: `const std::size_t arity` and `const bool check_for_registered_functions` are by-value. The top-level `const` is ignored in declarations and is a meaningless implementation detail to callers. It's inconsistent with modern C++ style but consistent with existing mCRL2 code, so this is informational only.
- Evidence type: code-inspection
- Suggested fix: No action needed (consistency with codebase style prevails).

## Proposed Tests
- [ ] `test_global_symbol_refcount`: After pool initialization, verify that `detail::g_as_int`'s underlying `_function_symbol` has the expected reference count (exactly the number of live `function_symbol` objects pointing to it). Failure signal: ref count mismatch indicating the leak.
- [ ] `test_concurrent_global_function_symbol_init` (TSan): Spawn two threads that simultaneously construct `global_function_symbol` instances; run under TSan. Failure signal: TSan data-race report on the `initialized` variable in `g_term_pool<true>()`.

## Open Questions
1. Is the +1 ref-count leak on built-in symbols intentional (to guarantee they are never GC'd)? If so, it should be documented. If not, the double-initialization should be eliminated.
2. Should `g_term_pool<true>()` be hardened for concurrent use, or should a `static_assert` / runtime check enforce it is only called before threads are spawned?
3. The `global_function_symbol` constructor calls `create_function_symbol` on the global pool directly (bypassing `g_thread_term_pool()`). Is there a lock protecting the pool's symbol table during concurrent `function_symbol` creation from other threads that have already registered via `g_thread_term_pool()`?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         8m 41s
Total session time:     9m 24s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         368.9k in, 13.2k out, 332.2k cached (Est. 1 Premium request)
 claude-haiku-4.5        1.0m in, 18.9k out, 976.7k cached (Est. 0 Premium requests)

```

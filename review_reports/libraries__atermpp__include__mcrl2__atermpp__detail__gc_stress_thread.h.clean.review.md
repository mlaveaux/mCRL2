---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/gc_stress_thread.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T11:57:42.472678+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/gc_stress_thread.h` (header) + implementation in `libraries/atermpp/source/gc_stress_thread.cpp` + call site in `global_aterm_pool.h`
- **Risk level:** Low (feature is disabled by default via `constexpr bool EnableGCStressThread = false`)
- **Overall verdict:** The header is clean and minimal. The implementation has a latent `std::terminate` path if `start_gc_stress_thread()` is called twice when the feature is enabled, and the call site has a data-race on the manual initialization flag. Both only manifest in the (default-off) stress-testing configuration.

## Findings (ordered by severity)

### [SEV-2] Double-call to `start_gc_stress_thread()` causes `std::terminate`
- **Location:** `libraries/atermpp/source/gc_stress_thread.cpp:36`
- **Why this is a problem:** Assigning a new `std::thread` to `g_gc_stress_thread` while it is still joinable (from a prior call) invokes `std::terminate()`. The function has no guard against re-entry.
- **Evidence type:** static-proof
- **Evidence:**
  - `std::thread::operator=` specification ([thread.thread.assign]): if `joinable()` is `true`, calls `std::terminate()` before the assignment.
  - If `start_gc_stress_thread()` is called twice without an intervening `stop`, and `EnableGCStressThread` is `true`, the second assignment to `g_gc_stress_thread` terminates the process.
  - Reproducer: set `EnableGCStressThread = true` in `aterm_configuration.h`, call `start_gc_stress_thread()` twice.
  - Run command: N/A (requires flipping a `constexpr` flag; no existing test exercises this path)
  - Observed result: `std::terminate()` called
  - Expected result: Idempotent or documented single-call precondition
- **Efficiency impact:** None
- **C++20-specific note:** Could use `std::jthread` which destructs cleanly (requests stop + joins), eliminating this class of bug entirely.
- **Suggested fix:** Either (a) guard with `if (g_gc_stress_thread.joinable()) return;` at the top of the enabled branch, or (b) switch to `std::jthread` + `std::stop_token` which makes stop/join automatic and eliminates the global atomic flag.

### [SEV-3] Call site data race on manual `initialized` flag
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:45–51`
- **Why this is a problem:** `static bool initialized` inside `g_term_pool<true>()` is read and written without synchronization. If two threads call `g_term_pool<true>()` concurrently before initialization completes, both may execute the body, double-constructing the pool and double-starting the stress thread (triggering SEV-2 above). C++ guarantees thread-safe *initialization* of block-scope statics, but subsequent reads/writes to an already-initialized `static bool` are ordinary data accesses.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Two threads enter `g_term_pool<true>()`, both read `initialized == false`, both execute the body.
  - Most likely trigger: early concurrent access to the term pool during library startup with `EnableGCStressThread = true`.
  - This is mitigated in practice because `g_term_pool<true>()` is typically called only during static initialization (single-threaded), but the invariant is undocumented.
  - Test plan: Enable `EnableGCStressThread`, spawn two threads that both call `g_term_pool<true>()` simultaneously, run under TSan.
  - Run command: `cmake --build build-tsan && ctest --test-dir build-tsan -R <new_test> --output-on-failure`
- **Efficiency impact:** None
- **C++20-specific note:** Replace with `static std::once_flag` + `std::call_once`, or use a properly guarded `static` object whose constructor does the work (leveraging C++ static init guarantees).
- **Suggested fix:** Use `std::call_once` or move initialization into a static local's constructor to get the language guarantee.

### [SEV-3] Missing `\pre`/`\post` contract annotations
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/gc_stress_thread.h:22–25`
- **Why this is a problem:** Per project convention (`.github/instructions/cpp.instructions.md`), new/changed public APIs require explicit `\pre`/`\post`. `start_gc_stress_thread()` has an implicit precondition (must not be called while a stress thread is already running) that is not documented.
- **Evidence type:** static-proof (convention violation)
- **Suggested fix:**
  ```cpp
  /// \brief Spawns a dedicated background thread that continuously triggers GC.
  /// \pre No GC stress thread is currently running.
  /// \post If EnableGCStressThread && GlobalThreadSafe, a background thread is active.
  void start_gc_stress_thread();
  
  /// \brief Stop the background GC stress thread if it is running.
  /// \post No GC stress thread is running.
  void stop_gc_stress_thread();
  ```

## Proposed Tests
- [ ] `test_gc_stress_double_start`: With `EnableGCStressThread = true`, call `start_gc_stress_thread()` twice without stopping. Currently causes `std::terminate`; after fix, should be idempotent. Signal: process survival vs. abort.
- [ ] `test_gc_stress_start_stop_cycle`: Start, stop, start, stop. Verify no TSan findings and no crash. Run under TSan.
- [ ] `test_concurrent_lazy_pool_init` (for call site): Two threads calling `g_term_pool<true>()` simultaneously under TSan. Signal: TSan data-race report on `initialized`.

## Open Questions
1. Is there a deliberate design assumption that `g_term_pool<true>()` is only ever called from a single thread (during static init)? If so, should it be enforced or at least asserted?
2. Has switching to `std::jthread` + `std::stop_token` been considered? It would eliminate the manual atomic flag, the `atexit` registration, and the double-start `std::terminate` risk in one change.
3. The `std::atexit` registration means `stop_gc_stress_thread()` runs during static destruction. Is there a guaranteed ordering that the thread-local `thread_aterm_pool` for the stress thread is destroyed *before* the global pool storage is reused? (Currently safe because the pool uses placement-new into a byte array that persists, but this is a fragile implicit invariant.)

## Copilot Raw Stdout

```text
● Explore: Find gc_stress_thread implementation
  ● Find GlobalThreadSafe references in atermpp
    $ grep -r "GlobalThreadSafe" /workspaces/mCRL2/libraries/atermpp --include="*.h" -n | head -20
    └ 17 lines...
  ● Find MCRL2_ENABLE_MULTITHREADING references
    $ grep -r "MCRL2_ENABLE_MULTITHREADING" /workspaces/mCRL2/libraries/atermpp --include="*.h" -n |
     head -10
    └ 3 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h lines 1-60
  └ 60 lines read

● Read libraries/atermpp/source/gc_stress_thread.cpp
  └ 58 lines read

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/gc_stress_thread.h` (header) + implementation in `libraries/atermpp/source/gc_stress_thread.cpp` + call site in `global_aterm_pool.h`
- **Risk level:** Low (feature is disabled by default via `constexpr bool EnableGCStressThread = false`)
- **Overall verdict:** The header is clean and minimal. The implementation has a latent `std::terminate` path if `start_gc_stress_thread()` is called twice when the feature is enabled, and the call site has a data-race on the manual initialization flag. Both only manifest in the (default-off) stress-testing configuration.

## Findings (ordered by severity)

### [SEV-2] Double-call to `start_gc_stress_thread()` causes `std::terminate`
- **Location:** `libraries/atermpp/source/gc_stress_thread.cpp:36`
- **Why this is a problem:** Assigning a new `std::thread` to `g_gc_stress_thread` while it is still joinable (from a prior call) invokes `std::terminate()`. The function has no guard against re-entry.
- **Evidence type:** static-proof
- **Evidence:**
  - `std::thread::operator=` specification ([thread.thread.assign]): if `joinable()` is `true`, calls `std::terminate()` before the assignment.
  - If `start_gc_stress_thread()` is called twice without an intervening `stop`, and `EnableGCStressThread` is `true`, the second assignment to `g_gc_stress_thread` terminates the process.
  - Reproducer: set `EnableGCStressThread = true` in `aterm_configuration.h`, call `start_gc_stress_thread()` twice.
  - Run command: N/A (requires flipping a `constexpr` flag; no existing test exercises this path)
  - Observed result: `std::terminate()` called
  - Expected result: Idempotent or documented single-call precondition
- **Efficiency impact:** None
- **C++20-specific note:** Could use `std::jthread` which destructs cleanly (requests stop + joins), eliminating this class of bug entirely.
- **Suggested fix:** Either (a) guard with `if (g_gc_stress_thread.joinable()) return;` at the top of the enabled branch, or (b) switch to `std::jthread` + `std::stop_token` which makes stop/join automatic and eliminates the global atomic flag.

### [SEV-3] Call site data race on manual `initialized` flag
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/global_aterm_pool.h:45–51`
- **Why this is a problem:** `static bool initialized` inside `g_term_pool<true>()` is read and written without synchronization. If two threads call `g_term_pool<true>()` concurrently before initialization completes, both may execute the body, double-constructing the pool and double-starting the stress thread (triggering SEV-2 above). C++ guarantees thread-safe *initialization* of block-scope statics, but subsequent reads/writes to an already-initialized `static bool` are ordinary data accesses.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Two threads enter `g_term_pool<true>()`, both read `initialized == false`, both execute the body.
  - Most likely trigger: early concurrent access to the term pool during library startup with `EnableGCStressThread = true`.
  - This is mitigated in practice because `g_term_pool<true>()` is typically called only during static initialization (single-threaded), but the invariant is undocumented.
  - Test plan: Enable `EnableGCStressThread`, spawn two threads that both call `g_term_pool<true>()` simultaneously, run under TSan.
  - Run command: `cmake --build build-tsan && ctest --test-dir build-tsan -R <new_test> --output-on-failure`
- **Efficiency impact:** None
- **C++20-specific note:** Replace with `static std::once_flag` + `std::call_once`, or use a properly guarded `static` object whose constructor does the work (leveraging C++ static init guarantees).
- **Suggested fix:** Use `std::call_once` or move initialization into a static local's constructor to get the language guarantee.

### [SEV-3] Missing `\pre`/`\post` contract annotations
- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/gc_stress_thread.h:22–25`
- **Why this is a problem:** Per project convention (`.github/instructions/cpp.instructions.md`), new/changed public APIs require explicit `\pre`/`\post`. `start_gc_stress_thread()` has an implicit precondition (must not be called while a stress thread is already running) that is not documented.
- **Evidence type:** static-proof (convention violation)
- **Suggested fix:**
  ```cpp
  /// \brief Spawns a dedicated background thread that continuously triggers GC.
  /// \pre No GC stress thread is currently running.
  /// \post If EnableGCStressThread && GlobalThreadSafe, a background thread is active.
  void start_gc_stress_thread();
  
  /// \brief Stop the background GC stress thread if it is running.
  /// \post No GC stress thread is running.
  void stop_gc_stress_thread();
  ```

## Proposed Tests
- [ ] `test_gc_stress_double_start`: With `EnableGCStressThread = true`, call `start_gc_stress_thread()` twice without stopping. Currently causes `std::terminate`; after fix, should be idempotent. Signal: process survival vs. abort.
- [ ] `test_gc_stress_start_stop_cycle`: Start, stop, start, stop. Verify no TSan findings and no crash. Run under TSan.
- [ ] `test_concurrent_lazy_pool_init` (for call site): Two threads calling `g_term_pool<true>()` simultaneously under TSan. Signal: TSan data-race report on `initialized`.

## Open Questions
1. Is there a deliberate design assumption that `g_term_pool<true>()` is only ever called from a single thread (during static init)? If so, should it be enforced or at least asserted?
2. Has switching to `std::jthread` + `std::stop_token` been considered? It would eliminate the manual atomic flag, the `atexit` registration, and the double-start `std::terminate` risk in one change.
3. The `std::atexit` registration means `stop_gc_stress_thread()` runs during static destruction. Is there a guaranteed ordering that the thread-local `thread_aterm_pool` for the stress thread is destroyed *before* the global pool storage is reused? (Currently safe because the pool uses placement-new into a byte array that persists, but this is a fragile implicit invariant.)


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         2m 30s
Total session time:     2m 50s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         119.6k in, 5.1k out, 75.6k cached (Est. 1 Premium request)
 claude-haiku-4.5        96.3k in, 4.3k out, 77.4k cached (Est. 0 Premium requests)

```

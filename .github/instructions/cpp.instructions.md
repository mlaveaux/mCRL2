---
applyTo: "{libraries,tools,tests}/**/*.{h,hpp,cpp}"
description: "Use for C++ implementation and review work in mCRL2 libraries, tools, and tests; enforces API contracts, C++20/23 policy, UB- and data-race-freedom, sanitizers, clang-format, and CMake/test expectations."
---

# mCRL2 C++ File Instructions

## What to optimize for
mCRL2 is a model checker: an incorrect result is the worst possible defect.
- Correctness first, then efficiency, then convenience.
- Clean, minimal public APIs with explicit contracts.
- Small, reviewable patches over broad code movement.

## API contracts (pre- and postconditions)
- Document every non-trivial public function with Doxygen comments: `/// \brief`, `\param`, `\returns`, and explicit `\pre` and `\post` clauses where behavior depends on them.
- Enforce cheap preconditions with `assert(...)`.
- Guard expensive validation with `#ifndef MCRL2_NO_SOUNDNESS_CHECKS` (see `libraries/core/include/mcrl2/core/detail/soundness_checks.h`), controlled by `MCRL2_ENABLE_SOUNDNESS_CHECKS`.
- Mark pure queries and must-use results `[[nodiscard]]`; use `noexcept` where guaranteed.

## Language and safety policy
- Baseline is C++20 (`CMAKE_CXX_STANDARD 20`). C++23 features may be used as soon as all minimum supported toolchains accept them (GCC 11, Clang 16, AppleClang 14, MSVC 19.31 / VS 2022 17.1); the project raises these requirements as features become widely available. Verify availability before use.
- Prefer modern constructs that bring safety close to Rust: RAII and ownership types over raw `new`/`delete`, `std::span`/`std::string_view` (mind dangling), `constexpr`, scoped enums, `template<...>` without space (modern style).
- SFINAE must always be avoided and replaced by concepts/`requires` clauses if possible: never introduce new `std::enable_if`-style constraints, and migrate existing SFINAE when touching code that uses it (prefer named concepts, e.g. in `mcrl2/data/concepts.h`).
- Use Clang lifetime analysis where APIs return references or views tied to argument lifetimes: annotate with `[[clang::lifetimebound]]` (guard via `__has_cpp_attribute` for portability) and keep `-Wdangling*` diagnostics clean.
- No undefined behavior: no signed overflow, out-of-bounds access, invalid casts, dangling references, uninitialized reads, or use of invalidated iterators. When in doubt, verify with UBSan.

## Concurrency
- All code must be data-race free. Multithreading is on by default (`MCRL2_ENABLE_MULTITHREADING`).
- Protect shared mutable state with the existing wrappers in `mcrl2/utilities/` (`mutex.h`, `shared_mutex.h`, thread-safe `indexed_set`, `thread_local.h`) so single-threaded builds compile the synchronisation away; do not add ad-hoc primitives.
- Document the thread-safety of every new public class or function in its doc comment.
- Concurrent code must pass ThreadSanitizer (`MCRL2_ENABLE_THREADSANITIZER=ON`); known-benign races are only suppressed via `cmake/thread_sanitizer.suppress`.

## Efficiency
- Justify algorithmic complexity on hot paths; avoid accidental O(n²) or worse.
- Avoid needless allocations and copies; use moves, `reserve`, and the existing shared-term machinery (atermpp) instead of ad-hoc caching.

## Style and formatting
- Format all new and touched code with `clang-format` using the repository `.clang-format`; do not reformat unrelated code.
- Follow local naming conventions in the surrounding code.
- Keep headers minimal; avoid unnecessary include growth.

## Guardrails
- Do not edit `3rd-party/` unless explicitly requested.
- Avoid changing public APIs unless the task explicitly requires it.
- Reuse existing abstractions in `libraries/*` before adding new ones; prefer deterministic behavior and explicit error handling.

## Validation expectations
- Fix bugs test-first: add a regression test that demonstrates the failure, observe it fail, then fix, then observe it pass (see the `mcrl2-add-random-or-regression-test` skill).
- Ensure the project still configures and builds with CMake; run relevant `ctest` subsets for impacted components.
- Run sanitizer builds for memory- or concurrency-relevant changes (see the `mcrl2-sanitizer-validation` skill): code must pass ASan+UBSan+LSan and TSan cleanly. Sanitizers and tests are the verification method for this codebase (the code itself is not formally model checked).
- If tests are unavailable for a changed behavior, document the gap.

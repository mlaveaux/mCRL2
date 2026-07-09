---
name: mcrl2-cpp20-modernization
description: "Active C++20 modernization priorities and preferred patterns for mCRL2. Use when writing new code or refactoring existing code to decide between legacy and modern idioms (SFINAE vs concepts, manual comparisons vs spaceship, pointer+size vs span, manual iterators vs ranges). Keywords: C++20, concepts, requires, enable_if, SFINAE, operator<=>, spaceship, std::span, std::ranges, modernization."
---

# C++20 Modernization

The codebase is compiled with `-std=c++20` and is actively modernizing. C++23 features are acceptable once all minimum supported compilers accept them (GCC 11, Clang 16, AppleClang 14, MSVC 19.31).

Current priorities:

1. **Concepts** — SFINAE must always be avoided and replaced by concepts/`requires` clauses whenever possible. Never introduce new SFINAE (`std::enable_if`, `std::void_t` tricks, tag dispatch for constraint purposes); when touching code that uses it, migrate it. The pattern `std::is_base_of_v<atermpp::aterm, T>` (typically inside `enable_if` guards) appears in ~60 headers and is being replaced with named concepts (e.g., `IsSubstitution` in `libraries/data/include/mcrl2/data/concepts.h`).
2. **Spaceship operator** — replacing manual comparison operators with `operator<=>`.
3. **`std::span`** — replacing `(const T*, size_t)` pairs and `const std::vector<T>&` where ownership isn't needed (watch lifetimes).
4. **`std::ranges`** — replacing `std::find_if(v.begin(), v.end(), pred)` with `std::ranges::find_if(v, pred)`, etc.

Also preferred: RAII and smart pointers, `[[nodiscard]]`, and `[[clang::lifetimebound]]` guarded with `__has_cpp_attribute`.

When writing new code or modifying existing code, prefer the C++20 approach. Update `scripts/code_generation/` templates when modernizing generated code so regeneration doesn't regress.

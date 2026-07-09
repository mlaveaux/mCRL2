---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/function_symbol.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T12:28:38.057481+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/function_symbol.h` (clean file review, not a diff)
- **Risk level:** Low — this is stable foundational infrastructure with sound design
- **Overall verdict:** Acceptable with minor modernization opportunities. No confirmed defects found. Several C++20 idiom improvements and minor robustness concerns noted.

## Findings (ordered by severity)

### [SEV-3] Missing `[[nodiscard]]` on pure query methods
- **Location:** function_symbol.h:51, 58, 65
- **Why this is a problem:** `defined()`, `name()`, and `arity()` are pure observers whose return values have no side effects. Discarding them is always a bug. Missing `[[nodiscard]]` allows silent misuse (e.g. `f.defined();` as a statement instead of `assert(f.defined())`).
- **Evidence type:** static-proof
- **Evidence:** These functions have no side effects; discarding the result is meaningless in every context.
- **C++20-specific note:** `[[nodiscard]]` is baseline C++17 and should be applied to all pure query functions per project style expectations.
- **Suggested fix:** Add `[[nodiscard]]` to `defined()`, `name()`, `arity()`, and all comparison operators.

---

### [SEV-3] Six manual comparison operators instead of `operator<=>`
- **Location:** function_symbol.h:73–117
- **Why this is a problem:** The underlying `shared_reference` already provides `operator<=>` returning `std::weak_ordering` (shared_reference.h:188). The six manually written relational operators (`==`, `!=`, `<`, `>`, `<=`, `>=`) are redundant boilerplate and a maintenance burden. They could be replaced by a defaulted `operator==` and a single delegating `operator<=>`.
- **Evidence type:** static-proof
- **Evidence:** `shared_reference<T>` declares `operator<=>` (shared_reference.h:188) and `operator==` (shared_reference.h:182). With C++20, declaring `auto operator<=>(const function_symbol&) const = default;` or delegating to `m_function_symbol <=> f.m_function_symbol` plus `friend bool operator==(…) = default;` synthesizes all six operators.
- **Efficiency impact:** None (both forms compile to identical code).
- **C++20-specific note:** Spaceship operator is baseline C++20 and supported by all minimum toolchains.
- **Suggested fix:**
  ```cpp
  bool operator==(const function_symbol& f) const = default;
  std::weak_ordering operator<=>(const function_symbol& f) const = default;
  ```
  Or if defaulting doesn't work due to the `detail::_function_symbol::ref` type:
  ```cpp
  [[nodiscard]] bool operator==(const function_symbol& f) const { return m_function_symbol == f.m_function_symbol; }
  [[nodiscard]] std::weak_ordering operator<=>(const function_symbol& f) const { return m_function_symbol <=> f.m_function_symbol; }
  ```

---

### [SEV-3] `name()` and `arity()` not marked `noexcept`
- **Location:** function_symbol.h:58–68
- **Why this is a problem:** The underlying `_function_symbol::name()` and `_function_symbol::arity()` are both `noexcept` (detail/function_symbol.h:39,44). The outer wrappers dereference a pointer (assert-guarded, no-throw) and call noexcept functions, so they are de facto noexcept. Failing to declare them as such prevents the compiler from optimizing call sites and breaks `noexcept` propagation in generic contexts (e.g. `std::is_nothrow_invocable`).
- **Evidence type:** static-proof
- **Evidence:** `shared_reference::operator->()` is `noexcept` (shared_reference.h:193). Called methods are `noexcept`. The entire path is non-throwing.
- **Suggested fix:** Add `noexcept` to `name()` and `arity()`.

---

### [SEV-3] `defined()` not marked `noexcept`
- **Location:** function_symbol.h:51
- **Why this is a problem:** `shared_reference::defined()` (shared_reference.h:176) is not marked `noexcept` in the underlying class, but it's a simple null check with no possible throw. The `function_symbol::defined()` wrapper should still be `noexcept` for documentation value, and ideally the underlying `shared_reference::defined()` should be fixed too.
- **Evidence type:** static-proof
- **Suggested fix:** Add `noexcept` to `function_symbol::defined()`.

---

### [SEV-3] Missing `\pre` documentation on `name()` and `arity()`
- **Location:** function_symbol.h:56–68
- **Why this is a problem:** Both methods dereference the internal shared pointer via `operator->()`, which has `assert(defined())`. The precondition that the function symbol must be defined (not default-constructed with null reference) is not documented. Callers must know this invariant.
- **Evidence type:** plausible-issue
- **Evidence:** If a default-constructed `function_symbol` is queried with `.name()` or `.arity()`, the `assert` in `shared_reference::operator->()` fires in debug mode and is UB (null dereference) in release mode.
- **Suggested fix:** Add `/// \pre defined()` to both `name()` and `arity()`.

---

### [SEV-3] `global_function_symbol` inherits publicly without virtual destructor
- **Location:** function_symbol.h:136–141
- **Why this is a problem:** `function_symbol` does not have a virtual destructor, yet `global_function_symbol` inherits from it publicly. Deleting a `global_function_symbol` through a `function_symbol*` pointer would be undefined behavior. However, `global_function_symbol` adds no data members, so there is no object slicing or resource leak risk in practice — it exists only to select a different constructor path.
- **Evidence type:** plausible-issue (low likelihood in practice)
- **Evidence:** No polymorphic deletion is expected in this codebase for these types; they are value-like. The risk is purely theoretical for future misuse.
- **Suggested fix:** Either mark `function_symbol` destructor as `protected` in a future refactor, or document that polymorphic deletion is not supported. Low priority.

## Proposed Tests
- [ ] **test_default_function_symbol_accessors:** Call `name()` and `arity()` on a default-constructed `function_symbol` under ASan to confirm the assertion fires (or document expected behavior). Run: `ctest --test-dir build -R function_symbol --output-on-failure`
- [ ] **test_spaceship_ordering_consistency:** If spaceship is adopted, verify that `<`, `>`, `<=`, `>=` produce results consistent with the current implementation for a set of function symbols with varied indices.

## Open Questions
1. Is `global_function_symbol` still necessary as a separate class, or could its registration behavior be moved to a factory function? The public inheritance without virtual destructor is a code smell, even if harmless today.
2. The `extern` globals `g_as_int`, `g_as_list`, `g_as_empty_list` — are these initialized before first use in all TU orderings, or is there a static-initialization-order fiasco risk? (Not within scope of this header alone, but worth investigating.)
3. The `friend struct std::hash<function_symbol>` is declared but the `std::hash` specialization is not visible in this header — is it always included transitively, or could a TU include only this header and get an incomplete `hash`?

## Copilot Raw Stdout

```text
● Read libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h
  └ 62 lines read

● Check clang-format config
  $ cd /workspaces/mCRL2 && cat .clang-format 2>/dev/null | head -30
  └ 31 lines...

● Read libraries/utilities/include/mcrl2/utilities/shared_reference.h
  └ 247 lines read

● Find function_symbol_generator class
  $ cd /workspaces/mCRL2 && grep -n "function_symbol_generator"
  libraries/atermpp/include/mcrl2/atermpp/*.h | head -20
  └ 15 lines...

● Check noexcept usage
  $ cd /workspaces/mCRL2 && grep -n "noexcept"
  libraries/atermpp/include/mcrl2/atermpp/function_symbol.h
  └ 6 lines...

● Find hash specialization
  $ cd /workspaces/mCRL2 && grep -n "hash.*function_symbol"
  libraries/atermpp/include/mcrl2/atermpp/*.h libraries/atermpp/source/*.cpp 2>/dev/null | head -20
  └ 2 lines...

● Check GlobalThreadSafe definition
  $ cd /workspaces/mCRL2 && grep -rn "GlobalThreadSafe"
  libraries/atermpp/include/mcrl2/atermpp/detail/aterm_configuration.h 2>/dev/null; grep -rn
  "GlobalThreadSafe" libraries/utilities/include/mcrl2/utilities/ 2>/dev/null | head -10
  └ 11 lines...

● Check accessor usage patterns
  $ cd /workspaces/mCRL2 && grep -n "name()\|arity()"
  libraries/atermpp/include/mcrl2/atermpp/function_symbol.h | head -10
  └ 5 lines...

● Check nodiscard usage
  $ cd /workspaces/mCRL2 && grep -n "nodiscard"
  libraries/atermpp/include/mcrl2/atermpp/function_symbol.h
  └ 1 line...

● Check for spaceship operator
  $ cd /workspaces/mCRL2 && grep -rn "spaceship\|operator<=>"
  libraries/atermpp/include/mcrl2/atermpp/function_symbol.h
  └ 1 line...

● Read libraries/utilities/include/mcrl2/utilities/configuration.h
  └ 25 lines read

● Check if name/arity called on default-constructed
  $ cd /workspaces/mCRL2 && grep -rn "function_symbol.*name()\|function_symbol.*arity()"
  libraries/atermpp/source/ tools/ | grep -v "Binary" | head -10
  └ 11 lines...

● Check defined() usage
  $ cd /workspaces/mCRL2 && grep -n "defined()"
  libraries/atermpp/include/mcrl2/atermpp/function_symbol.h
  └ 3 lines...

● Check accessor null safety
  $ cd /workspaces/mCRL2 && grep -rn "\.name()\|\.arity()"
  libraries/atermpp/include/mcrl2/atermpp/function_symbol.h
  └ 1 line...

● Check pointer dereference patterns
  $ cd /workspaces/mCRL2 && grep -rn "m_function_symbol->"
  libraries/atermpp/include/mcrl2/atermpp/function_symbol.h
  └ 3 lines...

● Check decrement behavior on zero refcount
  $ cd /workspaces/mCRL2 && grep -rn "decrement_reference_count"
  libraries/utilities/include/mcrl2/utilities/shared_reference.h | head -5
  └ 5 lines...

Now I have enough context to produce the review.

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/function_symbol.h` (clean file review, not a diff)
- **Risk level:** Low — this is stable foundational infrastructure with sound design
- **Overall verdict:** Acceptable with minor modernization opportunities. No confirmed defects found. Several C++20 idiom improvements and minor robustness concerns noted.

## Findings (ordered by severity)

### [SEV-3] Missing `[[nodiscard]]` on pure query methods
- **Location:** function_symbol.h:51, 58, 65
- **Why this is a problem:** `defined()`, `name()`, and `arity()` are pure observers whose return values have no side effects. Discarding them is always a bug. Missing `[[nodiscard]]` allows silent misuse (e.g. `f.defined();` as a statement instead of `assert(f.defined())`).
- **Evidence type:** static-proof
- **Evidence:** These functions have no side effects; discarding the result is meaningless in every context.
- **C++20-specific note:** `[[nodiscard]]` is baseline C++17 and should be applied to all pure query functions per project style expectations.
- **Suggested fix:** Add `[[nodiscard]]` to `defined()`, `name()`, `arity()`, and all comparison operators.

---

### [SEV-3] Six manual comparison operators instead of `operator<=>`
- **Location:** function_symbol.h:73–117
- **Why this is a problem:** The underlying `shared_reference` already provides `operator<=>` returning `std::weak_ordering` (shared_reference.h:188). The six manually written relational operators (`==`, `!=`, `<`, `>`, `<=`, `>=`) are redundant boilerplate and a maintenance burden. They could be replaced by a defaulted `operator==` and a single delegating `operator<=>`.
- **Evidence type:** static-proof
- **Evidence:** `shared_reference<T>` declares `operator<=>` (shared_reference.h:188) and `operator==` (shared_reference.h:182). With C++20, declaring `auto operator<=>(const function_symbol&) const = default;` or delegating to `m_function_symbol <=> f.m_function_symbol` plus `friend bool operator==(…) = default;` synthesizes all six operators.
- **Efficiency impact:** None (both forms compile to identical code).
- **C++20-specific note:** Spaceship operator is baseline C++20 and supported by all minimum toolchains.
- **Suggested fix:**
  ```cpp
  bool operator==(const function_symbol& f) const = default;
  std::weak_ordering operator<=>(const function_symbol& f) const = default;
  ```
  Or if defaulting doesn't work due to the `detail::_function_symbol::ref` type:
  ```cpp
  [[nodiscard]] bool operator==(const function_symbol& f) const { return m_function_symbol == f.m_function_symbol; }
  [[nodiscard]] std::weak_ordering operator<=>(const function_symbol& f) const { return m_function_symbol <=> f.m_function_symbol; }
  ```

---

### [SEV-3] `name()` and `arity()` not marked `noexcept`
- **Location:** function_symbol.h:58–68
- **Why this is a problem:** The underlying `_function_symbol::name()` and `_function_symbol::arity()` are both `noexcept` (detail/function_symbol.h:39,44). The outer wrappers dereference a pointer (assert-guarded, no-throw) and call noexcept functions, so they are de facto noexcept. Failing to declare them as such prevents the compiler from optimizing call sites and breaks `noexcept` propagation in generic contexts (e.g. `std::is_nothrow_invocable`).
- **Evidence type:** static-proof
- **Evidence:** `shared_reference::operator->()` is `noexcept` (shared_reference.h:193). Called methods are `noexcept`. The entire path is non-throwing.
- **Suggested fix:** Add `noexcept` to `name()` and `arity()`.

---

### [SEV-3] `defined()` not marked `noexcept`
- **Location:** function_symbol.h:51
- **Why this is a problem:** `shared_reference::defined()` (shared_reference.h:176) is not marked `noexcept` in the underlying class, but it's a simple null check with no possible throw. The `function_symbol::defined()` wrapper should still be `noexcept` for documentation value, and ideally the underlying `shared_reference::defined()` should be fixed too.
- **Evidence type:** static-proof
- **Suggested fix:** Add `noexcept` to `function_symbol::defined()`.

---

### [SEV-3] Missing `\pre` documentation on `name()` and `arity()`
- **Location:** function_symbol.h:56–68
- **Why this is a problem:** Both methods dereference the internal shared pointer via `operator->()`, which has `assert(defined())`. The precondition that the function symbol must be defined (not default-constructed with null reference) is not documented. Callers must know this invariant.
- **Evidence type:** plausible-issue
- **Evidence:** If a default-constructed `function_symbol` is queried with `.name()` or `.arity()`, the `assert` in `shared_reference::operator->()` fires in debug mode and is UB (null dereference) in release mode.
- **Suggested fix:** Add `/// \pre defined()` to both `name()` and `arity()`.

---

### [SEV-3] `global_function_symbol` inherits publicly without virtual destructor
- **Location:** function_symbol.h:136–141
- **Why this is a problem:** `function_symbol` does not have a virtual destructor, yet `global_function_symbol` inherits from it publicly. Deleting a `global_function_symbol` through a `function_symbol*` pointer would be undefined behavior. However, `global_function_symbol` adds no data members, so there is no object slicing or resource leak risk in practice — it exists only to select a different constructor path.
- **Evidence type:** plausible-issue (low likelihood in practice)
- **Evidence:** No polymorphic deletion is expected in this codebase for these types; they are value-like. The risk is purely theoretical for future misuse.
- **Suggested fix:** Either mark `function_symbol` destructor as `protected` in a future refactor, or document that polymorphic deletion is not supported. Low priority.

## Proposed Tests
- [ ] **test_default_function_symbol_accessors:** Call `name()` and `arity()` on a default-constructed `function_symbol` under ASan to confirm the assertion fires (or document expected behavior). Run: `ctest --test-dir build -R function_symbol --output-on-failure`
- [ ] **test_spaceship_ordering_consistency:** If spaceship is adopted, verify that `<`, `>`, `<=`, `>=` produce results consistent with the current implementation for a set of function symbols with varied indices.

## Open Questions
1. Is `global_function_symbol` still necessary as a separate class, or could its registration behavior be moved to a factory function? The public inheritance without virtual destructor is a code smell, even if harmless today.
2. The `extern` globals `g_as_int`, `g_as_list`, `g_as_empty_list` — are these initialized before first use in all TU orderings, or is there a static-initialization-order fiasco risk? (Not within scope of this header alone, but worth investigating.)
3. The `friend struct std::hash<function_symbol>` is declared but the `std::hash` specialization is not visible in this header — is it always included transitively, or could a TU include only this header and get an incomplete `hash`?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         1m 23s
Total session time:     1m 36s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         348.4k in, 3.7k out, 322.8k cached (Est. 1 Premium request)

```

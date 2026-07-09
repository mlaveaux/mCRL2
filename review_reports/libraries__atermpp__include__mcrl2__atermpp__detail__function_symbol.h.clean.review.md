---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T10:13:30.554566+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h` — a foundational internal class storing function symbol (name, arity) pairs for the aterm library.
- **Risk level:** Low (no correctness bugs; efficiency and hygiene concerns)
- **Overall verdict:** The class is simple and correct for its purpose. It has one concrete efficiency deficiency (missing move-constructor overload wastes a `std::string` allocation on every new symbol creation), a weak hash combination in the companion file, and several API contract/hygiene gaps per project standards.

## Findings (ordered by severity)

### [SEV-2] Constructor forces unnecessary string copy even when caller provides rvalue

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h:34`
- **Why this is a problem:** The only constructor signature is `_function_symbol(const std::string& name, std::size_t arity)`. The pool (`function_symbol_pool::create`) calls `m_symbol_set.emplace(std::move(name), arity)` (at `source/function_symbol_pool.cpp:77`), forwarding an rvalue `std::string&&` which the transparent hasher/equals path passes through without consuming. When the element doesn't exist and `emplace_impl` constructs the `_function_symbol`, the rvalue binds to `const std::string&` and is **copied** into `m_name`. This wastes one heap allocation per new function symbol.
- **Evidence type:** static-proof
- **Evidence:**
  - The constructor parameter `const std::string&` cannot invoke the move constructor of `m_name`. The `const` qualifier on `m_name` itself does *not* prevent move-initialization — only the parameter type does.
  - Reproducer: any path that creates a previously-unseen function symbol (`function_symbol("some_long_name_exceeding_SSO", 2)`) performs an avoidable allocation.
  - Run command: `perf stat -e cache-misses -- <any mCRL2 tool processing a large spec>` would show excess allocations, but this is provable by inspection.
  - Observed result: `m_name` is copy-constructed from the const reference.
  - Expected result: with a `_function_symbol(std::string&& name, std::size_t arity)` overload, `m_name` would be move-constructed (zero allocation for the string buffer).
- **Efficiency impact:** One heap allocation + memcpy per unique function symbol creation. For large specifications with thousands of unique symbols, this is measurable.
- **C++20-specific note:** A modern idiom is to take by value (`std::string name`) and `std::move` into the member, which handles both lvalue and rvalue callers with at most one move. Alternatively, provide both `const std::string&` and `std::string&&` overloads.
- **Suggested fix:**
  ```cpp
  _function_symbol(std::string name, std::size_t arity) :
     m_arity(arity),
     m_name(std::move(name))
  {}
  ```
  Or add a second constructor:
  ```cpp
  _function_symbol(std::string&& name, std::size_t arity) :
     m_arity(arity),
     m_name(std::move(name))
  {}
  ```

---

### [SEV-3] Weak hash combination (XOR with small integer) in companion hash

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_hash.h:41`
- **Why this is a problem:** `return (h ^ f.arity());` — XOR with arity (typically 0–10) barely perturbs the hash of the name. Function symbols with the same name but different arities (legal in mCRL2 process algebra) will cluster in the same hash buckets.
- **Evidence type:** plausible-issue
- **Evidence:**
  - For two symbols `("f", 2)` and `("f", 3)`, the hashes differ only in bits 0–1, causing bucket collisions in power-of-2 tables.
  - Reproducer: create a specification with many overloaded action names (same name, varying arities) and measure `function_symbol_pool` collision rate via `EnableHashtableMetrics`.
  - Expected result: a hash_combine pattern (e.g., `h ^= arity + 0x9e3779b9 + (h << 6) + (h >> 2)`) distributes uniformly.
- **Efficiency impact:** Degrades hash table performance from O(1) to O(k) for k same-name symbols; minor in typical use, potentially noticeable for heavily overloaded specs.
- **C++20-specific note:** No language feature concern; this is algorithmic.
- **Suggested fix:** Use a proper hash combine:
  ```cpp
  std::size_t h = string_hasher(f.name());
  h ^= f.arity() + 0x9e3779b9 + (h << 6) + (h >> 2);
  return h;
  ```

---

### [SEV-3] Non-standard include guard naming

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h:13`
- **Why this is a problem:** Guard `DETAIL_FUNCTION_SYMBOL_H` lacks the project-wide `MCRL2_` prefix. Other headers in the same directory use `MCRL2_ATERMPP_DETAIL_*` (e.g., `MCRL2_ATERMPP_DETAIL_FUNCTION_SYMBOL_HASH_H_`). The short name risks collision with third-party or system headers.
- **Evidence type:** static-proof (comparison with sibling headers)
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Rename to `MCRL2_ATERMPP_DETAIL_FUNCTION_SYMBOL_H`.

---

### [SEV-3] Missing `[[nodiscard]]` on pure query accessors

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h:39,44`
- **Why this is a problem:** `name()` and `arity()` are pure observers; discarding their results is always a logic error. Per project C++ instructions, pure queries should be `[[nodiscard]]`.
- **Evidence type:** static-proof (code inspection against project policy in `.github/instructions/cpp.instructions.md`)
- **Efficiency impact:** None (compile-time only).
- **C++20-specific note:** `[[nodiscard]]` is standard since C++17 and required by project conventions.
- **Suggested fix:**
  ```cpp
  [[nodiscard]] const std::string& name() const noexcept { ... }
  [[nodiscard]] std::size_t arity() const noexcept { ... }
  [[nodiscard]] bool operator==(const _function_symbol& f) const noexcept { ... }
  ```

---

### [SEV-3] Missing Doxygen documentation on public API

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h:34,39,44,49`
- **Why this is a problem:** Per project requirements (`.github/instructions/cpp.instructions.md`), public functions need `\brief`, `\param`, and `\pre`/`\post` where behavior depends on them. The constructor's implicit precondition (name is non-empty? arity has no max?) is undocumented.
- **Evidence type:** static-proof (policy gap)
- **Efficiency impact:** None.
- **Suggested fix:** Add `/// \brief` on each member and `\pre` on the constructor if there are constraints on inputs.

---

## Proposed Tests

- [ ] **perf_function_symbol_creation:** Benchmark creating N unique function symbols via `function_symbol(std::string&&, arity)` path vs. `const std::string&` path, measuring heap allocations. Purpose: quantify SEV-2 savings. Signal: reduced `malloc` count after adding move overload.
- [ ] **hash_collision_stress:** Insert 100 function symbols with the same name and arities 0–99 into a `function_symbol_pool`, then query `EnableHashtableMetrics` collision stats. Signal: collision rate > 50% confirms SEV-3 hash weakness.

## Open Questions

1. Is there a deliberate design reason for the `const` qualifier on `m_name` and `m_arity`? If immutability is the goal, it's achieved by the `noncopyable` + reference-counted pointer pattern already — the `const` members additionally prevent the class from being movable, which precludes any future refactoring that might need element relocation.
2. The pool's `unordered_set` is configured with `Resize = false` (last template param). Who manages resizing? If the set grows unbounded without rehashing, collision chains degrade performance independent of the hash quality issue.

## Copilot Raw Stdout

```text
● Explore: Explore function_symbol context
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Read libraries/utilities/include/mcrl2/utilities/shared_reference.h lines 1-80
  └ 80 lines read

● Read libraries/utilities/include/mcrl2/utilities/shared_reference.h lines 80-180
  └ 100 lines read

● Grep "_function_symbol\(" (/workspaces/mCRL2/libraries/atermpp)
  └ 108 lines found

● Read libraries/utilities/include/mcrl2/utilities/noncopyable.h
  └ 34 lines read

● Grep "hash.*_function_symbol|_function_symbol.*hash"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ 7 lines found

● Grep "new _function_symbol" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ No matches found

● Grep "function_symbol_pool" in *function_symbol_pool*
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ 24 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_hash.h
  └ 96 lines read

● Grep "emplace|construct.*_function_symbol" in *function_symbol_pool*
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ No matches found

● Grep "_function_symbol\(.*name.*arity" in *function_symbol_pool*
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ No matches found

● Grep "construct<_function_symbol|allocate.*_function_symbol|_function_symbol\(std::move"
(/workspaces/mCRL2/libraries/atermpp)
  └ 42 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_configuration.h
  └ 47 lines read

● Grep "function_symbol_pool::create\b" in *function_symbol_pool*"
(/workspaces/mCRL2/libraries/atermpp)
  └ No matches found

● Grep "function_symbol_pool::create\(" (/workspaces/mCRL2/libraries/atermpp)
  └ 55 lines found

● Read: Explore agent (Explore function_symbol context)
  └ Completed

Now I have enough context to produce the review. Let me verify one more thing about the unordered_set used by the pool:

● Grep "m_symbol_set"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_pool.h)
  └ 18 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_pool.h lines 77-86
  └ 9 lines read

● Grep "class unordered_set" (/workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities)
  └ 26 lines found

● Grep "emplace\(" (/workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities/unordered_set.h)
  └ 21 lines found

● Grep "emplace\(" in *unordered_set_implementation*
(/workspaces/mCRL2/libraries/utilities/include/mcrl2/utilities)
  └ 75 lines found

● Read libraries/utilities/include/mcrl2/utilities/detail/unordered_set_implementation.h lines
  100-145
  └ 46 lines read

I now have comprehensive context to produce the review. Here's my analysis:

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h` — a foundational internal class storing function symbol (name, arity) pairs for the aterm library.
- **Risk level:** Low (no correctness bugs; efficiency and hygiene concerns)
- **Overall verdict:** The class is simple and correct for its purpose. It has one concrete efficiency deficiency (missing move-constructor overload wastes a `std::string` allocation on every new symbol creation), a weak hash combination in the companion file, and several API contract/hygiene gaps per project standards.

## Findings (ordered by severity)

### [SEV-2] Constructor forces unnecessary string copy even when caller provides rvalue

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h:34`
- **Why this is a problem:** The only constructor signature is `_function_symbol(const std::string& name, std::size_t arity)`. The pool (`function_symbol_pool::create`) calls `m_symbol_set.emplace(std::move(name), arity)` (at `source/function_symbol_pool.cpp:77`), forwarding an rvalue `std::string&&` which the transparent hasher/equals path passes through without consuming. When the element doesn't exist and `emplace_impl` constructs the `_function_symbol`, the rvalue binds to `const std::string&` and is **copied** into `m_name`. This wastes one heap allocation per new function symbol.
- **Evidence type:** static-proof
- **Evidence:**
  - The constructor parameter `const std::string&` cannot invoke the move constructor of `m_name`. The `const` qualifier on `m_name` itself does *not* prevent move-initialization — only the parameter type does.
  - Reproducer: any path that creates a previously-unseen function symbol (`function_symbol("some_long_name_exceeding_SSO", 2)`) performs an avoidable allocation.
  - Run command: `perf stat -e cache-misses -- <any mCRL2 tool processing a large spec>` would show excess allocations, but this is provable by inspection.
  - Observed result: `m_name` is copy-constructed from the const reference.
  - Expected result: with a `_function_symbol(std::string&& name, std::size_t arity)` overload, `m_name` would be move-constructed (zero allocation for the string buffer).
- **Efficiency impact:** One heap allocation + memcpy per unique function symbol creation. For large specifications with thousands of unique symbols, this is measurable.
- **C++20-specific note:** A modern idiom is to take by value (`std::string name`) and `std::move` into the member, which handles both lvalue and rvalue callers with at most one move. Alternatively, provide both `const std::string&` and `std::string&&` overloads.
- **Suggested fix:**
  ```cpp
  _function_symbol(std::string name, std::size_t arity) :
     m_arity(arity),
     m_name(std::move(name))
  {}
  ```
  Or add a second constructor:
  ```cpp
  _function_symbol(std::string&& name, std::size_t arity) :
     m_arity(arity),
     m_name(std::move(name))
  {}
  ```

---

### [SEV-3] Weak hash combination (XOR with small integer) in companion hash

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_hash.h:41`
- **Why this is a problem:** `return (h ^ f.arity());` — XOR with arity (typically 0–10) barely perturbs the hash of the name. Function symbols with the same name but different arities (legal in mCRL2 process algebra) will cluster in the same hash buckets.
- **Evidence type:** plausible-issue
- **Evidence:**
  - For two symbols `("f", 2)` and `("f", 3)`, the hashes differ only in bits 0–1, causing bucket collisions in power-of-2 tables.
  - Reproducer: create a specification with many overloaded action names (same name, varying arities) and measure `function_symbol_pool` collision rate via `EnableHashtableMetrics`.
  - Expected result: a hash_combine pattern (e.g., `h ^= arity + 0x9e3779b9 + (h << 6) + (h >> 2)`) distributes uniformly.
- **Efficiency impact:** Degrades hash table performance from O(1) to O(k) for k same-name symbols; minor in typical use, potentially noticeable for heavily overloaded specs.
- **C++20-specific note:** No language feature concern; this is algorithmic.
- **Suggested fix:** Use a proper hash combine:
  ```cpp
  std::size_t h = string_hasher(f.name());
  h ^= f.arity() + 0x9e3779b9 + (h << 6) + (h >> 2);
  return h;
  ```

---

### [SEV-3] Non-standard include guard naming

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h:13`
- **Why this is a problem:** Guard `DETAIL_FUNCTION_SYMBOL_H` lacks the project-wide `MCRL2_` prefix. Other headers in the same directory use `MCRL2_ATERMPP_DETAIL_*` (e.g., `MCRL2_ATERMPP_DETAIL_FUNCTION_SYMBOL_HASH_H_`). The short name risks collision with third-party or system headers.
- **Evidence type:** static-proof (comparison with sibling headers)
- **Efficiency impact:** None.
- **C++20-specific note:** N/A.
- **Suggested fix:** Rename to `MCRL2_ATERMPP_DETAIL_FUNCTION_SYMBOL_H`.

---

### [SEV-3] Missing `[[nodiscard]]` on pure query accessors

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h:39,44`
- **Why this is a problem:** `name()` and `arity()` are pure observers; discarding their results is always a logic error. Per project C++ instructions, pure queries should be `[[nodiscard]]`.
- **Evidence type:** static-proof (code inspection against project policy in `.github/instructions/cpp.instructions.md`)
- **Efficiency impact:** None (compile-time only).
- **C++20-specific note:** `[[nodiscard]]` is standard since C++17 and required by project conventions.
- **Suggested fix:**
  ```cpp
  [[nodiscard]] const std::string& name() const noexcept { ... }
  [[nodiscard]] std::size_t arity() const noexcept { ... }
  [[nodiscard]] bool operator==(const _function_symbol& f) const noexcept { ... }
  ```

---

### [SEV-3] Missing Doxygen documentation on public API

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h:34,39,44,49`
- **Why this is a problem:** Per project requirements (`.github/instructions/cpp.instructions.md`), public functions need `\brief`, `\param`, and `\pre`/`\post` where behavior depends on them. The constructor's implicit precondition (name is non-empty? arity has no max?) is undocumented.
- **Evidence type:** static-proof (policy gap)
- **Efficiency impact:** None.
- **Suggested fix:** Add `/// \brief` on each member and `\pre` on the constructor if there are constraints on inputs.

---

## Proposed Tests

- [ ] **perf_function_symbol_creation:** Benchmark creating N unique function symbols via `function_symbol(std::string&&, arity)` path vs. `const std::string&` path, measuring heap allocations. Purpose: quantify SEV-2 savings. Signal: reduced `malloc` count after adding move overload.
- [ ] **hash_collision_stress:** Insert 100 function symbols with the same name and arities 0–99 into a `function_symbol_pool`, then query `EnableHashtableMetrics` collision stats. Signal: collision rate > 50% confirms SEV-3 hash weakness.

## Open Questions

1. Is there a deliberate design reason for the `const` qualifier on `m_name` and `m_arity`? If immutability is the goal, it's achieved by the `noncopyable` + reference-counted pointer pattern already — the `const` members additionally prevent the class from being movable, which precludes any future refactoring that might need element relocation.
2. The pool's `unordered_set` is configured with `Resize = false` (last template param). Who manages resizing? If the set grows unbounded without rehashing, collision chains degrade performance independent of the hash quality issue.


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         4m 19s
Total session time:     4m 4s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         754.0k in, 9.7k out, 703.0k cached (Est. 1 Premium request)
 claude-haiku-4.5        219.4k in, 5.0k out, 195.8k cached (Est. 0 Premium requests)

```

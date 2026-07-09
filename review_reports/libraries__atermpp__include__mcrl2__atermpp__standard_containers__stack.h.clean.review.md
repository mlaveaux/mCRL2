---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/standard_containers/stack.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T13:14:38.457969+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/stack.h` (full file, 168 lines)
- **Risk level:** Low-to-Medium — no UB or crash risk in current known usage, but API contract deviations and incorrect documentation
- **Overall verdict:** Functional but sloppy — the file is clearly copy-pasted from `deque.h` with incomplete adaptation. Multiple documentation errors, a standard-deviating `emplace()` return type, an unconditionally `noexcept` swap, and unnecessary iterator exposure. No confirmed runtime bugs in current usage (local stack in `aterm_io_binary.cpp`), but the interface quality is below what a public header should offer.

## Findings (ordered by severity)

### [SEV-2] `emplace()` returns `void` instead of `reference`
- **Location:** `stack.h:145-149`
- **Why this is a problem:** Since C++17, `std::stack::emplace` returns `decltype(auto)` (resolving to `reference`). The underlying `atermpp::deque::emplace_back()` already returns `reference`. This implementation discards that value, breaking the standard interface contract. Code like `auto& ref = stack.emplace(args...)` will fail to compile.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer: Any code attempting `auto& x = my_stack.emplace(value);`
  - Run command: Compile with any C++20 compiler
  - Observed result: compilation failure (void not convertible to reference)
  - Expected result: returns a reference to the emplaced element
- **Efficiency impact:** None
- **C++20-specific note:** C++17 changed the return type; this class fails to follow the updated contract.
- **Suggested fix:**
  ```cpp
  template< class... Args >
  decltype(auto) emplace( Args&&... args )
  {
    return m_container.emplace_back(std::forward<Args>(args)...);
  }
  ```

### [SEV-2] Unconditional `noexcept` on `swap` may cause `std::terminate`
- **Location:** `stack.h:156`
- **Why this is a problem:** `swap` is declared unconditionally `noexcept`, but `Container::swap` (or ADL `swap`) is not guaranteed to be `noexcept`. If the container's swap throws, `std::terminate` is called. While `atermpp::deque` likely has a non-throwing swap (inheriting from `std::deque`), a user-provided `Container` template argument could throw.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Instantiate `atermpp::stack<T, CustomContainer>` where `CustomContainer::swap` may throw
  - Test plan: Provide a container whose swap throws; verify `std::terminate` is triggered
  - Run command: `ctest --test-dir build -R stack_swap_noexcept --output-on-failure` (proposed test)
  - Expected result: conditional noexcept prevents `std::terminate`
- **Efficiency impact:** None
- **C++20-specific note:** Proper idiom is `noexcept(noexcept(m_container.swap(other.m_container)))` or use `std::is_nothrow_swappable_v`.
- **Suggested fix:**
  ```cpp
  void swap( stack& other ) noexcept(std::is_nothrow_swappable_v<Container>)
  {
    using std::swap; swap(m_container, other.m_container);
  }
  ```

### [SEV-2] Copy-paste documentation describes wrong class
- **Location:** `stack.h:10-16` (`\file` directive), `stack.h:30` (class `\brief`)
- **Why this is a problem:** The Doxygen `\file` says `mcrl2/data/standard_containers/deque.h` (wrong path AND wrong class name), and the `\brief` says "A deque class in which aterms can be stored." This makes the documentation actively misleading for anyone browsing API docs or using tooling that relies on Doxygen metadata.
- **Evidence type:** static-proof
- **Evidence:**
  - Observed result: `\file mcrl2/data/standard_containers/deque.h` and `\brief A deque class`
  - Expected result: `\file mcrl2/atermpp/standard_containers/stack.h` and `\brief A stack class`
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Replace all "deque" references with "stack" and fix the `\file` path.

### [SEV-3] Exposing `iterator`/`const_iterator` breaks stack adaptor abstraction
- **Location:** `stack.h:46-47`
- **Why this is a problem:** `std::stack` deliberately does not expose iterators — it's a LIFO adaptor. Publishing them invites users to traverse or randomly access the container, defeating the purpose of the restricted interface. The `mark()` method (the only operation needing iteration) delegates to the underlying container and does not need these typedefs at the stack level.
- **Evidence type:** static-proof
- **Evidence:**
  - `std::stack` exposes no `iterator` typedef per [container.adaptors.stack]
  - These typedefs serve no purpose in the current implementation
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Remove `iterator` and `const_iterator` type aliases.

### [SEV-3] Missing `[[nodiscard]]` on query methods
- **Location:** `stack.h:115-132` (`top()`, `empty()`, `size()`)
- **Why this is a problem:** Per the project instructions, pure queries should be `[[nodiscard]]`. Discarding the result of `empty()` or `top()` is almost certainly a bug.
- **Evidence type:** plausible-issue (missed-bug-prevention)
- **Evidence:**
  - Pattern: `stack.empty();` (forgot the `if`) silently discards the result
  - Note: no existing standard_containers in mCRL2 use `[[nodiscard]]`, so this matches current practice, but does not match stated project policy
- **Efficiency impact:** None
- **C++20-specific note:** `[[nodiscard]]` is zero-cost and prevents logic errors at compile time.
- **Suggested fix:** Add `[[nodiscard]]` to `top()`, `empty()`, and `size()`.

### [SEV-3] Unnecessary includes
- **Location:** `stack.h:22,24`
- **Why this is a problem:** `thread_aterm_pool.h` and `shared_mutex.h` are included but nothing from them is directly used. They are already transitively included by `deque.h`. Unnecessary includes increase compilation time and create spurious coupling.
- **Evidence type:** static-proof
- **Efficiency impact:** Marginal compile-time overhead
- **Suggested fix:** Remove both includes; they come transitively through `deque.h`.

### [SEV-3] Unconstrained template parameters (missed C++20 opportunity)
- **Location:** `stack.h:63,69,75,101`
- **Why this is a problem:** `InputIt`, `InputIterator`, and `Alloc` parameters are completely unconstrained. The codebase already uses concepts in `aterm_container.h`. Concepts would provide clear error messages on misuse.
- **Evidence type:** plausible-issue (poor diagnostics)
- **C++20-specific note:** Use `std::input_iterator` concept for `InputIt` and a suitable allocator concept for `Alloc`.
- **Suggested fix:** Add `requires` clauses, e.g.:
  ```cpp
  template<std::input_iterator InputIt>
  stack(InputIt first, InputIt last) : m_container(first, last) {}
  ```

### [SEV-3] Inconsistent template parameter naming
- **Location:** `stack.h:63` (`InputIt`) vs `stack.h:69` (`InputIterator`)
- **Why this is a problem:** Two constructors performing the same conceptual operation use different naming for their iterator parameter. This makes the code appear as two unrelated overloads rather than allocator-aware variants of the same operation.
- **Suggested fix:** Standardize on `InputIt` (matching STL convention).

### [SEV-3] Missing comparison operators
- **Location:** entire class
- **Why this is a problem:** `std::stack` provides `operator==`, `!=`, `<`, `<=`, `>`, `>=`, and (C++20) `<=>`. This class provides none. Users comparing stacks will get a compilation error.
- **Evidence type:** static-proof
- **C++20-specific note:** A single `friend auto operator<=>(const stack&, const stack&) = default;` would provide all six (if the container supports it), or delegate to `m_container` comparisons.
- **Suggested fix:** Add `operator==` and `operator<=>` delegating to `m_container`.

### [SEV-3] `mark()` parameter uses raw `std::stack<...>` instead of `term_mark_stack` typedef
- **Location:** `stack.h:161`
- **Why this is a problem:** `aterm_container.h` defines `using term_mark_stack = std::stack<std::reference_wrapper<detail::_aterm>>`. Using the raw type instead of the typedef reduces readability and risks divergence if the typedef changes.
- **Suggested fix:** `void mark(detail::term_mark_stack& todo) const`

## Proposed Tests
- [ ] `stack_emplace_returns_reference`: Verify `emplace()` returns a reference to the inserted element. Input: `stack.emplace(value)` assigned to `auto&`. Expected: compiles and reference is valid.
- [ ] `stack_swap_noexcept_conditional`: Verify swap with a potentially-throwing container doesn't unconditionally promise `noexcept`. Input: custom container with `swap() noexcept(false)`. Expected: `static_assert(!noexcept(stack.swap(other)))`.
- [ ] `stack_comparison_operators`: Verify equality and ordering. Input: two stacks with same/different elements. Expected: `==` and `<=>` work correctly.

## Open Questions
1. Is the exposure of `iterator`/`const_iterator` intentional for some usage not visible in the repository? If so, it should be documented.
2. Should `top()` and `pop()` assert `!empty()` as a cheap precondition, or is the cost of the check unacceptable on a hot path (the binary I/O serialization loop)?
3. Are there plans to support user-provided `Container` types beyond `atermpp::deque`? If not, the unconditional `noexcept` on swap is benign in practice (but still technically incorrect).

## Copilot Raw Stdout

```text
● Explore: Explore atermpp container patterns
  ● Check emplace return types in deque.h
    $ grep -A 3 "template.*emplace.*Args"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/deque.h
    └ 1 line...
  ● Check all emplace methods in deque
    $ grep -B 1 -A 3 "emplace"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/deque.h | head -30
    └ 27 lines...

● Read libraries/atermpp/include/mcrl2/atermpp/standard_containers/stack.h
  └ 168 lines read

● Explore: Check stack.h usage and threading
  ● Check stack.h top() method and const_reference
    $ grep -A 5 "class deque\|const_reference\|reference top()"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/stack.h
    └ 19 lines...
  ● Check if back() is protected in deque or inherited from std::deque
    $ grep -n "protected\|back()"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/standard_containers/deque.h
    └ 6 lines...

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/standard_containers/stack.h` (full file, 168 lines)
- **Risk level:** Low-to-Medium — no UB or crash risk in current known usage, but API contract deviations and incorrect documentation
- **Overall verdict:** Functional but sloppy — the file is clearly copy-pasted from `deque.h` with incomplete adaptation. Multiple documentation errors, a standard-deviating `emplace()` return type, an unconditionally `noexcept` swap, and unnecessary iterator exposure. No confirmed runtime bugs in current usage (local stack in `aterm_io_binary.cpp`), but the interface quality is below what a public header should offer.

## Findings (ordered by severity)

### [SEV-2] `emplace()` returns `void` instead of `reference`
- **Location:** `stack.h:145-149`
- **Why this is a problem:** Since C++17, `std::stack::emplace` returns `decltype(auto)` (resolving to `reference`). The underlying `atermpp::deque::emplace_back()` already returns `reference`. This implementation discards that value, breaking the standard interface contract. Code like `auto& ref = stack.emplace(args...)` will fail to compile.
- **Evidence type:** static-proof
- **Evidence:**
  - Reproducer: Any code attempting `auto& x = my_stack.emplace(value);`
  - Run command: Compile with any C++20 compiler
  - Observed result: compilation failure (void not convertible to reference)
  - Expected result: returns a reference to the emplaced element
- **Efficiency impact:** None
- **C++20-specific note:** C++17 changed the return type; this class fails to follow the updated contract.
- **Suggested fix:**
  ```cpp
  template< class... Args >
  decltype(auto) emplace( Args&&... args )
  {
    return m_container.emplace_back(std::forward<Args>(args)...);
  }
  ```

### [SEV-2] Unconditional `noexcept` on `swap` may cause `std::terminate`
- **Location:** `stack.h:156`
- **Why this is a problem:** `swap` is declared unconditionally `noexcept`, but `Container::swap` (or ADL `swap`) is not guaranteed to be `noexcept`. If the container's swap throws, `std::terminate` is called. While `atermpp::deque` likely has a non-throwing swap (inheriting from `std::deque`), a user-provided `Container` template argument could throw.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: Instantiate `atermpp::stack<T, CustomContainer>` where `CustomContainer::swap` may throw
  - Test plan: Provide a container whose swap throws; verify `std::terminate` is triggered
  - Run command: `ctest --test-dir build -R stack_swap_noexcept --output-on-failure` (proposed test)
  - Expected result: conditional noexcept prevents `std::terminate`
- **Efficiency impact:** None
- **C++20-specific note:** Proper idiom is `noexcept(noexcept(m_container.swap(other.m_container)))` or use `std::is_nothrow_swappable_v`.
- **Suggested fix:**
  ```cpp
  void swap( stack& other ) noexcept(std::is_nothrow_swappable_v<Container>)
  {
    using std::swap; swap(m_container, other.m_container);
  }
  ```

### [SEV-2] Copy-paste documentation describes wrong class
- **Location:** `stack.h:10-16` (`\file` directive), `stack.h:30` (class `\brief`)
- **Why this is a problem:** The Doxygen `\file` says `mcrl2/data/standard_containers/deque.h` (wrong path AND wrong class name), and the `\brief` says "A deque class in which aterms can be stored." This makes the documentation actively misleading for anyone browsing API docs or using tooling that relies on Doxygen metadata.
- **Evidence type:** static-proof
- **Evidence:**
  - Observed result: `\file mcrl2/data/standard_containers/deque.h` and `\brief A deque class`
  - Expected result: `\file mcrl2/atermpp/standard_containers/stack.h` and `\brief A stack class`
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Replace all "deque" references with "stack" and fix the `\file` path.

### [SEV-3] Exposing `iterator`/`const_iterator` breaks stack adaptor abstraction
- **Location:** `stack.h:46-47`
- **Why this is a problem:** `std::stack` deliberately does not expose iterators — it's a LIFO adaptor. Publishing them invites users to traverse or randomly access the container, defeating the purpose of the restricted interface. The `mark()` method (the only operation needing iteration) delegates to the underlying container and does not need these typedefs at the stack level.
- **Evidence type:** static-proof
- **Evidence:**
  - `std::stack` exposes no `iterator` typedef per [container.adaptors.stack]
  - These typedefs serve no purpose in the current implementation
- **Efficiency impact:** None
- **C++20-specific note:** N/A
- **Suggested fix:** Remove `iterator` and `const_iterator` type aliases.

### [SEV-3] Missing `[[nodiscard]]` on query methods
- **Location:** `stack.h:115-132` (`top()`, `empty()`, `size()`)
- **Why this is a problem:** Per the project instructions, pure queries should be `[[nodiscard]]`. Discarding the result of `empty()` or `top()` is almost certainly a bug.
- **Evidence type:** plausible-issue (missed-bug-prevention)
- **Evidence:**
  - Pattern: `stack.empty();` (forgot the `if`) silently discards the result
  - Note: no existing standard_containers in mCRL2 use `[[nodiscard]]`, so this matches current practice, but does not match stated project policy
- **Efficiency impact:** None
- **C++20-specific note:** `[[nodiscard]]` is zero-cost and prevents logic errors at compile time.
- **Suggested fix:** Add `[[nodiscard]]` to `top()`, `empty()`, and `size()`.

### [SEV-3] Unnecessary includes
- **Location:** `stack.h:22,24`
- **Why this is a problem:** `thread_aterm_pool.h` and `shared_mutex.h` are included but nothing from them is directly used. They are already transitively included by `deque.h`. Unnecessary includes increase compilation time and create spurious coupling.
- **Evidence type:** static-proof
- **Efficiency impact:** Marginal compile-time overhead
- **Suggested fix:** Remove both includes; they come transitively through `deque.h`.

### [SEV-3] Unconstrained template parameters (missed C++20 opportunity)
- **Location:** `stack.h:63,69,75,101`
- **Why this is a problem:** `InputIt`, `InputIterator`, and `Alloc` parameters are completely unconstrained. The codebase already uses concepts in `aterm_container.h`. Concepts would provide clear error messages on misuse.
- **Evidence type:** plausible-issue (poor diagnostics)
- **C++20-specific note:** Use `std::input_iterator` concept for `InputIt` and a suitable allocator concept for `Alloc`.
- **Suggested fix:** Add `requires` clauses, e.g.:
  ```cpp
  template<std::input_iterator InputIt>
  stack(InputIt first, InputIt last) : m_container(first, last) {}
  ```

### [SEV-3] Inconsistent template parameter naming
- **Location:** `stack.h:63` (`InputIt`) vs `stack.h:69` (`InputIterator`)
- **Why this is a problem:** Two constructors performing the same conceptual operation use different naming for their iterator parameter. This makes the code appear as two unrelated overloads rather than allocator-aware variants of the same operation.
- **Suggested fix:** Standardize on `InputIt` (matching STL convention).

### [SEV-3] Missing comparison operators
- **Location:** entire class
- **Why this is a problem:** `std::stack` provides `operator==`, `!=`, `<`, `<=`, `>`, `>=`, and (C++20) `<=>`. This class provides none. Users comparing stacks will get a compilation error.
- **Evidence type:** static-proof
- **C++20-specific note:** A single `friend auto operator<=>(const stack&, const stack&) = default;` would provide all six (if the container supports it), or delegate to `m_container` comparisons.
- **Suggested fix:** Add `operator==` and `operator<=>` delegating to `m_container`.

### [SEV-3] `mark()` parameter uses raw `std::stack<...>` instead of `term_mark_stack` typedef
- **Location:** `stack.h:161`
- **Why this is a problem:** `aterm_container.h` defines `using term_mark_stack = std::stack<std::reference_wrapper<detail::_aterm>>`. Using the raw type instead of the typedef reduces readability and risks divergence if the typedef changes.
- **Suggested fix:** `void mark(detail::term_mark_stack& todo) const`

## Proposed Tests
- [ ] `stack_emplace_returns_reference`: Verify `emplace()` returns a reference to the inserted element. Input: `stack.emplace(value)` assigned to `auto&`. Expected: compiles and reference is valid.
- [ ] `stack_swap_noexcept_conditional`: Verify swap with a potentially-throwing container doesn't unconditionally promise `noexcept`. Input: custom container with `swap() noexcept(false)`. Expected: `static_assert(!noexcept(stack.swap(other)))`.
- [ ] `stack_comparison_operators`: Verify equality and ordering. Input: two stacks with same/different elements. Expected: `==` and `<=>` work correctly.

## Open Questions
1. Is the exposure of `iterator`/`const_iterator` intentional for some usage not visible in the repository? If so, it should be documented.
2. Should `top()` and `pop()` assert `!empty()` as a cheap precondition, or is the cost of the check unacceptable on a hot path (the binary I/O serialization loop)?
3. Are there plans to support user-provided `Container` types beyond `atermpp::deque`? If not, the unconditional `noexcept` on swap is benign in practice (but still technically incorrect).


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         4m 2s
Total session time:     4m 23s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         167.9k in, 5.7k out, 143.9k cached (Est. 1 Premium request)
 claude-haiku-4.5        464.4k in, 8.9k out, 418.3k cached (Est. 0 Premium requests)

```

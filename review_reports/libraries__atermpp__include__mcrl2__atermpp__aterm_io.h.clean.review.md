---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/aterm_io.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-07T23:54:22.121424+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h` – stream I/O interfaces and container serialization for aterms
- **Risk level:** Low-Medium (no confirmed crash, but design gaps that can cause UB or resource exhaustion under specific conditions)
- **Overall verdict:** Functionally correct for its intended trust model but has several hardening gaps and a missing Rule-of-Five guard on a RAII type.

## Findings (ordered by severity)

### [SEV-2] `aterm_stream_state` missing deleted copy/move operations

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:103-116`
- **Why this is a problem:** The RAII class stores a reference to an `aterm_stream` and restores the transformer in its destructor. If accidentally copied (e.g., passed by value, stored in a container, or captured in a lambda), two destructors will restore the transformer—the second potentially corrupting the stream state. The class holds a reference member, making copy semantics nonsensical.
- **Evidence type:** static-proof
- **Evidence:**
  - Given current definition, `aterm_stream_state s1(stream); auto s2 = s1;` compiles without warning.
  - Both `s1.~aterm_stream_state()` and `s2.~aterm_stream_state()` will call `m_stream.set_transformer(m_transformer)`, the second restoring stale state after `s1` already restored.
  - Run command: N/A (static analysis—the compiler accepts the copy)
  - Observed result: Double-restore compiles silently
  - Expected result: Copy/move should be deleted; compiler rejects misuse
- **Efficiency impact:** None directly, but incorrect state restoration could cause wrong transformers to be applied.
- **C++20-specific note:** `= delete` for copy/move is the standard idiom for non-copyable RAII guards since C++11.
- **Suggested fix:**
  ```cpp
  class aterm_stream_state
  {
  public:
    aterm_stream_state(aterm_stream& stream) : m_stream(stream)
    {
      m_transformer = stream.get_transformer();
    }
    ~aterm_stream_state() { m_stream.set_transformer(m_transformer); }

    aterm_stream_state(const aterm_stream_state&) = delete;
    aterm_stream_state& operator=(const aterm_stream_state&) = delete;
    aterm_stream_state(aterm_stream_state&&) = delete;
    aterm_stream_state& operator=(aterm_stream_state&&) = delete;

  private:
    aterm_stream& m_stream;
    aterm_transformer* m_transformer;
  };
  ```

---

### [SEV-2] Unbounded loop in container deserialization on corrupted/untrusted input

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:140-152`
- **Why this is a problem:** `nof_elements.value()` returns `std::size_t`. A corrupted or adversarial binary stream can encode an arbitrary `std::size_t` value (up to ~18×10¹⁸). The loop will attempt that many read + insert operations with no upper bound, leading to memory exhaustion or near-infinite runtime.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: An aterm binary stream encodes an `aterm_int` with value `SIZE_MAX`. The container read operator loops `SIZE_MAX` times, calling `stream >> element` and inserting.
  - Triggering scenario: Processing a binary `.lps` or `.pbes` file that has been truncated or corrupted on disk.
  - Reproducer/test name: `container_deser_corrupted_count`
  - Run command: Craft a binary stream with an `aterm_int(SIZE_MAX)` followed by no further data; call `operator>>(stream, vector<aterm>)`.
  - Observed result: Unbounded iteration until stream read failure or OOM.
  - Expected result: Either throw an exception when count exceeds stream remaining bytes, or impose a configurable upper bound.
- **Efficiency impact:** Denial-of-service via resource exhaustion.
- **C++20-specific note:** N/A
- **Suggested fix:** Add a sanity check (e.g., maximum count derived from stream size or a hard limit), or rely on the stream's `get()` to throw on EOF and document that contract explicitly.

---

### [SEV-3] `read_list_from_string` / `read_int_from_string` / `read_appl_from_string` rely on assert-only type validation

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:189-216`
- **Why this is a problem:** `down_cast` performs a `reinterpret_cast` with only an `assert` for validation. In release builds (`NDEBUG` defined), passing a string that doesn't encode the expected type causes silent UB—the returned reference is to a wrongly-typed object.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: `read_int_from_string("f(x)")` in a release build: `down_cast<aterm_int>` reinterpret_casts an `aterm` (function application) to `aterm_int&`. Subsequent `value()` calls `reinterpret_cast<_aterm_int*>(m_term)->value()` on memory that is actually a `_aterm_appl`—reading garbage.
  - Reproducer/test name: `read_int_from_string_wrong_type_release`
  - Run command: `cmake --build build --config Release && ./build/libraries/atermpp/test/aterm_io_test --run_test=read_int_wrong_type` (hypothetical)
  - Observed result: UB—garbage integer or potential segfault
  - Expected result: Throw `std::runtime_error` or return a sentinel
- **Efficiency impact:** Negligible (one branch on the happy path).
- **C++20-specific note:** Could use `if consteval` or `[[assume]]` to make the check zero-cost in optimized builds, but a simple `if (!n.type_is_int()) throw ...` is negligible.
- **Suggested fix:** Replace `assert` with a runtime check that throws on type mismatch, or add a `checked_down_cast` variant used in these public-facing functions.

---

### [SEV-3] `std::distance(container.begin(), container.end())` where `std::ranges::size` suffices

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:126`
- **Why this is a problem:** For `std::ranges::sized_range` types (all standard containers, all mCRL2 wrapped containers), `std::ranges::size()` is O(1). `std::distance` with non-random-access iterators is O(n), doubling the traversal cost for forward-only ranges.
- **Evidence type:** static-proof
- **Evidence:**
  - The concept constraint is `std::ranges::range<T>` which includes non-sized ranges.
  - For a hypothetical `std::forward_list<aterm>` (sized-range is false), `std::distance` is unavoidable. But for all actual mCRL2 containers (`atermpp::vector`, `atermpp::deque`, `atermpp::unordered_set`), `size()` is O(1).
  - No runtime failure—purely an efficiency concern.
- **Efficiency impact:** O(n) extra traversal for sized containers; negligible for small containers, measurable for large ones.
- **C++20-specific note:** Use `if constexpr (std::ranges::sized_range<T>) { std::ranges::size(container); }` for compile-time dispatch.
- **Suggested fix:**
  ```cpp
  if constexpr (std::ranges::sized_range<T>)
    stream << aterm_int(std::ranges::size(container));
  else
    stream << aterm_int(std::ranges::distance(container));
  ```

---

### [SEV-3] `aterm_transformer` type alias limits transformer to plain function pointers

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:24,42`
- **Why this is a problem:** `using aterm_transformer = aterm(const aterm&)` is a function type. The stored `m_transformer` is a function pointer. This prevents using stateful lambdas, `std::function`, or functor objects as transformers. All current usages happen to be plain functions (`identity`, `remove_index_impl`, `add_index_impl`), but the interface is artificially restrictive for future extension.
- **Evidence type:** static-proof
- **Evidence:**
  - `stream.set_transformer([&ctx](const aterm& t){ return transform(ctx, t); });` would fail to compile because a capturing lambda cannot decay to a function pointer.
- **Efficiency impact:** No runtime cost, but forces users into global/static functions, potentially requiring less-efficient workarounds (global state, thread_local hacks).
- **C++20-specific note:** `std::move_only_function<aterm(const aterm&)>` (C++23) or `std::function` would enable stateful transformers without virtual dispatch overhead for the common case.
- **Suggested fix:** This is a design trade-off. If statefulness is intentionally excluded for performance, document it. Otherwise, consider `std::function<aterm(const aterm&)>` with a null check replacing the identity default.

---

### [SEV-3] Typo in documentation comment

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:101`
- **Why this is a problem:** `"Currently, onlt preserves"` → should be `"only"`.
- **Evidence type:** static-proof
- **Suggested fix:** `s/onlt/only/`

---

## Proposed Tests

- [ ] `aterm_stream_state_not_copyable`: Verify that `aterm_stream_state s2 = s1;` fails to compile. Use a `static_assert(!std::is_copy_constructible_v<aterm_stream_state>)` or a compilation-failure test.
- [ ] `container_deser_corrupted_count`: Construct a mock `aterm_istream` that yields `aterm_int(1000000)` then immediately signals EOF/default-aterm. Verify the container read operator handles this gracefully (exception or early termination) rather than looping indefinitely.
- [ ] `read_int_from_string_wrong_type`: Call `read_int_from_string("f(x)")` in a release-mode test. Expected: exception or well-defined error, not UB. (Currently fails—would require the suggested fix first.)
- [ ] `container_ser_size_consistency`: Round-trip a `std::vector<aterm>` of 10000 elements through `operator<<` / `operator>>` and verify the deserialized container equals the original.

## Open Questions

1. Is untrusted/corrupted binary input an intended threat model for aterm streams? If so, the unbounded-loop finding escalates to SEV-1.
2. Is the function-pointer-only transformer a deliberate performance decision, or an artifact of the original design predating `std::function`?
3. Should `read_*_from_string` functions be hardened for public API use, or are they considered internal with trusted inputs?

## Copilot Raw Stdout

```text
● Explore: Investigate aterm types
  ● Search for container deserialization in tests
    $ grep -rn "std::vector.*>>" /workspaces/mCRL2/libraries/atermpp/test/
    └ 1 line...
  ● Read libraries/atermpp/include/mcrl2/atermpp/aterm_core.h lines 33-80
    └ 48 lines read

● Explore: Check aterm_stream_state usage
  ● Display comprehensive review
    $ cat > /tmp/aterm_io_review.md << 'EOF'
    # aterm_io.h Code Review - Comprehensive Findings

    ## 1. Is `aterm_stream_state` Ever Copied?...
    └ 497 lines...
  ● Display final summary
    $ cat << 'EOF'

    ═══════════════════════════════════════════════════════════════════════════════
                              QUICK REFERENCE SUMMARY
    ═══════════════════════════════════════════════════════════════════════════════...
    └ 198 lines...

Now I have enough information to produce the review. Let me verify one more thing about the `std::inserter` behavior with specific container types used in mCRL2:

● Grep "aterm_stream_state"
  └ 46 lines found

● Grep "set_transformer\(identity\)|set_transformer\(nullptr\)"
  └ No matches found

● Read libraries/atermpp/include/mcrl2/atermpp/aterm_int.h lines 1-55
  └ 55 lines read

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h` – stream I/O interfaces and container serialization for aterms
- **Risk level:** Low-Medium (no confirmed crash, but design gaps that can cause UB or resource exhaustion under specific conditions)
- **Overall verdict:** Functionally correct for its intended trust model but has several hardening gaps and a missing Rule-of-Five guard on a RAII type.

## Findings (ordered by severity)

### [SEV-2] `aterm_stream_state` missing deleted copy/move operations

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:103-116`
- **Why this is a problem:** The RAII class stores a reference to an `aterm_stream` and restores the transformer in its destructor. If accidentally copied (e.g., passed by value, stored in a container, or captured in a lambda), two destructors will restore the transformer—the second potentially corrupting the stream state. The class holds a reference member, making copy semantics nonsensical.
- **Evidence type:** static-proof
- **Evidence:**
  - Given current definition, `aterm_stream_state s1(stream); auto s2 = s1;` compiles without warning.
  - Both `s1.~aterm_stream_state()` and `s2.~aterm_stream_state()` will call `m_stream.set_transformer(m_transformer)`, the second restoring stale state after `s1` already restored.
  - Run command: N/A (static analysis—the compiler accepts the copy)
  - Observed result: Double-restore compiles silently
  - Expected result: Copy/move should be deleted; compiler rejects misuse
- **Efficiency impact:** None directly, but incorrect state restoration could cause wrong transformers to be applied.
- **C++20-specific note:** `= delete` for copy/move is the standard idiom for non-copyable RAII guards since C++11.
- **Suggested fix:**
  ```cpp
  class aterm_stream_state
  {
  public:
    aterm_stream_state(aterm_stream& stream) : m_stream(stream)
    {
      m_transformer = stream.get_transformer();
    }
    ~aterm_stream_state() { m_stream.set_transformer(m_transformer); }

    aterm_stream_state(const aterm_stream_state&) = delete;
    aterm_stream_state& operator=(const aterm_stream_state&) = delete;
    aterm_stream_state(aterm_stream_state&&) = delete;
    aterm_stream_state& operator=(aterm_stream_state&&) = delete;

  private:
    aterm_stream& m_stream;
    aterm_transformer* m_transformer;
  };
  ```

---

### [SEV-2] Unbounded loop in container deserialization on corrupted/untrusted input

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:140-152`
- **Why this is a problem:** `nof_elements.value()` returns `std::size_t`. A corrupted or adversarial binary stream can encode an arbitrary `std::size_t` value (up to ~18×10¹⁸). The loop will attempt that many read + insert operations with no upper bound, leading to memory exhaustion or near-infinite runtime.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: An aterm binary stream encodes an `aterm_int` with value `SIZE_MAX`. The container read operator loops `SIZE_MAX` times, calling `stream >> element` and inserting.
  - Triggering scenario: Processing a binary `.lps` or `.pbes` file that has been truncated or corrupted on disk.
  - Reproducer/test name: `container_deser_corrupted_count`
  - Run command: Craft a binary stream with an `aterm_int(SIZE_MAX)` followed by no further data; call `operator>>(stream, vector<aterm>)`.
  - Observed result: Unbounded iteration until stream read failure or OOM.
  - Expected result: Either throw an exception when count exceeds stream remaining bytes, or impose a configurable upper bound.
- **Efficiency impact:** Denial-of-service via resource exhaustion.
- **C++20-specific note:** N/A
- **Suggested fix:** Add a sanity check (e.g., maximum count derived from stream size or a hard limit), or rely on the stream's `get()` to throw on EOF and document that contract explicitly.

---

### [SEV-3] `read_list_from_string` / `read_int_from_string` / `read_appl_from_string` rely on assert-only type validation

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:189-216`
- **Why this is a problem:** `down_cast` performs a `reinterpret_cast` with only an `assert` for validation. In release builds (`NDEBUG` defined), passing a string that doesn't encode the expected type causes silent UB—the returned reference is to a wrongly-typed object.
- **Evidence type:** plausible-issue
- **Evidence:**
  - Mechanism: `read_int_from_string("f(x)")` in a release build: `down_cast<aterm_int>` reinterpret_casts an `aterm` (function application) to `aterm_int&`. Subsequent `value()` calls `reinterpret_cast<_aterm_int*>(m_term)->value()` on memory that is actually a `_aterm_appl`—reading garbage.
  - Reproducer/test name: `read_int_from_string_wrong_type_release`
  - Run command: `cmake --build build --config Release && ./build/libraries/atermpp/test/aterm_io_test --run_test=read_int_wrong_type` (hypothetical)
  - Observed result: UB—garbage integer or potential segfault
  - Expected result: Throw `std::runtime_error` or return a sentinel
- **Efficiency impact:** Negligible (one branch on the happy path).
- **C++20-specific note:** Could use `if consteval` or `[[assume]]` to make the check zero-cost in optimized builds, but a simple `if (!n.type_is_int()) throw ...` is negligible.
- **Suggested fix:** Replace `assert` with a runtime check that throws on type mismatch, or add a `checked_down_cast` variant used in these public-facing functions.

---

### [SEV-3] `std::distance(container.begin(), container.end())` where `std::ranges::size` suffices

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:126`
- **Why this is a problem:** For `std::ranges::sized_range` types (all standard containers, all mCRL2 wrapped containers), `std::ranges::size()` is O(1). `std::distance` with non-random-access iterators is O(n), doubling the traversal cost for forward-only ranges.
- **Evidence type:** static-proof
- **Evidence:**
  - The concept constraint is `std::ranges::range<T>` which includes non-sized ranges.
  - For a hypothetical `std::forward_list<aterm>` (sized-range is false), `std::distance` is unavoidable. But for all actual mCRL2 containers (`atermpp::vector`, `atermpp::deque`, `atermpp::unordered_set`), `size()` is O(1).
  - No runtime failure—purely an efficiency concern.
- **Efficiency impact:** O(n) extra traversal for sized containers; negligible for small containers, measurable for large ones.
- **C++20-specific note:** Use `if constexpr (std::ranges::sized_range<T>) { std::ranges::size(container); }` for compile-time dispatch.
- **Suggested fix:**
  ```cpp
  if constexpr (std::ranges::sized_range<T>)
    stream << aterm_int(std::ranges::size(container));
  else
    stream << aterm_int(std::ranges::distance(container));
  ```

---

### [SEV-3] `aterm_transformer` type alias limits transformer to plain function pointers

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:24,42`
- **Why this is a problem:** `using aterm_transformer = aterm(const aterm&)` is a function type. The stored `m_transformer` is a function pointer. This prevents using stateful lambdas, `std::function`, or functor objects as transformers. All current usages happen to be plain functions (`identity`, `remove_index_impl`, `add_index_impl`), but the interface is artificially restrictive for future extension.
- **Evidence type:** static-proof
- **Evidence:**
  - `stream.set_transformer([&ctx](const aterm& t){ return transform(ctx, t); });` would fail to compile because a capturing lambda cannot decay to a function pointer.
- **Efficiency impact:** No runtime cost, but forces users into global/static functions, potentially requiring less-efficient workarounds (global state, thread_local hacks).
- **C++20-specific note:** `std::move_only_function<aterm(const aterm&)>` (C++23) or `std::function` would enable stateful transformers without virtual dispatch overhead for the common case.
- **Suggested fix:** This is a design trade-off. If statefulness is intentionally excluded for performance, document it. Otherwise, consider `std::function<aterm(const aterm&)>` with a null check replacing the identity default.

---

### [SEV-3] Typo in documentation comment

- **Location:** `libraries/atermpp/include/mcrl2/atermpp/aterm_io.h:101`
- **Why this is a problem:** `"Currently, onlt preserves"` → should be `"only"`.
- **Evidence type:** static-proof
- **Suggested fix:** `s/onlt/only/`

---

## Proposed Tests

- [ ] `aterm_stream_state_not_copyable`: Verify that `aterm_stream_state s2 = s1;` fails to compile. Use a `static_assert(!std::is_copy_constructible_v<aterm_stream_state>)` or a compilation-failure test.
- [ ] `container_deser_corrupted_count`: Construct a mock `aterm_istream` that yields `aterm_int(1000000)` then immediately signals EOF/default-aterm. Verify the container read operator handles this gracefully (exception or early termination) rather than looping indefinitely.
- [ ] `read_int_from_string_wrong_type`: Call `read_int_from_string("f(x)")` in a release-mode test. Expected: exception or well-defined error, not UB. (Currently fails—would require the suggested fix first.)
- [ ] `container_ser_size_consistency`: Round-trip a `std::vector<aterm>` of 10000 elements through `operator<<` / `operator>>` and verify the deserialized container equals the original.

## Open Questions

1. Is untrusted/corrupted binary input an intended threat model for aterm streams? If so, the unbounded-loop finding escalates to SEV-1.
2. Is the function-pointer-only transformer a deliberate performance decision, or an artifact of the original design predating `std::function`?
3. Should `read_*_from_string` functions be hardened for public API use, or are they considered internal with trusted inputs?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         6m 33s
Total session time:     7m 8s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         161.9k in, 8.2k out, 135.9k cached (Est. 1 Premium request)
 claude-haiku-4.5        1.0m in, 20.0k out, 974.2k cached (Est. 0 Premium requests)

```

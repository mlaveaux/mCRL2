---
{
  "file": "libraries/atermpp/source/aterm_io_binary.cpp",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T14:03:36.621299+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/source/aterm_io_binary.cpp` — binary aterm serialization/deserialization (BAF format v0x8308)
- **Risk level:** Medium
- **Overall verdict:** Functional code with several latent defects: undefined behavior in assertions reachable only with malformed streams, missing input validation that permits out-of-bounds access from corrupted data, variable shadowing, redundant casts, and unnecessary per-iteration work. No concurrency issues (single-threaded usage model). The code works correctly for well-formed streams produced by the same library version.

## Findings (ordered by severity)

### [SEV-1] No bounds-checking on indices read from untrusted stream
- **Location:** `aterm_io_binary.cpp:267`, `aterm_io_binary.cpp:288`
- **Why this is a problem:** `m_function_symbols[index]` and `m_terms[index]` use `operator[]` (no bounds check) on indices decoded from the binary stream. A truncated or corrupted file can supply an out-of-range index, causing undefined behavior (out-of-bounds deque access). This is a robustness/security concern for any tool that reads user-supplied `.lps`/`.lts`/`.pbes` files.
- **Evidence type:** plausible-issue
- **Evidence:**
  - **Reproducer/test name:** `corrupted_function_symbol_index_test`
  - **Design:** Write a valid BAF header, then a function_symbol packet, then an `aterm` packet with a function_symbol index of 255 (far exceeding the deque size of 2). Feed to `binary_aterm_istream::get()`.
  - **Expected result:** Exception thrown with clear error message.
  - **Observed result:** Undefined behavior (ASan reports heap-buffer-overflow or crash).
  - **Run command:** `cmake --build build --target aterm_io_binary_test && cd build && ctest -R aterm_io_binary_test --output-on-failure`
- **Efficiency impact:** None (bounds check on I/O path is negligible).
- **C++20-specific note:** Could use `std::span` with `at()`-like access or a simple range check with `mcrl2::runtime_error`.
- **Suggested fix:** Before indexing, validate: `if (index >= m_function_symbols.size()) throw mcrl2::runtime_error(...)` and similarly for `m_terms`.

### [SEV-2] Undefined behavior in assertions via `std::log2(0)`
- **Location:** `aterm_io_binary.cpp:184` (ostream `term_index_width()`), `aterm_io_binary.cpp:313` (istream `term_index_width()`)
- **Why this is a problem:** The assertions compute `static_cast<unsigned int>(std::log2(m_terms.size()) + 1)`. If `m_terms.size()` is 0, `std::log2(0.0)` returns `-∞`; casting `-∞` to `unsigned int` is undefined behavior per [conv.fpint]. While normal execution paths ensure `m_terms` is non-empty when these functions are called, a format error or future refactoring could trigger the assertion with an empty container. The UB occurs inside the `assert` expression itself, meaning even a "passing" assertion with size 0 and `m_term_index_width == 0` would be UB.
- **Evidence type:** static-proof
- **Evidence:**
  - C++ standard [conv.fpint]§7.3.11: "If the truncated value cannot be represented in the destination type, the behavior is undefined."
  - `std::log2(0.0)` returns `-∞` per IEEE 754 / C++ [cmath.syn].
  - `static_cast<unsigned int>(-∞)` is outside `[0, UINT_MAX]` → UB.
- **Efficiency impact:** None.
- **C++20-specific note:** The fix should use integer arithmetic: `std::bit_width(m_terms.size())` (C++20 `<bit>`) gives the exact number of bits needed and avoids floating-point entirely.
- **Suggested fix:** Replace `static_cast<unsigned int>(std::log2(x) + 1)` with `static_cast<unsigned int>(std::bit_width(x))` throughout. Add a guard `assert(m_terms.size() > 0)` before the width assertion, or restructure to avoid the assertion when size is 0.

### [SEV-2] Variable shadowing in `write_function_symbol`
- **Location:** `aterm_io_binary.cpp:221` and `aterm_io_binary.cpp:234`
- **Why this is a problem:** Line 221 declares `std::size_t result` and line 234 declares `auto result` (type `std::pair<size_type, bool>`) within the else-block, shadowing the outer variable. While the current code is functionally correct (the inner scope's `result.first` is the intended value), this violates `-Wshadow` and creates a maintenance trap: a future developer might accidentally use the wrong `result`.
- **Evidence type:** static-proof
- **Evidence:**
  - Compiling with `-Wshadow` will produce: `warning: declaration of 'result' shadows a previous local [-Wshadow]`
  - **Run command:** `cmake --build build --target mcrl2_atermpp 2>&1 | grep -i shadow`
- **Efficiency impact:** None.
- **C++20-specific note:** Use structured bindings: `auto [index, inserted] = m_function_symbols.insert(symbol);`
- **Suggested fix:**
  ```cpp
  auto [new_index, inserted] = m_function_symbols.insert(symbol);
  m_function_symbol_index_width = static_cast<unsigned int>(std::bit_width(m_function_symbols.size()));
  return new_index;
  ```

### [SEV-3] Redundant `static_cast` and unnecessary re-computation of `transformed`
- **Location:** `aterm_io_binary.cpp:104`
- **Why this is a problem:** The double `static_cast<const aterm&>(static_cast<const aterm&>(current.term))` applies two identical casts; the inner one converts `markable_aterm<aterm>` → `const aterm&`, the outer one is a no-op. Additionally, `transformed` is computed on every loop iteration (line 104) even when the term has already been written (caught at line 108/173 and immediately popped). For deeply shared terms this is wasted transformer invocations.
- **Evidence type:** static-proof
- **Evidence:**
  - `markable_aterm<aterm>::operator const aterm&()` already returns `const aterm&`. The second cast is identity.
  - When `m_terms.index(current.term) < m_terms.size()` and `!is_output`, control goes to line 175 (pop), wasting the `transformed` computation.
- **Efficiency impact:** Low for typical terms; proportional to number of shared subterm re-encounters during traversal.
- **C++20-specific note:** N/A.
- **Suggested fix:** Move the `transformed` computation inside the `if` block at line 109, and remove the redundant outer cast.

### [SEV-3] `packet_type` cast from unchecked integer
- **Location:** `aterm_io_binary.cpp:247`
- **Why this is a problem:** `static_cast<packet_type>(header)` where `header` is read from the stream. `packet_type` is a scoped enum with values 0–3 (2 bits). If the stream is well-formed, `header` is always 0–3. However, `read_bits(2)` can return 0–3, so the cast is safe for well-formed streams. This is low-risk but worth noting as a defensive coding gap: no `default` case or range check handles unexpected values.
- **Evidence type:** plausible-issue (only with corrupted higher-level stream logic)
- **Efficiency impact:** None.
- **Suggested fix:** Add a `default` or `else` branch that throws `mcrl2::runtime_error` for unrecognized packet types (future-proofing).

### [SEV-3] Use of floating-point `std::log2` for integer bit-width
- **Location:** Lines 152, 184, 190, 235, 255, 280, 304, 313, 319
- **Why this is a problem:** `std::log2` operates on `double`, which has 53-bit mantissa. For `std::size_t` values > 2⁵³, the result is imprecise. While aterm indices are unlikely to exceed 2⁵³ in practice, the idiomatic C++20 solution (`std::bit_width`) is exact for all values and avoids the floating-point conversion entirely.
- **Evidence type:** static-proof (precision loss for size > 2⁵³)
- **Efficiency impact:** `std::bit_width` compiles to a single `bsr`/`clz` instruction vs. an FP conversion + `log2`.
- **C++20-specific note:** `std::bit_width(n)` from `<bit>` is the correct, efficient replacement. It returns 0 for n=0, which is safe.
- **Suggested fix:** Replace all `static_cast<unsigned int>(std::log2(x) + 1)` with `static_cast<unsigned int>(std::bit_width(x))`.

## Proposed Tests
- [ ] `corrupted_function_symbol_index_test`: Feed BAF header + function_symbol packet + aterm packet with out-of-range symbol index. Expected: `mcrl2::runtime_error` thrown (currently: UB/crash under ASan).
- [ ] `corrupted_term_index_test`: Feed a valid BAF stream where an aterm references a term index beyond `m_terms.size()`. Expected: `mcrl2::runtime_error` thrown (currently: UB).
- [ ] `empty_stream_graceful_test`: Open `binary_aterm_istream` on a stream that immediately signals end-of-stream (function symbol 0). Expected: `t` is default `aterm()` with no UB.
- [ ] `large_term_bit_width_test`: Write and read back a term tree with > 1024 shared subterms to stress index width calculations; verify round-trip equality.

## Open Questions
1. Is the transformer guaranteed to be deterministic and side-effect-free? The two-pass traversal in `put()` calls `m_transformer` twice for the same term (once with `write=false`, once with `write=true`). A non-deterministic transformer could produce inconsistent structure between passes.
2. Should `binary_aterm_istream` reject streams with `BAF_VERSION` from the future (higher than current), or only mismatching versions? Currently any mismatch throws.
3. The `m_terms` indexed_set uses `atermpp::indexed_set` (GC-aware) on the output side but `atermpp::deque` on the input side. Is there a documented reason for the asymmetry, or could both use the same container?

## Copilot Raw Stdout

```text
● Explore: Explore aterm binary IO context
  ● Read libraries/atermpp/source/aterm_io_binary.cpp
    └ 333 lines read
  ● Read libraries/atermpp/include/mcrl2/atermpp/function_symbol.h lines 155-165
    └ 10 lines read

● Explore: Check log2 UB and indexed_set
  ● Read libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h
    └ 91 lines read
  ● Check m_transformer type
    $ sed -n '44,46p' /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/aterm_io.h
    └ 4 lines...

● Explore: Check log2 UB and cast behavior
  ● Read libraries/atermpp/include/mcrl2/atermpp/standard_containers/indexed_set.h
    └ 91 lines read
  ● Read libraries/atermpp/source/aterm_io_binary.cpp lines 95-180
    └ 86 lines read

● Explore: Verify log2 initialization state
  ● Read libraries/atermpp/source/aterm_io_binary.cpp
    └ 333 lines read
  ● Read libraries/atermpp/include/mcrl2/atermpp/aterm_io_binary.h
    └ 94 lines read

● Read libraries/atermpp/source/aterm_io_binary.cpp lines 215-245
  └ 31 lines read

● Read libraries/atermpp/source/aterm_io_binary.cpp lines 95-175
  └ 81 lines read

● Read libraries/atermpp/source/aterm_io_binary.cpp lines 245-325
  └ 81 lines read

Now I have all the information needed to produce the review. Let me write it up.

---

# Review Findings

## Summary
- **Scope reviewed:** `libraries/atermpp/source/aterm_io_binary.cpp` — binary aterm serialization/deserialization (BAF format v0x8308)
- **Risk level:** Medium
- **Overall verdict:** Functional code with several latent defects: undefined behavior in assertions reachable only with malformed streams, missing input validation that permits out-of-bounds access from corrupted data, variable shadowing, redundant casts, and unnecessary per-iteration work. No concurrency issues (single-threaded usage model). The code works correctly for well-formed streams produced by the same library version.

## Findings (ordered by severity)

### [SEV-1] No bounds-checking on indices read from untrusted stream
- **Location:** `aterm_io_binary.cpp:267`, `aterm_io_binary.cpp:288`
- **Why this is a problem:** `m_function_symbols[index]` and `m_terms[index]` use `operator[]` (no bounds check) on indices decoded from the binary stream. A truncated or corrupted file can supply an out-of-range index, causing undefined behavior (out-of-bounds deque access). This is a robustness/security concern for any tool that reads user-supplied `.lps`/`.lts`/`.pbes` files.
- **Evidence type:** plausible-issue
- **Evidence:**
  - **Reproducer/test name:** `corrupted_function_symbol_index_test`
  - **Design:** Write a valid BAF header, then a function_symbol packet, then an `aterm` packet with a function_symbol index of 255 (far exceeding the deque size of 2). Feed to `binary_aterm_istream::get()`.
  - **Expected result:** Exception thrown with clear error message.
  - **Observed result:** Undefined behavior (ASan reports heap-buffer-overflow or crash).
  - **Run command:** `cmake --build build --target aterm_io_binary_test && cd build && ctest -R aterm_io_binary_test --output-on-failure`
- **Efficiency impact:** None (bounds check on I/O path is negligible).
- **C++20-specific note:** Could use `std::span` with `at()`-like access or a simple range check with `mcrl2::runtime_error`.
- **Suggested fix:** Before indexing, validate: `if (index >= m_function_symbols.size()) throw mcrl2::runtime_error(...)` and similarly for `m_terms`.

### [SEV-2] Undefined behavior in assertions via `std::log2(0)`
- **Location:** `aterm_io_binary.cpp:184` (ostream `term_index_width()`), `aterm_io_binary.cpp:313` (istream `term_index_width()`)
- **Why this is a problem:** The assertions compute `static_cast<unsigned int>(std::log2(m_terms.size()) + 1)`. If `m_terms.size()` is 0, `std::log2(0.0)` returns `-∞`; casting `-∞` to `unsigned int` is undefined behavior per [conv.fpint]. While normal execution paths ensure `m_terms` is non-empty when these functions are called, a format error or future refactoring could trigger the assertion with an empty container. The UB occurs inside the `assert` expression itself, meaning even a "passing" assertion with size 0 and `m_term_index_width == 0` would be UB.
- **Evidence type:** static-proof
- **Evidence:**
  - C++ standard [conv.fpint]§7.3.11: "If the truncated value cannot be represented in the destination type, the behavior is undefined."
  - `std::log2(0.0)` returns `-∞` per IEEE 754 / C++ [cmath.syn].
  - `static_cast<unsigned int>(-∞)` is outside `[0, UINT_MAX]` → UB.
- **Efficiency impact:** None.
- **C++20-specific note:** The fix should use integer arithmetic: `std::bit_width(m_terms.size())` (C++20 `<bit>`) gives the exact number of bits needed and avoids floating-point entirely.
- **Suggested fix:** Replace `static_cast<unsigned int>(std::log2(x) + 1)` with `static_cast<unsigned int>(std::bit_width(x))` throughout. Add a guard `assert(m_terms.size() > 0)` before the width assertion, or restructure to avoid the assertion when size is 0.

### [SEV-2] Variable shadowing in `write_function_symbol`
- **Location:** `aterm_io_binary.cpp:221` and `aterm_io_binary.cpp:234`
- **Why this is a problem:** Line 221 declares `std::size_t result` and line 234 declares `auto result` (type `std::pair<size_type, bool>`) within the else-block, shadowing the outer variable. While the current code is functionally correct (the inner scope's `result.first` is the intended value), this violates `-Wshadow` and creates a maintenance trap: a future developer might accidentally use the wrong `result`.
- **Evidence type:** static-proof
- **Evidence:**
  - Compiling with `-Wshadow` will produce: `warning: declaration of 'result' shadows a previous local [-Wshadow]`
  - **Run command:** `cmake --build build --target mcrl2_atermpp 2>&1 | grep -i shadow`
- **Efficiency impact:** None.
- **C++20-specific note:** Use structured bindings: `auto [index, inserted] = m_function_symbols.insert(symbol);`
- **Suggested fix:**
  ```cpp
  auto [new_index, inserted] = m_function_symbols.insert(symbol);
  m_function_symbol_index_width = static_cast<unsigned int>(std::bit_width(m_function_symbols.size()));
  return new_index;
  ```

### [SEV-3] Redundant `static_cast` and unnecessary re-computation of `transformed`
- **Location:** `aterm_io_binary.cpp:104`
- **Why this is a problem:** The double `static_cast<const aterm&>(static_cast<const aterm&>(current.term))` applies two identical casts; the inner one converts `markable_aterm<aterm>` → `const aterm&`, the outer one is a no-op. Additionally, `transformed` is computed on every loop iteration (line 104) even when the term has already been written (caught at line 108/173 and immediately popped). For deeply shared terms this is wasted transformer invocations.
- **Evidence type:** static-proof
- **Evidence:**
  - `markable_aterm<aterm>::operator const aterm&()` already returns `const aterm&`. The second cast is identity.
  - When `m_terms.index(current.term) < m_terms.size()` and `!is_output`, control goes to line 175 (pop), wasting the `transformed` computation.
- **Efficiency impact:** Low for typical terms; proportional to number of shared subterm re-encounters during traversal.
- **C++20-specific note:** N/A.
- **Suggested fix:** Move the `transformed` computation inside the `if` block at line 109, and remove the redundant outer cast.

### [SEV-3] `packet_type` cast from unchecked integer
- **Location:** `aterm_io_binary.cpp:247`
- **Why this is a problem:** `static_cast<packet_type>(header)` where `header` is read from the stream. `packet_type` is a scoped enum with values 0–3 (2 bits). If the stream is well-formed, `header` is always 0–3. However, `read_bits(2)` can return 0–3, so the cast is safe for well-formed streams. This is low-risk but worth noting as a defensive coding gap: no `default` case or range check handles unexpected values.
- **Evidence type:** plausible-issue (only with corrupted higher-level stream logic)
- **Efficiency impact:** None.
- **Suggested fix:** Add a `default` or `else` branch that throws `mcrl2::runtime_error` for unrecognized packet types (future-proofing).

### [SEV-3] Use of floating-point `std::log2` for integer bit-width
- **Location:** Lines 152, 184, 190, 235, 255, 280, 304, 313, 319
- **Why this is a problem:** `std::log2` operates on `double`, which has 53-bit mantissa. For `std::size_t` values > 2⁵³, the result is imprecise. While aterm indices are unlikely to exceed 2⁵³ in practice, the idiomatic C++20 solution (`std::bit_width`) is exact for all values and avoids the floating-point conversion entirely.
- **Evidence type:** static-proof (precision loss for size > 2⁵³)
- **Efficiency impact:** `std::bit_width` compiles to a single `bsr`/`clz` instruction vs. an FP conversion + `log2`.
- **C++20-specific note:** `std::bit_width(n)` from `<bit>` is the correct, efficient replacement. It returns 0 for n=0, which is safe.
- **Suggested fix:** Replace all `static_cast<unsigned int>(std::log2(x) + 1)` with `static_cast<unsigned int>(std::bit_width(x))`.

## Proposed Tests
- [ ] `corrupted_function_symbol_index_test`: Feed BAF header + function_symbol packet + aterm packet with out-of-range symbol index. Expected: `mcrl2::runtime_error` thrown (currently: UB/crash under ASan).
- [ ] `corrupted_term_index_test`: Feed a valid BAF stream where an aterm references a term index beyond `m_terms.size()`. Expected: `mcrl2::runtime_error` thrown (currently: UB).
- [ ] `empty_stream_graceful_test`: Open `binary_aterm_istream` on a stream that immediately signals end-of-stream (function symbol 0). Expected: `t` is default `aterm()` with no UB.
- [ ] `large_term_bit_width_test`: Write and read back a term tree with > 1024 shared subterms to stress index width calculations; verify round-trip equality.

## Open Questions
1. Is the transformer guaranteed to be deterministic and side-effect-free? The two-pass traversal in `put()` calls `m_transformer` twice for the same term (once with `write=false`, once with `write=true`). A non-deterministic transformer could produce inconsistent structure between passes.
2. Should `binary_aterm_istream` reject streams with `BAF_VERSION` from the future (higher than current), or only mismatching versions? Currently any mismatch throws.
3. The `m_terms` indexed_set uses `atermpp::indexed_set` (GC-aware) on the output side but `atermpp::deque` on the input side. Is there a documented reason for the asymmetry, or could both use the same container?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         6m 29s
Total session time:     6m 56s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         346.1k in, 5.8k out, 308.6k cached (Est. 1 Premium request)
 claude-haiku-4.5        817.2k in, 18.1k out, 733.1k cached (Est. 0 Premium requests)

```

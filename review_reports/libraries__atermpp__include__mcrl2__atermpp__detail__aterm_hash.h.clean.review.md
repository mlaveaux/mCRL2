---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/aterm_hash.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T09:18:58.931355+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_hash.h` (full file, clean review)
- Risk level: **Medium-High** — contains strict aliasing UB and a latent correctness hazard for arity > 7
- Overall verdict: The file is functional for current usage (arity ≤ 7) but harbors undefined behavior via hardcoded type casts and a risky argument-comparison-before-symbol-check pattern that could cause out-of-bounds reads on hash collision.

## Findings (ordered by severity)

### [SEV-1] Hardcoded `_aterm_appl<8>` in `equal_args` and variadic `aterm_equals_finite<N>::operator()`
- Location: `aterm_hash.h:376, 383, 393`
- Why this is a problem: `equal_args()` unconditionally uses `_aterm_appl<8>` regardless of the enclosing `aterm_equals_finite<N>` template parameter. The cast at line 393 (`static_cast<const _aterm_appl<8>&>(term)`) violates strict aliasing when the actual object was constructed as `_aterm_appl<N>` with N ≠ 8. Since `function_application_storage<N>` is instantiated for N=1..7, every instantiation performs a type-punning cast. Additionally, if the pool were ever extended to N ≥ 9, `arg(I)` would access index I ≥ 8 on an `_aterm_appl<8>` — reading past the declared array.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Any instantiation of `aterm_equals_finite<N>::operator()(...args)` for N != 8; currently all N in {0,1,2,3,4,5,6,7} trigger it.
  - Run command: `cmake -B build-asan -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_ADDRESSSANITIZER=ON -DMCRL2_ENABLE_GUI_TOOLS=OFF && cmake --build build-asan && ctest --test-dir build-asan -R atermpp --output-on-failure`
  - Observed result: Strict aliasing UB (not flagged by ASan but is UB per [basic.lval]). On optimized builds compilers may miscompile.
  - Expected result: Cast should target `_aterm_appl<N>` matching the template parameter.
- Efficiency impact: None at runtime (works by accident due to identical memory layout).
- C++20-specific note: Violation of strict aliasing ([basic.lval] p11). Could use `std::start_lifetime_as` in C++23 or `std::launder`/careful placement to make this legal.
- Suggested fix:
  ```cpp
  template <std::size_t N, std::size_t I = 0, typename... Tp>
  bool equal_args(const _aterm_appl<N>&, const Tp&...)
    requires(I == sizeof...(Tp))
  { return true; }

  template <std::size_t N, std::size_t I = 0, typename... Tp>
  bool equal_args(const _aterm_appl<N>& term, const Tp&... t)
    requires(I < sizeof...(Tp))
  {
    return term.arg(I) == std::get<I>(std::forward_as_tuple(t...)) && equal_args<N, I+1>(term, t...);
  }

  template<std::size_t N>
  template<typename ...Args>
  bool aterm_equals_finite<N>::operator()(const _aterm& term, const function_symbol& symbol, const Args&... args) const noexcept
  {
    return term.function() == symbol && equal_args<N>(static_cast<const _aterm_appl<N>&>(term), args...);
  }
  ```

### [SEV-2] Equality check compares arguments before function symbol — potential out-of-bounds on hash collision
- Location: `aterm_hash.h:289-312`
- Why this is a problem: `aterm_equals<N>::operator()(first, second)` uses `first.function().arity()` (when N = DynamicNumberOfArguments) to determine how many arguments of `second` to read. On a hash collision where `second` has fewer arguments than `first`, this reads past `second`'s allocated arguments. The function symbol mismatch is only detected at line 311 *after* the out-of-bounds loop. The `arbitrary_function_application_storage` (arity ≥ 8) uses this comparator.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Construct two terms with different arities (≥ 8) whose hash collides; force a lookup triggering this equality check.
  - Run command: `ctest --test-dir build-asan -R atermpp --output-on-failure` (with ASan; requires triggering the collision)
  - Observed result: (plausible) ASan heap-buffer-overflow or use-of-uninitialized if term allocations don't have sufficient trailing padding.
  - Expected result: Function symbol should be compared first; early-return on mismatch.
- Efficiency impact: Even without UB, comparing arguments before the function symbol wastes cycles on collisions.
- Suggested fix: Move `if (first.function() != second.function()) return false;` before the argument loop.

### [SEV-2] `aterm_equals<N>::operator()(term, symbol, arguments*)` ignores N template parameter
- Location: `aterm_hash.h:320-333`
- Why this is a problem: This overload unconditionally uses `symbol.arity()` (line 324) instead of the N-based logic (`(N == DynamicNumberOfArguments) ? symbol.arity() : N`) used by all sibling overloads. When N is fixed, the enclosing hash table guarantees all terms have arity N. Using `symbol.arity()` is redundant but not harmful *if* the caller always passes a symbol with arity N. However, it's inconsistent and means a mismatched symbol could cause out-of-bounds reads on `term`.
- Evidence type: static-proof (inconsistency with other overloads; lines 175, 197, 217, 299, 343 all use the N-based conditional)
- Suggested fix: Use `const std::size_t arity = (N == DynamicNumberOfArguments) ? symbol.arity() : N;` for consistency and defense in depth.

### [SEV-3] `aterm_int_hasher::operator()(_aterm_int)` creates unnecessary temporary
- Location: `aterm_hash.h:270-273`
- Why this is a problem: `return aterm_int_hasher()(term.value());` constructs a temporary `aterm_int_hasher` to call the `size_t` overload. This is functionally correct but confusing and marginally wasteful.
- Evidence type: static-proof (code inspection)
- Efficiency impact: Negligible (optimizer likely eliminates), but confusing.
- Suggested fix: `return (*this)(term.value());` or directly inline the hash computation.

### [SEV-3] `static inline` on `xorshift` creates per-TU internal-linkage copies
- Location: `aterm_hash.h:277`
- Why this is a problem: `static` at namespace scope gives internal linkage, meaning each translation unit gets its own copy. The `inline` is redundant with `static`. For a tiny helper, `inline` alone (external linkage, single definition) is preferred per modern C++ convention.
- Evidence type: static-proof (language semantics)
- Efficiency impact: Negligible code bloat; prevents compiler from deduplicating across TUs.
- Suggested fix: Replace `static inline` with `inline` (or place in an anonymous namespace without `inline`).

### [SEV-3] Pointer hash `>> 4` discards information
- Location: `aterm_hash.h:43`
- Why this is a problem: Shifting `uintptr_t` right by 4 permanently discards the low 4 bits. While the comment explains the rationale (cache-friendly chaining), this produces poor distribution for open-addressing hash tables used elsewhere. The top bits (above virtual address space) are also always zero on 64-bit systems, effectively reducing hash entropy further.
- Evidence type: plausible-issue (documented performance trade-off in code comments)
- Efficiency impact: Potential clustering in hash tables that are not the primary chaining-based table.
- Suggested fix: No change if profiling confirms current behavior is optimal; otherwise consider a finalizer like `xorshift` on the result.

## Proposed Tests
- [ ] `aterm_hash_equal_args_arity_mismatch`: Instantiate `aterm_equals_finite<3>::operator()(term, symbol, a, b, c)` and verify correctness. Run under ASan to confirm no aliasing issues with the cast.
- [ ] `aterm_hash_dynamic_collision`: Create two terms with different arities (≥ 8) that hash to the same bucket in the dynamic pool. Verify equality comparison returns false without sanitizer violations.
- [ ] `aterm_hash_finite_all_arities`: For each N in {1..7}, verify `aterm_equals_finite<N>` variadic path returns correct results matching the array-based path for the same inputs.

## Open Questions
1. Is the `_aterm_appl<8>` hardcoding a leftover from when the pool had exactly 8 fixed-arity storages, or was it always intended as "large enough for any realistic arity"?
2. Was the arguments-before-function-symbol comparison order in `aterm_equals<DynamicNumberOfArguments>::operator()` validated by profiling, and does the dynamic pool's custom allocator guarantee sufficient trailing memory to make the out-of-bounds read safe in practice?
3. Is there a plan to migrate the flexible-array-member pattern to a more standard approach (e.g., `std::pmr` or a proper arena with explicit size tracking) to eliminate the structural reliance on UB?

## Copilot Raw Stdout

```text
● Explore: Explore aterm_hash dependencies
  └ Agent started in background with agent_id: agent-0. You can use read_agent to...

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_hash.h
  └ 409 lines read

● Explore: Find _term_appl typedef and ODR
  └ Agent started in background with agent_id: agent-1. You can use read_agent to...

● Grep "_term_appl" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ 14 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h lines 60-115
  └ 56 lines read

● Grep "class _aterm_appl|template.*_aterm_appl|struct _aterm_appl"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h)
  └ 3 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm.h lines 1-60
  └ 59 lines read

● Read: Explore agent (Explore aterm_hash dependencies)
  └ Completed

● Read: Explore agent (Find _term_appl typedef and ODR)
  └ Completed

● Grep "aterm_equals_finite.*operator\(\)|equal_args" in *.cpp (/workspaces/mCRL2/libraries/atermpp)
  └ No matches found

● Grep "aterm_equals_finite|aterm_hasher_finite" in *.h (/workspaces/mCRL2/libraries/atermpp)
  └ 8 lines found

● Grep "function_application_storage|aterm_equals_finite|aterm_hasher_finite"
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h)
  └ 11 lines found

● Grep "function_application_storage" in *.h (/workspaces/mCRL2/libraries/atermpp)
  └ 18 lines found

● Grep "aterm_pool_storage_implementation" in *.h
(/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail)
  └ 3 lines found

● Grep "\.emplace\(.*args\.\.\.|hasher_finite.*args|equals_finite.*args" in *.h
(/workspaces/mCRL2/libraries/atermpp)
  └ 2 lines found

● Grep "make_term_appl.*Args|variadic|args\.\.\." (/workspaces/mCRL2/libraries/atermpp/include/mcrl2
/atermpp/detail/aterm_pool_storage_implementation.h)
  └ 2 lines found

● Grep "make_term_appl\(.*result.*symbol.*args" in *.h (/workspaces/mCRL2/libraries/atermpp)
  └ No matches found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool_storage_implementation.h lines
  155-210
  └ 56 lines read

● Grep "emplace\(.*symbol.*arguments" (/workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/det
ail/aterm_pool_storage_implementation.h)
  └ 47 lines found

● Read libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h lines 25-40
  └ 16 lines read

Now I have all the context needed. Here's the review:

---

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/aterm_hash.h` (full file, clean review)
- Risk level: **Medium-High** — contains strict aliasing UB and a latent correctness hazard for arity > 7
- Overall verdict: The file is functional for current usage (arity ≤ 7) but harbors undefined behavior via hardcoded type casts and a risky argument-comparison-before-symbol-check pattern that could cause out-of-bounds reads on hash collision.

## Findings (ordered by severity)

### [SEV-1] Hardcoded `_aterm_appl<8>` in `equal_args` and variadic `aterm_equals_finite<N>::operator()`
- Location: `aterm_hash.h:376, 383, 393`
- Why this is a problem: `equal_args()` unconditionally uses `_aterm_appl<8>` regardless of the enclosing `aterm_equals_finite<N>` template parameter. The cast at line 393 (`static_cast<const _aterm_appl<8>&>(term)`) violates strict aliasing when the actual object was constructed as `_aterm_appl<N>` with N ≠ 8. Since `function_application_storage<N>` is instantiated for N=1..7, every instantiation performs a type-punning cast. Additionally, if the pool were ever extended to N ≥ 9, `arg(I)` would access index I ≥ 8 on an `_aterm_appl<8>` — reading past the declared array.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Any instantiation of `aterm_equals_finite<N>::operator()(...args)` for N != 8; currently all N in {0,1,2,3,4,5,6,7} trigger it.
  - Run command: `cmake -B build-asan -DCMAKE_BUILD_TYPE=Debug -DMCRL2_ENABLE_ADDRESSSANITIZER=ON -DMCRL2_ENABLE_GUI_TOOLS=OFF && cmake --build build-asan && ctest --test-dir build-asan -R atermpp --output-on-failure`
  - Observed result: Strict aliasing UB (not flagged by ASan but is UB per [basic.lval]). On optimized builds compilers may miscompile.
  - Expected result: Cast should target `_aterm_appl<N>` matching the template parameter.
- Efficiency impact: None at runtime (works by accident due to identical memory layout).
- C++20-specific note: Violation of strict aliasing ([basic.lval] p11). Could use `std::start_lifetime_as` in C++23 or `std::launder`/careful placement to make this legal.
- Suggested fix:
  ```cpp
  template <std::size_t N, std::size_t I = 0, typename... Tp>
  bool equal_args(const _aterm_appl<N>&, const Tp&...)
    requires(I == sizeof...(Tp))
  { return true; }

  template <std::size_t N, std::size_t I = 0, typename... Tp>
  bool equal_args(const _aterm_appl<N>& term, const Tp&... t)
    requires(I < sizeof...(Tp))
  {
    return term.arg(I) == std::get<I>(std::forward_as_tuple(t...)) && equal_args<N, I+1>(term, t...);
  }

  template<std::size_t N>
  template<typename ...Args>
  bool aterm_equals_finite<N>::operator()(const _aterm& term, const function_symbol& symbol, const Args&... args) const noexcept
  {
    return term.function() == symbol && equal_args<N>(static_cast<const _aterm_appl<N>&>(term), args...);
  }
  ```

### [SEV-2] Equality check compares arguments before function symbol — potential out-of-bounds on hash collision
- Location: `aterm_hash.h:289-312`
- Why this is a problem: `aterm_equals<N>::operator()(first, second)` uses `first.function().arity()` (when N = DynamicNumberOfArguments) to determine how many arguments of `second` to read. On a hash collision where `second` has fewer arguments than `first`, this reads past `second`'s allocated arguments. The function symbol mismatch is only detected at line 311 *after* the out-of-bounds loop. The `arbitrary_function_application_storage` (arity ≥ 8) uses this comparator.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Construct two terms with different arities (≥ 8) whose hash collides; force a lookup triggering this equality check.
  - Run command: `ctest --test-dir build-asan -R atermpp --output-on-failure` (with ASan; requires triggering the collision)
  - Observed result: (plausible) ASan heap-buffer-overflow or use-of-uninitialized if term allocations don't have sufficient trailing padding.
  - Expected result: Function symbol should be compared first; early-return on mismatch.
- Efficiency impact: Even without UB, comparing arguments before the function symbol wastes cycles on collisions.
- Suggested fix: Move `if (first.function() != second.function()) return false;` before the argument loop.

### [SEV-2] `aterm_equals<N>::operator()(term, symbol, arguments*)` ignores N template parameter
- Location: `aterm_hash.h:320-333`
- Why this is a problem: This overload unconditionally uses `symbol.arity()` (line 324) instead of the N-based logic (`(N == DynamicNumberOfArguments) ? symbol.arity() : N`) used by all sibling overloads. When N is fixed, the enclosing hash table guarantees all terms have arity N. Using `symbol.arity()` is redundant but not harmful *if* the caller always passes a symbol with arity N. However, it's inconsistent and means a mismatched symbol could cause out-of-bounds reads on `term`.
- Evidence type: static-proof (inconsistency with other overloads; lines 175, 197, 217, 299, 343 all use the N-based conditional)
- Suggested fix: Use `const std::size_t arity = (N == DynamicNumberOfArguments) ? symbol.arity() : N;` for consistency and defense in depth.

### [SEV-3] `aterm_int_hasher::operator()(_aterm_int)` creates unnecessary temporary
- Location: `aterm_hash.h:270-273`
- Why this is a problem: `return aterm_int_hasher()(term.value());` constructs a temporary `aterm_int_hasher` to call the `size_t` overload. This is functionally correct but confusing and marginally wasteful.
- Evidence type: static-proof (code inspection)
- Efficiency impact: Negligible (optimizer likely eliminates), but confusing.
- Suggested fix: `return (*this)(term.value());` or directly inline the hash computation.

### [SEV-3] `static inline` on `xorshift` creates per-TU internal-linkage copies
- Location: `aterm_hash.h:277`
- Why this is a problem: `static` at namespace scope gives internal linkage, meaning each translation unit gets its own copy. The `inline` is redundant with `static`. For a tiny helper, `inline` alone (external linkage, single definition) is preferred per modern C++ convention.
- Evidence type: static-proof (language semantics)
- Efficiency impact: Negligible code bloat; prevents compiler from deduplicating across TUs.
- Suggested fix: Replace `static inline` with `inline` (or place in an anonymous namespace without `inline`).

### [SEV-3] Pointer hash `>> 4` discards information
- Location: `aterm_hash.h:43`
- Why this is a problem: Shifting `uintptr_t` right by 4 permanently discards the low 4 bits. While the comment explains the rationale (cache-friendly chaining), this produces poor distribution for open-addressing hash tables used elsewhere. The top bits (above virtual address space) are also always zero on 64-bit systems, effectively reducing hash entropy further.
- Evidence type: plausible-issue (documented performance trade-off in code comments)
- Efficiency impact: Potential clustering in hash tables that are not the primary chaining-based table.
- Suggested fix: No change if profiling confirms current behavior is optimal; otherwise consider a finalizer like `xorshift` on the result.

## Proposed Tests
- [ ] `aterm_hash_equal_args_arity_mismatch`: Instantiate `aterm_equals_finite<3>::operator()(term, symbol, a, b, c)` and verify correctness. Run under ASan to confirm no aliasing issues with the cast.
- [ ] `aterm_hash_dynamic_collision`: Create two terms with different arities (≥ 8) that hash to the same bucket in the dynamic pool. Verify equality comparison returns false without sanitizer violations.
- [ ] `aterm_hash_finite_all_arities`: For each N in {1..7}, verify `aterm_equals_finite<N>` variadic path returns correct results matching the array-based path for the same inputs.

## Open Questions
1. Is the `_aterm_appl<8>` hardcoding a leftover from when the pool had exactly 8 fixed-arity storages, or was it always intended as "large enough for any realistic arity"?
2. Was the arguments-before-function-symbol comparison order in `aterm_equals<DynamicNumberOfArguments>::operator()` validated by profiling, and does the dynamic pool's custom allocator guarantee sufficient trailing memory to make the out-of-bounds read safe in practice?
3. Is there a plan to migrate the flexible-array-member pattern to a more standard approach (e.g., `std::pmr` or a proper arena with explicit size tracking) to eliminate the structural reliance on UB?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         6m 3s
Total session time:     5m 17s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         749.5k in, 13.2k out, 673.8k cached (Est. 1 Premium request)
 claude-haiku-4.5        324.7k in, 8.9k out, 268.9k cached (Est. 0 Premium requests)

```

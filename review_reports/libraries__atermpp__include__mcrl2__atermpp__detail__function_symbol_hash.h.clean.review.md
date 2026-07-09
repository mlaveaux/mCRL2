---
{
  "file": "libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_hash.h",
  "mode": "clean",
  "output_mode": "clean",
  "base_ref": null,
  "model": "claude-opus-4.6",
  "thinking": "on",
  "effort": "high",
  "generated_at_utc": "2026-07-09T10:17:37.609357+00:00",
  "dependency_count": 0,
  "copilot_return_code": 0,
  "session_limit_detected": false
}
---

# File Review

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_hash.h` — hash and equality functors for function symbol types
- Risk level: **Medium** — no immediate UB in normal builds, but portability issues, hash quality problems, and minor correctness concerns
- Overall verdict: **Functional but with quality gaps** — the hash functions work on common 64-bit Linux targets but have portability issues, use demonstrably weak combining, and ignore existing infrastructure (`hash_combine`) available in the same repository.

## Findings (ordered by severity)

### [SEV-2] Weak hash combining via plain XOR with small arity values
- Location: `function_symbol_hash.h:40` and `function_symbol_hash.h:77`
- Why this is a problem: `(h ^ f.arity())` is a poor hash combiner. Arity values are typically 0–5 in mCRL2 specifications. XOR with 0 is a no-op (all nullary symbols collide with their bare name hash). XOR with 1–5 only flips 1–3 low-order bits, causing systematic collision clusters. The codebase already provides `mcrl2::utilities::detail::hash_combine()` (boost-style golden-ratio mixer) in `libraries/utilities/include/mcrl2/utilities/hash_utility.h`.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Any symbol set with identical names but varying small arities (e.g., `f/0`, `f/1`, `f/2`, `f/3`)
  - Run command: N/A (hash quality analysis)
  - Observed result: For `hash("f") = H`, outputs are `H^0=H`, `H^1=H±1`, `H^2`, `H^3` — only low bits differ, all map to same bucket in power-of-2 tables unless the table is very small.
  - Expected result: Well-distributed hashes across the full `size_t` range.
- Efficiency impact: Increased collision rate in the function symbol pool (`unordered_set`), degrading O(1) lookup toward O(n) for workloads with many same-name/different-arity symbols.
- C++20-specific note: No language feature issue; this is an algorithm quality defect.
- Suggested fix: Replace `(h ^ f.arity())` with `mcrl2::utilities::detail::hash_combine(h, f.arity())` which is already available and used in `aterm_hash.h`.

### [SEV-2] Use of `std::uint64_t` instead of `std::uintptr_t` for pointer-to-integer cast
- Location: `function_symbol_hash.h:29`
- Why this is a problem: `reinterpret_cast<std::uint64_t>(ptr)` is implementation-defined when `sizeof(ptr) != sizeof(uint64_t)`. While this works on LP64 and LLP64 (Windows 64-bit), it is technically non-portable. The C++ standard provides `std::uintptr_t` specifically for this purpose (guaranteed to hold any pointer value). On ILP32 targets this would silently zero-extend, which works but is unidiomatic.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Compile on a 32-bit target
  - Run command: `cmake -S . -B build32 -DCMAKE_CXX_FLAGS=-m32 && cmake --build build32`
  - Observed result: Compiles with warnings about pointer size mismatch on strict compilers (`-Wpointer-to-int-cast` equivalent).
  - Expected result: Clean compilation on all supported targets.
- Efficiency impact: None on 64-bit.
- C++20-specific note: `std::uintptr_t` is the portable choice since C++11.
- Suggested fix: `return reinterpret_cast<std::uintptr_t>(f.m_function_symbol.get()) >> 5;`

### [SEV-3] Questionable `>> 5` shift assumption
- Location: `function_symbol_hash.h:26-29`
- Why this is a problem: The comment claims "Function symbols take 48 bytes" to justify a shift of 5 (dividing by 32). But: (1) `sizeof(_function_symbol)` is platform-dependent (std::string is 32 bytes on libstdc++/libc++ but 40 on some MSVC configurations, and the refcount field may be atomic); (2) the block allocator aligns to `alignof(_function_symbol)` not to 32 bytes; (3) if the size changes, the shift becomes wrong without compile-time protection. A `static_assert` on the size or using `alignof`-based shift would be more robust.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Build on platform where `sizeof(_function_symbol) != 48`
  - Run command: `static_assert(sizeof(atermpp::detail::_function_symbol) == 48)` in a test
  - Observed result: On MSVC with debug allocators or alternate STL implementations, size may differ.
  - Expected result: Hash should adapt to actual alignment/size or validate assumptions at compile time.
- Efficiency impact: If shift is too large, useful address bits are discarded and hash entropy is reduced.
- C++20-specific note: Could use `constexpr` shift computed from `alignof` for correctness.
- Suggested fix: Add `static_assert(sizeof(_function_symbol) >= 32, "...")` or compute shift as `constexpr int shift = std::countr_zero(alignof(_function_symbol));`.

### [SEV-3] Name shadowing in `function_symbol_hasher::operator()`
- Location: `function_symbol_hash.h:69`
- Why this is a problem: `const std::hash<_function_symbol> function_symbol_hasher;` declares a local variable with the same name as the enclosing struct. While valid C++, this shadows the struct name within the function body, making the code confusing and fragile for maintenance. A reader may mistake the local for a recursive call pattern.
- Evidence type: static-proof (code inspection)
- Evidence:
  - Reproducer/test name: N/A
  - Run command: Compile with `-Wshadow`
  - Observed result: Warning about shadowed name.
  - Expected result: Clean build under `-Wshadow`.
- Efficiency impact: None.
- C++20-specific note: N/A.
- Suggested fix: Rename local to `hasher` or simply `std::hash<_function_symbol>{}(symbol)`.

### [SEV-3] Missing `noexcept` on `std::hash` specializations
- Location: `function_symbol_hash.h:25` and `function_symbol_hash.h:38`
- Why this is a problem: The `std::hash<function_symbol>::operator()` and `std::hash<_function_symbol>::operator()` are not marked `noexcept`. Standard library implementations can optimize hash table operations when the hash functor is `noexcept`. The internal `function_symbol_hasher` correctly uses `noexcept`, but the `std::hash` specializations that it delegates to do not. Since `_function_symbol::name()` returns `const std::string&` (no allocation) and `std::hash<std::string>` is `noexcept` in all three major implementations, adding `noexcept` is safe and enables optimizations.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `static_assert(std::is_nothrow_invocable_v<std::hash<atermpp::function_symbol>, const atermpp::function_symbol&>)` — would fail.
  - Run command: N/A
  - Observed result: `is_nothrow_invocable` is false.
  - Expected result: Should be true for hash functions used in unordered containers.
- Efficiency impact: Minor; may prevent some container optimizations.
- C++20-specific note: Since C++17, containers may provide stronger exception guarantees when hash/equal are `noexcept`.
- Suggested fix: Add `noexcept` to both `operator()` overloads in the `std::hash` specializations.

### [SEV-3] Missing `[[nodiscard]]` on hash function operators
- Location: `function_symbol_hash.h:25, 38, 55-56`
- Why this is a problem: Hash function results are pure computations; discarding the return value is always a bug. `[[nodiscard]]` catches accidental misuse at compile time.
- Evidence type: static-proof (best practice)
- Efficiency impact: None.
- C++20-specific note: `[[nodiscard]]` is standard since C++17 and idiomatic in C++20.
- Suggested fix: Add `[[nodiscard]]` to all `operator()` overloads.

## Proposed Tests
- [ ] `function_symbol_hash_distribution_test`: Measure collision rate of `hash<_function_symbol>` for symbol set `{"f"/0, "f"/1, ..., "f"/255}`. Expected: collision rate < 5% in a 512-bucket table. Failure signal: collision rate > 50% with plain XOR, significantly better with `hash_combine`.
- [ ] `function_symbol_hash_static_assert`: Add `static_assert(sizeof(_function_symbol) >= 32)` to validate the `>> 5` assumption at compile time. Failure signal: compilation failure if size changes.
- [ ] `function_symbol_hash_noexcept_check`: `static_assert(std::is_nothrow_invocable_v<std::hash<function_symbol>, const function_symbol&>)` — verifies noexcept propagation. Currently fails.

## Open Questions
1. Is there a reason `hash_combine` from `mcrl2/utilities/hash_utility.h` was not used here, given it's used in the closely related `aterm_hash.h`? Was this an oversight or a deliberate performance tradeoff?
2. Has the 48-byte size assumption been validated on MSVC (where `sizeof(std::string)` can differ)? Are there CI builds that would catch a violation?
3. The `std::hash<function_symbol>` hashes by pointer identity (address of the pooled `_function_symbol`). This is correct only because symbols are interned. Is there documentation or an assertion ensuring that no `function_symbol` can ever refer to a non-interned `_function_symbol`?

## Copilot Raw Stdout

```text
● Explore: Gather context for hash review
  ● Search for hash function tests
    $ grep -r "std::hash.*function_symbol\|unordered_set.*function_symbol"
    /workspaces/mCRL2/libraries/atermpp/test/ --include="*.cpp" --include="*.h" 2>/dev/null | head
    -20
    └ 1 line...
  ● Check the ref type definition
    $ grep -A 5 "shared_reference<const _function_symbol>"
    /workspaces/mCRL2/libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h
    └ 7 lines...

✗ Check function symbol size assumptions
  $ cd /workspaces/mCRL2 && grep -n "sizeof\|alignof"
  libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol.h 2>/dev/null; echo "---";
  python3 -c "
  # Check _function_symbol size: shared_reference_counted has a reference count (size_t or
  atomic<size_t>), ...
  Permission denied and could not request permission from user

● Read libraries/utilities/include/mcrl2/utilities/shared_reference.h lines 1-50
  └ 50 lines read

● Check hash-related patterns in atermpp
  $ cd /workspaces/mCRL2 && grep -n "uintptr_t\|uint64_t\|reinterpret_cast"
  libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_hash.h; echo "---"; grep -rn
  "hash_combine\|hash_combine_cheap" libraries/atermpp/ | head -10
  └ 4 lines...

● Check shared_reference::get() return type
  $ cd /workspaces/mCRL2 && grep -n "get()"
  libraries/utilities/include/mcrl2/utilities/shared_reference.h | head -5; echo "---"; grep -n
  "class shared_reference" libraries/utilities/include/mcrl2/utilities/shared_reference.h | head -5
  └ 8 lines...

● Check pointer type used in shared_reference
  $ cd /workspaces/mCRL2 && grep -n "tagged_pointer"
  libraries/utilities/include/mcrl2/utilities/shared_reference.h | head -5; grep -n
  "T\*\|pointer\|get" libraries/utilities/include/mcrl2/utilities/tagged_pointer.h | head -20
  └ 23 lines...

● Read libraries/utilities/include/mcrl2/utilities/shared_reference.h lines 193-210
  └ 17 lines read

● Read libraries/utilities/include/mcrl2/utilities/tagged_pointer.h lines 63-110
  └ 48 lines read

● Read libraries/utilities/include/mcrl2/utilities/tagged_pointer.h lines 110-140
  └ 31 lines read

● Read libraries/utilities/include/mcrl2/utilities/tagged_pointer.h lines 140-160
  └ 21 lines read

Now I have all the context needed. Let me write the review.

# Review Findings

## Summary
- Scope reviewed: `libraries/atermpp/include/mcrl2/atermpp/detail/function_symbol_hash.h` — hash and equality functors for function symbol types
- Risk level: **Medium** — no immediate UB in normal builds, but portability issues, hash quality problems, and minor correctness concerns
- Overall verdict: **Functional but with quality gaps** — the hash functions work on common 64-bit Linux targets but have portability issues, use demonstrably weak combining, and ignore existing infrastructure (`hash_combine`) available in the same repository.

## Findings (ordered by severity)

### [SEV-2] Weak hash combining via plain XOR with small arity values
- Location: `function_symbol_hash.h:40` and `function_symbol_hash.h:77`
- Why this is a problem: `(h ^ f.arity())` is a poor hash combiner. Arity values are typically 0–5 in mCRL2 specifications. XOR with 0 is a no-op (all nullary symbols collide with their bare name hash). XOR with 1–5 only flips 1–3 low-order bits, causing systematic collision clusters. The codebase already provides `mcrl2::utilities::detail::hash_combine()` (boost-style golden-ratio mixer) in `libraries/utilities/include/mcrl2/utilities/hash_utility.h`.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Any symbol set with identical names but varying small arities (e.g., `f/0`, `f/1`, `f/2`, `f/3`)
  - Run command: N/A (hash quality analysis)
  - Observed result: For `hash("f") = H`, outputs are `H^0=H`, `H^1=H±1`, `H^2`, `H^3` — only low bits differ, all map to same bucket in power-of-2 tables unless the table is very small.
  - Expected result: Well-distributed hashes across the full `size_t` range.
- Efficiency impact: Increased collision rate in the function symbol pool (`unordered_set`), degrading O(1) lookup toward O(n) for workloads with many same-name/different-arity symbols.
- C++20-specific note: No language feature issue; this is an algorithm quality defect.
- Suggested fix: Replace `(h ^ f.arity())` with `mcrl2::utilities::detail::hash_combine(h, f.arity())` which is already available and used in `aterm_hash.h`.

### [SEV-2] Use of `std::uint64_t` instead of `std::uintptr_t` for pointer-to-integer cast
- Location: `function_symbol_hash.h:29`
- Why this is a problem: `reinterpret_cast<std::uint64_t>(ptr)` is implementation-defined when `sizeof(ptr) != sizeof(uint64_t)`. While this works on LP64 and LLP64 (Windows 64-bit), it is technically non-portable. The C++ standard provides `std::uintptr_t` specifically for this purpose (guaranteed to hold any pointer value). On ILP32 targets this would silently zero-extend, which works but is unidiomatic.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: Compile on a 32-bit target
  - Run command: `cmake -S . -B build32 -DCMAKE_CXX_FLAGS=-m32 && cmake --build build32`
  - Observed result: Compiles with warnings about pointer size mismatch on strict compilers (`-Wpointer-to-int-cast` equivalent).
  - Expected result: Clean compilation on all supported targets.
- Efficiency impact: None on 64-bit.
- C++20-specific note: `std::uintptr_t` is the portable choice since C++11.
- Suggested fix: `return reinterpret_cast<std::uintptr_t>(f.m_function_symbol.get()) >> 5;`

### [SEV-3] Questionable `>> 5` shift assumption
- Location: `function_symbol_hash.h:26-29`
- Why this is a problem: The comment claims "Function symbols take 48 bytes" to justify a shift of 5 (dividing by 32). But: (1) `sizeof(_function_symbol)` is platform-dependent (std::string is 32 bytes on libstdc++/libc++ but 40 on some MSVC configurations, and the refcount field may be atomic); (2) the block allocator aligns to `alignof(_function_symbol)` not to 32 bytes; (3) if the size changes, the shift becomes wrong without compile-time protection. A `static_assert` on the size or using `alignof`-based shift would be more robust.
- Evidence type: plausible-issue
- Evidence:
  - Reproducer/test name: Build on platform where `sizeof(_function_symbol) != 48`
  - Run command: `static_assert(sizeof(atermpp::detail::_function_symbol) == 48)` in a test
  - Observed result: On MSVC with debug allocators or alternate STL implementations, size may differ.
  - Expected result: Hash should adapt to actual alignment/size or validate assumptions at compile time.
- Efficiency impact: If shift is too large, useful address bits are discarded and hash entropy is reduced.
- C++20-specific note: Could use `constexpr` shift computed from `alignof` for correctness.
- Suggested fix: Add `static_assert(sizeof(_function_symbol) >= 32, "...")` or compute shift as `constexpr int shift = std::countr_zero(alignof(_function_symbol));`.

### [SEV-3] Name shadowing in `function_symbol_hasher::operator()`
- Location: `function_symbol_hash.h:69`
- Why this is a problem: `const std::hash<_function_symbol> function_symbol_hasher;` declares a local variable with the same name as the enclosing struct. While valid C++, this shadows the struct name within the function body, making the code confusing and fragile for maintenance. A reader may mistake the local for a recursive call pattern.
- Evidence type: static-proof (code inspection)
- Evidence:
  - Reproducer/test name: N/A
  - Run command: Compile with `-Wshadow`
  - Observed result: Warning about shadowed name.
  - Expected result: Clean build under `-Wshadow`.
- Efficiency impact: None.
- C++20-specific note: N/A.
- Suggested fix: Rename local to `hasher` or simply `std::hash<_function_symbol>{}(symbol)`.

### [SEV-3] Missing `noexcept` on `std::hash` specializations
- Location: `function_symbol_hash.h:25` and `function_symbol_hash.h:38`
- Why this is a problem: The `std::hash<function_symbol>::operator()` and `std::hash<_function_symbol>::operator()` are not marked `noexcept`. Standard library implementations can optimize hash table operations when the hash functor is `noexcept`. The internal `function_symbol_hasher` correctly uses `noexcept`, but the `std::hash` specializations that it delegates to do not. Since `_function_symbol::name()` returns `const std::string&` (no allocation) and `std::hash<std::string>` is `noexcept` in all three major implementations, adding `noexcept` is safe and enables optimizations.
- Evidence type: static-proof
- Evidence:
  - Reproducer/test name: `static_assert(std::is_nothrow_invocable_v<std::hash<atermpp::function_symbol>, const atermpp::function_symbol&>)` — would fail.
  - Run command: N/A
  - Observed result: `is_nothrow_invocable` is false.
  - Expected result: Should be true for hash functions used in unordered containers.
- Efficiency impact: Minor; may prevent some container optimizations.
- C++20-specific note: Since C++17, containers may provide stronger exception guarantees when hash/equal are `noexcept`.
- Suggested fix: Add `noexcept` to both `operator()` overloads in the `std::hash` specializations.

### [SEV-3] Missing `[[nodiscard]]` on hash function operators
- Location: `function_symbol_hash.h:25, 38, 55-56`
- Why this is a problem: Hash function results are pure computations; discarding the return value is always a bug. `[[nodiscard]]` catches accidental misuse at compile time.
- Evidence type: static-proof (best practice)
- Efficiency impact: None.
- C++20-specific note: `[[nodiscard]]` is standard since C++17 and idiomatic in C++20.
- Suggested fix: Add `[[nodiscard]]` to all `operator()` overloads.

## Proposed Tests
- [ ] `function_symbol_hash_distribution_test`: Measure collision rate of `hash<_function_symbol>` for symbol set `{"f"/0, "f"/1, ..., "f"/255}`. Expected: collision rate < 5% in a 512-bucket table. Failure signal: collision rate > 50% with plain XOR, significantly better with `hash_combine`.
- [ ] `function_symbol_hash_static_assert`: Add `static_assert(sizeof(_function_symbol) >= 32)` to validate the `>> 5` assumption at compile time. Failure signal: compilation failure if size changes.
- [ ] `function_symbol_hash_noexcept_check`: `static_assert(std::is_nothrow_invocable_v<std::hash<function_symbol>, const function_symbol&>)` — verifies noexcept propagation. Currently fails.

## Open Questions
1. Is there a reason `hash_combine` from `mcrl2/utilities/hash_utility.h` was not used here, given it's used in the closely related `aterm_hash.h`? Was this an oversight or a deliberate performance tradeoff?
2. Has the 48-byte size assumption been validated on MSVC (where `sizeof(std::string)` can differ)? Are there CI builds that would catch a violation?
3. The `std::hash<function_symbol>` hashes by pointer identity (address of the pooled `_function_symbol`). This is correct only because symbols are interned. Is there documentation or an assertion ensuring that no `function_symbol` can ever refer to a non-interned `_function_symbol`?


```

## Copilot Raw Stderr

```text

Total usage est:        1 Premium request
API time spent:         3m 47s
Total session time:     4m 2s
Total code changes:     +0 -0
Breakdown by AI model:
 claude-opus-4.6         308.5k in, 8.3k out, 282.3k cached (Est. 1 Premium request)
 claude-haiku-4.5        232.9k in, 5.2k out, 209.9k cached (Est. 0 Premium requests)

```

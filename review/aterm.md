# atermpp library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `atermpp` library is the foundation everything else stands on: it implements
*aterms* (maximally-shared, immutable, hash-consed terms), the reference-counting
garbage collector that reclaims them, the per-thread protection sets that keep live
terms reachable across a concurrent GC, the function-symbol pool, the binary/textual
term serialisation, and the GC-aware standard-container wrappers. It is the most
thread-safety- and performance-critical code in the toolset and is written by the
core developers, so the bar for a real defect here is high — several alarming-looking
constructs turned out to be deliberate and correct. The GC/pool internals, the
thread-protection machinery, the I/O, and the containers were reviewed file-by-file
with two assisting sweeps; the claims below were each checked against the source, and
a number of subagent "critical race" reports were **dismissed** (see the note at the
end of the section).

## 1. Thread-safety

### 1.1 Non-RAII manual locking in `index_traits` and the function-symbol pool

[index_traits.h](libraries/atermpp/include/mcrl2/atermpp/detail/index_traits.h#L88)

The variable/operator index assignment locks a per-type mutex by hand around map
operations that can throw:

```cpp
if constexpr (mcrl2::utilities::detail::GlobalThreadSafe) { variable_mutex<Variable, KeyType>().lock(); }
auto& m = variable_index_map<Variable, KeyType>();
auto i = m.find(x);
...
m[x] = value;                       // may rehash → may throw bad_alloc
...
if constexpr (mcrl2::utilities::detail::GlobalThreadSafe) { variable_mutex<Variable, KeyType>().unlock(); }
```

If any of the intervening `unordered_map`/stack operations throws (a rehash
allocation failure on `m[x] = value`, for instance), the mutex is never released and
every later index assignment of that type deadlocks. `erase`
([:122](libraries/atermpp/include/mcrl2/atermpp/detail/index_traits.h#L122)) has the
same shape, and `function_symbol_pool::create_helper`
([function_symbol_pool.cpp](libraries/atermpp/source/function_symbol_pool.cpp#L34))
locks `m_mutex` by hand and then calls `std::stoul` and map operations before a
manual unlock. These are the same non-RAII hazard flagged in the LTS, pbes and core
libraries; a `std::lock_guard`/`shared_guard` (which the rest of `atermpp` uses
almost everywhere) fixes them. They only apply to `GlobalThreadSafe` (multithreaded)
builds.

## 2. Robustness — deserialising untrusted input

### 2.1 Binary term loading indexes tables with unvalidated stream values

[aterm_io_binary.cpp](libraries/atermpp/source/aterm_io_binary.cpp#L268)

The binary reader looks up function symbols and argument terms with indices taken
directly from the input stream, using `operator[]`:

```cpp
function_symbol symbol = m_function_symbols[m_stream->read_bits(function_symbol_index_width())];
...
std::vector<aterm> arguments(symbol.arity());
for (std::size_t argument = 0; argument < symbol.arity(); ++argument)
{
  arguments[argument] = m_terms[m_stream->read_bits(term_index_width())];   // :287
}
```

A malformed or hostile `.lps`/`.pbes`/`.lts`/`.aterm` file can supply an index past
the end of `m_function_symbols`/`m_terms`, producing an out-of-bounds vector read.
Calibrated honestly this is a crash / denial-of-service (an OOB *read* of a
`std::vector`), not a memory-corruption/RCE primitive — but the mCRL2 tools do load
user-supplied files, so it should fail cleanly. Switching the two lookups to `.at()`
turns the corruption into a thrown `std::out_of_range` (which the surrounding loaders
already convert to a `runtime_error`). Relatedly, `symbol.arity()` is read from the
stream and used to size `std::vector<aterm> arguments(symbol.arity())`
([:285](libraries/atermpp/source/aterm_io_binary.cpp#L285)), so a hostile arity is a
large-allocation DoS.

### 2.2 Textual term loading: sign-cast and unbounded recursion

[aterm_io_text.cpp](libraries/atermpp/source/aterm_io_text.cpp#L200)

`parse_aterm_int` does `aterm_int(static_cast<std::size_t>(atol(number.data())))`: a
negative literal in the text format silently becomes `SIZE_MAX` rather than being
rejected, and `atol` performs no overflow detection (the 32-byte buffer is bounded,
so there is no overflow, only a wrong value). Separately, `parse_aterm` and
`parse_aterm_list` are mutually recursive with no depth bound, so a deeply nested
textual list overflows the stack. The textual format is rarely the untrusted-input
vector that the binary format is, so this is lower priority than 2.1.

## 3. Maintainability / style

- **Unreachable duplicated code.** `markable_aterm`'s assignment operators contain a
  dead second copy of their body, e.g.
  [aterm_container.h](libraries/atermpp/include/mcrl2/atermpp/detail/aterm_container.h#L125):
  ```cpp
  markable_aterm& operator=(const T& other) noexcept
  {
    m_t = other;
    return *this;
    m_t = other;   // unreachable
    return *this;  // unreachable
  }
  ```
  Harmless (the first `return` wins) and the same copy-paste artifact recurs in the
  sibling operators, but it is exactly what `-Wunreachable-code` exists to catch.
- **Acknowledged `reinterpret_cast` debt.** Several `aterm`↔`aterm_core` casts
  (e.g. [aterm_pool.h](libraries/atermpp/include/mcrl2/atermpp/detail/aterm_pool.h#L118)
  `empty_list()`) are layout-compatible by design and each already carries a `// TODO
  remove this reinterpret_cast` note; they are tech-debt to retire, not bugs.
- **Debug-only lock invariant.** `aterm_implementation.h` enforces the
  shared-lock-held precondition with `assert(detail::g_thread_term_pool().is_shared_locked())`,
  which disappears under `-DNDEBUG`; a violated locking contract then yields silent UB
  in release rather than a diagnostic.

## 4. Things that look good

- **GC marking is iterative, not recursive.** The mark phase uses an explicit
  `m_todo` work stack, so collecting a deeply nested term cannot overflow the stack —
  the kind of bound the recursive solvers elsewhere lack.
- **Careful cross-thread lifetime handoff.** `~thread_aterm_pool` transfers a worker's
  still-live variables and containers to the main-thread pool *before* unregistering,
  precisely so that static-local aterms first constructed on a worker thread are not
  reclaimed by a later GC — a subtle initialisation-order hazard handled deliberately
  and well-commented.
- **Containers cooperate with the collector.** The `vector`/`deque`/`unordered_set`/
  `unordered_map` wrappers acquire the thread pool's `shared_guard`/lock on mutating
  operations so their elements stay protected across a concurrent GC, and `shared_guard`
  RAII is used on the overwhelming majority of paths (the manual locks in 1.1 are the
  exception, not the rule).
- **I/O layer validates its own stream.** The underlying `bitstream` checks
  `good()`/`fail()` and throws, and the variable-size integer decoder has an explicit
  size bound and throws on a malformed encoding.

> Note: three subagent-reported "critical" issues were **checked and dismissed** as
> false positives: (a) the "data races from reading `std::atomic` members without an
> explicit `.load()`" — an implicit `std::atomic` conversion *is* a valid (seq-cst)
> atomic load, not a race; (b) the "TOCTOU on `g_main_thread_pool`" initialisation —
> it is a documented write-once performed by the main thread before any worker is
> spawned, and thread creation provides the necessary happens-before; and (c) the
> placement-`new`/`reinterpret_cast` "memory-safety" reports — these are the standard,
> intentional in-place term-construction and layout-compatible cast patterns of a
> hash-consing term library.

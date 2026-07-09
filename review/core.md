# core library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `core` library sits just above `atermpp`: it defines `identifier_string` (the
interned-string aterm used for every name in the toolset), the function-symbol
registry for the abstract syntax, the generic `builder`/`traverser` CRTP
infrastructure, the printing framework, and the wrapper around the third-party
*dparser* GLR parser that every `.mcrl2`/`.lps`/formula front end is built on. Most
of the type and visitor machinery is generated; the review focused file-by-file on
the hand-written glue — `dparser.{h,cpp}`, `identifier_string.h`, `load_aterm.h`,
`parser_utility.h`, `core.cpp`, and `detail/function_symbols.h`.

## 1. Thread-safety

### 1.1 Non-RAII mutex around the `DataAppl` function-symbol cache

[function_symbols.h](libraries/core/include/mcrl2/core/detail/function_symbols.h#L38)

`function_symbol_DataAppl_helper` grows a lazily-populated cache of `DataAppl`
function symbols under a manual lock:

```cpp
static std::mutex mutex;
i -= DataApplFixed;
mutex.lock();
do {
  function_symbols_DataAppl.push_back(std::make_unique<atermpp::function_symbol>(...));
} while (i >= function_symbols_DataAppl.size());
const atermpp::function_symbol& result = *function_symbols_DataAppl[i];
mutex.unlock();
```

The `std::make_unique`/`push_back` between `lock()` and `unlock()` can throw
`bad_alloc`; the mutex is then never released and every subsequent term
construction of a high-arity `DataAppl` deadlocks. This is the same non-RAII
locking hazard flagged in the LTS and pbes libraries — a `std::lock_guard` fixes it.
Two related accesses are documented-latent rather than active races: the
`i < DataApplFixed` fast path guards initialisation with a bare `static bool
initialised` (the comment concedes "this is not thread safe, but applications are
made during the initialisation process"), and
[gsIsDataAppl_no_check](libraries/core/include/mcrl2/core/detail/function_symbols.h#L90)
reads `*function_symbols_DataAppl[i]` **without** taking the mutex, which races the
vector's control block against a concurrent `push_back` reallocation. In practice
all needed arities are interned during single-threaded setup, so these are latent;
the non-RAII lock is the one to fix.

### 1.2 dparser error-message counters are unsynchronised global statics

[dparser.h](libraries/core/include/mcrl2/core/dparser.h#L36)

The diagnostic throttle is a process-global template static
(`dparser_error_message_count<T>::value`/`max_value`) mutated through free
`increment`/`reset`/`set_dparser_max_error_message_count` functions with no
synchronisation, and the `parser` constructor writes the global max on every
construction. Two parsers running on different threads would corrupt each other's
error counts. Parsing is assumed to be a non-concurrent setup activity, so this is
latent, but the shared mutable state is worth noting now that the toolset is
multithreaded by default.

## 2. Memory-safety

### 2.1 `parse_node`/`parser` own raw pointers but keep default copy semantics

[dparser.h](libraries/core/include/mcrl2/core/dparser.h#L78)

`parse_node` holds a raw `D_ParseNode* node` plus an optional `D_Parser* parser`,
and its destructor frees the node when `parser` is non-null:

```cpp
parse_node::~parse_node()
{
  if (parser) { free_D_ParseNode(parser, node); }
}
```

The struct has no user-declared copy/move operations, so it is freely copyable.
Copying the **root** node returned by `parser::parse` (the only `parse_node` with a
non-null `parser`) yields two wrappers that both free the same `D_ParseNode` —
a double-free. Child nodes from `child()`/`find_in_tree()` set `parser == nullptr`
so copying them is harmless, which is why this has not bitten in practice, but the
type is move-only in spirit and should `= delete` its copy operations (and add a
move) to make that explicit. The `parser` struct has the identical shape (it owns
`m_parser`, frees it in `~parser`, and keeps default copy), with the same latent
double-free.

## 3. Maintainability / style

- **Diagnostics bypass `mCRL2log`.** The debug helpers in `dparser.cpp` write
  directly to `std::cout`
  ([print_tree](libraries/core/source/dparser.cpp#L285),
  [announce](libraries/core/source/dparser.cpp#L302)) and `std::clog`
  ([parser_table::print](libraries/core/source/dparser.cpp#L192)); they are
  developer-only, but the rest of the toolset routes through `mCRL2log`.
  [print_tree](libraries/core/source/dparser.cpp#L283) also loops `for (i = 0; i <=
  child_count(); i++)`, reading one child past the end — benign only because
  `d_get_child` returns `nullptr` out of range and the `if (node)` guard absorbs it.
- **Reaching into dparser internals.** `dparser.cpp` `reinterpret_cast`s
  `D_Parser*` to the library-internal `Parser*`
  ([:434](libraries/core/source/dparser.cpp#L434),
  [:468](libraries/core/source/dparser.cpp#L468)) and reimplements the
  `D_ParseNode_to_PNode` pointer arithmetic
  ([:27](libraries/core/source/dparser.cpp#L27)) to recover error context. This is
  inherent to the third-party integration and documented, but it couples the code
  to dparser's private layout.
- **Placeholder briefs.** `dparser.h`, `dparser.cpp` and `add_binding.h` still carry
  `/// \brief add your file description here.`

## 4. Things that look good

- **Exemplary load error handling.**
  [load_aterm.h](libraries/core/include/mcrl2/core/load_aterm.h#L40) wraps both the
  binary and textual aterm readers in a `try`/`catch` that rethrows a
  `mcrl2::runtime_error` annotated with the format and source — the pattern the I/O
  reviews elsewhere asked for.
- **Thoughtful function-symbol cache.** The fixed `std::array<…, 100>` for common
  arities plus a `std::vector<std::unique_ptr<…>>` for the rest is a deliberate
  design: the `unique_ptr` indirection keeps returned `function_symbol` references
  stable across vector reallocation, with a comment explaining why a `deque` was
  rejected on performance grounds.
- **Carefully-documented disambiguation.** `ambiguity_fn` in `dparser.cpp` resolves
  the grammar's `RegFrm`/`ActFrm` and operator-precedence ambiguities with extensive
  comments and falls back to throwing a `runtime_error` listing the candidate parses
  rather than silently picking one.
- **Bounds-checked table access.** `parser_table::symbol_name` validates the symbol
  index and throws on overflow, and `parse` frees the partial tree and throws on a
  syntax error rather than returning a half-built node.
- **Clean core type.** `identifier_string` is a thin, correct alias of
  `atermpp::aterm_string` with a matching `std::hash` specialisation.

# symbolic library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `symbolic` library implements symbolic (LDD-based) state-space reachability on top
of the Sylvan decision-diagram package: the transition-relation `summand_group`s, the
bidirectional data-expression↔index maps, the LDD serialisation, and the reachability
driver. The whole library is gated behind `MCRL2_ENABLE_SYLVAN`, so it is an optional
component. It is compact and notably clean; it was read file-by-file.

## 1. Robustness — deserialising untrusted input

### 1.1 The LDD reader indexes its node table with unvalidated stream values

[ldd_stream.cpp](libraries/symbolic/source/ldd_stream.cpp#L181)

`binary_ldd_istream::get` rebuilds each LDD node from indices read straight from the
stream and looks them up with `operator[]`:

```cpp
std::size_t index = m_stream->read_bits(ldd_index_width());
return m_nodes[index];
...
std::size_t down_index  = m_stream->read_bits(ldd_index_width(true));
std::size_t right_index = m_stream->read_bits(ldd_index_width(true));
ldd down  = m_nodes[down_index];
ldd right = m_nodes[right_index];
```

`m_nodes` is a `std::vector` grown by `emplace_back`, so an out-of-range index in a
malformed or hostile `.ldd` file is an out-of-bounds read (crash / DoS, not a
corruption primitive — the same class and calibration as the `atermpp` binary reader,
finding 2.1). The fix is again `.at()`. To its credit the reader *does* validate the
`BLF_MAGIC` and `BLF_VERSION` header up front
([:162](libraries/symbolic/source/ldd_stream.cpp#L162),
[:168](libraries/symbolic/source/ldd_stream.cpp#L168)) and throws on a mismatch, so it
cleanly rejects non-LDD or version-incompatible files; only the per-node indices are
unchecked.

## 2. Things that look good

- **Iterative LDD traversal.** The serialiser's `node_iterator` walks the diagram
  bottom-up with an explicit `std::stack` rather than native recursion, so writing a
  very deep LDD cannot overflow the stack — the same discipline as `atermpp`'s
  iterative GC marker.
- **Header-validated I/O.** The magic/version handshake makes the binary format
  self-describing and rejects incompatible inputs with a clear `runtime_error`.
- **Clean core data structures.** `data_expression_index` is a tidy bidirectional
  `indexed_set` map (with a debug `assert` that inserted keys match the declared
  sort), and `summand_group` is a well-organised read/write transition relation whose
  bitset accesses are bounded by the parameter count.
- **Tidy throughout.** No placeholder briefs, no stray `std::cout`/`cerr`/`clog`, and
  no `assert(0)`. The one `reinterpret_cast<Context*>(context)` in
  [symbolic_reachability.h](libraries/symbolic/include/mcrl2/symbolic/symbolic_reachability.h#L125)
  is the idiomatic recovery of a typed context from Sylvan's C-style `void*` callback,
  which is unavoidable with that API — not a defect.

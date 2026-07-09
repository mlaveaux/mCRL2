# process library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `process` library holds the process-algebra AST (`process_expression` and its
~25 node types, `process_equation`, `process_specification`), the alphabet
operations and alphabet reduction (`block`/`hide`/`rename`/`comm`/`allow` pushing),
linearity/guardedness checks, equation cleanup (eliminate trivial/unused/duplicate
equations, SCC decomposition), parsing/printing/type-checking, and the generated
builder/traverser/replace infrastructure. It sits directly below `lps` (reviewed
above) and on top of `data`. Of ~16.5k lines, three files are almost entirely
generated ([traverser.h](libraries/process/include/mcrl2/process/traverser.h),
[builder.h](libraries/process/include/mcrl2/process/builder.h),
[process_expression.h](libraries/process/include/mcrl2/process/process_expression.h));
the hand-written algorithmic code is in the alphabet, elimination, linearity, and
I/O headers.

## 1. Correctness / robustness

### 1.1 Tarjan SCC is O(V·E) and recursive

[process_variable_strongly_connected_components.h](libraries/process/include/mcrl2/process/process_variable_strongly_connected_components.h#L57)

`tarjan_scc_algorithm::strongconnect` scans the **entire** edge list on every
vertex to find that vertex's successors:

```cpp
for (const edge& e: E)            // :66 — all edges, every call
{
  ...
  if (u != v) { continue; }       // skip edges not leaving v
  ...
}
```

That turns Tarjan's linear algorithm into `O(V·E)`; an adjacency list (build once,
index by source) restores `O(V+E)`. Separately, `strongconnect` recurses once per
tree edge, so the recursion depth equals the longest path in the process-dependency
graph — a deep or long chain of equations can blow the call stack. The classic
fix is an explicit work stack. The cache header even hints at the cost:
[pcrl_equation_cache.h](libraries/process/include/mcrl2/process/detail/pcrl_equation_cache.h#L23)
notes the alphabet mapping "can probably be computed more efficiently using the SCC
decomposition".

### 1.2 `remove_subsets` is order-dependent and incomplete

[multi_action_name.h](libraries/process/include/mcrl2/process/multi_action_name.h#L185)

```cpp
for (const multi_action_name& alpha: A)
  if (!includes(result, alpha))   // :190 — only checks already-inserted (smaller) elements
    result.insert(alpha);
```

The function promises to "remove elements of A that are a subset of another
element", but it only skips `alpha` when an **already-inserted** element is a
superset. `A` is a `std::set` iterated in ascending order, and a prefix-subset such
as `{a}` sorts *before* its superset `{a,b}`, so `{a}` is inserted first and never
removed. The result therefore keeps subsets the contract says to drop. It is not a
soundness bug — `remove_subsets` feeds `alphabet_operations::subsets`
([allow_set.h](libraries/process/include/mcrl2/process/allow_set.h)) where subsets
are implied anyway — but it defeats the intended compaction and leaves the
`allow_set` non-canonical, which can cause cache misses / missed fixpoint detection
in alphabet reduction. A correct version iterates largest-first or does the
pairwise check.

### 1.3 `eliminate_trivial_equations` can loop forever on a chain into a cycle

[eliminate_trivial_equations.h](libraries/process/include/mcrl2/process/eliminate_trivial_equations.h#L220)

`compute_chains` follows the `edges` map from each chain source:

```cpp
auto i = edges.find(source);
while (i != edges.end())          // :240
{
  const process_identifier& target = i->second;
  chain.push_back(target);
  i = edges.find(target);
}
```

The struct documents the precondition "there are no loops P1() = P2() = … = P1()",
but it is unguarded. A trivial-equation *tail* feeding a cycle
(`P0 = P1`, `P1 = P2`, `P2 = P1`) makes a `P0` chain walk `P1→P2→P1→…` forever
(those cycle nodes are not chain sources, so they are never pruned). A visited-set
guard would turn the hang into a no-op or a diagnostic.

### 1.4 `eliminate_unused_equations` dereferences a possibly-null equation pointer

[eliminate_unused_equations.h](libraries/process/include/mcrl2/process/eliminate_unused_equations.h#L113)

```cpp
new_equations.push_back(*equation_index[P]);
```

`equation_index` is a `std::map<process_identifier, process_equation*>`. If `P` is
reachable (from `init` or another rhs) but has **no** defining equation, `operator[]`
inserts a null `process_equation*` and the code dereferences it — a crash. This
only happens on an ill-formed/non-type-checked specification, but the elimination
runs without that guarantee; a `find` + explicit check would fail fast instead.

### 1.5 Unchecked lockstep iteration in the linearity check

[is_linear.h](libraries/process/include/mcrl2/process/is_linear.h#L56)

`detail::check_process_instance` walks the formal parameters and actual parameters
together but bounds only the first:

```cpp
for (; i != v.end(); ++i, ++j)    // :56 — j (actual params) is never bounded
  if (i->sort() != j->sort()) ...
```

If a `process_instance` carries fewer actual parameters than the equation's formal
parameters, `j` runs past `e.end()` (UB). The same class of precondition-trusting
lockstep loop was flagged in the LPS review (`multi_action.h`). Type-checked input
keeps the counts equal, but the check is also called directly on untyped fragments.

### 1.6 Unchecked file open in the process reader

[detail/process_io.h](libraries/process/include/mcrl2/process/detail/process_io.h#L31)

```cpp
std::ifstream from(input_filename, std::ifstream::in | std::ifstream::binary);
return process::parse_process_specification(from);
```

There is no `from.is_open()`/`fail()` check, so a missing input file is handed to
the parser as a failed stream (an obscure parse error rather than "cannot open
file"). The trailing `return result;` after the if/else is also dead code.

## 2. Performance

- **SCC O(V·E) + recursion** — see 1.1, the dominant performance defect.
- [remove_equations.h](libraries/process/include/mcrl2/process/remove_equations.h#L170):
  `std::min_element(i.begin(), i.end())` over a `std::set<iterator>` whose elements
  are already ordered by the same `operator<` — `*i.begin()` is the minimum, so the
  linear scan is redundant.
- [alphabet_operations.h](libraries/process/include/mcrl2/process/alphabet_operations.h#L159):
  `concat_includes` carries the comment `// TODO: Make the implementation more
  inefficient!` — a typo for *efficient*; it is an `O(|A1|·|A2|)` membership test
  that could short-circuit or index `A1`/`A2`.
- [parse.h](libraries/process/include/mcrl2/process/parse.h#L114): the
  parse-and-type-check process-expression entry point is self-flagged
  `N.B. Very inefficient!` (it builds a whole specification to check one
  expression).

## 3. Maintainability / style

- **Diagnostics bypass the logging framework.**
  [is_linear.h](libraries/process/include/mcrl2/process/is_linear.h) writes to
  `std::clog` *and* `std::cerr` and throws a raw
  `std::runtime_error("unexpected error in visit_seq")` instead of using
  `mCRL2log` and `mcrl2::runtime_error`. It is the only production header in the
  library that prints to `std::clog`/`std::cerr` (the rest are test files). This
  echoes the LPS `parelm` finding.
- **Stochastic / non-stochastic duplication.**
  [add_binding.h](libraries/process/include/mcrl2/process/add_binding.h#L78) carries
  `// TODO: get rid of this code duplication`, the same theme called out across the
  LPS transformation headers.
- **Argument shadowing in `apply_comm`.** In
  [alphabet_operations.h](libraries/process/include/mcrl2/process/alphabet_operations.h)
  the inner `for (const core::identifier_string& a : alpha)` shadows the outer
  `const core::identifier_string& a = c.name();`. The code is correct (the
  `beta.insert(a)` runs outside the inner loop) but the reuse of `a` is a reading
  hazard.
- **Dead member.** `replace_subterm_builder::ready` in
  [replace_subterm.h](libraries/process/include/mcrl2/process/replace_subterm.h) is
  declared but never set or read.
- **`inline` precedes the doc comment.** In
  [remove_equations.h](libraries/process/include/mcrl2/process/remove_equations.h)
  the free function is written `inline\n/// \brief …\nvoid remove_duplicate_equations(...)`,
  so the `\brief` sits between `inline` and the signature.
- **Placeholder file briefs.** Most headers still carry
  `/// \brief add your file description here.` (e.g.
  [process_specification.h](libraries/process/include/mcrl2/process/process_specification.h),
  [multi_action_name.h](libraries/process/include/mcrl2/process/multi_action_name.h),
  [alphabet_operations.h](libraries/process/include/mcrl2/process/alphabet_operations.h),
  the `eliminate_*` headers) — the same doc-stub pattern seen across the toolset.
- **C++20 modernisation.** `tarjan_scc_algorithm::vertex` hand-writes
  `operator==`/`operator!=` (defaultable); many raw iterator loops in the alphabet
  and elimination code are `std::ranges` / structured-binding candidates; and the
  `// TODO: check if this function is needed` on `enter(process_instance)` in
  [is_linear.h](libraries/process/include/mcrl2/process/is_linear.h#L185) should be
  resolved.

## 4. Things that look good

- **House style is followed here.** Unlike `utilities`, every header uses the
  canonical `MCRL2_PROCESS_*_H` guard (no trailing underscore) with the matching
  `#endif` comment, and the full copyright/Boost preamble is present.
- **aterm discipline is clean.** No `reinterpret_cast` and no `static_cast`/
  `dynamic_cast` on aterm types anywhere in the library; downcasts use
  `atermpp::down_cast` and the `static_cast<Derived&>(*this)` occurrences are the
  standard CRTP idiom, exactly as the aterm-patterns guidance prescribes.
- **Thoughtful memoisation.** `comm_inverse` is driven by a precomputed
  `comm_inverse_cache` keyed on action name, and pcrl-equation alphabets are cached
  ([pcrl_equation_cache.h](libraries/process/include/mcrl2/process/detail/pcrl_equation_cache.h),
  `push_block_cache`/`push_allow` caches) to avoid recomputation during alphabet
  reduction.
- **`process_specification` equality** is defined via canonical
  `process_specification_to_aterm` conversion, giving a single well-defined notion
  of equality consistent with the term representation.
- **Test breadth.** Nine test files cover alphabet reduction, type-checking,
  parsing, printing, find, replace, actions, and sort traversal.

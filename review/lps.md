# lps library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The LPS library implements Linear Process Specifications: the data structures for
linear processes, the state-space *explorer*, linearisation, and a large family of
LPS-to-LPS transformations (`constelm`, `parelm`, `sumelm`, `suminst`, `untime`,
`binary`, `lpsparunfold`, `confluence`, …). It sits on top of `data` and `process`
and is depended on by `pbes`, `lts`, and the symbolic tooling.

## 1. Correctness

### 1.6 `probabilistic_data_expression` hash vs. value-equality — *latent bug*

`include/mcrl2/lps/probabilistic_data_expression.h`

`std::hash<probabilistic_data_expression>` hashes the underlying aterm
syntactically (`:560-575`). In the **non**-`MCRL2_ENABLE_MACHINENUMBERS` build,
`operator==` is *value*-based — it rewrites `equal_to(*this, other)` (`:373-393`).
So `1/2` and `2/4` compare **equal** yet hash **differently**. The type is consumed
by the LTS layer (`stochastic_lts_builder.h`, `lts_lts.h`, `liblts_lts.cpp`). In
the `MCRL2_ENABLE_MACHINENUMBERS` build the hash *is* consistent — but only under
the invariant that every fraction is stored in gcd-reduced canonical form (see
1.7).

### 1.7 `probabilistic_data_expression`'s `data_expression` constructor bypasses normalisation

`:342-347`:
```cpp
explicit probabilistic_data_expression(const data::data_expression& d)
 : data::data_expression(d)
{ assert(d.sort()==data::sort_real::real_()); }
```
Unlike the `(enumerator, denominator)` constructors (which call
`remove_common_divisor`), this constructor stores `d` verbatim. A non-reduced
`creal(2,4)` can therefore enter the type; in `MCRL2_ENABLE_MACHINENUMBERS` mode
`operator==` is pure aterm comparison, so such a value compares unequal to a
canonical `1/2` — a wrong answer the reduced-form arithmetic path assumes can't
happen. Either normalise in the constructor or document the hard precondition.

### 1.8 `probabilistic_data_expression` relational operators may throw

In the non-machinenumbers build both `operator==` (`:393`) and `operator<`
(`:411`) `throw mcrl2::runtime_error` when the comparison fails to rewrite to
`true`/`false`. These operators are exactly what ordered/hashed containers invoke
internally; a comparator that throws mid-operation can leave a `std::set`/`std::map`
inconsistent. Guaranteeing canonical form so comparisons are total and
non-throwing removes the hazard.

### Linearisation (author-flagged uncertainties)

#### 1.11 `linearise.cpp` self-flagged "unlikely correct" spots

`source/linearise.cpp`
- `:5273` — `stochastic_distribution(variable_list(), sort_real::real_one())` is
  annotated `// TODO: UNLIKELY THAT THIS IS CORRECT.`
- `:7524` — a timed-constraint merge is left unsimplified with a TODO.

Author-acknowledged uncertainties in a 10k-line file; they deserve dedicated
regression tests around stochastic/timed linearisation.

#### 1.12 Unresolved algorithmic TODOs in `linearise_communication.h`

`:750` (`// TODO: Why don't we insert r here, as van Weerdenburg writes?`) and
`:770` (`// TODO: van Weerdenburg in his note only calculates S := S ∪ T.
Understand why that is not correct.`) question whether the implementation matches
van Weerdenburg's note. Worth resolving or documenting why any deviation is
intentional.

## 2. Thread-safety (state-space explorer)

### 2.1 DFS callbacks are explicitly *not* thread-safe

`include/mcrl2/lps/explorer_dfs.h`
- `:70`  `examine_transition(...); // TODO MAKE THREAD SAFE`
- `:100` `finish_state(0, s0); // TODO MAKE THREAD SAFE`
- `:190`, `:226` — same.

The DFS path passes hard-coded `thread_index = 0` / `number_of_threads = 1` and its
callbacks are not safe under the multi-threaded harness. This is documented only in
scattered TODO comments, not enforced in the API (e.g. an assert that
`number_of_threads == 1` for DFS).

### 2.2 Lock-discipline complexity in BFS

`include/mcrl2/lps/explorer_bfs.h` hand-manages `m_exclusive_state_access` with
conditional `lock()`/`unlock()` guarded by
`GlobalThreadSafe && number_of_threads > 1`. The lock is taken before the loop and
released inside it depending on `todo` emptiness — correct but fragile. A scoped
guard / clearer ownership model would reduce the risk of a future edit leaving the
mutex in the wrong state. Several `// TODO: join duplicate targets` indicate
incomplete optimisation.

### 2.3 Stray debugging artifact

`include/mcrl2/lps/explorer_bfs.h:50`:
```cpp
data::data_specification thread_data_specification = m_global_lpsspec.data(); /// XXXX Nodig??
```
A copy of the entire data specification is made per thread with a Dutch "is this
needed?" comment. Either it is required (document why) or it is an avoidable
per-thread copy.

## 3. Performance

### 3.1 `next_state` is a linear scan, called in quadratic context

`detail/lps_algorithm.h:85-94` implements `next_state(summand, v)` as a linear scan
over `summand.assignments()`. In `constelm.h` it is called inside the fixpoint loop
*per summand × per still-constant parameter*, making a single sweep roughly
`O(#summands · #params · #assignments)`. Building a per-summand `variable → rhs`
map once would cut a factor.

### 3.2 `constelm` re-rewrites the same expressions repeatedly

`constelm.h:190-191` evaluates `R(g_ij, sigma)` on line 190, then again on line 191
to form the log message / `z`. Line 190's result could be hoisted into a local and
reused. Rewriting is the dominant cost, so caching these is worthwhile.

### 3.3 `Invariant_Checker` copies the entire summand vector

`invariant_checker.h:245`: `f_summands = a_lps.process().action_summands();` stores
a **copy** of the summand vector in a checker that only reads it. Since `f_spec` is
already held by const reference, `f_summands` could be a reference.

### 3.4 `parelm1` is the quadratic default

`parelm.h::run(bool variant1 = true)` defaults to `parelm1()`, the worst-case
`O(#vars · #summands · #assignments)` fixpoint, while the graph-based `parelm2()`
(`boost::adjacency_list` + `reachable_nodes`) is asymptotically better. Revisit the
default, or at least document why parelm1 is preferred (it may avoid the Boost
graph construction overhead on small specs).

### 3.5 `find_equations` / `parameter_selection` rebuild work per call

- `lpsparunfoldlib.h:289-295` (`find_equations`) iterates **all** data-spec
  equations on every call; caching the result by head symbol would help.
- `detail/parameter_selection.h:66` constructs a fresh `std::regex` on every call.
  Making it a `static`/`thread_local static` avoids recompiling the pattern.

### 3.6 `index` counters typed `unsigned`/`int`

`constelm.h:145` uses `unsigned index = 0` and `parelm.h:146` uses `int index = 0`
while both feed a `std::map<data::variable, std::size_t>`. Use `std::size_t`.

### 3.7 `parelm` bypasses the logging framework

`parelm.h::report_results` writes verbose output to `std::clog` (`:94`, `:97`)
instead of `mCRL2log(log::verbose)`. Every other algorithm routes verbose output
through `mCRL2log`; this bypasses formatting, redirection, and log-level gating.

## 4. Maintainability / architecture

### 4.1 `linearise.cpp` is a 10k-line monolith

`source/linearise.cpp` is by far the largest file in the library. It mixes parsing
helpers, the communication/allow operators, timed and stochastic handling, and
post-processing (it even instantiates `constelm_algorithm` internally). The header
side has already been partially extracted (`linearise_communication.h`,
`linearise_allow_block.h`, `linearise_rename.h`, `linearise_utility.h`); the `.cpp`
would benefit from the same decomposition.

### 4.2 Heavy library-level dependencies

`libraries/lps/CMakeLists.txt` makes `mcrl2_lps` depend on `mcrl2_smt` **and**
`mcrl2_symbolic` for the whole library. Most LPS consumers need neither. Consider
moving the symbolic pieces (`lpsreach.h`, `symbolic_lts*.h`) and SMT usage behind
an optional sub-target so the core LPS data structures stay lightweight.

### 4.3 `check_well_typedness` invoked inside `assert(...)`

`source/lpsparunfoldlib.cpp:673`, `:836`; `include/mcrl2/lps/parelm.h:135`, `:191`
put an expensive whole-spec validation inside `assert`, so it only runs in debug
builds and gives no message on failure. This is the intended idiom, but be aware
the guarantee is *not* checked in production builds.

### 4.4 Recurring stochastic/non-stochastic code duplication

The most common maintainability theme: stochastic and non-stochastic variants are
copy-pasted rather than templated.
- `add_binding.h:114` — `// TODO: get rid of this code duplication`.
- `confluence.h:69`, `:77` — `// TODO: reuse this code`.
- `detail/linear_process_conversion_traverser.h:21` — `// TODO: join the stochastic
  and non-stochastic versions of the traversers`; `:272`/`:679` carry identical
  `// TODO: check if this is correct` bodies.
- `ltsmin.h:308` — `// TODO: get rid of this useless function`.

The library already uses `linear_process_base<ActionSummand>` and
`specification_base<...>` templating successfully — the same approach could fold
most of these.

### 4.5 Macro-divergent duplicated fraction arithmetic

`probabilistic_data_expression.h` carries two complete implementations of `+`, `-`,
`<`, `==` selected by `#ifdef MCRL2_ENABLE_MACHINENUMBERS`: ~250 lines of
hand-rolled multi-word fraction arithmetic vs. a few lines that defer to the data
rewriter. Only one compiles at a time, so the inactive branch can rot unnoticed.
Worth a dedicated unit test running the same fraction-arithmetic vectors regardless
of the macro, ideally extracting the arithmetic into one tested utility.

### 4.6 Unresolved stochastic-handling question in rewriters

`rewriters/one_point_condition_rewrite.h:132` carries
`// TODO: should the stochastic distribution be rewritten or not?` above
`update(lps::stochastic_action_summand&)`. The stochastic distribution handling
across the transformation algorithms deserves a single deliberate policy.

## 5. Minor / stylistic

- **Serialisation field order is non-obvious but consistent.** `source/lps_io.cpp`
  writes and reads `process_parameters`, `deadlock_summands`, `action_summands` in
  the same order; the on-disk order differs from the `linear_process_base` member
  declaration order. A comment would prevent a future "tidy-up" from
  desynchronising reader and writer.
- `action_rename.h:171` uses `static_cast<data::data_specification>(*i++)` where the
  codebase otherwise prefers explicit constructors / `down_cast`; `:131-135`
  `action_rename_rule::lhs()` returns `process::action` but the stored term may be a
  `ParamId(...)` the type checker rejects — a long-lived "temporary fix" worth
  tracking as real debt.
- `confluence_checker.h:907` — `f_intermediate` is sized `f_number_of_summands + 2`
  (over-allocation / in-band index sentinel); document why or use
  `std::optional<std::size_t>`.
- `lpsparunfoldlib.h:478` — `max_unfold_depth = 3` is a `constexpr` but its choice
  is unexplained; a one-line rationale would help.
- **C++20 modernisation:** many raw iterator loops remain (e.g. `multi_action.h`,
  `untime.h`, `linear_process.h`, `remove.h`, `is_well_typed.h`,
  `detail/specification_property_map.h`, `confluence_checker.h`) — mechanical
  `std::ranges`/structured-binding candidates. Comparison operators are still
  hand-written (`deadlock.h`, `multi_action.h`, `action_summand.h`,
  `stochastic_action_summand.h`); note `multi_action` cannot naively adopt the
  default spaceship because of its permutation-aware semantics (1.1).

## 6. Things that look good

- **C++20 concepts already partially adopted.** The rewriter/replace layer uses the
  `data::IsSubstitution` concept, and newer transformation headers use `requires`
  clauses instead of SFINAE.
- **Templating of the spec/process types** (`linear_process_base<ActionSummand>`,
  `specification_base<...>`) cleanly shares code between timed/stochastic and plain
  variants.
- **The explorer's callback design** (`discover_state`, `examine_transition`,
  `start_state`, `finish_state` as template parameters defaulting to
  `utilities::skip`, with `if constexpr (is_applicable<...>)`) is a nice
  zero-overhead extension point.
- **Caching in `lpsparunfold`** via `unfold_cache_element` keyed by sort, with an
  explicit comment warning that the case-function argument must be passed by value
  to avoid invalidation.
- **`suminst` error recovery** catches enumeration failures, rolls back partially
  added summands, and accounts correctly for "0 solutions ⇒ summand deleted".
- **The machine-number fraction path** divides out the gcd of the two denominators
  *before* multiplying to keep intermediate products inside one 64-bit word,
  falling back to arbitrary precision only on carry — a thoughtful fast/slow split.
- **`state`** is `atermpp::term_balanced_tree<data::data_expression>`, giving
  O(log n) structural sharing for the explorer's millions of states.
- **Test coverage breadth** is good: 37 test files, including dedicated
  linearisation tests (timed, stochastic, communication, allow/block, rename).

## Resolved since previous review (removed)

Verified fixed in the current source and dropped from this document: `1.3`
(`highway_todo_set` off-by-one, now `k < N`), `1.4` (`binary` now rewrites the
stochastic distribution via an `auto&` loop), `1.10` (`explorer_options` output
label), `4.7` (the dead `to_number_t` helper was removed), and the section-5 doc
stubs / file-tag / guard / preamble fixes. `1.9` (deadlock-summand substitution)
was investigated previously and confirmed **not** a bug.

# pres library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `pres` library is the real-valued sibling of `pbes`: it holds the PRES
(parameterised real equation system) data structures, the `lps2pres`/`lts2pres`
translations, the usual PBES-style transformations (`constelm`, `parelm`,
instantiation), and — uniquely — three RES solvers: an exact Gauss-elimination
solver (`ressolve_gauss_elimination.h`) and two numerical fixpoint-iteration
solvers (`ressolve_numerical.h`, `ressolve_numerical_directed.h`). Almost every
file is annotated "Based on pbes/…", so the library inherits both the sound
patterns *and* several defects of `pbes`. The solvers, rewriters, translations and
I/O were read; the generated builder/traverser/expression files were skimmed.

## 1. Correctness

### 1.1 `simplify_quantifiers_rewriter` drops the empty-domain / negation branch — same dead-code bug as pbes

[simplify_quantifiers_rewriter.h](libraries/pres/include/mcrl2/pres/rewriters/simplify_quantifiers_rewriter.h#L34)

The `infimum` and `supremum` handlers repeat the exact `if`/missing-`else`
structure flagged in the pbes review:

```cpp
if (variables.empty())   { result = true_(); }                 // :35
else if (is_minus(body)) { optimized_supremum(...); optimized_minus(result, result); }  // :39

if (is_and(body))        { ... result = ... }                  // :44  — bare `if`, not `else if`
else if (is_or(body))    { ... }
else                     { optimized_infimum(result, variables, body); }
```

The `supremum` handler is identical with `is_or` first
([:96](libraries/pres/include/mcrl2/pres/rewriters/simplify_quantifiers_rewriter.h#L96)).
Because the second chain always reassigns `result`, the empty-domain shortcut
(`true_()`/`false_()`) and the `minus`-distribution rewrite in the first chain are
**dead code**: when `body` is a `minus` it is neither `and` nor `or`, so control
falls to the final `else` and the just-computed minus-pushed value is discarded.
As in pbes this is a lost optimisation rather than a wrong answer, but the two
libraries should be fixed together.

### 1.2 `constelm`'s `sum` case is unimplemented and writes to `std::cerr`

[constelm.h](libraries/pres/include/mcrl2/pres/constelm.h#L362)

```cpp
void leave(const sum& )
{
std::cerr << "MUST STILL BE DONE\n";
}
```

The constant-elimination traverser's handler for the PRES-specific `sum` operator
is a stub: it prints an (unindented) message to `std::cerr` and does nothing. Any
PRES whose equations contain `sum` is silently mis-analysed by `prESconstelm`, and
the diagnostic bypasses `mCRL2log`. This is an incompleteness bug, not merely a
style issue.

## 2. Robustness

### 2.1 The numerical RES solvers have no iteration bound — divergence hangs

[ressolve_numerical.h](libraries/pres/include/mcrl2/pres/ressolve_numerical.h#L172)

```cpp
do {
  do { calculate_new_solution(base_equation_index, i); }
  while (!stable_solution_found(base_equation_index, i));
  apply_numerical_recursive_algorithm(i);
  calculate_new_solution(base_equation_index, i);
} while (!stable_solution_found(base_equation_index, i));
```

`stable_solution_found` returns `error <= pow(0.1, precision)`; there is **no
maximum-iteration cap**. A system that diverges, oscillates, or merely plateaus
just above the tolerance loops forever with no progress reporting beyond a debug
line. `ressolve_numerical_directed.h`
([~:255](libraries/pres/include/mcrl2/pres/ressolve_numerical_directed.h#L255))
has the identical unbounded structure. A configurable iteration limit (or a stall
detector) that aborts with a diagnostic would make the numerical solvers safe to
run on untrusted input.

### 2.2 `io.cpp` uses `reinterpret_cast` on an aterm (inherited from pbes)

[io.cpp](libraries/pres/source/io.cpp#L178)

```cpp
const data::function_symbol& y = reinterpret_cast<const data::function_symbol&>(x);
```

`add_index_impl` was copied verbatim from `pbes/source/io.cpp` and carries the same
`reinterpret_cast` on an `atermpp::aterm` that the pbes review flagged; the
guidelines require `atermpp::down_cast`.

### 2.3 Release-disabled float invariant in the directed solver

[ressolve_numerical_directed.h](libraries/pres/include/mcrl2/pres/ressolve_numerical_directed.h#L412)

```cpp
const double left_argument = data::sort_real::value<double>(pp.left());
assert(left_argument > 0);
```

The positivity of a `const_multiply` coefficient is checked only by `assert`, which
vanishes in release builds; if the invariant is ever violated the solver proceeds
with a non-positive multiplier instead of raising a `runtime_error`.

## 3. Performance

### 3.1 Numerical solver keys solutions by `identifier_string` in a hash map

[ressolve_numerical.h](libraries/pres/include/mcrl2/pres/ressolve_numerical.h#L139)

`calculate_new_solution` and `stable_solution_found` read and write
`m_new_solution`/`m_previous_solution` (both
`std::unordered_map<core::identifier_string, double>`) with `operator[]` inside the
innermost iteration loops, so every access hashes an identifier and may
default-insert. The Gauss solver already indexes equations positionally; a
`std::vector<double>` indexed by equation position would remove the hashing from
the hot loop.

## 4. Maintainability / style

- **Header-guard defects.**
  [is_res.h](libraries/pres/include/mcrl2/pres/is_res.h#L12) guards with
  `MCRL2_PRES_IS_BES_H` (copied from `is_bes.h`).
  [ressolve_numerical.h](libraries/pres/include/mcrl2/pres/ressolve_numerical.h#L12)
  opens `MCRL2_PRES_RESSOLVE_H` but closes `#endif // MCRL2_PRES_RESSOLVE_NUMERICAL_H`,
  and
  [ressolve_numerical_directed.h](libraries/pres/include/mcrl2/pres/ressolve_numerical_directed.h#L12)
  opens `MCRL2_PRES_RESSOLVE_DIRECTED_H` but closes `…_NUMERICAL_DIRECTED_H` — the
  open guard and its `#endif` comment disagree in both files (no actual collision,
  since the three macros differ, but neither matches its filename). The
  Gauss-elimination header gets this right.
- **`\file` / `\brief` copy-paste errors.**
  [ressolve_numerical_directed.h](libraries/pres/include/mcrl2/pres/ressolve_numerical_directed.h#L9)
  tags itself `\file mcrl2/pres/ressolve_numerical.h`;
  [enumerate_quantifiers_rewriter.h](libraries/pres/include/mcrl2/pres/rewriters/enumerate_quantifiers_rewriter.h#L9)
  tags `\file mcrl2/pbes/rewriters/…`; both numerical solvers describe themselves as
  "a gauss-elimination like algorithm" although they implement numerical iteration.
- **Leftover debug print.**
  [detail/parse.h](libraries/pres/include/mcrl2/pres/detail/parse.h#L124) contains
  `std::cerr << "EXPRESSION " << expression << "\n";` in the production parser.
- **Pervasive placeholder briefs.** Over twenty headers — `algorithms.h`,
  `print.h`, `replace.h`, `rewrite.h`, `typecheck.h`, `normalize_sorts.h`,
  `translate_user_notation.h`, `remove_equations.h`, `significant_variables.h`,
  `presinst_*`, `lts2pres.h`, `is_res.h`, … — still carry
  `/// \brief add your file description here.`
- **Unusual cast in `lts2pres`.**
  [lts2pres.h](libraries/pres/include/mcrl2/pres/lts2pres.h#L78) does
  `static_cast<const std::string&>(detail::mu_name(f))` on a `core::identifier_string`
  (an aterm wrapper) to build a string; the sibling lines pass the
  `identifier_string` straight into the `propositional_variable_instantiation`
  constructor, which is the cleaner idiom.

## 5. Things that look good

- **Cleanly single-threaded solvers.** Unlike `pbes`/`lts`, the RES solvers and
  instantiation use no manual `lock()`/`unlock()` and no shared mutable static
  state, so they carry none of the non-RAII-locking hazards flagged elsewhere.
- **Value caching in evaluation.** `detail::evaluate` memoises
  `data::sort_real::value<double>` conversions in a `value_cache`, and `run()`
  pre-rewrites every equation body once (`m_R(eq.formula())`) before iterating.
- **aterm discipline is otherwise good.** Casting uses `atermpp::down_cast`
  throughout (the `io.cpp` `reinterpret_cast` in 2.2 is the lone exception; the
  `static_cast<pres_expression>(x)` calls in the Gauss solver are safe *upcasts*
  from `infimum`/`supremum`/`sum`, not narrowing casts).
- **One header gets the hygiene right.** `ressolve_gauss_elimination.h` has a
  matching, correctly-named guard and a real `\brief`, showing the intended
  standard the numerical headers drifted from.

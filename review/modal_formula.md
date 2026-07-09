# modal_formula library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `modal_formula` library is the property front end: it defines action formulas,
regular formulas, and state formulas (the modal mu-calculus, recently extended with
*quantitative* operators — `minus`, `plus`, `const_multiply`, `infimum`,
`supremum`, `sum` — for real-valued model checking), and the preprocessing that
turns a user formula into the normalized, fixpoint-wrapped shape consumed by
`lps2pbes`/`lts2pbes` (and their `pres` counterparts). The bulk of the type
definitions, the builder, and the traverser are generated; the review concentrated
file-by-file on the hand-written algorithms — normalization, monotonicity, the
regular-formula translation, name-clash checking, and the nested-modal-operator
preprocessing.

## 2. Robustness

### 2.1 `regfrmtrans` returns an untranslated formula on an unexpected operator in release

[regfrmtrans.cpp](libraries/modal_formula/source/regfrmtrans.cpp#L230)

The recursive regular-formula translator ends its big type dispatch with

```cpp
else
{
  assert(0); //Unexpected state_formula.
}
return part;
```

In a release build the `assert` vanishes, so an unhandled state-formula operator is
returned **untranslated** instead of raising an error. It cannot be hit by today's
grammar, but the library has just grown a family of quantitative operators, and a
future addition that forgets to extend this dispatch would silently pass a formula
containing regular operators through to `lps2pbes`. A `throw mcrl2::runtime_error`
(as the sibling `is_monotonous`/`normalize` dispatchers already use for their
"unknown argument" case) would fail loudly instead.

## 3. Maintainability / style

- **Wrong `\file` path in `has_name_clashes.h`.**
  [has_name_clashes.h](libraries/modal_formula/include/mcrl2/modal_formula/has_name_clashes.h#L9)
  tags itself `\file mcrl2/modal_formula/detail/has_name_clashes.h`, but the file
  lives at `modal_formula/has_name_clashes.h` (and its guard correctly omits
  `DETAIL`) — a copy-paste leftover from when it presumably lived under `detail/`.
- **Pervasive placeholder briefs.** About twenty headers still carry
  `/// \brief add your file description here.` (with inconsistent `Add`/`add`
  capitalisation in `state_formula.h`, `normalize.h`, `state_formula_rename.h`),
  and both source files have an empty `/// \brief`
  ([modal_formula.cpp](libraries/modal_formula/source/modal_formula.cpp#L10),
  [regfrmtrans.cpp](libraries/modal_formula/source/regfrmtrans.cpp#L9)).
- **Inconsistent typed-view construction in `regfrmtrans`.** The translator builds
  typed views twice per node (`must(part).formula()` *and* `must(part).operand()`,
  `and_(part).left()` *and* `and_(part).right()`) where the `not_`/`minus` cases use
  a single `atermpp::down_cast`; harmless (aterm construction is cheap) but
  needlessly verbose.

## 4. Things that look good

- **Normalisation handles the quantitative dual correctly.**
  [normalize.h](libraries/modal_formula/include/mcrl2/modal_formula/normalize.h#L54)
  pushes negation through *both* the Boolean and the real-valued connectives
  (`not`↔`minus`, `and`↔`or`, `forall`↔`exists`, `supremum`↔`infimum`, `must`↔`may`,
  `mu`↔`nu` via `negate_variables`), and correctly refuses the cases that cannot be
  inverted — a negated free variable and the `plus` operator both raise a
  `runtime_error` rather than producing a wrong dual.
- **Complete monotonicity check.**
  [is_monotonous.h](libraries/modal_formula/include/mcrl2/modal_formula/is_monotonous.h#L25)
  tracks positive/negative variable polarity through every operator including the
  new quantitative ones and throws on an unknown argument.
- **Standard, capture-safe regular translation.** `regfrmtrans` is the textbook
  Kozen expansion (`[R*]φ → νX. φ ∧ [R]X`, `<R*>φ → μX. φ ∨ <R>X`) using an
  `xyz_identifier_generator` to mint fresh `X`, and terminates because each
  expansion strictly shrinks the regular sub-expression.
- **Clean name-clash checkers.** The nested-`mu`/`nu` and data-parameter clash
  checkers in `has_name_clashes.h` use tidy push/pop stack and set discipline.
- **Logging-clean, cast-clean.** No production header or source writes to
  `std::cerr`/`std::cout`/`std::clog` (only the tests do), `assert(0)` appears once,
  and casting uses `atermpp::down_cast` throughout — no `reinterpret_cast`/
  `static_cast` on aterms.

# data library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `data` library is the foundation of the whole toolset: abstract data types,
sorts, the data specification, the type checker and parser, the two term-rewriting
engines (the interpreting `jitty` rewriter and the C++-compiling `jittyc`
rewriter), the enumerator, and the BDD/SMT prover used for linear-arithmetic
reasoning. It is also the largest and most performance-critical library, with a lot
of hand-tuned low-level C++. The generated sort headers (`int.h`, `nat.h`, `set.h`,
…) and the generated class/builder/traverser code were skimmed; the review went
file-by-file through the rewrite engines, the prover, the enumerator, and the
linear-inequality machinery.

## 1. Thread-safety

### 1.1 File-scope mutable statics in the compiling rewriter's code generator

[jittyc.cpp](libraries/data/source/detail/rewrite/jittyc.cpp#L493)

The `jittyc` (compiling) rewriter keeps mutable state at namespace/function scope:

```cpp
static std::vector<std::size_t> treevars_usedcnt;          // :493
...
static std::size_t auxiliary_method_name_index = 0;        // :1710
```

`treevars_usedcnt` is rebuilt and indexed during match-tree construction, and
`auxiliary_method_name_index` numbers the generated C++ helper methods. Both are
shared across all `jittyc` instances in the process, so constructing two compiling
rewriters on different threads concurrently races (and a shared name counter could
even produce colliding symbol names in the generated code). In practice the
compiling rewriter is built once during setup rather than in the parallel hot path,
so this is latent rather than active — but it is undocumented shared state in a
library that is otherwise carefully cloned per thread.

### 1.2 Global enumerator iteration limit is unsynchronised

[enumerator_iteration_limit.h](libraries/data/include/mcrl2/data/detail/enumerator_iteration_limit.h#L30)

`set_enumerator_iteration_limit` writes a process-global static
(`enumerator_iteration_limit<std::size_t>::max_enumerator_iterations`) with no
synchronisation. It is a startup configuration knob (default 1000) so this is benign
by design, but it is genuinely global mutable state shared by every enumerator on
every thread.

## 2. Memory-safety

### 2.1 `BDD_Path_Eliminator` leaks its SMT solver (no destructor, no rule-of-three)

[bdd_path_eliminator.h](libraries/data/include/mcrl2/data/detail/prover/bdd_path_eliminator.h#L52)

```cpp
SMT_Solver* f_smt_solver;                               // :52  raw owning pointer
...
f_smt_solver = new mcrl2::data::detail::prover::cvc_smt_solver();   // :186
f_smt_solver = new mcrl2::data::detail::prover::z3_smt_solver();    // :195
```

The class allocates `f_smt_solver` with `new` in its constructor but defines **no
destructor** (and no copy control), so every `BDD_Path_Eliminator` leaks its solver,
and a copy would share/again-delete the same pointer. A `std::unique_ptr<SMT_Solver>`
would fix both the leak and the rule-of-three violation. This is in the experimental
BDD prover, so the impact is limited.

## 3. Robustness — the external SMT solver bridge

The SMT prover shells out to an external solver (`cvc3`/`z3`/`ario`).
[smt_lib_solver.cpp](libraries/data/source/detail/prover/smt_lib_solver.cpp#L40)
forks, wires up pipes, and `execlp`s the solver. To be clear about what is *not*
wrong: the SMT-LIB benchmark is fed to the child over **stdin**, and the arguments
to `execlp` are fixed literals, so there is **no command injection**. The real
issues are reliability ones:

- **Uninitialised buffer read.**
  [smt_lib_solver.cpp](libraries/data/source/detail/prover/smt_lib_solver.cpp#L80)
  declares `std::array<char, 64> output;` (uninitialised), checks only that
  `::read(pipe_stdout[0], output.data(), 8) > 0`, then does
  `strncmp(output.data(), "unsat", 5)` / `"unknown", 7` — comparing up to seven
  bytes when as few as one may have been read, i.e. reading uninitialised stack.
  (The child is the trusted solver, so this is a correctness bug, not a trust-
  boundary issue.)
- **Fd / zombie leak on the error path.** If `::write` of the benchmark fails the
  function throws while `pipe_stdout[0]`, `pipe_stderr[0]` and `pipe_stdin[0]` are
  still open and the child has not been `wait`ed, leaking descriptors and orphaning
  the process.
- **`usable()` throws then "returns".**
  [smt_lib_solver.cpp](libraries/data/source/detail/prover/smt_lib_solver.cpp#L131)
  has a `return false;` immediately after a `throw`, which is unreachable, and the
  function throws rather than returning a usability bool as its name implies.
- **PATH-based `execlp`.**
  [smt_lib_solver.h](libraries/data/include/mcrl2/data/detail/prover/smt_lib_solver.h#L874)
  runs whichever `z3`/`cvc3`/`ario` appears first on `PATH`; acceptable for a
  developer tool but worth noting.

## 4. Robustness — `assert` as control flow

[fourier_motzkin.h](libraries/data/include/mcrl2/data/fourier_motzkin.h#L127) ends
its positive/negative factor classification with a bare `assert(0)`:

```cpp
else if (is_negative(f,r)) { inequalities_with_negative_variable.push_back(e); }
else                       { assert(0); }
```

`is_positive`/`is_negative` rely on heuristic rewriting and are not guaranteed
exhaustive, so a factor the rewriter cannot classify falls through silently in a
release build. The compiling rewriter has the same pattern at
[jittyc.cpp](libraries/data/source/detail/rewrite/jittyc.cpp#L2299). A
`throw mcrl2::runtime_error` would fail loudly instead. (`negate` in
[linear_inequalities.h](libraries/data/include/mcrl2/data/linear_inequalities.h#L1413)
similarly returns a documented-unreachable default.)

## 5. Things that look good

- **Graceful rewrite-stack growth.**
  [rewrite_stack.h](libraries/data/include/mcrl2/data/detail/rewrite/rewrite_stack.h#L42)
  uses a thrown `recalculate_term_as_stack_is_too_small` to resize the rewrite
  vector and restart the current rewrite, so deeply-nested terms grow the heap-backed
  stack instead of overflowing the processor stack — and because the stack is an
  `atermpp::vector`, all intermediate terms are aterm-protected without per-entry
  bookkeeping.
- **Carefully-bounded stack allocation.** The hot matching path sizes its
  per-rule assignment buffer to `strat.number_of_variables()` and the matcher
  de-duplicates variables (each distinct variable is bound at most once), so the
  `MCRL2_SPECIFIC_STACK_ALLOCATOR` buffer cannot overflow — paired placement-new and
  explicit destruction keep it exception-correct.
- **Clean engine layering.** The interpreting (`jitty`) and compiling (`jittyc`)
  rewriters share a match-tree representation but are otherwise cleanly separated,
  and the per-thread `clone()` design keeps the rewriter state thread-local apart
  from the statics noted in 1.1.

> Note: the apparent unbounded write at
> [jitty.cpp](libraries/data/source/detail/rewrite/jitty.cpp#L376) (`assignments.size++`
> without a bound) is **not** a defect — the buffer is sized to the rule's variable
> count and the matcher binds each distinct variable only once; it is included here
> only to record that it was checked.

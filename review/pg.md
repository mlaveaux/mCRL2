# pg library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `pg` library solves parity games — the back end of the verification stack
(`pbespgsolve` turns a PBES into a parity game and applies one of several solvers:
the recursive/Zielonka solver, small-progress-measures with a family of lifting
strategies, priority promotion, and the decycle/deloop/component pre-solvers).
**It is imported external code** (copyright 2009–2013 University of Twente /
Michael Weber / Maks Verver) and therefore looks nothing like the rest of mCRL2:
`PascalCase` filenames and class names, `verti`/`edgei` integer typedefs, raw
`new[]`/`delete[]` ownership, `volatile`, `assert`-based validation and direct
`<iostream>` use instead of aterms, `snake_case`, RAII and `mCRL2log`. The review
therefore weighs two things differently from the other libraries: the wholesale
style deviation is *acceptable* for walled-off imported code, but the memory-safety
and untrusted-input handling still matter because the solver runs on user-supplied
games. The solvers, the graph/parity-game core, the I/O, the reference-counting and
abort machinery, and the lifting strategies were read.

## 1. Robustness

### 1.1 The recursive solver has unbounded recursion depth

[RecursiveSolver.cpp](libraries/pg/source/RecursiveSolver.cpp#L99)

`RecursiveSolver::solve` is the classic Zielonka algorithm with **two** recursive
calls per level (one on the minimum-priority subgame, one on the opponent's
attractor complement). Each level removes at least the lowest priority, so the
depth is `O(d)` in the number of distinct priorities — but for an adversarial input
`d` can be `O(V)`, and there is no depth guard (only an `aborted()` check at the top,
[:101](libraries/pg/source/RecursiveSolver.cpp#L101)). A game with a long priority
ladder overflows the stack. This is the same hazard flagged for the `process`
Tarjan SCC and the pbes solver; converting to an explicit work stack, or imposing a
depth limit that falls back to the iterative SPM solver, would make it safe on
untrusted games.

### 1.2 Untrusted PGSolver input is validated only by `assert`

[ParityGame_IO.cpp](libraries/pg/source/ParityGame_IO.cpp#L97)

```cpp
assert(prio_raw >= 0);
assert(prio_raw < 65536);
priority_t prio = prio_raw;
assert(player_raw == 0 || player_raw == 1);
...
assert(vertices[id] == invalid);   // :118
```

Stream *extraction* failures are handled (`while (is)`, `if (!(is >> …)) break;`),
but the **range** checks on a parsed `.pg` file are `assert`s, which vanish in
release builds. A malformed or hostile PGSolver file with a negative/huge priority,
a `player` field outside `{0,1}`, or a duplicate vertex id then silently produces a
corrupt `ParityGame` instead of a clean `runtime_error`. The
[FIXME at :115](libraries/pg/source/ParityGame_IO.cpp#L115) acknowledges that the
duplicate-vertex case (legal in the PGSolver format) is unsupported.

### 1.3 `read_raw` sizes allocations from an unvalidated binary header

[Graph.cpp](libraries/pg/source/Graph.cpp#L530)

```cpp
is.read(reinterpret_cast<char*>(&V), sizeof(V));
is.read(reinterpret_cast<char*>(&E), sizeof(E));
is.read(reinterpret_cast<char*>(&edge_dir), sizeof(edge_dir));
reset(V, E, edge_dir);            // allocates new verti[E], new edgei[V+1], …
```

`V`, `E` and `edge_dir` are read straight from the stream and handed to `reset`,
which performs the `new[]` allocations, with no sanity bound and no `is.good()`
check after the subsequent bulk reads. A corrupt/hostile raw graph file can request
an enormous allocation (DoS / `bad_alloc`) or, if the bulk reads hit EOF partway,
leave the successor/predecessor arrays partly uninitialised with no error raised.

## 2. Thread-safety

### 2.1 `Abortable` uses a `volatile` flag, not an atomic

[Abortable.h](libraries/pg/include/mcrl2/pg/Abortable.h#L24)

```cpp
static void abort_all() { global_abort_ = true; }
bool aborted() { return global_abort_; }
...
static volatile bool global_abort_;
```

`volatile` is not a synchronisation primitive in C++: a write in `abort_all()` and a
read in `aborted()` from different threads is a data race. It happens to work on the
usual platforms for a single bool, but the correct type is `std::atomic<bool>`.
Separately, the flag is a single process-global, so aborting one solver aborts
*every* `Abortable` in the process (including solvers on other threads) — a design
limitation worth noting now that the toolset is multithreaded by default.

### 2.2 `RefCounted` increments a non-atomic counter

[RefCounted.h](libraries/pg/include/mcrl2/pg/RefCounted.h#L34)

`ref()`/`deref()` mutate a `mutable std::size_t refs_` with plain `++`/`--`, so a
reference-counted object (lifting-strategy factory, etc.) shared across threads gets
a torn count and is freed early or leaked. `deref()` ends in `delete this`
([:43](libraries/pg/include/mcrl2/pg/RefCounted.h#L43)), which requires heap
allocation (the class documents this contract). In current use these objects are
per-thread, so this is latent rather than active, but it is the same non-atomic-
refcount pattern to avoid.

## 3. Correctness

### 3.1 Unsigned underflow in the linear lifting strategy

[LinearLiftingStrategy.cpp](libraries/pg/source/LinearLiftingStrategy.cpp#L55)

```cpp
vertex_ = vertex_ - failed_lifts_ - 1;
```

`vertex_` and `failed_lifts_` are unsigned `verti`. When `failed_lifts_ + 1 >
vertex_` this wraps to a huge value, so the alternating sweep restarts from a
nonsensical vertex index instead of the intended position. A `std::max`/saturating
guard would keep the index in range.

## 4. Maintainability / style

- **Whole-library convention deviation (by design).** The imported origin means the
  library breaks essentially every house rule — `PascalCase`, `verti`/`edgei`, raw
  `new[]`/`delete[]`, `volatile`, `assert`-validation, direct `<iostream>`. This is
  *acceptable* precisely because it is third-party code kept behind the
  `pbespgsolve` bridge; the recommendation is to keep it walled off (not to
  reformat it) while addressing the robustness items above that affect user input.
- **The one mCRL2-side bridge file lacks a brief.**
  [pbespgsolve.h](libraries/pg/include/mcrl2/pg/pbespgsolve.h#L10) carries the
  placeholder `/// \brief add your file description here.` — notable because it is
  the single file written in mCRL2 style; the pg-native headers use their own
  consistent `/*! \file X.h */` Doxygen style.
- **Acknowledged FIXME/TODO/HACK markers.**
  [DecycleSolver.cpp](libraries/pg/source/DecycleSolver.cpp#L69)
  (`has_succ` is not constant-time, making an inner loop `O(SCC²)`),
  [MaxMeasureLiftingStrategy.cpp](libraries/pg/source/MaxMeasureLiftingStrategy.cpp#L34)
  (heap operations swap more than necessary),
  [SmallProgressMeasures.cpp](libraries/pg/source/SmallProgressMeasures.cpp#L547)
  (`const_cast<StaticGraph&>` to mutate a `const` graph), and
  [PriorityPromotionSolver.cpp](libraries/pg/source/PriorityPromotionSolver.cpp#L242)
  (a `TODO` flagging an unresolved discrepancy between the paper and the
  implementation — a latent correctness question) are all worth tracking.

## 5. Things that look good

- **`StaticGraph` owns its arrays safely.** Despite holding raw `verti*`/`edgei*`
  members it `= delete`s both the copy constructor and copy assignment
  ([Graph.h](libraries/pg/include/mcrl2/pg/Graph.h#L267)) and provides a `noexcept`
  `swap`, so there is no shallow-copy double-free — the rule-of-three is handled
  deliberately rather than by accident.
- **Cooperative cancellation.** The `Abortable` mixin lets long solves be aborted;
  the recursive solver checks `aborted()` at entry and returns `false`, and callers
  test the returned strategy for emptiness before using it.
- **Clean architectural isolation.** The whole library is reachable only through
  `pbespgsolve`, the binary save/load paths are symmetric, the internal Doxygen
  `\file` style is consistent, and subgraph construction has an optional
  `MCRL2_ENABLE_MULTITHREADING` parallel path.

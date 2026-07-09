---
name: mcrl2-architecture
description: "High-level architecture of the mCRL2 toolset: library layout, the responsibilities of each library, and how tools are organized. Use when you need to know which library owns a concept, where to add reusable code, the standard library directory structure, or the meaning of the release/experimental/developer/deprecated tool tiers. Keywords: library layout, atermpp, data, lps, lts, pbes, pres, process, modal_formula, pg, smt, symbolic, utilities, gui, tools, namespace."
---

# mCRL2 Architecture

## Library Layout

Each library in `libraries/<name>/` follows:
```
<name>/
├── CMakeLists.txt
├── include/mcrl2/<name>/       # Public headers (namespace mcrl2::<name>)
│   └── detail/                 # Internal implementation (namespace ::detail)
├── source/                     # Compiled .cpp files
└── test/                       # Boost.Test unit tests (*_test.cpp)
```

## Key Libraries

| Library | Purpose |
|---------|---------|
| `atermpp` | Foundation: maximally-shared immutable terms (`aterm`, `aterm_list`) |
| `core` | Parser integration (dparser), identifiers |
| `data` | Data types, sorts, expressions, rewriters, traversers |
| `process` | Process algebra expressions |
| `lps` | Linear Process Specifications |
| `lts` | Labeled Transition Systems, bisimulation algorithms |
| `modal_formula` | State/action/regular formulas |
| `pbes` | Parameterized Boolean Equation Systems |
| `pres` | Parameterized Real Equation Systems |
| `pg` | Parity game solvers |
| `smt` | SMT solver integration |
| `symbolic` | LDD/Sylvan-based symbolic computation |
| `utilities` | Tool base classes, logging, CLI parsing, thread utilities |
| `gui` | Qt6 GUI support code |

## Tools

Tools live in `tools/{release,experimental,deprecated,developer}/`. Each tool is a small executable inheriting from the tool base classes in `utilities`. **Reusable code must go in a library, not a tool.**

| Tier | Meaning |
|------|---------|
| `release` | Stable, supported tools shipped to users |
| `experimental` | Research/in-progress tools (built with `MCRL2_ENABLE_EXPERIMENTAL`) |
| `developer` | Infrastructure/testing tools for maintainers (built with `MCRL2_ENABLE_DEVELOPER`) |
| `deprecated` | Legacy tools kept for compatibility (built with `MCRL2_ENABLE_DEPRECATED`) |

## Dependency Direction

The rough dependency stack (bottom builds up): `atermpp` → `core` / `utilities` → `data` → `process` → `lps` → `lts` / `modal_formula` → `pbes` / `pres` → `pg`. `smt` and `symbolic` plug into `data`/`lps`/`pbes` where solving or decision-diagram support is needed.

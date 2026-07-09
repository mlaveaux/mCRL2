---
name: mcrl2-code-style
description: "Coding conventions and style rules for mCRL2 C++ code. Use when writing or editing C++ source/headers, naming things, adding header guards, ordering includes, or adding the required file preamble. Keywords: clang-format, clang-tidy, snake_case, header guard, m_ prefix, file preamble, copyright, include order, namespace."
---

# mCRL2 Code Style

Enforced by `.clang-format` — 2-space indentation, 120-column limit, Allman-style braces, `PointerAlignment: Left`, `template<>` without a space. Run `clang-format` on touched code before committing; do not reformat unrelated code. Clang-tidy checks are in `.clang-tidy`; run `run-clang-tidy -p build/ -fix -format` to apply fixes.

## Conventions

- **Naming**: `snake_case` everywhere. Template parameters start uppercase. Macros are `UPPERCASE` with an `MCRL2_` prefix.
- **Private members**: `m_` prefix (e.g., `m_rewriter`).
- **Header guards**: Derived from file path under `include/`. Example: `include/mcrl2/data/variable.h` → `MCRL2_DATA_VARIABLE_H`. The `#endif` must have a trailing comment: `#endif // MCRL2_DATA_VARIABLE_H`.
- **Namespaces**: `mcrl2::<library>` for public API; internal details in `mcrl2::<library>::detail`. No `using namespace` in headers.
- **Includes**: Quotes for project headers (`"mcrl2/data/variable.h"`). Sorted case-sensitively per `.clang-format`. Keep includes minimal and consistent with neighboring code.
- **Standalone headers**: Every header must compile standalone (a `.cpp` with only `#include "header.h"` must compile).
- **Macros**: Avoid macros; prefer `constexpr` variables. When unavoidable, define macros in the build system (`cmake/ConfigureCompiler.cmake`), not in source code.
- **Contracts**: Public APIs document preconditions/postconditions with Doxygen `\pre`/`\post`, enforce cheap checks with `assert`, and put expensive checks behind `#ifndef MCRL2_NO_SOUNDNESS_CHECKS`.

## File Preamble

Every header and source file must start with:

```cpp
// Author(s): <NAMES>
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file <path/to/file.h>
/// \brief <Brief description>
```

## Commits

Bug-fix commits must reference the ticket in the commit message, formatted as `fixes #<number>`.

See `doc/sphinx/developer_manual/guidelines.rst` for the full coding guidelines.

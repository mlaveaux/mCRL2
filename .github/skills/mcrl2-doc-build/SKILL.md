---
name: mcrl2-doc-build
description: "Build and troubleshoot mCRL2 documentation. Use when asked to generate Sphinx/Doxygen docs, fix doc build issues, or validate documentation changes."
---

# mCRL2 Documentation Build Skill

Use this workflow for documentation tasks.

## 1. Prepare Python dependencies

```bash
pip install -r requirements.txt
```

A virtual environment is recommended when working locally.

## 2. Configure docs in CMake
Enable documentation and tune optional components as needed:
- `MCRL2_ENABLE_DOCUMENTATION=ON`
- `MCRL2_ENABLE_DOC_DOXYGEN`
- `MCRL2_ENABLE_DOC_PDFLATEX`
- `MCRL2_ENABLE_DOC_MANUAL`

## 3. Build docs
Build from the configured build directory:

```bash
cmake --build build --target doc
```

For quicker iterations:

```bash
cmake --build build --target fastdoc
```

## 4. Diagnose failures
When doc builds fail, report missing tools clearly (for example Doxygen or LaTeX packages) and suggest the smallest configuration change that unblocks progress.

## 5. Validate outputs
Confirm the expected artifacts were generated and point to their build output locations.

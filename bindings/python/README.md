# Building

This contains the Python bindings for mCRL2 implemented using `pybind11`. First of all, acquire the necessary submodules using the following command.

```
git submodule update --init --recursive
```

These bindings can either be generated as part of the CMake configuration by enabling `MCRL2_ENABLE_PYTHON_BINDINGS` and building the `mcrl2py` target.

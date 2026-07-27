# pyR2D2 packaging and build policy

This document records the supported build paths and the decisions made while
retiring the Fortran extension modules.

## Supported installation path

Use the PEP 517/setuptools build defined by `pyproject.toml` and `setup.py`:

```bash
python -m pip install .
python -m pip install -e .
```

The former package-local CMake/Makefile path has been removed so that there is
only one authoritative extension build configuration.

## Compiler policy

- C++20 is required.
- Linux builds enable OpenMP by default and therefore require a compiler and
  runtime providing `libgomp` or an equivalent OpenMP implementation.
- Portable builds do not use `-march=native`.
- A machine-local build may opt in with `PYR2D2_NATIVE=1`.
- Absolute build-environment RPATH entries are removed from the extension link
  command. A wheel intended for distribution must still be inspected in a
  clean build environment with a platform tool such as `auditwheel`.

## Python and optional dependencies

- Supported Python: 3.11 or newer.
- NumPy and SciPy are core dependencies.
- Zarr support: `python -m pip install '.[zarr]'`.
- Plotting environment: `python -m pip install '.[plot]'`.
- All optional features: `python -m pip install '.[all]'`.

Zarr is imported only when a Zarr operation is requested. Matplotlib is not
imported by the pyR2D2 package itself; the plotting extra is provided for the
examples and user analysis scripts.

## Distribution contents

Binary wheels contain Python modules, the compiled extension, and the JSON
files required at runtime. C++ sources and headers are kept in the source
distribution, but are not installed as wheel package data.

## Cleaning generated extensions

```bash
python setup.py clean --all
```

The custom clean command also removes in-place `cpp_util` extension binaries
for other Python ABIs, preventing stale `.so` or `.pyd` files from accumulating
in the source tree.

## Release checks

Before committing a packaging change:

1. Run `python -m pytest -q`.
2. Run `python setup.py clean --all`.
3. Build an sdist and wheel with `python -m build --no-isolation`.
4. Confirm the sdist includes every C++ header needed by `bindings.cpp`.
5. Confirm the wheel excludes C++ sources, headers, Makefiles, and retired
   modules.
6. Inspect the wheel extension for absolute RPATH entries and unexpected
   dynamic libraries.

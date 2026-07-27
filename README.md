# pyR2D2
Python code for analyzing results of R2D2 simulation.

## Installation

pyR2D2 requires Python 3.11 or newer and a C++20 compiler. On Linux the C++
extension is built with OpenMP support.

```bash
git clone git@github.com:hottahd/pyR2D2.git
cd pyR2D2
python -m pip install .
```

Install optional Zarr and plotting dependencies when needed:

```bash
python -m pip install '.[zarr]'
python -m pip install '.[plot]'
python -m pip install '.[all]'
```

Portable builds do not use `-march=native`. For a build that will only run on
the current machine, native CPU optimization can be requested explicitly:

```bash
PYR2D2_NATIVE=1 python -m pip install .
```

See [PACKAGING.md](PACKAGING.md) for build and wheel policies.

## Distribution

pyR2D2 is intended for R2D2 collaborators and is marked as a private package.
Do not upload it to PyPI or redistribute it outside the project rules.

## Quick start
You can generate `R2D2.Data` class instance as:

```python
import pyR2D2
datadir = '../run/d001/'
d = pyR2D2.Data(datadir)
```

## Documentation

https://hottahd.github.io/pyR2D2/master

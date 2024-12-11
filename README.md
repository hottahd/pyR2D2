# pyR2D2
Python code for analyzing results of R2D2 simulation.

## Installation

```bash
git clone git@github.com:hottahd/pyR2D2.git
cd pyR2D2
pip install .
```

## Quick start
You can generate `R2D2.Data` class instance as:

```python
import pyR2D2
datadir = '../run/d001/'
d = pyR2D2.Data(datadir)
```

## Utilities by Fortran

`pyR2D2` provides several utilities written in Fortran language.

You need to `make` at `pyR2D2/` directory to use these utilities.

## Documentation

https://hottahd.github.io/pyR2D2/master

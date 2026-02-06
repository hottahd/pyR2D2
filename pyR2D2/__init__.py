"""
This is a package for dealing with R2D2 simulation data.

Example
-------
A :class:`pyR2D2.Data` instance can be created as follows:

.. code-block:: python

    import pyR2D2
    datadir = '../run/d001/data/'
    d = pyR2D2.Data(datadir)

"""

from . import cpp_util, fortran_util, util, write
from .color import color
from .constant import constant
from .cpp_util import EOS
from .data import Data
from .data_io import zarr_util
from .data_io.parameters import Parameters
from .data_io.read import (
    FullData,
    ModelS,
    MPIRegion,
    OnTheFly,
    OpticalDepth,
    RestrictedData,
    Slice,
    XSelect,
    ZSelect,
)
from .sync.sync import Sync

__all__ = [
    "Data",
    "Parameters",
    "XSelect",
    "ZSelect",
    "MPIRegion",
    "FullData",
    "RestrictedData",
    "OpticalDepth",
    "OnTheFly",
    "Slice",
    "ModelS",
    "Sync",
    "color",
    "constant",
    "write",
    "util",
    "fortran_util",
    "cpp_util",
]

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"

# from setuptools_scm import get_version

# try:
#     __version__ = '.'.join(get_version(root="..", relative_to=__file__).split('.')[:3])
# except Exception as e:
#     __version__ = "unknown"

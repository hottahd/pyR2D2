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
from .data import Data
from .data_io.parameters import Parameters
from .data_io.read import XSelect, ZSelect, MPIRegion, FullData, RestrictedData, OpticalDepth, OnTheFly, Slice, ModelS
from .sync.sync import Sync
from .color import color
from .constant import constant
from . import write
from . import util
from . import fortran_util

__all__ = ['Data',
           'Parameters',
           'XSelect',
           'ZSelect',
           'MPIRegion',
           'FullData',
           'RestrictedData',
           'OpticalDepth',
           'OnTheFly',
           'Slice',
           'ModelS',
           'Sync',
           'sync',
           'color',
           'constant',
           'write',
           'util',
           ]

try:
    from ._version import version as __version__
except ImportError:
    __version__ = 'unknown'

# from setuptools_scm import get_version

# try:
#     __version__ = '.'.join(get_version(root="..", relative_to=__file__).split('.')[:3])
# except Exception as e:
#     __version__ = "unknown"

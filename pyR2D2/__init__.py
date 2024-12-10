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
from .read.read import Read
from .sync.sync import Sync
from .color import color
from .constant import constant
from . import write
from . import util
from . import fortran_util

__all__ = ['Data','Read','Sync','color','constant','write','util']

from setuptools_scm import get_version

try:
    __version__ = '.'.join(get_version(root="..", relative_to=__file__).split('.')[:3])
except Exception as e:
    __version__ = "unknown"

"""
This is a package for dealing with R2D2 simulation data.

Example
-------
A :class:`R2D2.Data` instance can be created as follows:

.. code-block:: python

    import R2D2
    datadir = '../run/d001/data/'
    d = R2D2.Data(datadir)
        
"""
from .data import Data
from .read.read import Read
from .sync.sync import Sync
from .color import color
from .constant import constant
from . import write
from . import util

__all__ = ['Data','Read','Sync','color','constant','write','util']

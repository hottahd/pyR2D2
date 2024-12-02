"""
R2D2 package

Deals with R2D2 data.

Example
-------
A R2D2.Data instance can be created as follows:

.. code-block:: python

    import R2D2
    datadir = '../run/d001/data/'
    d = R2D2.Data(datadir)
        
"""
from .data import Data
from .read import Read
from .sync import Sync
from .color import Color
from .constant import Constant
from . import write
from . import util
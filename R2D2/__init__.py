"""
R2D2 package

Deals with R2D2 data.

Example
-------
A R2D2_data instance can be created as follows:

.. code-block:: python

    import R2D2
    datadir = '../run/d001/data/'
    d = R2D2.R2D2_data(datadir)
        
"""
from .r2d2_data import R2D2_data
from .r2d2_read import R2D2_read
from .r2d2_sync import R2D2_sync
from .color import Color
from .constant import Constant
from . import write
from . import util
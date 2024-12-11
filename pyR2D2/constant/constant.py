"""
Provides Sun-related constants

.. table:: Provided constants
    :widths: 25 50 25
    :align: left

    ======== ========================== =========
    Name     description                value    
    ======== ========================== =========
    MSUN     Solar mass [g]             1.988e33
    RSUN     Solar radius [cm]          6.957e10
    LSUN     Solar luminosity [erg/s]   3.828e33
    ASUN     Solar age [yr]             4.570e9
    ======== ========================== =========

Example
-------

.. code-block:: python

    import pyR2D2
    print(pyR2D2.constant.lsun)

"""
# Solar parameters
MSUN = 1.988e33 # g
RSUN = 6.957e10 # cm
LSUN = 3.828e33 # erg/s
ASUN = 4.570e9 # yr
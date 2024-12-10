"""
Provides Sun-related constants

.. table:: Provided constants
    :widths: 25 50 25
    :align: left

    ======== ========================== =========
    Name     description                value    
    ======== ========================== =========
    mass     Solar mass [g]             1.988e33
    rsun     Solar radius [cm]          6.957e10
    lsun     Solar luminosity [erg/s]   3.828e33
    asun     Solar age [yr]             4.570e9
    ======== ========================== =========

Example
-------

.. code-block:: python

    import pyR2D2
    print(pyR2D2.constant.lsun)

"""
# Solar parameters
msun = 1.988e33 # g
rsun = 6.957e10 # cm
lsun = 3.828e33 # erg/s
asun = 4.570e9 # yr
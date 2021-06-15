Heat Index
======================================

The Heat index, takes into account air temperature and relative humidity, \
to be an indication of how hot it feels.

The method of this index uses a multiple regression and Apparent Temperature's
calculation of relative humidity.

How To Use
======================================
You need 2m temperature in kelvin and optional relative humidity
as water vapour pressure, because this can be calculated from 2m temperature.

.. code-block:: python
    calculate_heat_index(2m temperature,relative humidity)


Interpret the Output
======================================

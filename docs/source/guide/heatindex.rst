Heat Index
======================================

The Heat index, takes into account air temperature and relative humidity, \
to be an indication of how hot it feels.

The method of this index uses a multiple regression and Apparent Temperature's
calculation of relative humidity.

Here we have two methods: Heat Index Simplified and Heat Index Adjusted.
As is set out by https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml \

We carry out the calculation in Fahrenheit and then convert to Celsius in keeping \
with the other thermal indices in this library.

How To Use
------------------
You need 2m temperature in Kelvin and, optionally, relative humidity
such as saturation vapour pressure (because this can be directly calculated from 2m temperature)

.. code-block:: python

    calculate_heat_index_simplified(2m_temperature,relative_humidity)

use numpy array

You need 2m temperature and 2m dew point temperature in Kelvin.

.. code-block:: python

    calculate_heat_index_adjusted(2m_temperature,2m_dew_point_temperature)


Interpret the Output
----------------------
Here is a suggested way for you to interpret heat index outputs. However, it is by no means the only way to go about defining thermal stress.
Heat Index and Apparent Temperature have the same thresholds.

.. csv-table:: Heat Index Thresholds
    :file: atthresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1

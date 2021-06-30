Heat Index
======================================

The Heat index, takes into account air temperature and relative humidity, \
to be an indication of how hot it feels.

The method of this index uses a multiple regression and Apparent Temperature's
calculation of relative humidity.

Here we have two methods Heat Index Simplified and Heat Index Adjusted.
As is set out by https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml \

We carry out the calculation in fahrenheit and then convert to celcius in keeping \
with the other thermal indices in this library.

How To Use
======================================
You need 2m temperature in kelvin and optional relative humidity
as water vapour pressure, because this can be calculated from 2m temperature.

.. code-block:: python
    calculate_heat_index(2m temperature,relative humidity)


Interpret the Output
======================================
Here is a suggested way for you to interpret heat index outputs, it is by no means the only way to go about defining thermal stress.
Heat Index and Apparent Temperature have the same thresholds.

.. csv-table:: Heat Index Thresholds
    :file: "atthresholds.csv"
    :header-rows: 1
    :class: longtable
    :widths: 1 1

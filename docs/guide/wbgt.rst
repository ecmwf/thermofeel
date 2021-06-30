Wet Bulb Globe Temperature
======================================
This is one of the most used heat indexes and was developed in the US by the Army and Marines.
Traditionally it is calculated using natural wet-bulb temperature, globe temperature and dry bulb temperature.

Here we present a contemporary WBGT method from *De Dear* calculated using Mean Radiant Temperature and one of WBGT
approximations from the Australian Bureau of Meteorology.

More Information: https://www.sciencedirect.com/science/article/abs/pii/S0378778817335971?via%3Dihub

How To Use
======================================

**Wet Bulb Globe Temperature Simple/ Approximation**

You need 2m Temperature in Kelvin

.. code-block:: python

   calculate_wbgts(2m temperature)

**Wet Bulb Globe Temperature**

**This method is not tested for Windows**

You need 2m Temperature in Kelvin, Mean Radiant Temperature in Kelvin and 10 meter height wind speed.

.. code-block:: python

   calculate_wbgt(2m temperature, mean radiant temperature, wind speed)

Interpret the Output
======================================

Here is a suggested way for you to interpret wet bulb globe temperature outputs, it is by no means the only way to go about defining thermal stress.
These are based upon the Wet Bulb Globe Temperature and as such might have a different accuracy for the approximation.

.. csv-table:: WBGT Thresholds
    :file: wbgtthresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1

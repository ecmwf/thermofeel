Wet Bulb Globe Temperature
======================================
The wet bulb globe temperature is one of the most used heat indexes and was developed in the US by the Army and Marines.
Traditionally it is calculated using natural wet-bulb temperature, globe temperature and dry bulb temperature.

Here we present a contemporary WBGT method that uses globe temperature from *De Dear* calculated using Mean Radiant Temperature and WBGT from Stull et al. 2011 and one of WBGT
approximations from the Australian Bureau of Meteorology.

More Information: https://www.sciencedirect.com/science/article/abs/pii/S0378778817335971?via%3Dihub

How To Use
---------------

**Wet Bulb Globe Temperature Simple**

You need 2m Temperature in Kelvin. It returns the wet bulb globe temperature in Celsius.

.. code-block:: python

   calculate_wbgts(2m_temperature)

**Wet Bulb Temperature**

You need 2m temperature in Celsius and relative humidity percent. It returns wet bulb temperature in Celsius.

.. code-block:: python

    calculate_wbt(2m_temperature, relative_humidity_percent)

**Globe Temperature**

You need 2m temperature in Kelvin, mean radiant temperature and 10m wind speed

.. code-block:: python

    calculate_gbt(2m_temperature, mean_radiant_temperature, 10m_wind_speed)

**Wet Bulb Globe Temperature**

**This method is not tested for Windows**

You need 2m temperature in Kelvin, mean radiant temperature in Kelvin and 10 m wind speed in m/s and  2m dew point temperature in Celsius. It returns wet bulb globe temperature in Celsius

.. code-block:: python

   calculate_wbgt(2m_temperature, mean_radiant_temperature, 10m_wind_speed)


Interpret the Output
---------------------

Here is a suggested way for you to interpret wet bulb globe temperature outputs. However, it is by no means the only way to go to classify thermal stress.
These are based upon the wet bulb globe temperature and as such might have a different accuracy for the approximation.

.. csv-table:: WBGT Thresholds
    :file: wbgtthresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1

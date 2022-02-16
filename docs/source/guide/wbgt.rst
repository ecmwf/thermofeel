Wet Bulb Globe Temperature
======================================
The wet bulb globe temperature is one of the most used heat indexes and was developed in the US by the Army and Marines *Minard (1961)*.
Traditionally it is calculated using natural wet-bulb temperature, globe temperature and dry bulb temperature.

Here we present a contemporary WBGT method that uses globe temperature from *De Dear (1987)* calculated using Mean Radiant Temperature and one of WBGT
approximations from the Australian Bureau of Meteorology.

More Information:

- https://www.sciencedirect.com/science/article/abs/pii/S0378778817335971?via%3Dihub

- De Dear, R. (1987). Ping-pong globe thermometers for mean radiant temperatures. Heating and Ventilation Engineer and Journal of Air Conditioning, 60, 10–11. Retrieved from https://ci.nii.ac.jp/naid/10030966825

- Minard, D. (1961). Prevention of heat casualties in Marine Corps recruits. Period of 1955-60, with comparative incidence rates and climatic heat stresses in other training categories. Military Medicine, 126(4), 261–272. https://doi.org/10.1093/milmed/126.4.261

How To Use
---------------

**Wet Bulb Globe Temperature Simple**

You need 2m temperature in Kelvin.

It returns the wet bulb globe temperature in Celsius.

.. code-block:: python

   calculate_wbgts(2m_temperature)

**Wet Bulb Temperature**

You need 2m temperature in Celsius and relative humidity percent.

It returns the wet bulb temperature in Celsius.

.. code-block:: python

    calculate_wbt(2m_temperature, relative_humidity_percent)

**Globe Temperature**

You need 2m temperature and mean radiant temperature in Kelvin and 10m wind speed in m/s.

It returns the bulb globle temperature in Celsius. 

.. code-block:: python

    calculate_bgt(2m_temperature, mean_radiant_temperature, 10m_wind_speed)

**Wet Bulb Globe Temperature**

**This method is not tested for Windows**

You need 2m temperature and mean radiant temperature in Kelvin and 10 m wind speed in m/s.

It returns the wet bulb globe temperature in Celsius.

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

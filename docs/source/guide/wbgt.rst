Wet Bulb Globe Temperature
======================================
The Wet Bulb Globe Temperature (WBGT) represents the thermal environment to which an individual is exposed and 
its value gives a first approximation of the heat stress on a person. 
The WBGT is an International Organization for Standardization (ISO) screening method to establish 
the presence or absence of heat stress. 

Traditionally it is calculated using natural wet-bulb temperature, globe temperature and dry bulb temperature.

Here we present a contemporary WBGT method that uses globe temperature from *De Dear (1987)* calculated using 
the mean radiant temperature and one of WBGT approximations from the Australian Bureau of Meteorology.

More information:
* Australian Government Bureau of Meteorology. Thermal Comfort observations - About the approximation to the WBGT used by the Bureau of Meteorology. https://www.bom.gov.au/info/thermal_stress/#approximation
* De Dear, R. Ping-pong globe thermometers for mean radiant temperatures. H and V Engineer 60(681), 10–11 (1988)
* Guo, H., Teitelbaum, E., Houchois, N., Bozlar, M., & Meggers, F. Revisiting the use of globe thermometers to estimate radiant temperature in studies of heating and ventilation. Energy Build 180, 83-94 (2018). https://doi.org/10.1016/j.enbuild.2018.08.029
* Minard, D. Prevention of heat casualties in Marine Corps recruits. Period of 1955-60, with comparative incidence rates and climatic heat stresses in other training categories. Mil Med 126(4), 261–272 (1961). https://doi.org/10.1093/milmed/126.4.261
* Stull, R. Wet-bulb temperature from relative humidity and air temperature. J Appl Meteorol Climatol 50(11), 2267-2269 (2011). https://doi.org/10.1175/JAMC-D-11-0143.1

How To Use
---------------

**Wet Bulb Globe Temperature Simple**

You need 2m temperature in Kelvin and relative humidity as a percentage.

It returns the WBGT in Kelvin.

.. code-block:: python

   calculate_wbgt_simple(2m_temperature, relative_humidity_percent)

**Wet Bulb Temperature**

You need 2m temperature in Kelvin and relative humidity as a percentage.

It returns the wet bulb temperature in Kelvin.

.. code-block:: python

    calculate_wbt(2m_temperature, relative_humidity_percent)

**Globe Temperature**

You need 2m temperature and mean radiant temperature in Kelvin and 10m wind speed in meters per second.

It returns the bulb globe temperature in Kelvin. 

.. code-block:: python

    calculate_bgt(2m_temperature, mean_radiant_temperature, 10m_wind_speed)

**Wet Bulb Globe Temperature**

You need 2m air temperature, 2m dew point temperature and mean radiant temperature in Kelvin and 10 m wind speed in meters per second.

It returns the WBGT in Kelvin.

.. code-block:: python

   calculate_wbgt(2m_temperature, mean_radiant_temperature, 10m_wind_speed, 2m_dew_point_temperature)


Interpret the Output
---------------------

The Wet Bulb Globe Temperature is described in terms of the heat-related risk to human health. Scales are generally tailored 
to the geographical area of interest. One scale that has been adopted at the global scale is:

.. csv-table:: WBGT Thresholds
    :file: wbgtthresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1

Normal Effective Temperature (NET)
======================================

Normal effective temperature (NET) links effective temperature, which indicates \
the effects on comfort through air temperature and relative humidity, \
and an organism’s thermoregulatory capacity.

More Information: https://www.sciencedirect.com/topics/engineering/effective-temperature

How To Use
-----------------
You need 2m temperature  and 2m dew point temperature in Kelvin and 10 m wind speed in m/s.
The wind speed in this method is converted to 2m wind speed as an approximation of 1.2m wind speed.

It returns the normal effective temperature in Celsius.

.. code-block:: python

    calculate_net_effective_temperature(2m_temperature, 10m_wind_speed, 2m_dew_point_temperature)

Interpret the Output
------------------------
Here is a suggested way for you to interpret Net Effective Temperature outputs. However, it is by no means the only way to go about defining thermal stress.

.. csv-table:: NET Thresholds
    :file: netthresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1

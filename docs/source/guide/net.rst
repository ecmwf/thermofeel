Normal Effective Temperature
======================================

The Normal Effective Temperature (NET) is defined as the temperature felt by a human body and 
can indicate thermal exchange between the human body and the environment.

More information: Li, P. W., & Chan, S. T. Application of a weather stress index for alerting the public to stressful weather in Hong Kong. Meteorol Appl 7(4), 369-375 (2000). https://doi.org/10.1017/S1350482700001602

How To Use
-----------------
You need 2m air temperature in Kelvin, 10 m wind speed in meters per second and relative humidity as a percentage.

It returns the NET in Kelvin.

.. code-block:: python

    calculate_net_effective_temperature(2m_temperature, 10m_wind_speed, relative_humidity_percent)

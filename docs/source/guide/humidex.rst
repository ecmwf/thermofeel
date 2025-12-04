Humidex
======================================

The Humidex is defined as the temperature (in °C) the human body perceives in hot, humid weather.

More information: Blazejczyk, K., Epstein, Y., Jendritzky, G. et al. Comparison of UTCI to selected thermal indices. Int J Biometeorol 56, 515–535 (2012). https://doi.org/10.1007/s00484-011-0453-2

How To Use
------------------
You need 2m air temperature and 2m dew point temperature in Kelvin.

It returns the Humidex in Kelvin. 

.. code-block:: python

   calculate_humidex(2m_temperature, 2m_dew_point_temperature)


Interpret the Output
-----------------------
The Humidex is described in terms of comfort.

.. csv-table:: Humidex Thresholds
    :file: humidexthresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1

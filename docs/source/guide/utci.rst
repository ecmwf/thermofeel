Universal Thermal Climate Index
======================================
The Universal Thermal Climate Index (UTCI) is a measure of the thermal stress on the human body from outdoor conditions,
defined as the equivalent air temperature of a reference environment that would cause the same physiological response
as the actual conditions. 

More information: Jendritzky, G., de Dear, R. & Havenith, G. UTCI—Why another thermal index?. Int J Biometeorol 56, 421–428 (2012). https://doi.org/10.1007/s00484-011-0513-7


How To Use
----------------

You need 2m air temperature and mean radiant temperature in Kelvin, 2m dew point temperature in Kelvin or 
water vapour pressure in hPa, and 10m wind speed in m/s.
Please use with numpy arrays.

The UTCI is returned in Kelvin.

.. code-block:: python

    rh_pc = calculate_relative_humidity_percent(2m_temperature, 2m_dew_point_temperature)
    ehPa = calculate_saturation_vapour_pressure(2m_temperature) * rh_pc / 100.0
    utci = calculate_utci(2m_temperature, 10m_wind_speed, mean_radiant_temperature, 2m_dew_point_temperature=None, ehPa=None)


Interpret the Output
-------------------------
The UTCI is classified into a 10-category scale. Each category, defined by a specific range of UTCI values, corresponds 
to a well-defined set of human physiological responses to the outdoor environment.

.. csv-table:: UTCI Thresholds
    :file: utcithresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1

Relative Humidity and Saturation Vapour Pressure
================================================

There are two methods to compute Relative Humidity used in *thermofeel* and these are used together in some \
index calculations, such as in the UTCI.

Relative Humidity Percent: https://www.theweatherprediction.com/habyhints/186/
Saturation Vapour Pressure over Water: http://www.thunderscientific.com/tech_info/reflibrary/its90formulas.pdf

How To Use
--------------

**Saturation Vapour Pressure over Water**

The input is 2m temperature in Kelvin. The output is in hPa. 

.. code-block:: python

   calculate_saturation_vapour_pressure(2m_temperature)

**Relative Humidity Percent**

The inputs are 2m temperature and dew point temperature in Kelvin. The output is in %. 

.. code-block:: python

   calculate_relative_humidity_percent(2m_temperature,2m_dew_point_temperature)




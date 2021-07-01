Relative Humidity and Water Vapour Pressure
======================================

There are two types of Relative Humidity used in *thermofeel* and sometimes such as in the UTCI these are \
be used together.

Relative Humidity Percent:https://www.theweatherprediction.com/habyhints/186/
Saturation/Water Vapour Pressure: http://www.thunderscientific.com/tech_info/reflibrary/its90formulas.pdf

How To Use
======================================

**saturation/water vapour pressure**
You will only need 2m temperature

.. code-block:: python
    calculate_saturation_vapour_pressure(2m temperature)

***Relative Humidity Percent**

.. code-block:: python
    calculate_relative_humidity_percent(2m temperature,dew point temperature)




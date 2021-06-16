Apparent Temperature
======================================
Was developed by *Steadman,1984* to describe the thermal comfort/resistance of an average adult walking
when they are exposed to certain combination of temperatures, relative humidity's and wind speed.

The method used here follows the Apparent Temperature method used by the *Australian Bureau of Meteorology*
and *Blazejczyk et al., 2011 Int J Biometeorol*

Apparent temperature is similar to the Heat Index method also available to calculate in *thermofeel*

How To Use
======================================
You need 2m temperature in kelvin, wind speed at 10 meters
and optional relative humidity as water vapour pressure,
because this can be calculated from 2m temperature.

The wind speed in this method is converted to 2 meters as
an approximation of 1.2 meter wind speed.

.. code-block:: python
    calculate_heat_index(2m temperature,relative humidity,wind speed)


Interpret the Output
======================================

Here is a suggested way for you to interpret apparent temperature outputs, it is by no means the only way to go about defining thermal stress.

.. csv-table:: Apparent Temperature Thresholds
   :file: "atthresolds.csv"
   :widths: 30, 70
   :header-rows: 1


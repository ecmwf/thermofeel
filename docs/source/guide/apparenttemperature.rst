Apparent Temperature
======================================
Was developed by *Steadman,(1984)* to describe the thermal comfort/resistance of an average adult walking
when they are exposed to certain combination of temperatures, relative humidity's and wind speed.

The method used here follows the Apparent Temperature method used by the *Australian Bureau of Meteorology*
and *Blazejczyk et al. (2011)*

Apparent temperature is similar to the Heat Index methods also available to calculate in *thermofeel*

More Information: 

- https://journals.ametsoc.org/view/journals/apme/23/12/1520-0450_1984_023_1674_ausoat_2_0_co_2.xml

- https://link.springer.com/article/10.1007/s00484-011-0453-2 


How To Use
-----------
You need 2m temperature in Kelvin, 10m wind speed in m/s and, optionally, relative humidity such as water vapour pressure 
(because this can be calculated from 2m temperature). The wind speed in this method is converted to 
2 meters as an approximation of 1.2 meter wind speed.

It returns the apparent temperature in Kelvin

.. code-block:: python

   calculate_apparent_temperature(2m_temperature, 10m_wind_speed, relative_humidity)
    

Interpret the Output
--------------------

Here is a suggested way for you to interpret apparent temperature outputs, it is by no means the only way to go about defining thermal stress.

.. csv-table:: Apparent Temperature
    :file: atthresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1



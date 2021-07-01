Net Effective Temperature (NET)
======================================

Net effective temperature (NET) is also known as normal effective temperature. It links effective temperature which indicates \
the effects on comfort through air temperature and relative humidity \
and an organism’s thermoregulatory capacity.

More Information:https://www.sciencedirect.com/topics/engineering/effective-temperature

How To Use
======================================
You need 2m temperature in kelvin, wind speed at 10 meters height and 2m dew point temperature.

The wind speed in this method is converted to 2 meters as
an approximation of 1.2 meter wind speed.

.. code-block:: python
    calculate_net_effective_temperature(2m temperature,wind speed, 2m dew point temperature)

Interpret the Output
======================================
Here is a suggested way for you to interpret Net Effective Temperature outputs, it is by no means the only way to go about defining thermal stress.

.. csv-table:: NET Thresholds
    :file: netthresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1
Humidex
======================================

Humidex, aims to describe how hot and humid weather feels to an average person, \
and was developed by the Canadian meteorological service.


How To Use
======================================
You need 2m temperature in kelvin and dew point temperature in kelvin.

.. code-block:: python
    calculate_humidex(2m temperature,dew point temperature)


Interpret the Output
======================================

.. csv-table:: Humidex Thresholds
   :file: humidexthresolds.csv
   :widths: 30, 70
   :header-rows: 1

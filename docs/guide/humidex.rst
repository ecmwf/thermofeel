Humidex
======================================

Humidex, aims to describe how hot and humid weather feels to an average person, \
and was developed by the Canadian meteorological service.
More information: http://www.csgnetwork.com/canhumidexcalc.html



How To Use
======================================
You need 2m temperature in kelvin and dew point temperature in kelvin.

.. code-block:: python
   calculate_humidex(2m temperature,dew point temperature)


Interpret the Output
======================================

.. csv-table:: Humidex Thresholds
    :file: humidexthresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1
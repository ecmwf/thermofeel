Wind Chill
======================================
Wind Chill is an indication of cold thermal conditions.
Originally developed by Siple and Passel (1945)
And updated by the Canadian Meteorological Service in 2001.

More information: http://www.ec.gc.ca/meteo-weather/default.asp?lang=n&n=5FBF816A-1#wc6

How To Use
======================================

You will need 2m temperature in kelvin and 10 meter height wind speed.

.. code-block:: python
    calculate_wind_chill(2m temperature, wind speed)

Interpret the Output
======================================
.. csv-table:: Wind Chill Thresholds
    :file: windchillthresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1
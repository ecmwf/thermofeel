Universal Thermal Climate Index (UTCI)
======================================
The universal thermal climate index (UTCI), is a bioclimatological model of an average human body’s response
to different thermal conditions, where the subject is
not acclimatised to the climate and is outdoors, doing minimal work.

More Information: https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/gdj3.102


How To Use
======================================
You will need 2m temperature in kelvin, mean radiant temperature in kelvin,
relative humidity (calculated as shown below,using dew point temperature) and 10 meter height wind speed.

.. code-block:: python
    rh_pc = tfc.calculate_relative_humidity_percent(self.t2m, self.td)
    ehPa = tfc.calculate_saturation_vapour_pressure(self.t2m) * rh_pc / 100.0
    calculate_utci(2m temperature,mean radiant temperature,
    ehPa, wind speed)


Interpret the Output
======================================

.. csv-table:: UTCI Thresholds
    :file: utcithresholds.csv
    :header-rows: 1
    :class: longtable
    :widths: 1 1
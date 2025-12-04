Apparent Temperature
======================================
The Apparent Temperature (AT) is defined as the temperature giving the same discomfort as under the current ambient 
temperature and humidity.  

the AT is based on a mathematical model of human body heat balance for an adult, walking outdoors in the shade.   

Absolute humidity conforming with a dew point of 14°C is chosen as a reference.  If ambient humidity is:
* higher than the reference humidity level, the apparent temperature will be higher than the ambient temperature.
* lower than the reference humidity level, the apparent temperature will be lower than the ambient temperature.

More information:
* Steadman, R. G. A universal scale of apparent temperature.  J Appl Meteorol Climatol 23(12), 1674-1687 (1984). https://doi.org/10.1175/1520-0450(1984)023%3C1674:AUSOAT%3E2.0.CO;2
* Australian Government Bureau of Meteorology. Thermal Comfort observations - About the formula for the apparent temperature. https://www.bom.gov.au/info/thermal_stress/#atapproximation

How To Use
-----------
You need 2m air temperature in Kelvin, 10m wind speed in meters per second and relative humidity as a percentage.

It returns the apparent temperature in Kelvin

.. code-block:: python

   calculate_apparent_temperature(2m_temperature, 10m_wind_speed, relative_humidity_percent)
    
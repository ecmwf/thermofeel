Mean Radiant Temperature
======================================
The Mean Radiant Temperature (MRT) is defined as that uniform temperature of a fictive black-body radiation enclosure 
(emission coefficient ε = 1) which would result in the same net radiation energy exchange with the subject as the actual, 
more complex radiation environment.
MRT is a key component of calculating indices such as the Wet Bulb Globe Temperature and Universal Thermal Climate Index.

More information: Di Napoli, C., Hogan, R.J. & Pappenberger, F. Mean radiant temperature from global-scale numerical weather prediction models. Int J Biometeorol 64, 1233–1245 (2020). https://doi.org/10.1007/s00484-020-01900-5

How To Use
------------

There are two different methods to calculate MRT in *thermofeel*.

**Mean Radiant Temperature**

It requires the following variables:
* ssrd is the surface solar radiation downwards
* ssr is the surface net solar radiation
* dsrp is the direct radiation from the Sun
* strd is the surface thermal radiation downwards
* fdir is the total sky direct solar radiation at surface
* strr is the surface net thermal radiation
* cossza is cosine of solar zenith angle

All the radiation variables are in :math:`W/{m}^{-2}` and the cosine of solar zenit angle is dimensionless. 
Please use numpy arrays.

It returns the mean radiant temperature in Kelvin.

.. code-block:: python

  calculate_mean_radiant_temperature(ssrd, ssr, dsrp, strd, fdir, strr, cossza)

**Mean Radiant Temperature from Globe Temperature**

2m air temperature and 2m wet bulb temperature are expressed in Kelvin and 10m wind speed in meters per second. 
Please use numpy arrays.

It returns the mean global radiant temperature in Kelvin.

.. code-block:: python

  calculate_mrt_from_bgt(2m_temperature, 2m_wet_bulb_globe_temperature, 10m_wind_speed)
  
  

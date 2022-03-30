Mean Radiant Temperature
======================================
Mean radiance temperature is the incidence of radiation on the body from all directions.
This is a key component of calculating Wet Bulb Globe Temperature and Universal Thermal Climate Index.

More information: https://link.springer.com/article/10.1007/s00484-020-01900-5

How To Use
------------

There are two different methods that will return slightly different values. This is to be anticipated.

**Mean Radiant Temperature**

Requires variables:
#. ssrd is surface solar radiation downwards [J/m^-2]
#. ssr is surface net solar radiation [J/m^-2]
#. dsrp is direct radiation from the Sun [J/m^-2]
#. strd is Surface thermal radiation downwards [J/m^-2]
#. fdir is Total sky direct solar radiation at surface [J/m^-2]
#. strr is Surface net thermal radiation [J/m^-2]
#. cossza is cosine of solar zenith angle

All the radiation variables are in :math:`J/{m}^{-2}` and the cosine of solar zenit angle in degrees. Please use numpy arrays.

It returns the mean radiant temperature in Kelvin.

.. code-block:: python

  calculate_mean_radiant_temperature(ssrd, ssr, dsrp, strd, fdir, strr, cossza)

**Mean Radiant Temperature from Globe Temperature**

2m temperature and 2m wet bulb temperaature are expressed in Kelvin and 10m wind speed in m/s. Please use numpy arrays.

It returns the mean global radiant temperature in Kelvin.

.. code-block:: python

  calculate_mrt_from_bgt(2m_temperature,10m_wind_speed,2m_wet_bulb_globe_temperature)
  
  

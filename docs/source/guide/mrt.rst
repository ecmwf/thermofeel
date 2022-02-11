Mean Radiant Temperature
======================================
Mean radiance temperature is the incidence of radiation on the body from all directions.
This is a key component of calculating Wet Bulb Globe Temperature and Universal Thermal Climate Index.

More information: https://link.springer.com/article/10.1007/s00484-020-01900-5

How To Use
------------

There are two different methods that will return slightly different values. This is to be anticipated.

**Mean Radiant Temperature**

All the radiation variables are in :math:`J/{m}^{-2}` and the cosine of solar zenit angle in degrees. Please use numpy arrays

.. code-block:: python

    calculate_mean_radiant_temperature(surface_solar_radiation_downwards,surface_net_solar_radiation,
    Total_sky_direct_solar_radiation_at_surface, Surface_thermal_radiation_downwards,Surface_net_thermal_radiation,
    cosine_of_solar_zenith_angle)

**Mean Radiant Temperature from Globe Temperature**

2m temperature and 2m wet bulb temperaature are expressed in Kelvin and 10m wind speed in m/s. 

.. code-block:: python

  calculate_mrt_from_bgt(2m_temperature,10m_wind_speed,2m_wet_bulb_globe_temperature)
  
  

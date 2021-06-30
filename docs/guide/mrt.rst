Mean Radiant Temperature
======================================
Is the incidence of radiation on the body from all directions.
It is a key component of calculating Wet Bulb Globe Temperature and Universal Thermal Climate Index.
More information: https://link.springer.com/article/10.1007/s00484-020-01900-5

How To Use
======================================

These methods will return slightly different values. This is to be anticipated.

**Mean Radiant Temperature**

.. code-block:: python

    calculate_mean_radiant_temperature(surface-solar-radiation-downwards,surface-net-solar-radiation,
    Total-sky-direct-solar-radiation-at-surface, Surface-thermal-radiation-downwards,Surface-net-thermal-radiation,
    cosine-of-solar-zenith-angle)

**Mean Radiant Temperature from Wet Bulb Globe Temperature**

.. code-block:: python
  calculate_mrt_from_wbgt(2m temperature,wind speed,wet bulb globe temperature)
  
  

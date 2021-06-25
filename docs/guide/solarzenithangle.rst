Cosine of the Solar Zenith Angle
======================================
is the cosine of the angle between the Sun and directly overhead during daylight hours
it is a key component of calculating mean radiant temperature from radiation.

Here 2 methods are presented one takes the Cosine of the Solar Zenith Angle Instantly and one is designed
to calculate the Solar Zenith Angle over a forecasting system time step.

They do have a difference in values but once this is used to calculate mean radiant temperature then UTCI or WBGT
the difference in the heat index output is in the mean about 0.004 degrees Celsius different.

How To Use
======================================
** Cosine Solar Zenith Angle Instantaneous**

.. code-block:: python
    calculate_solar_zenith_angle()


** Cosine Solar Zenith Angle Integrated**

.. code-block:: python
    calculate_solar_zenith_angle()


Interpret the Output
======================================
The maximum value is 0.9 the values should range between 0 and 0.9.
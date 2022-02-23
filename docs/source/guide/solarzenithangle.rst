Cosine of the Solar Zenith Angle
======================================
The Cosine of the Solar Zenith Angle is the cosine of the angle between the Sun and directly overhead
during daylight hours. It is a key component of calculating mean radiant temperature from radiation.

Here 2 methods are presented: one takes the Cosine of the Solar Zenith Angle Instantly and one is designed
to calculate the Solar Zenith Angle over a forecasting system time step.

There is a difference in output values of both methods, but once these values are used to calculate the
radiant heat temperature in either UTCI or WBGT, the difference in their heat index outputs is, in average, 0.004 degrees Celsius.

more information: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GL066868

How To Use
----------------

**Cosine Solar Zenith Angle Instantaneous**

Inputs are h (hour), latitude (in degrees), longitude (in degrees), y (year), m (month) and d (day).

It returns the average of cosine of the solar zenith angle in degrees.

.. code-block:: python

    calculate_solar_zenith_angle(h, lat, lon, y, m, d)


**Cosine Solar Zenith Angle Integrated**

Inputs are latitude (in degrees), longitude (in degrees), y (year), m (month), d (day), h (hour), tbegin (offset in 
hours from forecast time to begin of time interval for integration) and tend (offset in hours from forecast time
to end of time interval for integration). 

It returns the average of cosine of the solar zenith angle during interval in degrees.

.. code-block:: python

    calculate_cos_solar_zenith_angle_integrated(lat, lon, y, m, d, h, tbegin, tend)


Interpret the Output
-------------------------

The maximum value is 0.9 the values should range between 0 and 0.9.

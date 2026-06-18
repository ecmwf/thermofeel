Heat Force (Liljegren WBGT)
======================================
Heat force (Dutch: *hittekracht*) is a 0--10 scale that communicates
environmental heat stress to the public, introduced operationally by KNMI (the
Royal Netherlands Meteorological Institute) from June 2026. It is analogous to
the familiar wind-force and UV-index scales.

Heat force is derived from the Wet Bulb Globe Temperature (WBGT) computed with
the physically based **Liljegren method** (Liljegren et al., 2008), which is
regarded as the "gold standard" for deriving WBGT from standard meteorological
variables. The Liljegren method solves the steady-state energy balance of the
natural wet-bulb thermometer and the black globe by iteration, and combines them
with the air temperature as ``WBGT = 0.7*Tnw + 0.2*Tg + 0.1*Ta``.

More information:
* Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S., & Sharp, R. Modeling the Wet Bulb Globe Temperature Using Standard Meteorological Measurements. J Occup Environ Hyg 5(10), 645-655 (2008). https://doi.org/10.1080/15459620802310770
* Kong, Q. & Huber, M. Explicit Calculations of Wet-Bulb Globe Temperature Compared With Approximations and Why It Matters for Labor Productivity. Earth's Future 10 (2022). https://doi.org/10.1029/2021EF002334
* KNMI Technical Report TR-26-04. From Wet Bulb Globe Temperature (WBGT) to Heat Force (2026).

How To Use
---------------

**Wet Bulb Globe Temperature (Liljegren)**

You need 2m air temperature in Kelvin, relative humidity as a percentage,
surface air pressure in hectopascals, 10m wind speed in metres per second,
instantaneous downward shortwave radiation at the surface (SSRD) in W/m\ :sup:`2`,
the fraction of that radiation that is direct beam (0--1), and the cosine of the
solar zenith angle.

It returns the WBGT in Kelvin.

.. code-block:: python

   calculate_wbgt_liljegren(
       2m_temperature, relative_humidity_percent, pressure_hPa,
       10m_wind_speed, ssrd, fdir_fraction, cos_solar_zenith_angle
   )

The KNMI operational guards are applied internally: wind below 0.62 m/s at 10 m
is raised to that floor (then scaled to the 2 m sensor height); the direct-beam
fraction is clamped to [0, 0.9] and set to 0 when the sun is at or below 89.5
degrees zenith. The cosine of the solar zenith angle can be obtained from the
`earthkit-meteo <https://github.com/ecmwf/earthkit-meteo>`_ library. ``NaN`` is
returned where the iteration does not converge.

By default the 10 m wind is converted to 2 m with the KNMI/Liljegren
stability-dependent profile (``wind_scaling="liljegren"``); pass
``wind_scaling="brode"`` to use the generic :func:`scale_windspeed` log profile
instead. The Liljegren conversion is also available on its own:

.. code-block:: python

   calculate_wind_speed_2m_liljegren(10m_wind_speed, cos_solar_zenith_angle, ssrd)

**Heat Force**

You need the WBGT in Kelvin (for example from ``calculate_wbgt_liljegren``).

It returns the heat force as a whole number from 0 to 10.

.. code-block:: python

   calculate_heat_force(wbgt)


Interpret the Output
---------------------

Heat force uses fixed 2 °C bands of WBGT. Heat force 0 means no heat load; heat
force 10 is extremely rare in the current Dutch climate. Heat force describes
only the *environmental* heat load -- the actual health risk also depends on an
individual's "heat fitness" (health, activity level, clothing, acclimatisation).

.. csv-table:: Heat-force scale
    :header: "Heat force", "WBGT (°C)"
    :widths: 1 1

    "0", "< 14"
    "1", "14 -- 16"
    "2", "16 -- 18"
    "3", "18 -- 20"
    "4", "20 -- 22"
    "5", "22 -- 24"
    "6", "24 -- 26"
    "7", "26 -- 28"
    "8", "28 -- 30"
    "9", "30 -- 32"
    "10", "≥ 32"

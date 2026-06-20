# Overview of thermofeel

*thermofeel* is a thermal-indices library which allows these indices to be
integrated in numerical weather prediction systems. It was developed and is
maintained by ECMWF (European Centre for Medium-Range Weather Forecasts).

*thermofeel* calculates the following thermal indices:

- Apparent Temperature
- Heat Index
- Humidex
- Normal Effective Temperature
- Universal Thermal Climate Index
- Wet Bulb Globe Temperature
- Wet Bulb Globe Temperature (Liljegren method)
- Heat Force (KNMI 0–10 heat-stress scale)
- Excess Heat Factor and Excess Cold Factor
- Wind Chill

In support of the above, it also calculates:

- Globe Temperature
- Mean Radiant Temperature
- Mean Radiant Temperature from Globe Temperature
- Relative Humidity Percentage
- Saturation Vapour Pressure
- Wet Bulb Temperature

## Calling convention

Every function is **vectorised over NumPy** and works elementwise on arrays of
any shape. Inputs and outputs are **SI** (temperatures in Kelvin, wind in m/s at
10 m, relative humidity in %, vapour pressure in hPa); the excess heat/cold
factors are the documented unit-agnostic exception.

Pass **NumPy arrays**, and **wrap a scalar in an array**, e.g.:

```python
import numpy as np
import thermofeel as tf

utci = tf.calculate_utci(
    t2_k=np.array([300.0]), va=np.array([3.0]), mrt=np.array([310.0]),
    td_k=np.array([290.0]),
)
```

Most functions also accept a bare Python scalar, but those that branch
internally (e.g. `calculate_heat_index_simplified`, `approximate_dsrp`,
`calculate_saturation_vapour_pressure_multiphase`) require array input — so
wrapping scalars in an array is the reliable convention. Inputs are not clamped
to each index's validity range (see the per-function docstrings); masking
out-of-range points is the caller's responsibility.

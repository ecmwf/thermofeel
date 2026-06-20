# Wind Chill

The Wind Chill is defined as the air temperature of an equivalent environment
that, under calm wind conditions, would entail the same skin surface heat loss to
the environment as in the actual, windy, environment. It takes into account the
assumptions of convective and radiative heat loss described in modern heat
transfer theory, and assumes no impact from the sun.

More information: Blazejczyk, K., Epstein, Y., Jendritzky, G. et al. Comparison
of UTCI to selected thermal indices. Int J Biometeorol 56, 515–535 (2012).
<https://doi.org/10.1007/s00484-011-0453-2>

## How to use

You will need 2 m air temperature in Kelvin and 10 m wind speed in metres per
second.

It returns the wind chill in Kelvin.

```python
calculate_wind_chill(2m_temperature, 10m_wind_speed)
```

## Interpret the output

The wind chill is described in terms of the risk incurred by human skin, based on
the rate of heat loss caused by exposure to wind and low temperatures.

| Wind Chill (°C) | Risk |
| --- | --- |
| 0 to −10 | Low risk |
| −10 to −28 | Moderate risk |
| −28 to −40 | High risk |
| −40 to −48 | Very high risk |
| −48 to −55 | Severe risk |
| ≤ −55 | Extreme risk |

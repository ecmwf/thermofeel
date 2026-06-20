# Universal Thermal Climate Index

The Universal Thermal Climate Index (UTCI) is a measure of the thermal stress on
the human body from outdoor conditions, defined as the equivalent air
temperature of a reference environment that would cause the same physiological
response as the actual conditions.

More information: Jendritzky, G., de Dear, R. & Havenith, G. UTCI—Why another
thermal index?. Int J Biometeorol 56, 421–428 (2012).
<https://doi.org/10.1007/s00484-011-0513-7>

## How to use

You need 2 m air temperature and mean radiant temperature in Kelvin, 2 m dew
point temperature in Kelvin or water vapour pressure in hPa, and 10 m wind speed
in m/s. Please use with NumPy arrays.

The UTCI is returned in Kelvin.

```python
rh_pc = calculate_relative_humidity_percent(2m_temperature, 2m_dew_point_temperature)
ehPa = calculate_saturation_vapour_pressure(2m_temperature) * rh_pc / 100.0
utci = calculate_utci(2m_temperature, 10m_wind_speed, mean_radiant_temperature, td_k=None, ehPa=ehPa)
```

## Interpret the output

The UTCI is classified into a 10-category scale. Each category, defined by a
specific range of UTCI values, corresponds to a well-defined set of human
physiological responses to the outdoor environment.

| Universal Thermal Climate Index (°C) | Thermal Stress |
| --- | --- |
| > 46 | Extreme heat stress |
| 38 to 46 | Very strong heat stress |
| 32 to 38 | Strong heat stress |
| 26 to 32 | Moderate heat stress |
| 9 to 26 | No thermal stress |
| 0 to 9 | Slight cold stress |
| −13 to 0 | Moderate cold stress |
| −27 to −13 | Strong cold stress |
| −40 to −27 | Very strong cold stress |
| < −40 | Extreme cold stress |

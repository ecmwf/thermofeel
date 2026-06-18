# Heat Index

The Heat Index (HI) is defined as the temperature the human body perceives in
shady conditions when perspiration is limited due to increased relative
humidity.

*thermofeel* computes the HI with two methods: Heat Index Simplified and Heat
Index Adjusted.

More information:

- Blazejczyk, K., Epstein, Y., Jendritzky, G. et al. Comparison of UTCI to selected thermal indices. Int J Biometeorol 56, 515–535 (2012). <https://doi.org/10.1007/s00484-011-0453-2>
- NOAA/National Weather Service, Weather Prediction Center. The Heat Index Equation. <https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml>

## How to use

### Simplified Heat Index

You need 2 m air temperature in Kelvin and relative humidity as a percentage.
Please use NumPy arrays.

It returns the HI in Kelvin.

```python
calculate_heat_index_simplified(2m_temperature, relative_humidity_percent)
```

### Adjusted Heat Index

You need 2 m air temperature and 2 m dew point temperature in Kelvin.
Please use NumPy arrays.

It returns the HI in Kelvin.

```python
calculate_heat_index_adjusted(2m_temperature, 2m_dew_point_temperature)
```

## Interpret the output

The HI is described in terms of danger of heat-related illnesses (NOAA/NWS
classification):

| Heat Index (°C) | Risk |
| --- | --- |
| 27 – 32 | Caution |
| 32 – 39 | Extreme caution |
| 39 – 51 | Danger |
| ≥ 51 | Extreme danger |

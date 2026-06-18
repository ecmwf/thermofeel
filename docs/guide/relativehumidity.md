# Relative Humidity and Saturation Vapour Pressure

There are two methods to compute relative humidity in *thermofeel* and these are
used together in some index calculations, such as in the UTCI.

## How to use

### Saturation Vapour Pressure over Water

The input is 2 m air temperature in Kelvin. The output is in hPa.

```python
calculate_saturation_vapour_pressure(2m_temperature)
```

### Relative Humidity Percent

The inputs are 2 m temperature and dew point temperature in Kelvin. The output is
a percentage.

```python
calculate_relative_humidity_percent(2m_temperature, 2m_dew_point_temperature)
```

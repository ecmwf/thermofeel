# Humidex

The Humidex is defined as the temperature (in °C) the human body perceives in
hot, humid weather.

More information: Blazejczyk, K., Epstein, Y., Jendritzky, G. et al. Comparison
of UTCI to selected thermal indices. Int J Biometeorol 56, 515–535 (2012).
<https://doi.org/10.1007/s00484-011-0453-2>

## How to use

You need 2 m air temperature and 2 m dew point temperature in Kelvin.

It returns the Humidex in Kelvin.

```python
calculate_humidex(2m_temperature, 2m_dew_point_temperature)
```

## Interpret the output

The Humidex is described in terms of comfort:

| Humidex (°C) | Discomfort levels |
| --- | --- |
| 20 – 30 | Little discomfort |
| 30 – 40 | Some discomfort |
| 40 – 46 | Great discomfort, avoid exertion |
| ≥ 54 | Dangerous, possible heat stroke |

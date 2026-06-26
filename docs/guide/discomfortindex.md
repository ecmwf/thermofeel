# Discomfort Index

The Discomfort Index (DI), also called Thom's Discomfort Index or the
Temperature-Humidity Index, estimates human heat discomfort from air temperature
and relative humidity alone. It was introduced by Thom (1959) — originally from
dry-bulb and wet-bulb temperatures — and is implemented here in the Celsius and
relative-humidity form of Giles et al. (1990). It remains widely used as a simple
heat-stress indicator.

More information:

Thom, E.C. (1959). The Discomfort Index. Weatherwise, 12(2), 57-61.
<https://doi.org/10.1080/00431672.1959.9926960> (origin of the index)

Giles, B.D., Balafoutis, C. & Maheras, P. (1990). Too hot for comfort: the
heatwaves in Greece in 1987 and 1988. International Journal of Biometeorology,
34(2), 98-104. <https://doi.org/10.1007/BF01093455> (the temperature /
relative-humidity formulation used here)

Epstein, Y. & Moran, D.S. (2006). Thermal comfort and the heat stress indices.
Industrial Health, 44(3), 388-398. <https://doi.org/10.2486/indhealth.44.388>
(review and discomfort categories)

## How to use

You need 2 m air temperature in Kelvin and relative humidity in percent. It
returns the Discomfort Index in Kelvin.

```python
calculate_discomfort_index(2m_temperature, relative_humidity)
```

Relative humidity can be obtained from temperature and dew point with
`calculate_relative_humidity_percent`.

## Interpret the output

The index uses the relative-humidity form
DI = T - 0.55 (1 - 0.01 RH)(T - 14.5), with T in °C and RH in %. At 100% relative
humidity the index equals the air temperature; in drier air it falls below it,
reflecting the reduced heat stress of low humidity. This is Thom's index on the
Celsius scale; the same index also appears in a Fahrenheit relative-humidity form
(with 58 °F in place of 14.5 °C) and in Thom's original dry-bulb/wet-bulb form.

The commonly used thermal-sensation categories (Giles et al. 1990; Epstein &
Moran 2006) are:

| DI (°C) | Thermal sensation |
| --- | --- |
| < 21 | No discomfort |
| 21 - 24 | Less than 50% of the population feels discomfort |
| 24 - 27 | More than 50% of the population feels discomfort |
| 27 - 29 | Most of the population feels discomfort |
| 29 - 32 | Everyone feels severe stress |
| ≥ 32 | State of medical emergency |

The output is not clamped to this range: out-of-range inputs return the raw
index and the caller is responsible for any masking. Convert the returned Kelvin
value back to °C (e.g. with `kelvin_to_celsius`) to look it up against the table.

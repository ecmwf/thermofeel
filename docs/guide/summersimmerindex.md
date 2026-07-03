# Summer Simmer Index

The Summer Simmer Index (SSI) estimates warm-season human heat discomfort from
air temperature and relative humidity alone. It is computed internally in
Fahrenheit and returned in Kelvin. The inner term of the formula is Thom's
Temperature-Humidity Index in Fahrenheit, so the SSI is an affine transform of
that index — the Fahrenheit sibling of the
[Discomfort Index](discomfortindex.md).

More information:

Pepi, J.W. (1987). The Summer Simmer Index. Weatherwise, 40(3), 143-145.
<https://doi.org/10.1080/00431672.1987.9933356>

!!! warning "Provenance caveat"

    The 1987 article is not openly available, so the equation implemented here is
    reproduced from **secondary sources**. It is the **common 1987 closed form**

        SSI = 1.98 (Tf - (0.55 - 0.0055 RH)(Tf - 58)) - 56.83

    with `Tf` the air temperature in °F and `RH` in % — an affine transform of
    Thom's Fahrenheit Temperature-Humidity Index. It is **not** the author's
    later, tabulated *New Summer Simmer Index*, which is a different (steeper)
    relationship. The interpretation bands below are likewise the commonly
    reproduced comfort scale for this common form, not the New SSI categories.

## How to use

You need 2 m air temperature in Kelvin and relative humidity in percent. It
returns the Summer Simmer Index in Kelvin.

```python
calculate_summer_simmer_index(2m_temperature, relative_humidity_percent)
```

Relative humidity can be obtained from temperature and dew point with
`calculate_relative_humidity_percent`.

## Interpret the output

The SSI is conventionally read on the Fahrenheit scale on which it is defined.
The commonly reproduced comfort bands for the common form are:

| SSI (°F) | Perceived comfort |
| --- | --- |
| < 70 | Cool |
| 70 – 77 | Comfortable |
| 77 – 83 | Slightly uncomfortable |
| 83 – 91 | Increasing discomfort; caution |
| 91 – 100 | Significant discomfort; heat fatigue with prolonged exposure |
| 100 – 112 | Dangerous; heat exhaustion and heat cramps likely |
| ≥ 112 | Extreme danger; heatstroke imminent |

These thresholds come from secondary sources and describe the common form only;
they are **not** the tabulated New Summer Simmer Index bands.

The output is not clamped to this range: out-of-range inputs return the raw index
and the caller is responsible for any masking. Convert the returned Kelvin value
to °F (e.g. with `kelvin_to_fahrenheit`) to look it up against the table.

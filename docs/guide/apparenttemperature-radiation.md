# Apparent Temperature (radiation)

The radiation-inclusive Apparent Temperature (AT) extends Steadman's apparent
temperature with a term for the net radiation absorbed by the body, so that
solar and thermal loading raise the perceived temperature above the shade value.

It is the operational form published by the Australian Bureau of Meteorology:

`AT = Ta + 0.348·e − 0.70·va + 0.70·q/(va + 10) − 4.25`  (Ta in °C)

where `e` is the ambient water-vapour pressure in hPa, `va` is the 10 m wind
speed and `q` is the net radiation absorbed per unit area of body surface.

More information:

- Steadman, R. G. Norms of apparent temperature in Australia. Aust. Met. Mag. 43, 1–16 (1994). <https://doi.org/10.1071/es94001>
- Australian Government Bureau of Meteorology. Thermal Comfort observations — About the formula for the apparent temperature. <https://www.bom.gov.au/info/thermal_stress/#atapproximation>

## How to use

You need 2 m air temperature in Kelvin, 10 m wind speed in metres per second,
relative humidity as a percentage and the net radiation absorbed per unit
body-surface area `q` in W m⁻².

It returns the apparent temperature in Kelvin.

```python
calculate_apparent_temperature_radiation(2m_temperature, 10m_wind_speed, relative_humidity_percent, q)
```

## Interpret the output

Relative to the non-radiation `calculate_apparent_temperature`, this form adds
the radiation term `0.70·q/(va + 10)`, so a positive radiation load makes the
apparent temperature warmer.

`q` is **caller-supplied**: it is the net radiation absorbed per unit area of
body surface (W m⁻²), analogous to how `cossza` is supplied to other functions.
It is **not** an NWP surface radiation flux and **not** the mean radiant
temperature — you must derive and supply `q` yourself for your application.

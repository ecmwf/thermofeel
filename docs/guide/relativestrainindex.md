# Relative Strain Index

The Relative Strain Index (RSI) is a dimensionless warm-season heat-stress index.
It expresses the strain a warm, humid environment places on a young, healthy
adult from just two variables: air temperature and the ambient water-vapour
pressure of the air.

Index origin: Lee, D.H.K. & Henschel, A. (1966) *Effects of physiological and
clinical factors on response to heat*, Annals of the New York Academy of
Sciences 134(2):743–749. <https://doi.org/10.1111/j.1749-6632.1966.tb43059.x>

Implemented closed form and assessment bands: Asghari, M. et al. (2020)
*Feasibility of Relative Strain Index (RSI) for the Assessment of Heat Stress in
Outdoor Environments: Case Study in Three Different Climates of Iran*, The Open
Ecology Journal 13:11–18. <https://doi.org/10.2174/1874213002013010011>

## How to use

You need 2 m air temperature in Kelvin and relative humidity in percent. The
ambient water-vapour pressure `e` (in hPa) is derived internally from these two
inputs with `calculate_nonsaturation_vapour_pressure`, so no separate vapour
input is required.

It returns the Relative Strain Index, which is dimensionless.

```python
calculate_relative_strain_index(2m_temperature, relative_humidity_percent)
```

The implemented form is `RSI = (Ta − 21) / (58 − e)`, with `Ta` the air
temperature in °C and `e` the ambient water-vapour pressure in hPa. RSI is a
summer index meant for air temperatures up to about 35 °C, and it rises with
both temperature and humidity (through the vapour-pressure term `e`). At lower
air temperatures it stays near the comfort level across most of the humidity
range, though very humid conditions can still lift it past the comfort threshold.

## Interpret the output

The index is read against five strain levels (Asghari et al. 2020, Table 2):

| Relative Strain Index | Interpretation |
| --- | --- |
| < 0.15 | Climate comfort |
| 0.15 – 0.25 | Discomfort for sensitive people (the elderly and children) |
| 0.25 – 0.35 | Discomfort for all people |
| 0.35 – 0.45 | Risk of exposure to excessive heat (heat stroke) for 50% of people or more |
| ≥ 0.45 | Risk of hyperthermia for all people |

## Provenance and variants

The index originates with Lee & Henschel (1966); the specific hectopascal closed
form implemented here, together with the five-level assessment table above, is
taken from Asghari et al. (2020), which states the formula with units and was
read directly from the open-access source. The same thresholds
(0.15 / 0.25 / 0.35 / 0.45) recur in the Romanian bioclimatology literature
(Ciulache 2006; Ionac & Ciulache 2008). Some secondary sources attribute the
level table to "Błażejczyk 2011", but no Błażejczyk 2011 publication tabulating
RSI could be found under that citation — the 2011 *Miscellanea Geographica* paper
cited defines the unrelated Bioclimatic Contrast Index — so this guide cites the
verified Asghari (2020) tabulation instead.

A second, different closed form appears in the literature,
`(10.7 + 0.74·(Ta − 35)) / (44 − Pa)`, where `Pa` is expressed in another unit
(most likely mmHg). It could not be verified against the primary text and is
**not** implemented here; only the peer-reviewed hPa form is provided.

Note the domain edge: as `e` approaches 58 hPa (near-saturation around 35–36 °C)
the denominator `58 − e` shrinks, so the index grows without bound (±infinity),
and for `e` above 58 hPa (very hot and near-saturated, beyond the ~35 °C validity
range) the denominator turns negative and RSI becomes negative — a spurious value
for a heat-strain index. These out-of-range values are returned as-is and are
deliberately not clamped, so mask or guard the result if your inputs can reach
that regime.

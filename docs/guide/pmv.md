# Predicted Mean Vote (PMV) and Predicted Percentage of Dissatisfied (PPD)

The Predicted Mean Vote (PMV) is Fanger's steady-state thermal-comfort index: it
predicts the mean thermal-sensation vote of a large group of people on the
7-point ASHRAE scale, from a full human heat balance. The Predicted Percentage
of Dissatisfied (PPD) maps a PMV value onto the percentage of occupants likely
to be thermally dissatisfied.

More information: ISO 7730:2005, *Ergonomics of the thermal environment —
Analytical determination and interpretation of thermal comfort using calculation
of the PMV and PPD indices and local thermal comfort criteria* (Annex D).
Originating work: Fanger, P.O. (1970) *Thermal Comfort: Analysis and
Applications in Environmental Engineering*, McGraw-Hill.

## How to use

You need the 2 m air temperature and the mean radiant temperature in Kelvin, the
relative air velocity at the body in m/s, and one humidity input — either the
relative humidity in percent (`rh`) or the water-vapour pressure in hPa
(`vapour_pressure_hpa`). You may also set the metabolic rate (`met`), the
clothing insulation (`clo`), and any external work (`wme`).

```python
pmv = calculate_pmv(
    2m_temperature,          # K
    mean_radiant_temperature,  # K
    relative_air_velocity,   # m/s, at the body (see note below)
    rh=relative_humidity,    # %  (or vapour_pressure_hpa=...)
    met=1.2,                 # met (ISO 8996 sedentary/office activity)
    clo=0.5,                 # clo (light indoor clothing)
    wme=0.0,                 # met (external work)
)

calculate_ppd(pmv)
```

`calculate_pmv` returns the dimensionless PMV; `calculate_ppd` returns the PPD in
percent. The clothing-surface temperature is solved by the ISO 7730 Annex D
fixed-point iteration over the whole array at once; any element that does not
converge within the iteration cap is returned as `NaN`.

### Important: `var` is body-level velocity, not the 10 m wind

`var` is the **relative air velocity at the body** — the still-air speed the
occupant experiences plus any air movement induced by their own activity. It is
**not** the 10 m meteorological wind speed. Do not pass a model wind field
directly; PMV is an indoor/occupant comfort index.

### `met` and `clo` are parameters, with cited default conventions

PMV is only defined for a stated activity and clothing level, so `met` and `clo`
are explicit inputs. The defaults are conventions rather than constants of the
standard: `met=1.2` corresponds to sedentary/office activity (ISO 8996), and
`clo=0.5` to light indoor clothing. Set them to match your scenario.

## Interpret the output

PMV is reported on the ASHRAE 7-point thermal-sensation scale:

| PMV | Thermal sensation |
| --- | --- |
| +3 | Hot |
| +2 | Warm |
| +1 | Slightly warm |
| 0 | Neutral |
| −1 | Slightly cool |
| −2 | Cool |
| −3 | Cold |

PPD is smallest (5%) at thermal neutrality (PMV = 0) and rises symmetrically as
PMV departs from zero in either direction.

ISO 7730 defines three categories of thermal environment from the PMV/PPD pair:

| Category | Predicted Percentage of Dissatisfied (PPD) | Predicted Mean Vote (PMV) |
| --- | --- | --- |
| A | < 6 % | −0.2 < PMV < +0.2 |
| B | < 10 % | −0.5 < PMV < +0.5 |
| C | < 15 % | −0.7 < PMV < +0.7 |

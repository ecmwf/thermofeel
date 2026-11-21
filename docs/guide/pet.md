# Physiological Equivalent Temperature (PET)

The Physiological Equivalent Temperature is the outdoor companion to UTCI and
one of the most widely used indices in human biometeorology. It is defined as
**the air temperature of a standard indoor reference environment at which the
human heat balance closes with the same core and skin temperature as in the
outdoor environment being assessed**, solved from the MEMI model (Munich
Energy-balance Model for Individuals): a steady-state core/skin two-node model
plus a clothing node.

Because activity and clothing are fixed by the definition, PET depends only on
the meteorological inputs. That is what makes it a *climatic* index rather than
a behavioural one, and what lets it be read against a layperson's indoor
experience.

More information: Höppe, P. (1999) *The physiological equivalent temperature — a
universal index for the biometeorological assessment of the thermal
environment*, International Journal of Biometeorology 43:71–75,
[doi:10.1007/s004840050118](https://doi.org/10.1007/s004840050118). For the
full equation set, the corrections applied here and a critique of the model,
see Walther, E. & Goestchel, Q. (2018) *The P.E.T. comfort index: Questioning
the model*, Building and Environment 137:1–10,
[doi:10.1016/j.buildenv.2018.03.054](https://doi.org/10.1016/j.buildenv.2018.03.054).

## How to use

You need the 2 m air temperature and the mean radiant temperature in Kelvin, the
air velocity **at the body** in m/s, and one humidity input — either the
relative humidity in percent (`rh`) or the water-vapour pressure in hPa
(`vapour_pressure_hpa`).

```python
pet = calculate_pet(
    2m_temperature,            # K
    mean_radiant_temperature,  # K
    air_velocity_at_body,      # m/s (see note below)
    rh=relative_humidity,      # %  (or vapour_pressure_hpa=...)
)
```

The result is the PET in Kelvin. Everything else has a default that reproduces
Höppe's reference subject and the PET reference activity and clothing:

| Parameter | Default | Meaning |
|---|---|---|
| `m_act_w` | `80.0` | activity metabolism [W], **whole body** |
| `clo` | `0.9` | clothing insulation [clo] |
| `height`, `mass`, `age`, `sex` | 1.80 m, 75 kg, 35, `"male"` | reference subject |
| `p_atm_hpa` | `1013.25` | atmospheric pressure [hPa] |
| `f_eff` | `0.725` | effective radiative area fraction |

### The reference environment, and the identity it implies

PET is defined against this indoor reference environment:

| Quantity | Value |
|---|---|
| mean radiant temperature | equal to the air temperature |
| air velocity | 0.1 m/s |
| water vapour pressure | 12 hPa |
| clothing | 0.9 clo |
| activity metabolism | 80 W (whole body) |

Feeding that environment to the function therefore returns the air temperature
itself:

```python
calculate_pet(t, t, 0.1, vapour_pressure_hpa=12.0)   # == t
```

The test suite pins this identity, because it simultaneously fixes the
reference environment, the clothing convention and the activity-unit
convention.

### Important: `va` is body-level velocity, not the 10 m wind

`va` is the air velocity **at the body**, not the 10 m meteorological wind
speed. Do not pass a model wind field directly — scale it first with
[`scale_windspeed`](relativehumidity.md), as for PMV.

### `m_act_w` is whole-body watts, not met units

`m_act_w` is the activity metabolism in **watts for the whole body**, matching
the source model, whose reference value of 80 W defines PET. This is *not* the
ISO met unit used by [`calculate_pmv`](pmv.md) (1 met = 58.15 W m⁻²), and it is
also not the convention used by `pythermalcomfort.pet_steady`, whose `met`
argument is effectively activity watts ÷ 58.2.

## Interpretation

PET is a **comparator, not an absolute measure of thermal strain**. Höppe is
explicit that a person at PET = 20 °C is cold in swimming trunks and sweating in
a coat. Read it against the published thermal-stress classes below (Matzarakis &
Mayer 1996; Matzarakis, Mayer & Iziomon 1999,
[doi:10.1007/s004840050119](https://doi.org/10.1007/s004840050119)) — but note
the caveat that follows.

| PET (°C) | Thermal perception | Physiological stress |
|---|---|---|
| < 4 | very cold | extreme cold stress |
| 4–8 | cold | strong cold stress |
| 8–13 | cool | moderate cold stress |
| 13–18 | slightly cool | slight cold stress |
| 18–23 | comfortable | no thermal stress |
| 23–29 | slightly warm | slight heat stress |
| 29–35 | warm | moderate heat stress |
| 35–41 | hot | strong heat stress |
| > 41 | very hot | extreme heat stress |

!!! warning "The stress classes carry a model-variant uncertainty"
    These classes were derived with the **original Höppe/VDI** skin-diffusion
    model. thermofeel implements the **Woodcock** clothing-aware evaporative
    resistance, which is the state of the art and what the other maintained
    implementations use, but which moves PET by −7 to +2.6 K relative to the
    original. Treat the class boundaries as carrying that uncertainty, and state
    which variant produced any operational threshold.

## Model variant and known divergences

The published sources disagree in several places. thermofeel follows the
VDI/Walther reference implementations, because those are what the published PET
tables and stress classes were generated with. Each choice is recorded in the
function docstring; the two that matter most:

- **The clothing temperature is frozen** at its actual-environment value while
  the reference air temperature is solved. Both papers' prose says it is
  adjusted alongside the air temperature, but every released implementation
  freezes it — and re-solving it moves PET by up to 24 K in hot, high
  radiant-load conditions, precisely the regime that matters for heat-warning
  work.
- **Evaporation uses the Woodcock resistance** (see the warning above). The
  original Höppe constant skin-diffusion resistance is deliberately *not*
  implemented, because it could not be verified against a primary source.

Accuracy against Höppe's own published Table 1 is ~1 K rms. The two source
papers already disagree with each other by up to 1.3 K on those same cases, so
do not report PET to a precision the model does not support.

## Performance

The published implementations hand the 3×3 non-linear system to a per-point
`scipy.optimize.fsolve`, which does not vectorise. thermofeel exploits the fact
that the system is triangularisable — the core-node balance contains no
clothing temperature, and its forcing term contains no body temperature at all —
which collapses the problem to two nested one-dimensional *monotone* root finds,
each solved by bisection over the whole array at once.

The result is numpy-only, with no SciPy dependency, and is roughly two orders of
magnitude faster than a per-point solver, which is what makes PET tractable over
global forecast grids and ensembles.

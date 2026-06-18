# Architecture

thermofeel is a small, pure-Python library of human thermal comfort indices
built on NumPy. This document explains how the pieces fit together.

For the design rationale and the SI-unit contract, see [DESIGN.md](DESIGN.md).
For why the library exists, see [MOTIVATION.md](MOTIVATION.md). For release
history, see [../ChangeLog.rst](../ChangeLog.rst).

## High-Level Structure

```
                ┌─────────────────────────────────────────┐
                │             User / NWP pipeline           │
                │   (scalars or NumPy lat/lon/step grids)   │
                └─────────────────────┬─────────────────────┘
                                      │  import thermofeel
                        ┌─────────────▼─────────────┐
                        │      thermofeel/           │
                        │  __init__.py  (public API) │
                        └─────────────┬─────────────┘
                                      │  from .thermofeel import *
              ┌───────────────────────▼───────────────────────┐
              │              thermofeel/thermofeel.py           │
              │                                                 │
              │  comfort indices   supporting quantities        │
              │  ───────────────   ─────────────────────        │
              │  utci              relative_humidity_percent     │
              │  apparent_temp     saturation_vapour_pressure    │
              │  heat_index_*      nonsaturation_vapour_pressure  │
              │  humidex           mean_radiant_temperature       │
              │  normal_eff_temp   bgt / mrt_from_bgt             │
              │  wbgt / wbgt_simple wbt                           │
              │  wind_chill        dew_point / scale_windspeed    │
              │                    approximate_dsrp               │
              └───────────────────────┬───────────────────────┘
                                      │  uses
                        ┌─────────────▼─────────────┐
                        │     thermofeel/helpers.py  │
                        │  SI unit converters         │
                        │  (C/K/F interconversion)    │
                        └─────────────┬─────────────┘
                                      │
                              ┌───────▼───────┐
                              │     numpy      │
                              └───────────────┘

   thermofeel/experimental_wbgt.py — iterative Liljegren WBGT (NOT exported,
   not part of the public API; depends on the public functions above)
```

## Package Layout

| Path | Role |
|------|------|
| `thermofeel/__init__.py` | Public surface: `from .thermofeel import *` and the canonical `__version__` (single source of truth; `pyproject.toml` reads it dynamically) |
| `thermofeel/thermofeel.py` | All index and supporting-quantity functions. The substance of the library |
| `thermofeel/helpers.py` | Pure SI unit converters used by the index functions |
| `thermofeel/experimental_wbgt.py` | Iterative Liljegren WBGT/globe/wet-bulb solver — experimental, unexported, not covered by the public test suite |
| `tests/` | pytest suites + stored CSV reference outputs |
| `examples/` | Runnable usage examples (incl. ECMWF GRIB/eccodes pipeline examples) |
| `docs/` | Sphinx documentation (Read the Docs), one guide page per index |

## Data Flow Between Functions

Indices are composed from supporting quantities. The important dependency
chains (each arrow = "is computed from"):

```
relative_humidity_percent ← t2_k, td_k
saturation_vapour_pressure ← t2_k
nonsaturation_vapour_pressure ← t2_k, rh

UTCI            ← t2_k, va, mrt, (td_k or ehPa)
                  └ water vapour pressure ← saturation_vapour_pressure × rh

WBGT            ← bgt(t2_k, mrt, va) + wbt(t2_k, rh) + t2_k
WBGT (simple)   ← t2_k + nonsaturation_vapour_pressure
BGT             ← t2_k, mrt, scale_windspeed(va, 1.1)
MRT             ← ssrd, ssr, dsrp, strd, fdir, strr, cossza
MRT from BGT    ← t2_k, bgt_k, scale_windspeed(va, 1.1)

Apparent Temp   ← t2_k, va, nonsaturation_vapour_pressure
Normal Eff Temp ← t2_k, scale_windspeed(va, 1.2), rh
Wind Chill      ← t2_k, va (converted to km/h internally)
Heat Index      ← t2_k, rh        (simplified)
Heat Index adj. ← t2_k, td_k → rh (NOAA regression, in °F)
Humidex         ← t2_k, td_k
```

The two cross-cutting helpers that almost everything routes through:
- **`scale_windspeed(va, h)`** — converts the public 10 m wind to the height a
  given formula assumes (1.1 m globe, 1.2 m NET).
- **the unit converters in `helpers.py`** — every formula that is published in
  °C/°F/km/h converts at entry and exit so the public boundary stays SI.

## Solar Geometry Boundary

`calculate_mean_radiant_temperature` needs the cosine of the solar zenith angle
(`cossza`) and direct radiation. These are **inputs**, supplied by the caller
(e.g. `earthkit-meteo.solar.calculate_cos_solar_zenith_angle`). `approximate_dsrp`
is a convenience to derive direct radiation from `fdir` and `cossza` when the
dataset lacks it — explicitly documented as approximate near low sun angles.

## Distribution

- **PyPI:** `pip install thermofeel` (pure-Python wheel/sdist).
- **CI:** ECMWF `downstream-ci` (`.github/workflows/ci.yml`) runs the test suite
  and the `python_qa` lint gate; `cd.yml` publishes to PyPI on a bare
  `MAJOR.MINOR.MICRO` git tag.
- **Docs:** Read the Docs builds `docs/` (Sphinx).

## Testing Topology

Two complementary suites (see [TEST.md](TEST.md)):
- `tests/test_thermofeel.py` — array-level regression against stored
  `tests/*.csv` reference outputs driven by `thermofeel_testcases.csv`.
- `tests/test_scalars.py` — pointwise checks against individually verified
  reference values.

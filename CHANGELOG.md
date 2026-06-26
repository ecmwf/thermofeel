# ChangeLog

## 2.3.0

- Added `calculate_discomfort_index`: Thom's Discomfort Index (Temperature-
  Humidity Index) from 2 m temperature and relative humidity, returning Kelvin.
  Implements the Celsius/relative-humidity formulation of Thom's index given by
  Giles et al. (1990, https://doi.org/10.1007/BF01093455),
  DI = T - 0.55 (1 - 0.01 RH)(T - 14.5); the index originates with Thom (1959,
  https://doi.org/10.1080/00431672.1959.9926960). Validated against analytic
  reference values (including the RH=100% -> DI=Ta and T=14.5 °C -> DI=14.5 °C
  identities) and pinned by a regression CSV. Adds a guide page. See also the
  review Epstein & Moran (2006, https://doi.org/10.2486/indhealth.44.388), which
  uses a different wet-bulb formulation.

## 2.2.0

- Added inline type hints across the public API and a `py.typed` marker
  (PEP 561), so downstream type-checkers pick up thermofeel's signatures
- Exported `fahrenheit_to_celsius` at the top level, completing the C/K/F
  converter set
- Fixed silent integer-dtype truncation in
  `calculate_saturation_vapour_pressure_multiphase` and `approximate_dsrp`: an
  integer-typed input array no longer truncates the float result
- Added the missing citation for `calculate_relative_humidity_percent`
  (Magnus-Tetens; Murray 1967) and documented validity ranges in the UTCI, wet
  bulb, simple WBGT, NET, apparent-temperature and humidex docstrings; documented
  the array calling convention and that relative humidity is not clamped
- `calculate_bgt` (and therefore the Stull `calculate_wbgt`) now returns the
  mean radiant temperature at exactly zero wind — the analytic calm-air limit
  (no convection ⇒ globe at radiative equilibrium) — instead of `NaN`. Values for
  any positive wind are unchanged
- Ran a first numerical-robustness pass (`ROBUSTNESS.md`): audited the
  `log`/`sqrt`/`power`/division domains and added `NaN`-propagation and
  error-contract tests
- The published wheel now ships only the `thermofeel` package (`examples/`,
  `scripts/` and `tests/` are no longer packaged)
- Raised the minimum supported Python to 3.10 (`requires-python = ">=3.10"`),
  refreshed the trove classifiers and added Python 3.14
- Repaired the `examples/compute-thermal-indices.py` example to the 2.x API:
  solar geometry now obtained from `earthkit-meteo`
  (`solar.cos_solar_zenith_angle_integrated`), `calculate_utci`/`calculate_wbgt`/
  `calculate_bgt`/`calculate_normal_effective_temperature` return Kelvin directly,
  renamed/removed functions and the obsolete timing decorators dropped
- Removed dead, overwritten code in `scale_windspeed` (kept the computed log-law
  coefficient) and corrected the `approximate_dsrp` threshold comment
- Added the `thermofeel.excess_heat` submodule: Excess Heat Factor (EXHF),
  Excess Cold Factor (EXCF), heatwave severity, and the supporting daily mean
  temperature, significance and acclimatisation indices, after Nairn & Fawcett
  (2014). The factors are also exposed at the top level as
  `calculate_excess_heat_factor` / `calculate_excess_cold_factor`. The functions
  implement the per-day formulas; the temporal aggregations are left to upstream
  tooling (e.g. earthkit-transforms)
- Added `calculate_wbgt_liljegren`: physically based Wet Bulb Globe Temperature
  using the Liljegren et al. (2008) method, validated against Liljegren's
  reference implementation (github.com/mdljts/wbgt)
- Added `calculate_wind_speed_2m_liljegren`: the KNMI/Liljegren
  stability-dependent 10 m-to-2 m wind-speed profile, used by default in
  `calculate_wbgt_liljegren` (selectable via its `wind_scaling` argument)
- Added `calculate_heat_force`: the KNMI 0–10 heat-force (hittekracht) scale
  derived from WBGT (KNMI Technical Report TR-26-04)
- Removed the broken, unexported `experimental_wbgt` module, now superseded by
  `calculate_wbgt_liljegren`
- Fixed the `calculate_bgt` argument order in the scalar test suite

## 2.1.7

- Change project to single branch main

## 2.1.6

- Updated humidex unit tests

## 2.1.5

- Fix descriptions
- Fix variable names
- Updated humidex calculation and reference

## 2.1.4

- Fix setup of project
- Removed unused dependency on earthkit-meteo

## 2.1.3

- Fix documentation

## 2.1.2

- Fix calculation of adjusted Heat Index (function `calculate_heat_index_adjusted`)

## 2.1.1

- Small bug fix in the UTCI computation

## 2.1.0

- Fix computation of normal effective temperature
- Fix documentation links to new repository URL

## 2.0.0

**Standardisation**

- I/O variables converted into International System of Units (SI) (e.g., K)
- Docstrings standardised in their descriptions (e.g., inputs as float arrays) and references
- Variable names standardised across function declarations (i.e., t2m vs t2k vs t_k → now t2_k)
- Use of temperature converter functions
- Input variables are limited to those explicitly needed in the function formula (see e.g., `calculate_normal_effective_temperature` function)

**New Functions**

- `scale_windspeed` function to scale 10 m wind speed to height h with h < 10 m
- `calculate_nonsaturation_vapour_pressure` to calculate vapour pressure at a given relative humidity and air temperature
- `calculate_wbgt_simple` function renamed
- `calculate_normal_effective_temperature` function renamed

**Improvements**

- thermofeel library docstring lists computed variables in alphabetical order
- the cosine of the solar zenith angle now computed via the [earthkit-meteo](https://github.com/ecmwf/earthkit-meteo) library and must be provided as input
- `calculate_bgt` function is calculated via a 4× faster formula
- `calculate_saturation_vapour_pressure_multiphase` formulas replaced with those used in the IFS
- changeable threshold in `approximate_dsrp` function (set to 0.1 by default)
- invalidity outside input variables range specified in `calculate_wind_chill` function docstring

**Bug Fixes**

- fixed `approximate_dsrp` to avoid fdir being overwritten with dsrp when calculating MRT
- fixed `calculate_wbgt_simple` constant value and vapour pressure calculated from non-saturated formula
- fixed `calculate_bgt` wind speed at 1.1 m in input
- fixed `calculate_mrt_from_bgt` wind speed at 1.1 m in input
- fixed `calculate_normal_effective_temperature` wind speed at 1.2 m in input
- fixed `calculate_apparent_temperature` wind speed at 10 m in input; vapour pressure calculated from non-saturated formula
- fixed `calculate_wind_chill` wind speed in km/h and operation symbols in main formula
- fixed `calculate_heat_index_simplified` a wrong sign and missing constant value in hiarray; hi set to 2 m air temperature when the latter is below 20 °C
- fixed `fahrenheit_to_kelvin` converter function

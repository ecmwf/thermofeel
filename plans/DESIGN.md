# Design: thermofeel — Human Thermal Comfort Index Library

Repo: ecmwf/thermofeel

> For **why** thermofeel exists and what problem it solves, see `MOTIVATION.md`.
> For the **package/module layout**, see `ARCHITECTURE.md`.
> For **release history**, see `../CHANGELOG.md`.

## Design Premises

1. **A library of pure functions, not a framework.** Each index is a standalone
   function `calculate_<index>(...)` that takes physical inputs and returns the
   index. No classes, no global state, no configuration objects. This keeps the
   API trivially composable and easy to audit against the literature.

2. **One SI unit contract for the whole library.** This is the central design
   decision of the 2.0 redesign and must be upheld:
   - Temperatures: **Kelvin** (`*_k`)
   - Wind speed: **metres per second** at 10 m (`va`)
   - Relative humidity: **percent** (`rh`)
   - Vapour pressure: **hectopascal** (hPa == mbar)
   Internally a formula may convert to °C/°F or km/h (via `helpers.py`), but it
   MUST convert back so the public boundary is always SI. Callers never have to
   think about units.

   *Documented exception:* the Excess Heat / Excess Cold Factors
   (`thermofeel.excess_heat`) are **unit-agnostic** — they operate on temperature
   *differences* and accept K or °C, returning the input unit *squared* (e.g.
   K²). This is intrinsic to those indices (they combine temperature anomalies),
   not a contract violation; it is called out in their docstrings and guide page.

3. **Vectorised over NumPy.** Every function must work elementwise on NumPy
   arrays of any shape, and equivalently on Python scalars wrapped in arrays.
   No Python-level loops over grid points. Masking and conditional branches use
   `np.where` / boolean indexing, never per-element `if`.

4. **Inputs limited to what the formula needs.** A function declares exactly the
   variables its formula consumes (e.g. `calculate_normal_effective_temperature`
   takes `t2_k, va, rh` and nothing else). Derived quantities are computed by
   calling the dedicated helper (e.g. relative humidity, vapour pressure) so
   there is one implementation of each physical relationship.

5. **Every formula is cited.** The docstring of each function names the
   reference (author, year) and a DOI/URL. A formula with no citation is not
   ready to merge.

6. **Solar geometry is an input, not a responsibility.** The cosine of the solar
   zenith angle (`cossza`) and direct radiation are provided by the caller. The
   2.0 release removed the in-library solar computation; use `earthkit-meteo`
   upstream. This keeps thermofeel's scope to thermal indices.

## Public API Shape

All functions are exported at the top level via `from .thermofeel import *`
(re-exported in `thermofeel/__init__.py`), so users write `thermofeel.calculate_utci(...)`.

Three families:

- **Comfort indices** — `calculate_utci`, `calculate_apparent_temperature`,
  `calculate_heat_index_simplified`, `calculate_heat_index_adjusted`,
  `calculate_humidex`, `calculate_normal_effective_temperature`,
  `calculate_wbgt`, `calculate_wbgt_simple`, `calculate_wbgt_liljegren`,
  `calculate_heat_force`, `calculate_excess_heat_factor`,
  `calculate_excess_cold_factor`, `calculate_wind_chill`.
- **Supporting physical quantities** — `calculate_mean_radiant_temperature`,
  `calculate_bgt`, `calculate_mrt_from_bgt`, `calculate_relative_humidity_percent`,
  `calculate_saturation_vapour_pressure`,
  `calculate_saturation_vapour_pressure_multiphase`,
  `calculate_nonsaturation_vapour_pressure`, `calculate_wbt`,
  `calculate_dew_point_from_relative_humidity`, `scale_windspeed`,
  `calculate_wind_speed_2m_liljegren`, `approximate_dsrp`.
- **Unit converters** (`helpers.py`) — `celsius_to_kelvin`, `kelvin_to_celsius`,
  `kelvin_to_fahrenheit`, `fahrenheit_to_celsius`, `fahrenheit_to_kelvin`.

`calculate_excess_heat_factor` / `calculate_excess_cold_factor` are thin
top-level wrappers over the `thermofeel.excess_heat` submodule, which also
provides the supporting per-day functions (`daily_mean_temperature`,
`significance_index`, `acclimatisation_index`, `heatwave_severity`); those stay
namespaced rather than re-exported at the top level.

## Key Design Decisions

### Wind speed is always supplied at 10 m

The public contract takes 10 m wind speed (`va`). Formulas that need wind at a
different height (1.1 m for the globe thermometer, 1.2 m for NET) call
`scale_windspeed(va, h)` internally. Callers never pre-scale. This was a
recurring source of bugs before 2.0 and is now centralised.

The Liljegren WBGT needs the 2 m wind. It defaults to the KNMI/Liljegren
stability-dependent power-law profile (`calculate_wind_speed_2m_liljegren`,
which derives a Pasquill-Gifford stability class from sun, radiation and wind),
matching KNMI's operational pipeline; `calculate_wbgt_liljegren(...,
wind_scaling="brode")` selects the generic `scale_windspeed` log profile
instead.

### Vapour pressure: saturated vs. non-saturated

Two distinct relationships exist and must not be confused:
- `calculate_saturation_vapour_pressure` (Hardy 1998) — for UTCI's water vapour
  pressure derived from RH.
- `calculate_nonsaturation_vapour_pressure` (Bureau of Meteorology) — for
  Apparent Temperature and simple WBGT, which want the actual (non-saturated)
  vapour pressure at a given RH.
Picking the wrong one silently biases the index. The choice per index is part of
the formula's specification and is covered by tests.

### Globe Temperature via a closed-form root

`calculate_bgt` solves the globe-temperature energy balance with a closed-form
quartic root rather than iteration (faster and deterministic over large grids).
The root is real-valued across the documented input domain for non-zero wind; at
exactly zero wind speed it returns `NaN` (documented in the docstring and in
`ROBUSTNESS.md`).

### Three WBGT implementations

The library offers three WBGT routines for different needs:
- `calculate_wbgt_simple` — single empirical regression (ACSM), temperature + RH.
- `calculate_wbgt` — globe temperature from the closed-form `calculate_bgt` plus
  the Stull wet bulb; needs `mrt` and `td`.
- `calculate_wbgt_liljegren` — the physically based Liljegren (2008) method that
  solves the globe and natural-wet-bulb energy balances by fixed-point
  iteration. This is the "gold standard" used operationally by KNMI and is the
  basis of `calculate_heat_force`. Its implementation (physical constants,
  property functions, the two energy-balance solvers, the stability-based 2 m
  wind profile) lives in the internal **`thermofeel/liljegren.py`** submodule,
  transcribed from Liljegren's reference C code and validated bit-for-bit
  against it; only the thin public wrappers stay in `thermofeel.py`. NaN is
  returned where the iteration does not converge.

### Output stays in SI even when the formula is empirical

Heat Index and Wind Chill are defined in °F / °C and km/h in their source
papers. The implementation converts at the boundaries and returns Kelvin, with
the validity range documented in the docstring (e.g. Wind Chill is only valid
for −50 °C…+5 °C and 5…80 km/h). Out-of-range inputs are not clamped — the
caller is responsible for masking, and the docstring states this.

### Numerical reproducibility is a feature

Stored CSV reference outputs (`tests/*.csv`) pin the numerical results. A change
that moves these values is either a deliberate formula correction (update the
reference, document in `CHANGELOG.md`, cite the reason) or a bug. The test
suite is the guard; see `TEST.md`.

## Dependencies

- **numpy** — the only runtime dependency; provides the vectorised maths.
- **pytest** — test-only (optional `[test]` extra).
- Formerly depended on `numba` (removed in 2.0 — the closed-form `bgt` and pure
  NumPy made it unnecessary) and on in-library solar-angle code (removed in 2.0,
  now an input).

## Non-Goals

- Not a plotting or I/O library (examples may use `eccodes`/`matplotlib`, but the
  core does not).
- Not a solar-position calculator (caller supplies `cossza`).
- Not a units framework — the SI contract is a convention enforced by review and
  tests, deliberately not a `pint`-style runtime unit system (keeps the core
  dependency-light and fast).

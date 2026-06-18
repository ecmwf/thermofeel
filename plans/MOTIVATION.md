# Motivation: Why thermofeel Exists

## The Problem

Human thermal comfort and heat/cold stress are not measured directly by weather
models — they are *derived* from the basic forecast variables (temperature,
humidity, wind, radiation) through published thermal-index formulas. Those
formulas are scattered across decades of biometeorology literature, each with
its own unit conventions, input assumptions, and validity ranges.

Anyone wanting to turn a weather forecast into a heat-stress product hits the
same recurring friction:

1. **The formulas are easy to get subtly wrong.** Wind speed at 10 m vs. at the
   height of a globe thermometer; vapour pressure from a saturated vs. a
   non-saturated formula; °C vs. K vs. °F; km/h vs. m/s. A transcription slip
   produces a plausible-looking but wrong index.
2. **Implementations are not reusable.** A formula buried in a one-off script,
   in non-SI units, with no tests and no reference citation, cannot be trusted or
   reused by the next pipeline.
3. **They need to run at forecast scale.** A thermal index must be computed over
   global gridded fields (millions of points, many forecast steps), so the
   implementation has to be vectorised, not a per-point Python loop.

## What we are building

A **library of human thermal comfort indices** that is:

1. **Correct and cited** — every index implements a specific peer-reviewed
   formula, named and linked in its docstring.
2. **Consistent in units** — a single SI-based input/output contract across the
   whole library (see `DESIGN.md`): temperatures in Kelvin, wind in m/s,
   humidity in %, pressures in hPa.
3. **Vectorised** — built on NumPy so the same call works on a scalar, a 1-D
   series, or a full lat/lon/step grid, with no Python-level loops.
4. **Minimal-dependency** — the core depends only on NumPy, so it drops cleanly
   into an operational NWP pipeline.

It is not a weather model and not a plotting tool — it is the trusted middle
layer that turns forecast variables into thermal-stress variables.

## Currently computed

Thermal comfort indices:

- Universal Thermal Climate Index (UTCI)
- Apparent Temperature
- Heat Index (simplified and adjusted)
- Humidex
- Normal Effective Temperature
- Wet Bulb Globe Temperature (full and simple)
- Wind Chill

Supporting physical quantities:

- Globe Temperature
- Mean Radiant Temperature (and from Globe Temperature)
- Relative Humidity percentage
- Saturation / non-saturation vapour pressure
- Wet Bulb Temperature
- Dew point from relative humidity
- Wind-speed scaling between heights

## Target users

- Numerical weather prediction and climate centres producing heat/cold-stress
  forecast products (the original ECMWF use case).
- Researchers in biometeorology, public health, and climate-impact studies who
  need a reproducible, citable implementation of standard indices.
- Anyone with their own meteorological data who wants comparable thermal indices
  without re-deriving the formulas.

The goal is transparency and reproducibility: a published, tested, referenced
implementation that users can adopt directly or audit against the literature.

## Constraints

- SI units in and out, consistently (the 2.0 standardisation).
- Vectorised over NumPy arrays; no mandatory dependency beyond NumPy.
- Each formula traceable to a published reference.
- Solar-geometry inputs (cosine of solar zenith angle) are provided by the
  caller (e.g. via `earthkit-meteo`), not recomputed here — thermofeel stays
  focused on the thermal indices themselves.

## Success criteria

1. **Correctness:** every index reproduces its reference formula and known worked
   examples within a documented tolerance.
2. **Consistency:** all functions honour the SI input/output contract.
3. **Scalability:** the same function runs unchanged on scalars and on global
   gridded arrays.
4. **Reproducibility:** results are stable across releases unless a formula is
   deliberately corrected, and every such correction is recorded in
   `ChangeLog.rst`.

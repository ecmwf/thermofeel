# Numerical Robustness & Failure Modes

> **Status:** first hardening pass done (2.2.0). The `np.log`/`np.sqrt`/
> `np.power`/division domains were audited, the guards confirmed, and
> `NaN`/`Inf` behaviour pinned by `tests/test_robustness.py`. Findings are in §5.
> This document records how thermofeel behaves under bad, edge, or out-of-domain
> input: the failure modes to guard against, plus the checklist and findings log.

## 0. Why robustness is the risk surface

thermofeel is a **pure-Python numerical library** whose only runtime dependency
is NumPy. Its core does not parse external byte streams, manage native memory or
use FFI, execute or deserialise anything, or touch the network or filesystem.
The entire risk surface is therefore **numerical**: does each formula return a
correct, finite, well-defined value across the range of inputs it can be handed,
and does it behave predictably at the edges?

(The `examples/` scripts read GRIB/NetCDF via external tools, but examples are
not part of the distributed library surface.)

## 1. Scope

- **In scope:** `thermofeel/thermofeel.py`, `thermofeel/helpers.py`,
  `thermofeel/__init__.py`.
- **Out of scope:** `examples/`, `docs/`, `tests/`.

## 2. Failure modes

thermofeel runs inside trusted pipelines on trusted arrays, so the realistic
failures come from **bad or out-of-domain input data**.

| # | Class | Concrete vector in this library |
|---|-------|---------------------------------|
| A | Silent wrong answer from out-of-range input | Wind Chill / Heat Index evaluated outside their documented validity windows return plausible but invalid numbers |
| B | Domain error producing `NaN`/`Inf` | `np.log`/`np.sqrt`/`np.power` on non-positive or negative arguments (e.g. RH = 0, negative discriminant in the `bgt` closed-form root) |
| C | Division by near-zero | `approximate_dsrp` dividing by `cossza` as it approaches zero (already guarded by a `threshold`); RH denominators |
| D | Unit-contract violation | A caller passing °C where Kelvin is expected — produces a wrong index with no error |

## 3. Robustness invariants

The unit contract (SI in/out), the "validity ranges are documented, not clamped"
rule, and determinism are design guarantees — see `DESIGN.md`. The
robustness-specific guarantees this document tracks, layered on top, are:

1. **No crash on legitimate array input.** Functions operate elementwise and must
   not raise on finite, correctly-shaped arrays; domain edge cases yield
   `NaN`/`Inf` consistent with NumPy semantics rather than exceptions — except
   where an explicit `ValueError` is the correct contract (e.g. UTCI given
   neither `ehPa` nor `td_k`).
2. **`NaN`/`Inf` behaviour is defined, not incidental.** Where an input drives a
   formula out of its mathematical domain, the resulting `NaN`/`Inf` is a
   documented, tested outcome rather than a surprise.

## 4. Robustness checklist (first hardening pass)

- [x] Audit every `np.log`, `np.sqrt`, `np.power`, and division for the input
      domain that produces `NaN`/`Inf`; document the behaviour or guard it where
      a guard is clearly correct (do not mask real out-of-domain use). → R-1…R-6.
- [x] Confirm `approximate_dsrp`'s `threshold` guard fully prevents
      divide-by-zero across the input range, and that the documented large-error
      regime near low sun angles is called out. → R-5 (guard confirmed).
- [x] Confirm the `calculate_bgt` closed-form root is real-valued (no negative
      discriminant) across the documented input domain, including zero/low wind.
      → R-1 (real for `va > 0`; `NaN` at exactly `va == 0`, now documented).
- [x] Add explicit tests for `NaN`/`Inf` propagation behaviour so it is a
      defined contract rather than incidental. → `tests/test_robustness.py`.

## 5. Findings log

> Each entry: ID, severity, surface, description, status, mitigation, regression
> test.

- **R-1 — MEDIUM — `calculate_bgt` (closed-form globe temperature) — FIXED
  (2.2.0).** At exactly zero wind (`va == 0`) the convective term `d` vanishes,
  so the closed form is a `0/0` indeterminate (previously `NaN` with a NumPy sqrt
  warning). The analytic limit as `va -> 0` is the mean radiant temperature (no
  convection ⇒ globe at radiative equilibrium, `bgt -> mrt`), so `mrt` is now
  returned at `va == 0` and the warnings are suppressed via `np.errstate`. This
  also makes the only internal consumer, `calculate_wbgt`, finite at `va == 0`.
  Invalid negative wind still yields `NaN`. No reference data changed (the test
  inputs have `va >= 0.02`; all `va > 0` outputs are byte-identical). **Tests:**
  `test_bgt_zero_wind_returns_mrt`, `test_wbgt_zero_wind_is_finite`,
  `test_bgt_negative_wind_is_nan`.
- **R-2 — LOW — all public indices, `NaN` propagation.** A `NaN` input yields a
  `NaN` output element-wise; no function raises on finite, correctly-shaped
  arrays. **Status:** confirmed and pinned. **Test:**
  `test_nan_temperature_propagates`.
- **R-3 — LOW — `liljegren.solve_globe`/`solve_wetbulb` non-convergence.** Each
  element that does not converge within `MAX_ITER` (e.g. a `NaN` input, which
  never satisfies the tolerance) returns `NaN`; converging neighbours are
  unaffected (per-element). **Status:** by design, confirmed. **Test:**
  `test_wbgt_liljegren_nan_element_propagates`.
- **R-4 — LOW (by design) — `np.log` at RH = 0.**
  `calculate_dew_point_from_relative_humidity` and Liljegren `dew_point` evaluate
  `np.log` of a vapour pressure that is 0 when RH = 0 → `-inf`. Documented domain
  is RH > 0. **Status:** accepted; no guard (masking would hide real bad input).
- **R-5 — guarded, no action — divide-by-(near-)zero.** `approximate_dsrp`
  divides by `cossza` only where `cossza > threshold` (default 0.1); the Liljegren
  solvers guard `cza → 0` (`cza_safe`), `tan(sza)`, and the `0 * inf` direct-beam
  term. **Status:** guards reviewed and confirmed correct.
- **R-6 — LOW — `liljegren.solve_wetbulb` `(pair - ewick)` denominator.** Safe
  for the meteorological domain (`ewick` ≪ `pair`); would only approach zero near
  water's boiling point, far outside valid wet-bulb inputs. **Status:** accepted.

## 6. Severity definitions

- **HIGH** — a crash/exception on legitimate finite array input, or a silent
  wrong answer for in-domain input.
- **MEDIUM** — undocumented `NaN`/`Inf` for plausible inputs; missing validity
  documentation.
- **LOW** — hardening, hygiene, documentation gaps.

## 7. Recommended follow-ups

- Run the §4 checklist and record results in §5.
- Add `hypothesis` property tests for physical-sanity invariants (e.g. UTCI
  increases with temperature at fixed other inputs; indices stay finite over the
  documented input domain). See `IDEAS.md`.

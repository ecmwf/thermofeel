# Numerical Robustness & Failure Modes

> **Status:** scaffold — no formal robustness pass has been run yet.
> This document records how thermofeel behaves under bad, edge, or out-of-domain
> input: the failure modes to guard against, plus a checklist and findings log
> for a first hardening pass.

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
- **Out of scope:** `examples/`, `docs/`, `tests/`. `thermofeel/experimental_wbgt.py`
  is experimental and unexported — see `TODO.md`.

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

- [ ] Audit every `np.log`, `np.sqrt`, `np.power`, and division for the input
      domain that produces `NaN`/`Inf`; document the behaviour or guard it where
      a guard is clearly correct (do not mask real out-of-domain use).
- [ ] Confirm `approximate_dsrp`'s `threshold` guard fully prevents
      divide-by-zero across the input range, and that the documented large-error
      regime near low sun angles is called out.
- [ ] Confirm the `calculate_bgt` closed-form root is real-valued (no negative
      discriminant) across the documented input domain, including zero/low wind.
- [ ] Add explicit tests for `NaN`/`Inf` propagation behaviour so it is a
      defined contract rather than incidental.

## 5. Findings log

> Each entry: ID, severity, surface, description, status, mitigation, regression
> test. None recorded yet — populate as the checklist above is worked through.

_(no robustness pass run yet)_

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

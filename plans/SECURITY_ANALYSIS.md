# Security & Robustness Analysis

> **Status:** scaffold — no formal audit pass has been run yet.
> This document records the (deliberately light) threat model for thermofeel and
> the input-robustness checklist that stands in for a heavyweight security audit.

## 0. Why this is light

thermofeel is a **pure-Python numerical library** whose only runtime dependency
is NumPy. It does **not**:

- parse untrusted byte streams or files (no wire format, no deserialisation),
- contain `unsafe`/native memory management or hand-written FFI,
- execute, `eval`, `pickle`, or fetch anything,
- hold credentials, open sockets, or touch a filesystem in its core.

So the memory-corruption / decompression-bomb / FFI-contract threat classes that
dominate a binary-format library do not apply. The realistic concerns are
**numerical robustness** and **supply chain**, addressed below. (The `examples/`
scripts use `eccodes`/`matplotlib`/`numpy` I/O, but examples are not part of the
distributed library surface.)

## 1. Scope

- **Core library:** `thermofeel/thermofeel.py`, `thermofeel/helpers.py`,
  `thermofeel/__init__.py`.
- **Out of scope:** `examples/` (illustrative, may read GRIB/NetCDF via external
  tools), `docs/`, `tests/`. `thermofeel/experimental_wbgt.py` is experimental
  and unexported — see `TODO.md`.

## 2. Threat / failure model

The "adversary" here is realistically **bad or out-of-domain input data**, not a
malicious actor — thermofeel runs inside trusted pipelines on trusted arrays.

| # | Class | Concrete vector in this library |
|---|-------|---------------------------------|
| A | Silent wrong answer from out-of-range input | Wind Chill / Heat Index evaluated outside their documented validity windows return plausible but invalid numbers |
| B | Domain error producing `NaN`/`Inf` | `np.log`/`np.sqrt`/`np.power` on non-positive or negative arguments (e.g. RH = 0, negative discriminant in the `bgt` closed-form root) |
| C | Division by near-zero | `approximate_dsrp` dividing by `cossza` as it approaches zero (already guarded by a `threshold`); RH denominators |
| D | Unit-contract violation | A caller passing °C where Kelvin is expected — produces a wrong index with no error |
| E | Supply chain | a CVE in NumPy or the build/test toolchain |

## 3. Robustness invariants (the contract we aim to keep)

1. **No silent unit surprises.** Inputs and outputs are SI; every function
   documents its units. (Enforced by review + `DESIGN.md`, not at runtime.)
2. **Documented validity ranges.** Any index defined only over a sub-domain
   states that range in its docstring; the library does **not** clamp — masking
   is the caller's responsibility, and that responsibility is documented.
3. **No crash on array input.** Functions operate elementwise and should not
   raise on legitimately-shaped finite arrays; domain edge cases yield
   `NaN`/`Inf` consistent with NumPy semantics rather than exceptions, except
   where an explicit `ValueError` is the correct contract (e.g. UTCI with neither
   `ehPa` nor `td_k`).
4. **Determinism.** Same inputs → same outputs across runs and platforms (pinned
   by the `tests/*.csv` regression baselines).

## 4. Robustness checklist (to work through in a first pass)

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
- [ ] Wire `pip-audit` (or equivalent) into CI to catch dependency CVEs.

## 5. Findings log

> Each entry: ID, severity, surface, description, status, mitigation, regression
> test. None recorded yet — populate as the checklist above is worked through.

_(no audit pass run yet)_

## 6. Severity definitions

- **HIGH** — a crash/exception on legitimate finite array input, or a silent
  wrong answer for in-domain input.
- **MEDIUM** — undocumented `NaN`/`Inf` for plausible inputs; missing validity
  documentation; a dependency CVE reachable from the core.
- **LOW** — defence-in-depth, hygiene, documentation gaps.

## 7. Recommended follow-ups

- Run the §4 checklist and record results in §5.
- Add `hypothesis` property tests for physical-sanity invariants (see `IDEAS.md`).
- Add `pip-audit` to the CI lint/QA stage.

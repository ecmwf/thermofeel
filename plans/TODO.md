# Features Decided to Implement

Accepted features and fixes that are planned but not yet done. Each entry should
carry enough intent that an implementer knows the aim and how to verify it. Code
agents are encouraged to ask questions to get the design right before coding.

For speculative, not-yet-accepted ideas, see `IDEAS.md`.

## Features

- [ ] **Output-unit mode** (the long-standing units-of-output request)
    Add a mode to control the units of the returned values. Today every index
    returns SI (Kelvin); some users want °C directly. Design questions to settle
    before implementing:
    - Per-call argument (`units="K"|"C"`) vs. a module-level setting? A per-call
      argument keeps functions pure and is preferred over global state.
    - Does it apply only to temperature-like outputs, or also to vapour pressure
      / wind? Probably temperature-like only; document explicitly.
    - The SI-in/SI-out contract in `DESIGN.md` stays the default; the mode is an
      opt-in convenience at the output boundary only.
    Verify: round-trip tests that `units="C"` equals `kelvin_to_celsius(units="K")`
    for every affected index.

- [ ] **Discomfort Index** (new thermal index)
    Add a `calculate_*` function for a Discomfort Index, promoted from the
    "More thermal indices" idea in `IDEAS.md`. The exact variant, its
    peer-reviewed reference, the input variables it consumes, and any
    discomfort/validity bands are to be settled with the maintainer before
    coding (there are several discomfort indices — e.g. Thom's
    Temperature-Humidity Index). It must follow the library contract: a pure,
    vectorised NumPy function with SI inputs/outputs (temperature in Kelvin,
    relative humidity in %), converting at the boundary if the source formula is
    in °C, and citing the reference in its docstring.
    Verify: pointwise scalar checks against published worked examples
    (`tests/test_scalars.py`) plus an array regression case with a stored
    reference CSV (`tests/test_thermofeel.py`); a guide page under `docs/guide/`
    and a runnable `examples/` snippet; add it to the index list in `README.md`,
    the `thermofeel/thermofeel.py` module docstring and `docs/guide/overview.md`
    (kept in step); and a `CHANGELOG.md` entry.

## Robustness

- The first numerical-robustness hardening pass is **done** (see
  `ROBUSTNESS.md` §5 and `tests/test_robustness.py`), including the R-1 calm-air
  fix for `calculate_bgt`. Possible follow-up, not yet accepted: `hypothesis`
  property tests (see `IDEAS.md`).

## Notes

- This file (`plans/TODO.md`) is the single canonical accepted backlog. It
  supersedes the former repo-root `TODO.rst` and `TODO.md`, whose contents have
  been folded in here.

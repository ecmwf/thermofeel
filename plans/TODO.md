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

## Robustness

- [ ] **Work through the numerical-robustness checklist in `ROBUSTNESS.md`.**
    The first hardening pass: audit `np.log`/`np.sqrt`/`np.power`/division for
    `NaN`/`Inf`-producing domains, confirm the `approximate_dsrp` divide-by-zero
    guard, confirm the `calculate_bgt` root stays real-valued, and add explicit
    `NaN`/`Inf` propagation tests. Record findings in `ROBUSTNESS.md`.

## Notes

- This file (`plans/TODO.md`) is the single canonical accepted backlog. It
  supersedes the former repo-root `TODO.rst` and `TODO.md`, whose contents have
  been folded in here.

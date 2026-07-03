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

## Indices — the 2.4.0 batch programme

Accepted batch of new thermal indices, promoted from the IDEAS assessment. The
full programme spec (provenance, formulas, validation gates, agent topology,
open decisions) is in **`plans/NEW_INDICES.md`** — that is the source of truth;
this is the backlog pointer. Delivered as a single `2.4.0` PR, fanned out one
agent per index, each proven by its own validation gates (G1–G4 in
`NEW_INDICES.md` §3).

- [ ] **Apparent Temperature (radiation form)** — `calculate_apparent_temperature_radiation`
      (Steadman 1994 / BoM). Gate status: ready. Spec §5.1.
- [ ] **Relative Strain Index** — `calculate_relative_strain_index`
      (Lee & Henschel; peer-reviewed hPa form). Gate status: ready. Spec §5.2.
- [ ] **Summer Simmer Index** — `calculate_summer_simmer_index` (Pepi 1987).
      Gate status: conditional on provenance (analytic THI gate). Spec §5.3.
- [ ] **PMV / PPD** — `calculate_pmv`, `calculate_ppd` (ISO 7730:2005). Gate
      status: ready; richest validation table. Spec §5.4.
  PET (`calculate_pet`) was assessed for this batch but is **deferred to
  `IDEAS.md`** (not citable from open sources; no published reference values) —
  see `NEW_INDICES.md` §5.5/§6.

## Robustness

- The first numerical-robustness hardening pass is **done** (see
  `ROBUSTNESS.md` §5 and `tests/test_robustness.py`), including the R-1 calm-air
  fix for `calculate_bgt`. Possible follow-up, not yet accepted: `hypothesis`
  property tests (see `IDEAS.md`).

## Notes

- This file (`plans/TODO.md`) is the single canonical accepted backlog. It
  supersedes the former repo-root `TODO.rst` and `TODO.md`, whose contents have
  been folded in here.

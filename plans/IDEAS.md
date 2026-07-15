# Ideas

> [!Warning]
> This file is a dump of ideas. It contains speculative thoughts, some generated
> by AI agents, that may not reflect the vision of the maintainers. Nothing here
> is accepted. Do not implement until an idea is promoted to `TODO.md`.

Speculative directions for possible future work, grouped loosely. Ideas that
mature and get decided move up to `TODO.md`.

## Indices and physics

- [ ] **More thermal indices.** Speculative candidates that need a citable
  equation set + reference values before they can leave IDEAS. **Standard
  Effective Temperature (SET)** — Gagge two-node model; heavy, and overlaps PMV.

  *Already implemented in 2.3.0 (see `CHANGELOG.md` and the per-index guide pages
  under `docs/guide/`): Apparent Temperature radiation form, Relative Strain
  Index, Summer Simmer Index, PMV/PPD. WBGT via the full Liljegren model is
  `calculate_wbgt_liljegren`; the Discomfort Index is `calculate_discomfort_index`.*

- [ ] **Physiological Equivalent Temperature (PET).** Höppe (1999),
  DOI `10.1007/s004840050118`. The outdoor "gold-standard" companion to UTCI:
  PET is the air temperature of a reference indoor setting (MRT = Ta, still air
  v = 0.1 m/s, vapour pressure = 12 hPa) in which a reference person reaches the
  same core and skin temperature as in the actual environment, solved from the
  MEMI human energy-balance model (Mayer & Höppe 1987). Assessed for the 2.3.0
  batch and **deferred**; the full details are retained below so it can be picked
  up directly.
  - **Inputs (SI):** `t2_k` (K), `mrt_k` (K), 10 m wind (m/s), relative humidity
    (%), plus a standardised reference person (Höppe uses ~35-yr male, 1.75 m,
    75 kg, ~80 W work metabolism, 0.9 clo — verify each value against a primary
    source). Output is an equivalent temperature in °C → return **Kelvin**.
    Proposed signature:
    `calculate_pet(t2_k, mrt_k, va, rh, *, pressure_hpa=1013.25, met=…, clo=0.9, …)`.
  - **Why it is deferred (the blocker):** the full MEMI energy-balance equations
    are **not visible in open primary sources**, and **no published verbatim
    `(inputs → PET)` validation rows** were found during research. The only
    available implementations are `pythermalcomfort.pet_steady` (MIT, but SciPy
    `fsolve` — not numpy-only) and Ladybug (**AGPL — must not be copied**).
    Shipping now would mean transcribing equations we cannot cite and validating
    only against another implementation's output, which fails `DESIGN.md` §5
    ("a formula with no citation is not ready to merge") and this file's own gate.
  - **Unblock condition:** obtain a fully-cited equation set — **Walther &
    Goestchel (2018)** `10.1016/j.buildenv.2018.03.054` ("first exhaustive
    explanation of PET"; the corrected steady-state model) and/or **VDI 3787
    Part 2** — via ECMWF institutional access, **and** a set of published
    reference values (VDI worked examples / Höppe 1999 or Matzarakis et al. 1999
    tables). Then implement **independently** from those equations: numpy-only,
    bounded deterministic iteration (mirror the `liljegren.py` solver),
    non-convergence → `NaN`; no SciPy, no AGPL code.
  - **Categories (docs):** the PET thermal-perception / physiological-stress
    bands are from Matzarakis & Mayer (1996).
  - **Pitfalls:** PET vs. mPET (modified PET) are different indices; RayMan and
    `pythermalcomfort` are documented to disagree numerically — pin the *cited*
    equation set, not a tool's output. `pythermalcomfort.pet_steady` may be used
    as a development cross-check oracle only, never a runtime dependency.
- [ ] **Validity-range helpers.** Several indices are only defined over a range
  (Wind Chill, Heat Index). A small companion that returns a validity mask for a
  given index + inputs would let pipelines blank out-of-range points cleanly
  instead of each caller re-deriving the bounds.
- [ ] **Vapour-pressure formula audit.** Document/justify the choice of each
  empirical vapour-pressure relation (Hardy, Bureau of Meteorology, IFS
  multiphase) and whether a single consistent one could serve more indices
  without changing results meaningfully.

## API & ergonomics

- [ ] **Uniform scalar input support.** Today pure-arithmetic functions accept
  bare Python scalars while the internally-branching ones
  (`calculate_heat_index_simplified`, `approximate_dsrp`,
  `calculate_saturation_vapour_pressure_multiphase`) require arrays. Coercing
  inputs at the boundary (`np.atleast_1d(np.asarray(x, dtype=float))`) would let
  every function accept scalars uniformly — decide whether to return a scalar or
  a 1-element array for scalar input. (The array contract is documented for now.)
- [ ] **Optional unit-aware boundary.** Beyond the accepted output-unit mode in
  `TODO.md`, explore an optional `pint`-style typed-quantity wrapper as a
  *separate* opt-in layer (never in the dependency-light core) for users who
  want unit safety.
- [ ] **`xarray` / `dask` friendliness.** The functions are already elementwise
  NumPy, so they likely work under `xarray.apply_ufunc` / dask arrays today.
  Verify, document the pattern, and add an example — without taking a hard
  dependency.

## Performance

- [ ] **Benchmark suite.** A small reproducible benchmark over a global-grid-sized
  array per index, to catch performance regressions and to justify (or retire)
  the closed-form `bgt` vs. iterative trade-off quantitatively.
- [ ] **Optional acceleration.** If a hot path appears, evaluate an optional
  accelerated backend (e.g. `numba`/`numexpr`) strictly as an extra — the 2.0
  redesign deliberately removed the mandatory `numba` dependency, so any return
  must stay opt-in.

## Quality & tooling

- [ ] **Property-based tests** (`hypothesis`): monotonicity and physical-sanity
  properties (e.g. UTCI increases with temperature at fixed other inputs; indices
  stay finite over the documented input domain).
- [ ] **Cross-implementation parity.** Compare outputs against an independent
  reference implementation of each index (where one exists) as a one-off
  validation artefact, separate from the regression CSVs.

## Documentation & outreach

- [ ] **Worked-example notebook per index** in `docs/`, showing inputs, the
  formula, the citation, and a small plot of the index over a parameter sweep.
- [ ] **A "choosing an index" guide** — a short decision page on which thermal
  index suits which question (heat stress vs. cold stress, outdoor vs. indoor,
  data availability), since users often pick the wrong one.

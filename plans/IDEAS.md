# Ideas

> [!Warning]
> This file is a dump of ideas. It contains speculative thoughts, some generated
> by AI agents, that may not reflect the vision of the maintainers. Nothing here
> is accepted. Do not implement until an idea is promoted to `TODO.md`.

Speculative directions for possible future work, grouped loosely. Ideas that
mature and get decided move up to `TODO.md`.

## Indices and physics

- [ ] **More thermal indices.** Speculative candidates that need a citable
  equation set + reference values before they can leave IDEAS:
  - **Standard Effective Temperature (SET)** — Gagge two-node model; heavy,
    overlaps PMV.
  - **Physiological Equivalent Temperature (PET)** — Höppe 1999. Assessed for the
    2.3.0 batch but **deferred**: the full MEMI equations are not openly citable
    and no published `(inputs → PET)` validation rows were found; reference
    implementations are SciPy-based or AGPL. Revisit only with a fully-cited
    equation set (Walther & Goestchel 2018 `10.1016/j.buildenv.2018.03.054` /
    VDI 3787-2) **and** published reference values; never reverse-engineer the
    AGPL code. Retained spec: `plans/NEW_INDICES.md` §5.5.

  *Promoted out of IDEAS into the 2.3.0 batch (`plans/NEW_INDICES.md`): Apparent
  Temperature radiation form, Relative Strain Index, Summer Simmer Index,
  PMV/PPD.* (WBGT via the full Liljegren model is implemented as
  `calculate_wbgt_liljegren`; the Discomfort Index as `calculate_discomfort_index`.)
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

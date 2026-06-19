# Ideas

> [!Warning]
> This file is a dump of ideas. It contains speculative thoughts, some generated
> by AI agents, that may not reflect the vision of the maintainers. Nothing here
> is accepted. Do not implement until an idea is promoted to `TODO.md`.

Speculative directions for possible future work, grouped loosely. Ideas that
mature and get decided move up to `TODO.md`.

## Indices and physics

- [ ] **More thermal indices.** Candidates frequently requested in
  biometeorology: Physiological Equivalent Temperature (PET), Standard Effective
  Temperature (SET), Discomfort Index (Thom). Each needs a citation and reference
  values before it can leave IDEAS. (WBGT via the full Liljegren model is now
  implemented as `calculate_wbgt_liljegren`.)
- [ ] **Validity-range helpers.** Several indices are only defined over a range
  (Wind Chill, Heat Index). A small companion that returns a validity mask for a
  given index + inputs would let pipelines blank out-of-range points cleanly
  instead of each caller re-deriving the bounds.
- [ ] **Vapour-pressure formula audit.** Document/justify the choice of each
  empirical vapour-pressure relation (Hardy, Bureau of Meteorology, IFS
  multiphase) and whether a single consistent one could serve more indices
  without changing results meaningfully.

## API & ergonomics

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

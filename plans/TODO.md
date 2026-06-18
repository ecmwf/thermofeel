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

## Correctness fixes (from codebase review)

- [ ] **Repair the stale example `examples/compute-thermal-indices.py`.**
    It calls functions that no longer exist in the 2.x API and would crash:
    - `thermofeel.calculate.heat_index_adjusted` — wrong namespace
      (should be `thermofeel.calculate_heat_index_adjusted`)
    - `calculate_cos_solar_zenith_angle_integrated` — removed in 2.0; solar
      geometry is now a caller input (use `earthkit-meteo`)
    - `calculate_net_effective_temperature` — renamed to
      `calculate_normal_effective_temperature`
    - `kelvin_to_celcius` — misspelling of `kelvin_to_celsius`
    - `thermofeel.func_timers` / the `@thermofeel.timer` decorator — no longer exist
    - old UTCI kwargs `va_ms=`, `mrt_k=`, `e_hPa=` — now `va=`, `mrt=`, `ehPa=`
    Bring it in line with the current API or clearly mark it legacy. An example
    that calls a removed function is a bug, not documentation (see `AGENTS.md`).

- [ ] **Fix the `calculate_bgt` argument order in `tests/test_scalars.py::test_bgt`.**
    It is called as `calculate_bgt(t2_k, va, mrt)` but the signature is
    `calculate_bgt(t2_k, mrt, va)`. The test passes only because the chosen
    values are near-symmetric; the "negative va treated as 0" comment does not
    reflect the actual call. Correct the order and make the comment match.

## Code cleanliness

- [ ] **Remove dead code in `scale_windspeed`** (`thermofeel/thermofeel.py`).
    `c = 1 / np.log10(10 / 0.01)` is immediately overwritten by the hardcoded
    `c = 0.333333333333`. Keep one — preferably the computed expression with a
    comment — and drop the misleading first line.

- [ ] **Fix the `approximate_dsrp` comment.** It says "for cossza <= 0.01" while
    the `threshold` parameter defaults to `0.1`. Make the comment reference the
    parameter so the two cannot drift.

- [ ] **Sync the docs version.** `docs/source/conf.py` has `release = "v1"`,
    stale and out of step with `thermofeel/__init__.py` `__version__`. Either set
    it from the package version (single source of truth) or add a check. Also
    review the Python-version trove classifiers in `pyproject.toml` (claim
    3.6–3.13; the library runs on 3.14).

## Experimental code

- [ ] **Decide the fate of `thermofeel/experimental_wbgt.py` (Liljegren).**
    Its iterative solvers use `np.where(np.abs(...))` as a convergence mask;
    `np.abs(x)` returns values, not a boolean condition, so `np.where` selects
    all non-zero elements every iteration — almost certainly not the intended
    convergence test. Either:
    - finish it properly (correct the convergence loop, add reference-value tests
      and a citation, then export it through the public API), or
    - keep it clearly marked experimental and out of the public surface.
    Do not leave it in a half-working, untested middle state.

## Notes

- This file (`plans/TODO.md`) is the single canonical accepted backlog. It
  supersedes the former repo-root `TODO.rst` and `TODO.md`, whose contents have
  been folded in here.

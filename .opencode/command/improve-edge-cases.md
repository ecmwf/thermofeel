---
description: Systematic numerical edge-case audit and patching across the library
agent: build
---

# Improve Edge Case Handling

Perform a systematic edge-case audit across the library. thermofeel is a pure
numerical library, so the risk surface is numerical: does each formula return a
correct, finite, well-defined value across the range of inputs it can be handed,
and does it behave predictably at the edges? Use `plans/ROBUSTNESS.md` (findings
R-1…R-8) and `tests/test_robustness.py` as the baseline — extend them, don't
duplicate them.

## What to Look For

### Input-value edge cases
- Empty arrays, single-element arrays, 0-D scalars wrapped in arrays, very large
  arrays
- `NaN`, `+Inf`, `-Inf`, negative zero, subnormals — confirm they propagate
  elementwise to `NaN`/`Inf` rather than raising
- Boundary humidity (`rh = 0` → `np.log` domain), boundary/zero/negative wind
  (the `calculate_bgt` closed-form root; the calm-air `bgt -> mrt` limit),
  `cossza → 0` (the `approximate_dsrp` threshold guard and the Liljegren
  `cza_safe` / `tan(sza)` guards)
- Integer-typed input arrays — must not silently truncate a float result (see
  R-7); preserve dtype with `dtype=float` where results are fractional

### Domain / validity-range cases
- Each index that is only defined over a range (Wind Chill −50…+5 °C and
  5…80 km/h; the Heat-Index 20 °C gate and the adjusted Heat-Index cold / low-
  humidity / high-humidity branches; UTCI polynomial ranges; Stull wet-bulb
  RH 5…99 %). Inputs are **not** clamped — confirm the documented behaviour
  (out-of-range in → out-of-range out) and that it is stated in the docstring
- The `calculate_utci` `ValueError` when neither `ehPa` nor `td_k` is supplied
- Liljegren non-convergence returning `NaN` (the C reference's −9999 sentinel)
  and `calculate_heat_force` propagating that `NaN` rather than a spurious band

### Calling-convention cases
- Functions that branch internally with `np.where` / boolean indexing
  (`calculate_heat_index_simplified`, `approximate_dsrp`,
  `calculate_saturation_vapour_pressure_multiphase`) require array input —
  confirm the array contract is documented and behaves consistently

## Process

1. Scan each module and identify unhandled or under-tested edge cases
2. For ambiguous behaviour, ask the user to clarify the intended semantics
   (mask vs. propagate vs. raise) before changing anything — do NOT add a guard
   that masks genuine out-of-domain input
3. Add handling only where a guard is clearly correct (e.g. the analytic calm-air
   limit), preserving the "validity ranges are documented, not clamped" rule
4. Write tests for each edge case in `tests/test_robustness.py` (and the scalar
   suite where a specific value matters)
5. Document the behaviour: update the function docstring, `docs/`, and the
   `plans/ROBUSTNESS.md` findings log (ID, severity, surface, status, the
   regression test that pins it)
6. Run `make all` to verify no regressions
7. Summarize findings and changes (grouped, with the ROBUSTNESS.md IDs)

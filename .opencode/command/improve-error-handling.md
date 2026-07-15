---
description: Audit the error / NaN-Inf contracts across the library and tighten them up
agent: build
---

# Improve Error Handling

Perform an error-handling audit across the library. thermofeel runs inside
trusted pipelines on trusted NumPy arrays; it does not parse byte streams, use
FFI, or touch the network/filesystem in its core. The "error model" is therefore
about two things: the few explicit `ValueError` contracts, and the defined
`NaN`/`Inf` behaviour at domain edges. See `plans/ROBUSTNESS.md` §3 for the
invariants.

## Checks

### Exception contracts
- `ValueError` is raised where it is the correct contract, with a clear,
  actionable message — e.g. `calculate_utci` given neither `ehPa` nor `td_k`;
  `calculate_wbgt_liljegren` given an unknown `wind_scaling`. Confirm each is
  tested in `tests/test_robustness.py`
- Functions do NOT raise on finite, correctly-shaped array input — domain edge
  cases yield `NaN`/`Inf` per NumPy semantics, not exceptions
- No bare `except:` and no overly broad exception swallowing anywhere in
  `thermofeel/`

### NaN / Inf as a defined contract
- A `NaN` input yields a `NaN` output elementwise (never a spurious finite
  value, never a raise); converging neighbours in an array are unaffected
- Where an input drives a formula out of its mathematical domain (`rh = 0` →
  `log`; non-convergent Liljegren iteration; the zero-wind `bgt` indeterminate),
  the resulting `NaN`/`Inf` is the **documented, tested** outcome — not incidental
- Guards that exist to avoid spurious warnings (`np.errstate`, `cza_safe`, the
  `approximate_dsrp` threshold) must not mask genuinely bad input (e.g. negative
  wind still yields `NaN`)

### Documentation
- Each public function's docstring states its units, validity range, and any
  `NaN`/`Inf` / `ValueError` behaviour
- `plans/ROBUSTNESS.md` records each finding (ID, severity, surface, status,
  mitigation, regression test); update it when behaviour changes
- `docs/` describes the calling convention (array input; not clamped) so callers
  know masking out-of-range points is their responsibility

## Process

1. Scan all of `thermofeel/` for raise sites, bare excepts, and error-swallowing
   patterns; and for `log`/`sqrt`/`power`/division domains that can go `NaN`/`Inf`
2. For each finding: classify (wrong/missing exception, undocumented `NaN`/`Inf`,
   poor message, masks real bad input)
3. Fix each issue — raise the right error with a clear message, or document and
   test the `NaN`/`Inf` outcome; never add a guard that hides genuine
   out-of-domain input
4. Run `make all` to verify no regressions
5. Update `plans/ROBUSTNESS.md` and the affected docstrings/`docs/`
6. Summarize all changes made

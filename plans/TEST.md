# Test Plan

Repo: ecmwf/thermofeel

> This document describes the *shape* of the test suite — what is tested, where,
> and why. It deliberately avoids frozen counts (they drift with every release).
> Run `make test` (or `pytest tests/`) for the current source of truth.

## Coverage shape

| Component | Where tests live | Shape |
|-----------|------------------|-------|
| Index + supporting functions (array) | `tests/test_thermofeel.py` | Regression against stored `tests/*.csv` reference outputs, driven by `tests/thermofeel_testcases.csv`; uses `numpy.testing.assert_array_almost_equal` |
| Index + supporting functions (scalar) | `tests/test_scalars.py` | Pointwise `pytest.approx` checks against individually verified reference values, including documented edge cases |
| Unit converters (`helpers.py`) | exercised indirectly via the above + `celsius_to_kelvin`-wrapped scalar inputs | Round-trip correctness through every index path |

The two suites are complementary: the array suite guards *numerical stability
across a realistic spread of inputs*; the scalar suite guards *individual values
against the literature / known answers* and is where edge-case reasoning lives.

## Reference-data contract

`tests/*.csv` are committed reference outputs (`rh.csv`, `es.csv`, `utci.csv`,
`wbgt.csv`, …). They are the regression baseline.

- A test failure that moves these numbers means one of two things:
  1. **A bug** — fix the code, not the reference.
  2. **A deliberate formula correction** — regenerate the affected CSV (each test
     has a commented `np.savetxt(...)` line that produced it), record the change
     and its scientific justification in `../ChangeLog.rst`, and cite the source.
- Never weaken the `decimal=` tolerance to make a failing assertion pass without
  understanding why it drifted.

## Key interactions to verify

- Every index reproduces its reference formula within tolerance for the standard
  test-case spread.
- Unit boundaries: inputs in SI (Kelvin, m/s, %, hPa) → outputs in SI. Converter
  round-trips (`celsius_to_kelvin`/`kelvin_to_celsius`, the Fahrenheit pair) are
  exact.
- Composition paths: `bgt` → `wbgt`, `mrt` → `utci`, `bgt` → `mrt_from_bgt`,
  `saturation_vapour_pressure` → UTCI water-vapour-pressure path.
- `scale_windspeed` produces the right height adjustment used inside `bgt`
  (1.1 m), `mrt_from_bgt` (1.1 m), and `normal_effective_temperature` (1.2 m).
- Saturated vs. non-saturated vapour-pressure selection per index (a wrong choice
  is silently plausible — assert the actual value).
- UTCI accepts either `ehPa` or `td_k`, and raises `ValueError` when neither is
  supplied.

## Edge cases

- **Heat Index simplified**: result equals 2 m air temperature below the 20 °C
  gate; regression formula applies above it.
- **Heat Index adjusted**: the cold branch (≤ 40 °F), the low-humidity/high-temp
  and high-humidity/moderate-temp adjustment regions, and the warm-mask switch
  between the simple and regression formulas.
- **Wind Chill**: documented validity window (−50 °C…+5 °C, 5…80 km/h) — verify
  reference values from independent sources (e.g. the Environment Canada /
  Wikipedia worked examples) at the edges; the formula is *not* clamped, so the
  test documents that out-of-range inputs return out-of-range numbers.
- **Saturation vapour pressure multiphase**: liquid-water vs. ice branch
  selection by `phase`.
- **`approximate_dsrp`**: behaviour near low solar zenith angles (the documented
  large-error regime) and that `fdir` is not overwritten.
- **Negative / zero wind speed** into `bgt`/`wbgt` (numerical behaviour of the
  closed-form root).

## What is intentionally NOT covered here

- `thermofeel/experimental_wbgt.py` (Liljegren iterative solver) is experimental,
  unexported, and not part of the regression suite. If it is promoted to the
  public API, it must come with its own round-trip + reference-value tests and a
  citation-backed expected output before merging (see `TODO.md`).

## Running

```bash
make test            # via uv-managed .venv
# or directly:
pytest tests/ -v
```

CI runs the same suite through ECMWF `downstream-ci` on every push/PR, plus the
`python_qa` lint gate (`black` + `isort` + `flake8`).

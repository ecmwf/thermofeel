---
description: Measure coverage and fill gaps until the package is at >=95%
agent: build
---

# Improve Code Coverage

Perform a comprehensive code coverage analysis and fill gaps to reach 95%+ coverage.

## Process

### 1. Measure Current Coverage

```bash
make test       # runs pytest with --cov=thermofeel --cov-report=term-missing
```

Or directly, to iterate on one area:
```bash
.venv/bin/python -m pytest tests/ --cov=thermofeel --cov-report=term-missing
```

Coverage config lives in `pyproject.toml` (`[tool.coverage]`, branch coverage
on); CI uploads to Codecov (`codecov.yml`).

### 2. Identify Gaps

For each module below 95% coverage (`thermofeel/thermofeel.py`,
`thermofeel/liljegren.py`, `thermofeel/excess_heat.py`, `thermofeel/helpers.py`):
- Read the source to understand which code paths and branches are untested
- Categorize each gap as:
  - **(a) Needs new tests** — testable code paths with no coverage
  - **(b) Domain/edge path** — a `NaN`/`Inf` or out-of-domain branch (the
    closed-form `bgt` calm-air limit, Liljegren non-convergence, a `ValueError`
    contract); test it via `tests/test_robustness.py`
  - **(c) Dead code** — unreachable code that should be removed

### 3. Prioritize by Impact

Focus on the modules with the largest absolute gaps first — typically the
index functions in `thermofeel.py` and the Liljegren solvers in `liljegren.py`.

### 4. Write Tests

- Add tests to the EXISTING suites, matching their patterns and style:
  - `tests/test_thermofeel.py` — array-level regression against the stored
    `tests/*.csv` reference outputs
  - `tests/test_scalars.py` — pointwise checks against individually verified
    reference values (this is where edge-case reasoning lives)
  - `tests/test_robustness.py` — `NaN`/`Inf` propagation and error contracts
  - `tests/test_excess_heat.py` — the excess-heat submodule
- Each test should target a specific uncovered code path or branch
- Test error paths and domain edges, not just the happy path
- Validate any new expected value against the literature / a worked example and
  record its provenance in the test — never just pin whatever the code emits

### 5. Remove Dead Code

If coverage analysis reveals unreachable code:
- Verify it's truly unreachable (not just untested)
- **Remove it** rather than adding `# pragma: no cover`. Per AGENTS.md,
  suppression to dodge coverage is prohibited; the only accepted pragma is on the
  `unittest.main()` entrypoints.

### 6. Verify

- Run `make all` to confirm all tests pass and lint is clean
- Re-measure coverage to confirm improvement
- Report before/after comparison by module

## Target

Aim for at least **95% line + branch coverage** across `thermofeel/` (the suite
has historically reached 100%). Acceptable exceptions are limited to the
`if __name__ == "__main__": unittest.main()` entrypoints.

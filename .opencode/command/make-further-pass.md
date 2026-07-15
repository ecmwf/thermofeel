---
description: "Run a strictness-graded quality review pass — usage: /make-further-pass <pass-number>"
agent: build
---

# Further Pass Review

Perform a quality review pass over the codebase or the specified area.

**This is pass number $ARGUMENTS (default: 2 if not specified).**

Track the pass number and increase strictness with each successive pass.
thermofeel is a small, **published**, pure-Python (numpy-only) numerical
library; reviews focus on scientific correctness, the SI-unit contract,
numerical robustness, and keeping the duplicated sources of truth in step — not
on systems-level concerns.

## Pass 1-2 (Foundation)
- Simplification opportunities — reduce complexity, remove redundancy
- Naming quality — variables, functions, parameters, modules
- Comments and docstring quality — accurate, helpful, not redundant; each
  public function still names its scientific reference (author, year + DOI/URL)
- Identify duplicate sources of truth introduced or touched by this change
  (the same value or generated content kept in two places without a guard)
- Running the required formatter/lint/tests (`make fmt`, then `make all`)

## Pass 3-4 (Hardening)
Everything from Pass 1-2, PLUS:
- Scan for edge-cases and logical regressions (see `/improve-edge-cases`)
- Vectorisation: no Python-level loops over array elements; conditional logic
  uses `np.where` / boolean indexing, never per-element `if`
- SI-unit boundary: any formula that converts to °C/°F/km-h internally converts
  back so the public boundary stays SI (Kelvin, m/s, %, hPa); the documented
  exception (excess heat/cold factors) is intentional
- All documentation up-to-date with changes (`docs/`, docstrings, `README.md`)
- `NaN`/`Inf` behaviour at domain edges is documented and intentional, not
  incidental (cross-check `plans/ROBUSTNESS.md`)
- Public API surface is minimal and clean; no dead code, no unused imports, no
  stale TODOs
- NO suppression smell: no `# noqa`, `# type: ignore`, `# fmt: off`,
  `@pytest.mark.skip`, or `# pragma: no cover` added to dodge a real problem
  (the existing `unittest.main()` pragma is the one accepted case)

- Duplicate source-of-truth drift scan. Each hit must be inspected; the greps
  are not auto-fail rules — most matches are legitimate, but every match is a
  candidate for a value that has drifted from its canonical source.

  ```bash
  # The package version lives ONLY in thermofeel/__init__.py (__version__);
  # pyproject.toml reads it dynamically. Any other literal version is suspect.
  git grep -nE '[0-9]+\.[0-9]+\.[0-9]+' -- pyproject.toml thermofeel/ \
    | grep -viE 'python_requires|requires-python|target-version|py3[0-9]|>=|CY[0-9]'
  make version            # must match thermofeel/__init__.py
  make check              # imports the package and resolves the version

  # The list of computed indices must stay in step across these three places:
  git grep -nE 'Excess Heat|Wind Chill|Heat Force|Liljegren|Universal Thermal' \
    -- README.md thermofeel/thermofeel.py docs/guide/overview.md
  ```

  For the index list: `README.md`, the `thermofeel/thermofeel.py` module
  docstring, and `docs/guide/overview.md` must agree (AGENTS.md "Derived
  artefacts have one source of truth"). For the version: only
  `thermofeel/__init__.py` carries it.

- Reference-data freshness. The committed `tests/*.csv` are the regression
  baseline. If this change moved any of them, confirm it is a deliberate formula
  correction (regenerate via the test's commented `np.savetxt(...)` line, record
  it in `CHANGELOG.md` with the scientific justification and citation) — never a
  silently re-baselined bug. Never weaken a `decimal=` tolerance to pass.

## Pass 5+ (Polish)
Everything from Pass 3-4, PLUS with ZERO TOLERANCE:
- Every public function has a complete docstring with units and a citation
- Every error path and `ValueError` contract is tested
- Every documented edge case has a test (scalar AND array suites where it makes
  sense — `tests/test_scalars.py`, `tests/test_thermofeel.py`,
  `tests/test_robustness.py`)
- Code reads like well-written prose — a newcomer could follow each formula back
  to its reference
- Performance: no unnecessary array copies or temporaries on the global-grid
  hot paths; dtype is preserved (no silent integer truncation)
- Consistency: similar patterns handled the same way everywhere (wind scaling
  via `scale_windspeed`, unit conversion via `helpers.py`)
- Fresh-environment check for tooling/packaging changes. If this change touches
  `pyproject.toml`, the `Makefile`, `requirements.txt`, or the wheel contents,
  rebuild from a clean env (`make clean && make venv && make all`, and
  `make build` to inspect `dist/`) — packaging regressions only show up on a
  fresh build.

## Process

1. State which pass number this is and what strictness level applies
2. Scan the codebase (or specified files/area)
3. List all findings grouped by category
4. Fix each finding, verifying the fix passes tests
5. Run the formatter and linter: `make fmt` then `make lint`
6. Run the full gate: `make all` (and `make docs` if docs changed)
7. Summarize what was changed and what the next pass should focus on

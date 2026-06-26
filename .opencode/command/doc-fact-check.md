---
description: Run docs code examples and verify prose claims against the codebase
agent: build
---

# Documentation Fact-Check

Run code examples in the docs and verify prose claims against the actual codebase.

## Arguments

`$ARGUMENTS` — optional file paths relative to `docs/` to narrow the scope.

If no arguments given, check all `.md` files under `docs/`. Accepts individual
files (`guide/utci.md`) or directories (`guide/`).

Examples: `/doc-fact-check guide/utci.md`, `/doc-fact-check guide/`, `/doc-fact-check`

## Setup

```bash
make venv                                   # uv-managed .venv with the package + test deps
.venv/bin/python -c "import thermofeel; print('OK')"
```

If any setup step fails, report it and continue with what works. Blocks that
fail due to environment issues (missing tool, unbuilt dependency) are **not**
doc bugs — report them separately as setup problems.

Note: `docs/api.md` is generated from the source docstrings by **mkdocstrings**,
so signatures there cannot drift from the code — but the hand-written prose and
examples in `docs/index.md`, `docs/examples.md`, and `docs/guide/*.md` can.

## 1. Run Code Examples

For each fenced code block in the selected docs:

### Classify each block
- **Runnable** — self-contained Python: has its own imports, no `...` ellipsis
  placeholders, does not need user-specific data files
- **Continuation** — follows a runnable block on the same page and reuses its
  variables; concatenate with its predecessor
- **Skip** — partial snippets, signatures-only, output samples (lines starting
  with `>>>` output or `$`), unlabeled code blocks (no language tag),
  install/setup commands (`pip install`, `make ...`), or blocks needing
  user-specific files (GRIB/NetCDF inputs, `*.grib`, `*.nc`)

### Execute runnable blocks

**Python:** write to a temp `.py` file, run `.venv/bin/python <file>`. 30s
timeout per block. Most examples import `numpy` and `thermofeel` and pass NumPy
arrays — remember several functions require array (not scalar) input.

**Bash:** only commands that can run without user-specific data. Run and check
exit code 0.

**Never** invent missing imports or variables. If a block isn't self-contained, skip it.

Record pass/fail/skip for every block.

## 2. Check Claims Against Code

Read each doc file and find verifiable factual claims — these are inline code
spans (`` `calculate_utci` ``), numbers, units, citations, and table cells, not
prose paragraphs. For each one, check the source.

**What to check:**
- **API names / signatures** — function names, parameter names, defaults
  mentioned in prose → verify they exist and match in `thermofeel/thermofeel.py`,
  `thermofeel/helpers.py`, `thermofeel/excess_heat.py` (e.g. `threshold=0.1`,
  `wind_scaling="liljegren"`, `clip=False`)
- **Units** — every documented input/output unit must match the SI contract the
  function actually implements (Kelvin, m/s at 10 m, %, hPa); flag any place a
  guide states the wrong unit
- **Scientific citations** — author/year and DOI/URL in a guide page must match
  the reference in the corresponding function docstring
- **Validity ranges** — documented windows (e.g. Wind Chill −50…+5 °C and
  5…80 km/h, the Heat-Index 20 °C gate, UTCI polynomial ranges) must match the
  docstring and the code's branch conditions
- **The index list** — the computed-indices list must agree across `README.md`,
  the `thermofeel/thermofeel.py` module docstring, and `docs/guide/overview.md`
  (the "keep in step" rule in AGENTS.md)
- **Numerical claims** — any worked-example value quoted in the docs → recompute
  it and confirm

**What NOT to check:** subjective claims ("accurate", "physically based"),
high-level explanatory prose that makes no testable assertion.

## 3. Report

For each finding:
```
[ERROR|STALE|DRIFT] file.md:LINE — Summary
  Docs say: <what the docs claim>
  Code says: <what the source actually shows>
```

- **ERROR** — code example fails, API doesn't exist, wrong signature/unit
- **STALE** — outdated value, renamed API, removed function, wrong citation
- **DRIFT** — minor mismatch (index list out of step, slightly different default)

End with a summary: total files checked, blocks run/passed/failed/skipped,
number of claim issues found.

For ERROR findings: propose a fix but do **not** apply it — present for user review.

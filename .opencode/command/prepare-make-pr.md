---
description: Run pre-flight checks then commit, push and open the PR
agent: build
---

# Prepare and Create Pull Request

Final preparation check and PR creation workflow.

## Pre-flight Checks

Run the full gate. If it fails, STOP and report the failure — do not open the PR.

```bash
make all          # check + test + lint (pure Python: import-check, pytest, ruff)
```

`make all` is the single source of truth for the pre-commit gate; it wraps the
underlying `pytest` and `ruff` commands (run `make help` for the individual
targets). Do NOT run a narrower subset in its place. Run `make fmt` first to
auto-apply lint fixes and formatting.

If the change touches `docs/` or any docstring, also build the docs:

```bash
make docs         # mkdocs build --strict
```

## PR Creation

If all checks pass:

1. Review what files to include — stage only source, test, doc, and config files.
2. Do NOT stage: `.venv/`, `dist/`, `build/`, `site/`, `*.egg-info/`,
   `__pycache__/`, `.coverage`, `coverage.xml`, or any local data
   (`*.grib`, `*.npz`). `.github/` and `.opencode/` ARE tracked.
3. If not already on a feature branch, create one off `main` using the
   convention `<type>/<kebab-summary>` (`feat`, `fix`, `docs`, `chore`,
   `refactor`, `test`, `ci`, `perf`, `build`).
4. Commit with a clear message matching the repo style.
5. Record every user-facing change in `CHANGELOG.md` under the current
   release section before opening the PR.
6. Push to origin.
7. Create a pull request against `ecmwf/thermofeel` `main` with:
   - Title: concise summary of the change
   - Body: what changed and why; for any formula change, the reference and
     rationale; how it was verified (tests/examples)
   - Link to any related issues
8. Report the PR URL.

Only commit, push, or open a PR when the user has asked you to.

---
description: "Cut a new release — usage: /make-release <X.Y.Z>"
agent: build
---

# Make Release

Create a new release of thermofeel (a pure-Python package published to PyPI).

## Arguments

Provide the version as the first argument, e.g. `/make-release 2.3.0`.
If no version is given, default to a MINOR bump of the current
`thermofeel/__init__.py` `__version__`.

> **Canonical tooling:** `make release` drives `scripts/release.py`; the version
> contract is in `AGENTS.md` "Version control". This command automates that
> procedure — keep the two in sync.

`__version__` in `thermofeel/__init__.py` is the **single source of truth**;
`pyproject.toml` reads it dynamically. Tags are **bare** `MAJOR.MINOR.MICRO`
(NO leading `v`). thermofeel is a **published** library — follow SemVer: MINOR
for new features, MICRO for fixes/docs, and **NEVER bump MAJOR unless the user
explicitly says so**.

## Pre-release Checks

Run these in order. If ANY step fails, STOP and prompt the user.

### 1. Clean, pushed working tree
```bash
git status                      # must be clean — all changes committed
git fetch origin
git rev-parse HEAD origin/main  # HEAD must equal origin/main (pushed, in sync)
```
`scripts/release.py` blocks unless you are on `main`, the tree is clean, and
`HEAD` matches `origin/main` — releases are cut from a pushed `main`.

### 2. CHANGELOG entry (curated by hand)
The `## X.Y.Z` section must already exist in `CHANGELOG.md` before releasing —
it is never auto-generated. Add it by hand, matching the existing format (a
`## X.Y.Z` heading followed by flat bullets; this repo does NOT use
`[Unreleased]` or Keep-a-Changelog Added/Changed/Fixed sections). Summarise the
user-facing changes since the last tag (`git log <last-tag>..HEAD --oneline`),
and for any formula change cite the reference and rationale.

### 3. Build and inspect the artefacts
```bash
make build      # builds sdist + wheel into dist/ and runs `twine check`
```
Confirm the wheel ships only the `thermofeel` package (no `tests/`, `examples/`,
`scripts/`).

### 4. Full gate
```bash
make all        # check + test + lint (the local gate)
```
`make release ... CONFIRM=1` re-runs `make all` itself before tagging (unless
`SKIP_GATE=1`), but run it yourself first so failures surface early.

## Release Process

### 1. Plan the release (no changes made)
```bash
make release X.Y.Z
```
This validates everything (SemVer relation, tag not already present, CHANGELOG
section, clean+synced `main`) and prints the plan — it makes NO changes.

### 2. Apply: bump, commit, tag
```bash
make release X.Y.Z CONFIRM=1
```
With `CONFIRM=1` the tooling runs the gate, bumps `__version__` forward to
`X.Y.Z` if needed (NEVER backward), commits `Release X.Y.Z`, and creates the
annotated bare tag. It **NEVER pushes** — it prints the exact `git push`
commands. (A MAJOR bump additionally requires `ALLOW_MAJOR=1`, which you only set
when the user has explicitly approved a MAJOR release.)

### 3. Publish (push the tag)
Pushing the tag triggers `.github/workflows/cd.yml`, which builds and publishes
to PyPI. **Only push with explicit user approval.**
```bash
git push origin main        # if the version-bump commit was made
git push origin X.Y.Z       # pushing the bare tag triggers the cd workflow -> PyPI
```

### 4. Verify on PyPI
Wait for the `cd` workflow to finish and for PyPI to propagate (~1-2 min), then
smoke-test in a throwaway environment:
```bash
pip install thermofeel==X.Y.Z
python -c "import thermofeel; print(thermofeel.__version__)"
```

**IMPORTANT:** Version tags are bare semver (e.g. `2.3.0`), NEVER prefixed with `v`.
**IMPORTANT:** NEVER bump MAJOR unless the user explicitly says so.
**IMPORTANT:** NEVER push to remote without explicit user approval.

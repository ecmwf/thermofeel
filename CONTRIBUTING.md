# Contributing to *thermofeel*

Please report [bug reports](https://github.com/ecmwf-projects/thermofeel/issues)
or [pull requests](https://github.com/ecmwf-projects/thermofeel/issues) on
[GitHub](https://github.com/ecmwf-projects/thermofeel).

We want your feedback, please e-mail: <user-services@ecmwf.int>

The package is installed from PyPI with:

```bash
pip install thermofeel
```

How you could use thermofeel:

- Test out the code in this library with data from different models and in
  different structures
- Compare the methods with other outputs
- Plot outputs, examples are provided in the examples directory
- Go further than our examples and let us know how it goes

## Development setup

This project is pure Python (runtime dependency: `numpy`) and uses a top-level
`Makefile` as the single entry point. Run `make help` to see every target.

```bash
make venv     # create a local .venv (uv) with the package + test deps
make all      # the full gate: check + test + lint
make fmt      # apply ruff lint fixes + formatting
make docs     # build the MkDocs documentation site
```

## Releasing

Releases are published to PyPI by the `cd` workflow when a bare
`MAJOR.MINOR.MICRO` tag (no leading `v`) is pushed. The version is single-sourced
in `thermofeel/__init__.py` (`__version__`); `pyproject.toml` reads it
dynamically.

```bash
make build              # build + metadata-check the sdist/wheel in dist/ (no upload)
make release X.Y.Z      # plan only: validate everything and show what would happen
make release X.Y.Z CONFIRM=1   # apply: bump __version__ forward if needed, commit, tag
```

`make release` only ever moves the version **forward** (never backward), requires
the `## X.Y.Z` `CHANGELOG.md` section to exist, must run on a clean, pushed `main`,
and **never pushes** — it prints the `git push` commands for you to run. A MAJOR
bump additionally requires `ALLOW_MAJOR=1`.

# (C) Copyright 2024- ECMWF and individual contributors.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# Top-level Makefile — the single entry point for building, testing, and
# linting thermofeel. Run `make help` to see every target.
#
# thermofeel is pure Python (runtime dependency: numpy). The environment is
# managed by `uv`; targets wrap the raw tool commands so there is one source of
# truth (see CONTRIBUTING.md). The QA tool is ruff (linter + formatter, the
# single replacement for flake8 + isort + black); its config lives in
# pyproject.toml under [tool.ruff].

.PHONY: help all check test lint fmt docs docs-serve version build release clean venv

# ── Tooling ─────────────────────────────────────────────────────────────────

UV      ?= uv
VENV    ?= .venv
PYTHON  ?= $(VENV)/bin/python

# ruff (linter + formatter) runs in an ephemeral uv environment (`--no-project`
# so uv does not build thermofeel just to run a linter). Config lives in
# pyproject.toml under [tool.ruff].
RUFF    ?= $(UV) run --no-project --with ruff ruff

# Code that the QA tools operate on (the package, its tests, and helper scripts).
PY_SRC  ?= thermofeel tests scripts

# ── Defaults ──────────────────────────────────────────────────────────────

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}'

all: check test lint ## Run all checks, tests, and lints (the gate; run before committing)

# ── Environment ─────────────────────────────────────────────────────────────

venv: ## Create the local .venv (uv) and install the package + test deps
	@test -d $(VENV) || $(UV) venv $(VENV)
	$(UV) pip install --python $(PYTHON) -e ".[test]"

# ── Check ─────────────────────────────────────────────────────────────────

check: venv ## Import-check the package and verify the version resolves
	$(PYTHON) -c "import thermofeel; print('thermofeel', thermofeel.__version__, 'imports OK')"

# ── Test ──────────────────────────────────────────────────────────────────

test: venv ## Run the pytest suite with coverage (array + scalar tests)
	$(PYTHON) -m pytest tests/ -v --cov=thermofeel --cov-report=term-missing

# ── Lint / format ───────────────────────────────────────────────────────────

lint: ## Lint and check formatting with ruff (check-only)
	$(RUFF) check $(PY_SRC)
	$(RUFF) format --check --diff $(PY_SRC)

fmt: ## Apply ruff lint fixes (incl. import sorting) and formatting in place
	$(RUFF) check --fix $(PY_SRC)
	$(RUFF) format $(PY_SRC)

# ── Docs ────────────────────────────────────────────────────────────────────
# MkDocs (Material) + mkdocstrings. mkdocstrings introspects the installed
# package, so these run in the project env (thermofeel importable) with the doc
# tooling added on top from docs/requirements.txt.

docs: ## Build the documentation site (MkDocs) into site/
	$(UV) run --with-requirements docs/requirements.txt mkdocs build --strict

docs-serve: ## Serve the documentation locally with live reload
	$(UV) run --with-requirements docs/requirements.txt mkdocs serve

# ── Versioning ────────────────────────────────────────────────────────────
# Single source of truth: thermofeel/__init__.py __version__ (pyproject reads it
# dynamically). To release, bump it there, add a CHANGELOG.md entry, then tag
# the bare MAJOR.MINOR.MICRO (no leading 'v').

version: ## Print the canonical version from thermofeel/__init__.py
	@$(UV) run --no-project python -c "import re,pathlib; print(re.search(r'__version__\s*=\s*[\"\047]([^\"\047]+)', pathlib.Path('thermofeel/__init__.py').read_text()).group(1))"

# ── Build / Release ───────────────────────────────────────────────────────
# `make build` produces the dist/ artefacts locally and checks their metadata
# (no upload). `make release X.Y.Z` validates the release and, with CONFIRM=1,
# bumps __version__ forward to X.Y.Z (never backward), commits, and creates the
# bare tag — it NEVER pushes. Publishing happens in CI (.github/workflows/cd.yml)
# when the tag is pushed. See scripts/release.py and AGENTS.md for the protocol.

build: ## Build the sdist + wheel into dist/ and check metadata (no upload)
	rm -rf dist/
	$(UV) build
	$(UV) run --no-project --with twine twine check dist/*

release: ## Validate + tag a release: make release X.Y.Z [CONFIRM=1] (never pushes)
	@CONFIRM="$(CONFIRM)" ALLOW_MAJOR="$(ALLOW_MAJOR)" ALLOW_BRANCH="$(ALLOW_BRANCH)" SKIP_GATE="$(SKIP_GATE)" \
		$(UV) run --no-project python scripts/release.py "$(or $(VERSION),$(filter-out release,$(MAKECMDGOALS)))"

# Allow a bare version goal (`make release 2.3.0`) to be passed through to the
# release target instead of being treated as its own target. Non-version typos
# still error rather than silently no-op.
%:
	@case "$@" in \
	  [0-9]*.[0-9]*.[0-9]*) : ;; \
	  *) printf "make: *** No rule to make target '%s'.  Stop.\n" "$@" >&2; exit 2 ;; \
	esac

# ── Cleanup ───────────────────────────────────────────────────────────────

clean: ## Remove the venv, build artefacts, and caches
	rm -rf $(VENV)
	rm -rf build/ dist/ *.egg-info thermofeel.egg-info
	rm -rf site/
	rm -rf .pytest_cache
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

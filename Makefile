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
# truth (see CONTRIBUTING.rst). The QA tools mirror what CI's python_qa step
# enforces: black + isort (profile=black) + flake8 (config in tox.ini).

.PHONY: help all check test lint fmt docs version clean \
        venv black-check isort-check flake8

# ── Tooling ─────────────────────────────────────────────────────────────────

UV      ?= uv
VENV    ?= .venv
PYTHON  ?= $(VENV)/bin/python

# Formatting / lint tools run in ephemeral uv environments (`--no-project` so uv
# does not build thermofeel just to run a linter). Pinned config lives in tox.ini.
BLACK   ?= $(UV) run --no-project --with black black
ISORT   ?= $(UV) run --no-project --with isort isort
FLAKE8  ?= $(UV) run --no-project --with flake8 flake8

# Code that the QA tools operate on (the package + its tests).
PY_SRC  ?= thermofeel tests

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

test: venv ## Run the pytest suite (array regression + scalar pointwise tests)
	$(PYTHON) -m pytest tests/ -v

# ── Lint / format ───────────────────────────────────────────────────────────

flake8: ## Run flake8 (config in tox.ini)
	$(FLAKE8) $(PY_SRC)

isort-check: ## Check import ordering (isort, profile=black)
	$(ISORT) --check-only --diff $(PY_SRC)

black-check: ## Check code formatting (black)
	$(BLACK) --check --diff $(PY_SRC)

lint: flake8 isort-check black-check ## Run all lints (flake8 + isort + black, check-only)

fmt: ## Apply isort + black formatting in place
	$(ISORT) $(PY_SRC)
	$(BLACK) $(PY_SRC)

# ── Docs ────────────────────────────────────────────────────────────────────

docs: ## Build the Sphinx documentation into docs/build/html
	$(UV) run --no-project --with-requirements docs/requirements.txt \
		sphinx-build -b html docs/source docs/build/html

# ── Versioning ────────────────────────────────────────────────────────────
# Single source of truth: thermofeel/__init__.py __version__ (pyproject reads it
# dynamically). To release, bump it there, add a ChangeLog.rst entry, then tag
# the bare MAJOR.MINOR.MICRO (no leading 'v').

version: ## Print the canonical version from thermofeel/__init__.py
	@$(UV) run --no-project python -c "import re,pathlib; print(re.search(r'__version__\s*=\s*[\"\047]([^\"\047]+)', pathlib.Path('thermofeel/__init__.py').read_text()).group(1))"

# ── Cleanup ───────────────────────────────────────────────────────────────

clean: ## Remove the venv, build artefacts, and caches
	rm -rf $(VENV)
	rm -rf build/ dist/ *.egg-info thermofeel.egg-info
	rm -rf docs/build/
	rm -rf .pytest_cache
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

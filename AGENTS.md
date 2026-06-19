# Claude and Other Agents

This file (`AGENTS.md`) is the canonical agent-instructions document.
`CLAUDE.md` is a symlink to it for cross-tool compatibility — edit `AGENTS.md`
and refer to it by that name.

# Guidelines

- If an LSP / symbol-navigation tool is available, prefer it over Grep/Read for
  code navigation — use it to find definitions, references, and workspace symbols.
  Fall back to Grep/Read when no such tool is exposed.

- CRITICAL: NEVER suppress warnings, lint errors, or test failures with
  annotations (`# noqa`, `# type: ignore`, `# fmt: off`, `@pytest.mark.skip`,
  `pytestmark = pytest.mark.xfail`, etc.) unless the suppression itself is the
  correct semantic choice. If a lint fires, fix the underlying code. If a test
  fails, fix the bug. If a warning appears, resolve the root cause. Quick
  workarounds that hide problems are strictly prohibited.

- CRITICAL: Always prefer proper solutions over quick fixes.
  When facing a problem, invest the effort to understand the root cause and fix
  it correctly rather than applying a workaround. Specifically:
  - Do NOT skip Python versions or configurations to avoid fixing a failure.
  - Do NOT weaken numerical tolerances in tests to make them pass — if a result
    drifts, understand *why* before touching the reference value.
  - Do NOT remove or weaken tests to make CI pass.
  - Do NOT add TODO/FIXME/HACK comments as a substitute for doing the work now.
  - If a proper fix requires changing multiple modules, do it.
  - If a proper fix requires understanding unfamiliar code, read it first.
  - If you are unsure whether a fix is proper, ask before proceeding.

- IMPORTANT: thermofeel is a PUBLISHED library with downstream users (PyPI,
  used inside NWP pipelines). Unlike a pre-release project, API and behaviour
  changes have real consequences:
  - Treat the public function signatures in `thermofeel/thermofeel.py` and the
    converters in `thermofeel/helpers.py` as a contract.
  - Follow Semantic Versioning (see *Version control*). A breaking API change is
    a MAJOR bump and must be called out in `CHANGELOG.md` with a migration note
    (as was done for the 1.x → 2.0 standardisation).
  - When a thermal-index formula changes, the change is observable to every user
    — document the reference and the rationale in the docstring and `CHANGELOG.md`.

- IMPORTANT: Scientific correctness is the whole point of this library.
  - Every index function cites a peer-reviewed reference in its docstring. When
    you add or change a formula, cite the source and keep units explicit (SI:
    Kelvin, m/s, %, hPa). See `plans/DESIGN.md` for the unit conventions.
  - Validate new/changed formulas against published worked examples or an
    independent reference implementation, not just against a self-generated
    number. Record the expected value and its provenance in the test.

- IMPORTANT: when planning and before you do any work:
  - ALWAYS mention how you would verify and validate that work is correct
  - include TDD tests in your plan (this repo has both array-based regression
    tests and pointwise scalar tests — add to both where it makes sense)
  - take a behaviour-driven approach
  - you are very much ENCOURAGED to ask questions to get the design correct
  - ALWAYS seek clarifications to sort out ambiguities
  - ALWAYS provide a summary of the Design and implementation Plan

- IMPORTANT: when you build code and new features:
  - ALWAYS document those features in `docs/` (MkDocs)
  - Remember to add examples (see below)

- IMPORTANT:
  - when you commit your work, make sure it passes all checks, tests and lints —
    by running `make all`

- IMPORTANT: Derived artefacts have one source of truth.
  If a value is conceptually one piece of data but appears in many places, pick
  one canonical location and either generate the others from it or add a
  consistency check. Concretely for thermofeel:
  - The package version lives ONLY in `thermofeel/__init__.py` (`__version__`).
    `pyproject.toml` reads it dynamically (`version = { attr = ... }`) — never
    hard-code a version into `pyproject.toml`.
  - The list of computed indices appears in `README.md`, the
    `thermofeel/thermofeel.py` module docstring, and `docs/`. Keep them in step.

- IMPORTANT: NEVER use process-ephemeral references in code, comments,
  docstrings, commit messages, or planning documents. The code outlives the
  workflow that produced it; references to that workflow age into noise.

  Banned vocabulary (each becomes meaningless once the PR merges):
  - Workflow ordinals: "Phase 2", "Pass 5", "Round 1", "Step 3"
  - Issue-tracker references: "PR-1", "sub-task 4", "issue #94"
  - Review-feedback severity buckets: "Critical #1", "High #2"
  - History phrases inside code: "before this PR", "after review feedback",
    "fixed in commit abc1234"

  Replace each with a name that describes what the thing **is** or **does**.
  Git history records chronology; comment and commit text must record substance.

# Design & Purpose

- `README.md` — short entry-level overview (also the GitHub landing page)
- `plans/MOTIVATION.md` — why thermofeel exists and what we're building (long form)
- `plans/DESIGN.md` — design rationale, the SI-unit I/O contract, key decisions
- `plans/ARCHITECTURE.md` — package/module map and the data-flow between indices
- `plans/TODO.md` — features and fixes decided to implement (accepted backlog)
- `plans/IDEAS.md` — speculative ideas not yet decided
- `plans/TEST.md` — test plan and coverage shape
- `plans/ROBUSTNESS.md` — numerical failure modes and the input-robustness
  checklist (the risk surface for a pure numerical library)
- `CONTRIBUTING.md` — contributor setup and workflow
- `CHANGELOG.md` — release history (the single backward-looking record)

Follow `plans/DESIGN.md` principles in all code.

# Build / lint / test (required before marking done)

This project is pure Python (depends only on `numpy`). The top-level `Makefile`
is the single entry point — run `make help` to list every target. The full gate
before marking work done is:

```bash
make all          # = make check test lint  (the gate; run before committing)
```

Useful focused targets (run `make help` for the complete list):

| Target | What it does |
|--------|--------------|
| `make venv` | Create the local `.venv` via `uv` and install the package + test deps |
| `make check` | Import-check the package and verify the version is resolvable |
| `make test` | Run the pytest suite (array regression + scalar pointwise tests) |
| `make lint` | `ruff check` + `ruff format --check` |
| `make fmt` | Apply `ruff check --fix` + `ruff format` in place |
| `make docs` | Build the MkDocs documentation site (`mkdocs build --strict`) |
| `make version` | Print the canonical version from `thermofeel/__init__.py` |
| `make clean` | Remove `.venv`, build artefacts, caches |

- The toolchain is `ruff` — a single linter + formatter replacing `flake8` +
  `isort` + `black`. Its config lives in `pyproject.toml` under `[tool.ruff]`.
  CI runs the same `ruff check` + `ruff format --check` in the `qa` job. Run
  `make fmt` before `make lint` to auto-fix lint and formatting.
- The Makefile uses `uv` to manage the virtual environment. The targets wrap the
  raw tool commands so there is one source of truth — see the *Building* and
  *Testing* sections of `CONTRIBUTING.md`.

# Version control

- Git project at github.com/ecmwf/thermofeel; single `main` branch.
- IMPORTANT:
  - versions are tagged using Semantic Versioning 'MAJOR.MINOR.MICRO'
  - NEVER bump MAJOR unless the user says so (MAJOR = breaking API/behaviour change)
  - Increment MINOR for new features, MICRO for bugfixes and documentation updates
- NEVER prepend git tags or releases with 'v' (the CD workflow triggers on tags
  matching `[0-9]+.[0-9]+.[0-9]+`).
- NOTE: SINGLE SOURCE OF TRUTH FOR VERSION — `thermofeel/__init__.py`
  `__version__` is canonical. `pyproject.toml` derives it dynamically. On a
  release: update `__version__`, add the `CHANGELOG.md` entry, commit, then tag.
- REMEMBER on releases:
  - check all is committed and pushed upstream, otherwise STOP and warn the user
  - bump `__version__` in `thermofeel/__init__.py`
  - add the release section to `CHANGELOG.md`
  - git tag with the bare version (no `v`) and push; the `cd` workflow publishes
    to PyPI

# Tracking Work Done

Record every user-facing change in `CHANGELOG.md` as it merges — this is the
single backward-looking record. Keep entries concise and reference-bearing for
formula changes. Durable design decisions go in `plans/DESIGN.md`; behaviour and
usage go in `docs/`.

# Documentation

Create and maintain documentation under `docs/` (MkDocs + Material with
mkdocstrings, published on Read the Docs).
- Easy to follow by an average technical person, with well-separated topics
- One guide page per index under `docs/guide/`; the API reference (`docs/api.md`)
  is generated from the source docstrings via mkdocstrings, so keep docstrings
  accurate and well-formed
- Always state the input/output units and cite the scientific reference
- Add worked examples when a formula is hard to follow, especially edge cases
  (e.g. the validity range of Wind Chill, the temperature gating in Heat Index)

# Examples

Maintain runnable usage examples under `examples/`.
- Show the most common cases: computing each index from typical inputs
- Show how the indices compose (e.g. MRT → UTCI, BGT → WBGT)
- Keep examples in step with the current API — an example that calls a removed
  or renamed function is a bug, not documentation

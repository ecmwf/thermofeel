---
description: Onboard the agent to the thermofeel project by reading the docs and surveying the code
agent: plan
---

# Onboard to the Project

Build a working mental model of this repository so you can act competently
on it. Do NOT rely on any knowledge baked into this command — derive
everything from the repository's own files as they exist right now.

## Arguments

`$ARGUMENTS` — optional focus area (e.g. `utci`, `wbgt`, `liljegren`,
`excess-heat`, `robustness`, `docs`). If given, go deeper on that area during
the code survey (Step 3) while still doing the full doc pass. If empty, do a
balanced survey across the whole library.

## 1. Read the project-level docs

Read these first, in order. They are the source of truth for *intent*:

- The root markdown files: `README.md`, `AGENTS.md` (canonical agent
  instructions; `CLAUDE.md` is a symlink to it), `CONTRIBUTING.md`,
  `CHANGELOG.md` (read the most recent release section(s) for current
  direction — the changelog is hand-curated `## X.Y.Z` flat bullets, NOT
  Keep-a-Changelog and with no `[Unreleased]` section). Also skim
  `ACKNOWLEDGEMENTS.md`.
- Every `plans/*.md` file. Read each one fully. They cover motivation, the
  design rationale + SI-unit contract, the architecture / data-flow between
  indices, the test plan, the numerical-robustness analysis, the accepted-work
  backlog, and speculative ideas.

If a root markdown file or a `plans/*.md` file exists that is not named above,
read it too — discover them with a glob rather than assuming a fixed list.
Treat `plans/IDEAS.md` as speculative/unreviewed (it says so itself), not as
committed direction.

## 2. Map the repository structure

Without reading every file, establish the layout:

- The Python package `thermofeel/` and its modules (`thermofeel.py` — the
  public API; `liljegren.py` — internal WBGT solvers; `excess_heat.py`;
  `helpers.py` — SI unit converters; `__init__.py` — the public surface and the
  canonical `__version__`).
- Where tests, examples, docs, and release tooling live (`tests/`, `examples/`,
  `docs/`, `scripts/`).
- The build/test entry point: start from the `Makefile` — run nothing, just
  read the target graph (`make help`). Note that the package depends only on
  `numpy` at runtime and is managed with `uv`.

## 3. Survey the code (representative, not exhaustive)

Read enough source to understand *how the project achieves its purpose* — the
public surface and the core data flow — NOT every file. Aim for:

- The public `calculate_*` API in `thermofeel/thermofeel.py` and the SI-unit
  input/output contract every function honours (Kelvin, m/s at 10 m, %, hPa).
- How indices are composed from supporting quantities — follow a couple of the
  dependency chains in `plans/ARCHITECTURE.md` end to end (e.g.
  `saturation_vapour_pressure` → UTCI; `bgt` → `wbgt`; `mrt` → UTCI).
- The three WBGT implementations and how `calculate_wbgt_liljegren` delegates to
  the internal `thermofeel/liljegren.py` energy-balance solvers and 2 m wind
  profile.
- The unit-converter boundary in `thermofeel/helpers.py` and how each empirical
  formula converts in/out so the public boundary stays SI.
- The error / NaN-Inf model and the numerical invariants the code enforces
  (cross-check what `plans/ROBUSTNESS.md` claims against what the code does).

Prefer LSP / symbol navigation where available; fall back to Grep/Read.
When the docs and the code disagree, trust the code and note the drift.

## 4. Report back

Produce a concise onboarding summary (do not change any files):

1. **Purpose** — what the library is and the problem it solves, in 2-3 lines.
2. **Architecture** — the package modules and how they relate (a short list or
   a small diagram).
3. **Core data flow** — one index traced end to end, from SI inputs through any
   supporting quantities to the SI result.
4. **Conventions that constrain edits** — the non-negotiable rules an agent must
   respect here (pulled from `AGENTS.md` / `plans/DESIGN.md`): the SI in/out
   contract, the version single-source-of-truth (`thermofeel/__init__.py`),
   SemVer for a *published* library (never bump MAJOR unbidden), every formula
   cited, vectorised over NumPy (no Python loops), no suppressed lints/tests,
   the `make all` gate, the commit/PR workflow.
5. **How to build, test, and lint** — the canonical commands (the gate to run
   before marking work done): `make all`.
6. **Current direction & open work** — what's active vs. speculative, from the
   latest `CHANGELOG.md` section and `plans/TODO.md` (accepted) vs.
   `plans/IDEAS.md` (speculative).
7. **Anything stale or contradictory** you noticed between docs and code (e.g.
   the index list must stay in step across `README.md`, the
   `thermofeel/thermofeel.py` module docstring, and `docs/`).

Keep it tight and skimmable. The goal is a shared mental model, not a copy of
the docs.

# Programme: New Thermal Indices (target 2.4.0)

Batch delivery of new indices promoted from the IDEAS assessment. Sprint 0
(external research: provenance + validation) is **done**; the per-index findings
are captured below. Implementation is **fanned out one agent per index** and
**overlaid in a single PR** (target `2.4.0`). Each index must pass its own
validation gates (§3) independently; **any index that cannot clear its gates is
dropped from the batch without blocking the others** and returned to `IDEAS.md`
with a written reason.

> This is a planning document (like the rest of `plans/`). It is the single
> source of truth for the programme. Behaviour/usage still lands in `docs/`;
> durable rationale in `plans/DESIGN.md`; the backlog pointer in `plans/TODO.md`.

## 1. Scope & decisions

- **Version:** one MINOR bump to `2.4.0` for the whole batch (new features;
  never MAJOR). A single `## 2.4.0` `CHANGELOG.md` section; `__version__` bumped
  once, by the orchestrator, at integration.
- **Single PR**, on `feat/new-indices` (branched off `feat/discomfort-index`,
  i.e. stacked on PR #54). GitHub retargets to `main` once #53/#54 merge.
- **SI contract (`DESIGN.md` §2) is non-negotiable:** inputs are Kelvin (`*_k`),
  `rh` in %, vapour pressure in hPa, wind in m/s; temperature-like outputs return
  Kelvin. Dimensionless indices (PMV, PPD, RSI) return dimensionless; PPD is a %.
- **Every formula cited** (author, year, DOI) in the docstring. Reference values
  in tests are **pinned literals** traceable to a published source or an
  independent MIT-licensed oracle; the oracle is **never** a runtime dependency
  (`TEST.md`). Core stays **numpy-only**.
- **No index ships on trust.** Each implementation agent must *re-derive/verify*
  the reference values against the cited formula before pinning them (Sprint 0
  already found one dossier oracle value that does not reconcile with its stated
  formula — see §5.1).

## 2. Indices in this batch

| # | Index | Function(s) | Inputs (SI) | Formula source (DOI) | Feasibility | Gate status |
|---|-------|-------------|-------------|----------------------|-------------|-------------|
| 1 | Apparent Temperature (radiation) | `calculate_apparent_temperature_radiation` | `t2_k, va, rh, q` | Steadman 1994 `10.1071/es94001`; BoM operational page | Easy | **PASS** |
| 2 | Relative Strain Index | `calculate_relative_strain_index` | `t2_k, rh` | Lee & Henschel 1966 `10.1111/j.1749-6632.1966.tb43059.x`; form via Asghari 2020 `10.2174/1874213002013010011` | Easy | **PASS** (variant note) |
| 3 | Summer Simmer Index | `calculate_summer_simmer_index` | `t2_k, rh` | Pepi 1987 `10.1080/00431672.1987.9933356` | Easy | **CONDITIONAL** (secondary-source formula; analytic gate) |
| 4 | PMV / PPD | `calculate_pmv`, `calculate_ppd` | `t2_k, mrt_k, var, rh, met, clo` | ISO 7730:2005; Fanger 1970 | Medium (bounded iteration) | **PASS** |
| 5 | PET | `calculate_pet` | `t2_k, mrt_k, va, rh (+person)` | Höppe 1999 `10.1007/s004840050118`; Walther & Goestchel 2018 `10.1016/j.buildenv.2018.03.054`; VDI 3787-2 | Hard | **BLOCKED** — see §5.5 |

## 3. Validation gates (every index must pass all four)

- **G1 — Provenance.** DOI-backed citation; the exact formula/algorithm
  transcribed from a *fetched* source (not memory). Any constant whose source
  could not be read is flagged and resolved before merge.
- **G2 — Analytic.** Hand-derived identities / limits encoded as tests
  (e.g. SSI = affine transform of Thom THI; `PMV=0 ⇒ PPD=5%`; AT-radiation with
  `q=0` differs from the non-radiation form only by its own constants).
- **G3 — External reference.** ≥5 pinned `(inputs → output)` rows from a
  published table or an independent MIT oracle, in `tests/test_scalars.py`, each
  re-verified against the cited formula.
- **G4 — Regression + gate.** A 50-row array CSV over
  `tests/thermofeel_testcases.csv` (drift guard) in `tests/test_thermofeel.py`;
  `make all` green (100% coverage) and `make docs --strict` clean.

## 4. Delivery model & agent topology

**Integration branch:** `feat/new-indices`. **One implementation agent per
index**, each in its **own git worktree** on a sub-branch `impl/<index>`, each
with its own `uv` venv (`make venv`). Research is done; agents implement + test +
document their single index only.

### 4.1 Conflict-avoidance protocol (critical for parallelism)

Each agent edits **only**:
1. its function(s) in `thermofeel/thermofeel.py`, inserted at its **pre-assigned
   anchor** (§4.2) — additive, distinct regions;
2. its **own** test additions: one method in `tests/test_scalars.py`, one in
   `tests/test_thermofeel.py`, and one `setUp` load line; its reference CSV
   `tests/<index>.csv` (force-added — `.gitignore` has `*.csv`);
3. its **own new** `docs/guide/<index>.md`.

Agents **must NOT touch** the shared "list" surfaces — these are reconciled
**once, by the orchestrator**, at integration:
`README.md` index list · the `thermofeel/thermofeel.py` module-docstring list ·
`docs/guide/overview.md` · `mkdocs.yml` nav · `CHANGELOG.md` ·
`thermofeel/__init__.py` (`__version__`). This is the DI "keep-in-step" rule and
the single-source-of-truth version rule; centralising these avoids five agents
colliding on the same lines.

### 4.2 Pre-assigned insertion anchors (distinct → mostly auto-merge)

| Index | `thermofeel.py` anchor | Test-method neighbours |
|-------|------------------------|------------------------|
| SSI | after `calculate_discomfort_index` | after `test_discomfort_index` |
| RSI | after `calculate_normal_effective_temperature` | after `test_normal_effective_temperature` |
| AT-radiation | after `calculate_apparent_temperature` | after `test_apparent_temperature` |
| PMV / PPD | end of module, after `calculate_heat_index_adjusted` | end of test class |
| PET | after the PMV/PPD block (if unblocked) | end of test class |

### 4.3 Sprints (indices run in parallel; grouped for review sanity)

- **Sprint 0 — Research** ✅ done (this document).
- **Sprint 1 — Easy closed-forms (parallel ×3):** AT-radiation, RSI, SSI.
- **Sprint 2 — PMV / PPD** (iterative; larger; has the richest validation table).
- **Sprint 3 — PET** — **only if the provenance gate is unblocked** (§5.5);
  otherwise deferred to `IDEAS.md` with a written reason.
- **Sprint 4 — Integration:** orchestrator merges sub-branches, reconciles the
  shared lists, does the single `2.4.0` bump + CHANGELOG, runs `make all` +
  `make docs`, opens the single PR, drives `/copilot-review-loop`.

## 5. Per-index specifications

### 5.1 Apparent Temperature — radiation form

- **Provenance (G1 PASS).** Steadman, R.G. (1994) *Norms of apparent temperature
  in Australia*, Aust. Met. Mag. 43:1–16, DOI `10.1071/es94001`; operational
  formula reproduced on the BoM page <http://www.bom.gov.au/info/thermal_stress/>.
- **Formula (°C internally):**
  `AT = Ta + 0.348·e − 0.70·ws + 0.70·q/(ws + 10) − 4.25`,
  with `e` = water-vapour pressure (hPa), `ws` = 10 m wind (m/s), `q` = **net
  radiation absorbed per unit body-surface area (W·m⁻²)**. `e` uses BoM's
  `e = (rh/100)·6.105·exp(17.27·Ta/(237.7+Ta))` — constants identical to
  thermofeel's `calculate_nonsaturation_vapour_pressure` (reuse it).
- **Signature:** `calculate_apparent_temperature_radiation(t2_k, va, rh, q)` →
  Kelvin. Separate function (not an optional arg on the existing AT) so the `q`
  requirement and the different constants are explicit.
- **`q` is a caller-supplied input** (like `cossza`): it is body-absorbed net
  radiation, **not** an NWP surface flux and **not** MRT. Document this; do not
  fabricate a mapping.
- **Validation:**
  - G2: with `q=0` the radiation form is *not* equal to the non-radiation form
    (different constants: 0.348 vs 0.33, −4.25 vs −4.00) — encode that.
  - G3 anchor (VERIFIED consistent with the BoM formula in Sprint 0):
    `Ta=23 °C, rh=70 %, ws=1 m/s, q=50 → AT≈28.1 °C` (validation-data, MIT).
  - ⚠️ The dossier's `Ta=25, rh=30, v=0.1, q=100 → 25.3 °C` value **does NOT
    reconcile** with the BoM formula (hand-calc ≈ 30.9 °C). **Do not pin it**
    until the discrepancy is explained; derive the remaining rows analytically
    from the cited formula and document them as exact evaluations.
- **Oracle:** `pythermalcomfort.at`, `jsthermalcomfort` (both MIT) — cross-check
  only, once the formula convention is pinned.

### 5.2 Relative Strain Index

- **Provenance (G1 PASS, variant note).** Origin Lee, D.H.K. & Henschel, A.
  (1963/1966) DOI `10.1111/j.1749-6632.1966.tb43059.x`; the implementable
  closed form is stated with units by Asghari et al. (2020),
  DOI `10.2174/1874213002013010011`.
- **Formula (implement the hPa form):** `RSI = (Ta − 21) / (58 − e)`, `Ta` in °C,
  `e` = ambient water-vapour pressure in **hPa** (from
  `calculate_nonsaturation_vapour_pressure(t2_k, rh)`). Dimensionless.
  - **Variant caveat:** the alternative literature form
    `(10.7 + 0.74·(Ta−35))/(44 − Pa)` uses `Pa` in a different unit (likely
    mmHg) and is **[UNVERIFIED]** from primary text — do **not** implement it;
    document its existence only.
- **Signature:** `calculate_relative_strain_index(t2_k, rh)` → dimensionless.
- **Validation:**
  - G3: Asghari et al. Table 4 rows (15-yr summer monthly means, rounded), e.g.
    `Arak Jul: Ta=28.12, rh=22.95 → 0.151`; `Bandar Abbas Jul: 34.27, 66.5 → 0.60`;
    `Sari Aug: 27.54, 73.7 → 0.21` (source URL in dossier). Re-verify each against
    the formula + thermofeel's vapour pressure; pin at the tolerance the rounding
    supports.
  - Bands (docs only): Błażejczyk 2011 `10.2478/v10288-012-0004-7`.

### 5.3 Summer Simmer Index

- **Provenance (G1 CONDITIONAL).** Pepi, J.W. (1987) *The Summer Simmer Index*,
  Weatherwise 40(3):143–145, DOI `10.1080/00431672.1987.9933356` (bibliographic,
  DOI-backed). The **equation** is only available via secondary sources — ship
  the **common 1987 closed form** and state the provenance limitation in the
  docstring + guide. Do **not** claim to implement the later "New SSI" table
  (it does not match the closed form).
- **Formula (°F internally):** `SSI_F = 1.98·[Tf − (0.55 − 0.0055·rh)(Tf − 58)] − 56.83`.
- **Analytic gate (G2, strong):** the inner bracket is Thom's THI in Fahrenheit,
  so `SSI_F = 1.98·THI_F − 56.83`. In °C this is Thom's DI with the reference
  constant `58 °F = 14.444…°C` (thermofeel's `calculate_discomfort_index` uses
  14.5 °C, so SSI is a near-exact affine image of DI — cross-check to that, not
  bit-exact). Encode the affine identity as a test.
- **Signature:** `calculate_summer_simmer_index(t2_k, rh)` → Kelvin.
- **Validation:** G3 is the weak point — published tables are the author's *New*
  SSI (mismatch). Lean on the **analytic identity** (G2) plus a small number of
  exact closed-form evaluations; document that authoritative independent point
  tables are scarce. **If the maintainer judges provenance insufficient, this
  index is dropped to `IDEAS.md`.**

### 5.4 PMV / PPD (ISO 7730:2005 baseline)

- **Provenance (G1 PASS).** ISO 7730:2005 (equations, Annex D program +
  validation table); Fanger 1970. Corroborated by CBE `pythermalcomfort`
  (MIT), CBE Comfort Tool, Ladybug.
- **Formula:** full ISO set — `p_a = rh·10·exp(16.6536 − 4030.183/(t_a+235))` Pa;
  `I_cl = 0.155·clo`; `M = 58.15·met`; `f_cl` piecewise; `h_c = max(2.38·|t_cl−t_a|^0.25,
  12.1·√v_ar)`; the `t_cl` fixed-point (init + `EPS=0.00015`, fail if `N>150`);
  `PMV = (0.303·e^{−0.036M}+0.028)·L`; `PPD = 100 − 95·exp(−0.03353·PMV⁴ −
  0.2179·PMV²)`. (Full term list in the Sprint-0 dossier; transcribe from the ISO
  PDF, not memory.)
- **Signatures:**
  `calculate_pmv(t2_k, mrt_k, var, rh=None, vapour_pressure_hpa=None, met=1.2, clo=0.5, wme=0.0)`
  → dimensionless; `calculate_ppd(pmv)` → %. Require exactly one of
  `rh`/`vapour_pressure_hpa`. Vectorise the iteration (whole-array fixed point,
  like the Liljegren solver); non-convergence → `NaN` (documented).
- **Key semantics:** `var` is **relative air velocity at the body**, *not* 10 m
  wind — call this out in the docstring (do not silently pass `va`).
- **Defaults:** `met=1.2` (ISO 8996 sedentary), `clo=0.5` (light) — cited as
  conventions, not ISO universals; maintainer confirms (§6).
- **Validation (G3, strong):** ISO 7730 Annex D Table D.1 — 13 rows captured in
  Sprint 0, e.g. `(t_a=22, t_r=22, v=0.1, rh=60, met=1.2, clo=0.5) → PMV=−0.75,
  PPD=17`; `(27,27,0.1,60,1.2,0.5) → 0.77, 17`; `(27,27,0.3,60,1.2,0.5) → 0.44,
  9`; … plus `(19,18,0.1,40,1.2,1.0) → −0.70, 15.3`. Pin all; also G2
  `PMV=0 ⇒ PPD=5`. Oracle: `pythermalcomfort.pmv_ppd_iso`.

### 5.5 PET — Physiological Equivalent Temperature  ⚠️ BLOCKED

- **Definition (G1 partial):** Höppe 1999 `10.1007/s004840050118` defines PET
  (air temperature of a reference indoor setting giving the same core+skin state)
  — **but the full MEMI energy-balance equations are not visible in open
  primary sources**, and **no published verbatim `(inputs → PET)` validation
  rows were found**. Existing implementations are `pythermalcomfort` (MIT, but
  SciPy `fsolve`) and Ladybug (**AGPL — must not be copied**).
- **Why blocked:** shipping PET now would mean transcribing equations we cannot
  fully cite and validating only against another implementation's output — which
  violates `DESIGN.md` §5 ("a formula with no citation is not ready to merge")
  and the IDEAS gate ("needs a citation and reference values").
- **Unblock condition (maintainer decision, §6):** obtain a fully-cited equation
  set — **Walther & Goestchel 2018** (`10.1016/j.buildenv.2018.03.054`, "first
  exhaustive explanation of PET") and/or **VDI 3787 Part 2** — via ECMWF
  institutional access, *and* a set of published reference values (VDI examples /
  Höppe 1999 tables). If obtained → Sprint 3 proceeds with those as the cited
  source + G3 data. If not → **PET stays in `IDEAS.md`** (with this reason) and
  the 2.4.0 batch ships the other four.
- **If unblocked, target:** the Walther & Goestchel corrected steady PET,
  implemented independently from the published equations (numpy-only, bounded
  deterministic iteration; no SciPy, no AGPL code), signature
  `calculate_pet(t2_k, mrt_k, va, rh, *, pressure_hpa=1013.25, met=…, clo=0.9, …)`
  → Kelvin, with the VDI reference person/environment documented and cited.

## 6. Open decisions for the maintainer

1. **PET disposition** (biggest): (a) attempt now **only if** ECMWF can supply
   Walther & Goestchel 2018 / VDI 3787-2 full text + reference values, else
   (b) defer PET to `IDEAS.md` and ship the other four as 2.4.0. Never
   reverse-engineer the AGPL implementation.
2. **SSI provenance:** ship the common Pepi-1987 closed form with the affine-THI
   analytic gate + a provenance caveat, or defer to `IDEAS.md`?
3. **PMV parameters:** accept `met`/`clo`/`var` as explicit args with cited
   defaults (`met=1.2`, `clo=0.5`, `var` = body-level velocity), acknowledging
   this broadens the "standard met inputs" premise for physiological indices?
4. **Version/PR:** one `2.4.0` batch in a single PR stacked on #54 (confirm), or
   split?

## 7. Risk register

| Risk | Sev | Mitigation |
|------|-----|------------|
| PET equations not fully citable; no published validation | High | Provenance gate §5.5; defer unless sources obtained |
| SSI formula from secondary sources; New-SSI table mismatch | Med | Analytic THI gate; ship common form w/ caveat or defer |
| RSI has two literature forms w/ different units | Med | Implement only the peer-reviewed hPa form; document the other |
| Dossier oracle values not always matching the formula (AT-rad 25.3) | Med | Mandatory re-derivation gate (G1/G2) before pinning |
| 5 agents editing shared list files | Med | Anchor protocol §4.1–4.2; orchestrator-only shared edits |
| PMV `var` mistaken for 10 m wind | Low | Explicit docstring; distinct arg name `var` |
| Multi-worktree merge friction | Low | Distinct anchors; sequential integration onto `feat/new-indices` |

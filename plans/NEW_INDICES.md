# New Indices: DIRINT and PET

> **Status:** both **accepted** for implementation (decided 2026-07-31), promoted
> out of `IDEAS.md`. This file carries the design detail; `TODO.md` remains the
> canonical list of accepted work and points here.

Two additions are planned:

| Item | What | Blocker status | Where the value is |
|---|---|---|---|
| **DIRINT** (Perez et al. 1992) | A correction on top of the existing `approximate_fdir_disc` | **None** — primary source is open access (verified) | Only with the stability index, on <= 1 h data (measured: −11.4 % RMSE) |
| **PET** (Höppe 1999) | Physiological Equivalent Temperature, the outdoor companion to UTCI | Citable equation set (see §2) | Vectorisation: pythermalcomfort is ~2700x slower than thermofeel's comparable solver (measured) |

---

# 1. DIRINT — Perez et al. (1992)

## 1.1 What it is

DIRINT is **not a new separation model**: it is a *multiplicative correction to
DISC*, driven by a 4-dimensional binned lookup table. thermofeel already ships
the DISC core (`approximate_fdir_disc`), so the addition is small and shares
almost all of its work:

```
kt, airmass, dni_disc  = <existing DISC core>
kt'                    = kt / (1.031 exp(-1.4 / (0.9 + 9.4/AM)) + 0.1)     [eqn 1]
delta_kt'              = 0.5 (|kt'(t) - kt'(t+1)| + |kt'(t) - kt'(t-1)|)   [eqn 2, 3]
w                      = exp(0.07 Td[degC] - 0.075)                        [eqn 4]
dni_dirint             = dni_disc * C[bin(kt'), bin(zen), bin(dkt'), bin(w)]  [eqn 5]
fdir                   = dni_dirint * cossza      (clipped to [0, ssrd])
```

Reference: **Perez, R.R., Ineichen, P., Seals, R.D., Maxwell, E.L. & Zelenka, A.
(1992). Dynamic global-to-direct irradiance conversion models. *ASHRAE
Transactions* 98(1), 354–369.** Open access (published version):
<https://archive-ouverte.unige.ch/unige:38583>

## 1.2 Evidence base — what DIRINT actually buys us

Measured before committing, on the shipped SURFRAD validation harness
(3 stations x 4 seasonal months of 2023, BSRN-style QC, hourly, strict common
subset **N = 3362**; pvlib used as the reference implementation):

| Variant | RMSE (W m⁻²) | vs DISC |
|---|---:|---:|
| DISC (shipped) | 43.2 | — |
| DIRINT, no stability, no dew point | 42.4 | −1.9 % |
| DIRINT **+ dew point** | 48.2 | **+11.7 % (worse)** |
| DIRINT **+ stability (delta_kt')** | 38.3 | **−11.4 %** |
| DIRINT + stability + dew point | 42.9 | −0.7 % |

Three conclusions drive the design below:

1. **All of the benefit is in the stability index.** Without it DIRINT is
   statistically indistinguishable from DISC (0.8 W m⁻² on 43). The API must
   therefore make `delta_kt'` a first-class input, not an afterthought.
2. **The dew-point term degrades accuracy on our data**, consistently at every
   station (bon +3.1, dra +5.8, gwn +5.0 W m⁻²). It must be **opt-in and
   documented**, not on by default. (Caveat: dew point was derived from measured
   T/RH; it round-trips through thermofeel's independent Magnus variant to
   0.23 %, and the degradation is large and universal, so the finding is robust
   — but a re-check with a true measured/model dew point belongs in validation.)
3. `delta_kt'` needs the *neighbouring time steps* of the same grid point, and
   Perez/pvlib guidance is that it must **not** be used when spacing is
   >= 1.5 h. ECMWF open data is 3 h/6 h stepped, so **the open-data path gets no
   benefit**; the gain is for <= 1-hourly data (MARS/operational, reanalysis).
   This must be stated prominently in the docstring and guide page.

## 1.3 Prerequisite (step 1): verify the coefficient table against the paper

`DESIGN.md` §5 — "a formula with no citation is not ready to merge" — applies to
the 1260-entry table exactly as it applies to a formula. **pvlib's table is
transcribed from MATLAB code, not from the paper** (its own source comment says
so), so we must not simply copy it.

- [ ] Download the published version from the UNIGE archive (link above).
- [ ] Verify the **1260 coefficients** (6 x 6 x 7 x 5), **all four sets of bin
      edges**, and **eqns 1–5** against the printed tables.
- [ ] Transcribe into thermofeel **from the paper**; use pvlib only as a
      cross-check oracle (same pattern as Erbs/DISC and PMV).
- [ ] Record the verification in the validation README (what was checked
      against which table/page).

**If the coefficients turn out not to be printed in the paper, stop and
escalate** — that would put DIRINT in the same position as PET (§2), and we do
not ship an uncitable table.

## 1.4 Public API

```python
def approximate_fdir_dirint(
    ssrd, cossza, doy, *,
    pressure_hpa=1013.25,
    td_k=None,              # dew point [K] -> Perez w (eqn 4); OFF by default (§1.2)
    kt_prime_prev=None,     # kt' of the previous step, same grid (eqn 2)
    kt_prime_next=None,     # kt' of the next step, same grid
    solar_constant=SOLAR_CONSTANT,
    min_cossza=0.065,
    max_airmass=12.0,
) -> np.ndarray
```

plus a small public helper so callers can produce the neighbours in one line
each:

```python
def clearness_index_zenith_independent(
    ssrd, cossza, doy, *, pressure_hpa=1013.25,
    solar_constant=SOLAR_CONSTANT, min_cossza=0.065, max_airmass=12.0,
) -> np.ndarray     # kt'  (eqn 1)
```

Rationale:

- Takes `doy`, **not** a `times`/`DatetimeIndex` — consistent with
  `approximate_fdir_disc`, and pandas is not a thermofeel dependency.
- Neighbours are passed as **kt'**, not as raw `ssrd_prev`/`cossza_prev`/
  `ssrd_next`/`cossza_next`: kt' is the quantity the model actually needs, it
  halves the argument count, and it matches thermofeel's established
  "caller supplies the derived input" style (`cossza`, `q` in AT-radiation).
- Omitting both neighbours selects the model's own **"stability unavailable"**
  bin; omitting `td_k` selects its **"dew point unavailable"** bin. Both are
  designed-in states of the table, not error paths.
- One-sided differencing when only one neighbour is given (eqn 3).
- Returns direct **horizontal** radiation clipped to `[0, ssrd]`, NaN-propagating
  and night-zeroed — identical contract to `approximate_fdir_erbs`/`_disc`.

## 1.5 Specification (bin edges — to be confirmed against the paper in step 1)

Table shape `(6, 6, 7, 5)`, indexed `[kt'_bin, zenith_bin, delta_kt'_bin, w_bin]`.

| Variable | Bins | Edges |
|---|---:|---|
| kt' | 6 | 0.24, 0.40, 0.56, 0.70, 0.80, (<= 1) |
| zenith [deg] | 6 | 25, 40, 55, 70, 80, (>= 80) |
| delta_kt' | 7 | 0.015, 0.035, 0.07, 0.15, 0.30, (<= 1), **+ "unavailable"** |
| w [cm] | 5 | 1, 2, 3, (>= 3), **+ "unavailable"** |

## 1.6 Performance design

This is the explicit ask, so it is specified as targets plus concrete technique.

**Measured baseline** (N = 1 038 240, a 0.25° global field, same machine):

| | time |
|---|---:|
| `thermofeel.approximate_fdir_disc` (today) | **0.29 s** |
| `pvlib.irradiance.disc` | 0.44 s |
| `pvlib.irradiance.dirint` (no stability, no dew) | **0.54 s** |
| `pvlib.irradiance.dirint` (+ dew point) | 0.59 s |

**Targets:** `approximate_fdir_dirint` <= **1.4x** `approximate_fdir_disc`
(i.e. ~0.40 s at this size), and comfortably faster than pvlib's 0.54 s. No
more than two extra full-length temporaries beyond the DISC path.

Techniques:

1. **No pandas anywhere.** pvlib requires a `DatetimeIndex` purely for its Series
   bookkeeping; for a field library that is pure overhead.
2. **Share the DISC core.** Refactor an internal `_disc_core(...) -> (dni, kt,
   airmass)` used by *both* `approximate_fdir_disc` and DIRINT, so nothing is
   computed twice. The existing public DISC behaviour must be bit-identical
   after the refactor (regression CSV + validation re-run prove it).
3. **Branch-free binning via `np.searchsorted`** on four tiny constant edge
   arrays (<= 7 entries), replacing pvlib's ~24 boolean-mask assignments — each
   of which is a full pass over N. Four passes instead of ~24.
4. **Fused flat gather.** Flatten the table to 1260 entries and compute a single
   index `idx = ((kb*6 + zb)*7 + db)*5 + wb`, then one `TABLE.take(idx)` —
   replacing 4-array fancy indexing with one gather.
5. **Sentinels as real bins, not NaN branches.** "No dew point" and "no
   stability" are dedicated bins, so the default path (neither supplied) folds
   to *constant* bin indices with no masking work at all.
6. **Inherit the existing edge contract** (finite mask, night zeroing, clip to
   `[0, ssrd]`) from the Erbs/DISC implementations rather than re-deriving it.
7. **Measure, then decide** on daytime-mask compression (~50 % of a global field
   is night, so gathering on a compressed subset may pay). Prototype it; keep it
   only if it beats the simple path. Do not assume.

## 1.7 Validation plan

Extends the existing campaign rather than starting a new one.

- **Tier 0 — provenance:** the §1.3 table/bin/equation verification against the
  paper, recorded in the README.
- **Tier 1 — implementation fidelity:** element-wise vs `pvlib.irradiance.dirint`
  on the real SURFRAD inputs, constants matched, to machine precision — the
  pattern already used for Erbs and DISC. Must cover all four availability
  combinations (with/without stability x with/without dew point).
- **Tier 2 — ground truth:** re-run the SURFRAD harness in
  `validation/fdir/`, reporting the five-variant table of §1.2 as the shipped
  result, with pre-registered gates. Add the per-station dew-point effect.
- **Tier 3 — global field:** extend `validate_fdir_global.py` with the same
  exact bounds/night gates; note that the open-data path exercises the
  "stability unavailable" bin (§1.2 conclusion 3).
- Unit tests: bin-edge boundaries (values exactly on each edge), sentinel bins,
  NaN propagation, one-sided differencing, night zeroing, `[0, ssrd]` clipping,
  integer-dtype inputs. Regression CSV as for the other estimators.
- Coverage must stay at 100 %.

## 1.8 Documentation

- `docs/guide/approximations.md`: add DIRINT alongside Erbs/DISC, **with the
  measured table of §1.2** and an explicit "when does this help?" statement
  (only <= 1 h stepped data with neighbours; not ECMWF open data).
- Docstring: state the dew-point finding and the >= 1.5 h contraindication.
- `docs/api.md`, `CHANGELOG.md`, and the `--approximate-fdir` option of
  `examples/compute-thermal-indices.py` (add `dirint`; it must degrade
  gracefully to the no-stability bin when neighbouring steps are unavailable).
- Fold the new function into `plans/ROBUSTNESS.md` §5's contract.

## 1.9 Risks and open questions

| | Risk | Mitigation |
|---|---|---|
| R1 | Coefficients not printed in the paper | Step 1 checks this *first*; if absent, stop and escalate (do not ship an uncitable table) |
| R2 | Dew-point term degrades accuracy (measured) | Opt-in, off by default, documented; re-check with a true dew point in validation |
| R3 | Stability index invalid at >= 1.5 h spacing, so open data gains nothing | Document prominently; position DIRINT as for <= 1-hourly data |
| R4 | First **data table** in a library of pure formulas | Module-level constant with a provenance comment naming the paper table; no data files, no new deps |
| R5 | Users adopt DIRINT expecting improvement and get DISC-equivalent results | "When does this help?" section in the guide, and the measured table in the docs |

Open questions for the maintainers:

- Should `approximate_fdir_dirint` accept `ssrd_prev`/`cossza_prev` (+`_next`)
  as an ergonomic alternative to `kt_prime_prev`/`_next`? Current proposal
  favours kt' + helper; the alternative is friendlier but adds four arguments.
- Should the example gain a `--approximate-fdir dirint` path that fetches
  neighbouring steps automatically, or stay single-step (and therefore always
  hit the "stability unavailable" bin)?

## 1.10 Definition of done

- [ ] Coefficients, bin edges and eqns 1–5 verified against Perez et al. (1992)
- [ ] `approximate_fdir_dirint` + `clearness_index_zenith_independent`, numpy-only
- [ ] DISC core shared; existing DISC output bit-identical
- [ ] pvlib parity to machine precision across all four availability modes
- [ ] Performance targets of §1.6 met and recorded in the benchmark suite (§3)
- [ ] SURFRAD + global validation legs updated; `validation/fdir/README.md`
      carries the measured comparison
- [ ] Guide page, `api.md`, `CHANGELOG.md`, `ROBUSTNESS.md` updated
- [ ] `make all` green at 100 % coverage

---

# 2. PET — Höppe (1999)

**Plan pending:** an external analysis is being supplied and will be folded in
here. Recorded now so the context is not lost:

- **Motivation (measured):** `pythermalcomfort.pet_steady` is a pure Python loop
  — 2.4 ms/point, i.e. **~41 minutes for one 0.25° global field**. thermofeel's
  comparable *iterative* solver (`calculate_pmv`) runs at 0.89 µs/point, a
  ~2700x gap. Note this is specific to `pet_steady`: `set_tmp` in the same
  library **is** vectorised (~15 s for the same grid), so the slowness is an
  implementation choice (SciPy `fsolve` per point), not a property of the model.
- **Known blocker (from `IDEAS.md`):** the MEMI equations are not in open
  primary sources and no published `(inputs -> PET)` reference rows were found.
  Höppe (1999) `10.1007/s004840050118` and Walther & Goestchel (2018)
  `10.1016/j.buildenv.2018.03.054` are **both closed access** (verified);
  Ladybug's implementation is AGPL and must not be copied.
- **What the incoming analysis must supply** to unblock: (a) a fully-cited
  equation set, (b) published reference values to validate against, (c) an
  unambiguous statement of *which* PET variant (PET vs mPET; RayMan and
  pythermalcomfort are documented to disagree).
- **Implementation template, once unblocked:** the PMV/PPD work in 2.3.0 —
  independent implementation from the standard, vectorised bounded iteration
  mirroring `liljegren.py`, non-convergence -> NaN, validated against published
  rows with pythermalcomfort as a secondary oracle only (never a runtime
  dependency).

---

# 3. Supporting: benchmark suite

Accepted alongside the above (previously `IDEAS.md` "Performance"). Needed to
hold DIRINT and PET to their performance targets rather than asserting them.

- A small reproducible harness timing each public index over a global-grid-sized
  array, reporting time and µs/point, runnable on demand (not in the CI gate —
  timings are machine-dependent).
- Seeded synthetic inputs; report per-index cost and flag regressions against a
  recorded baseline.
- Seed values already measured on the reference machine: `calculate_pmv`
  0.89 µs/pt, `calculate_utci` 1.15 µs/pt, `approximate_fdir_disc` 0.28 µs/pt
  (0.29 s at N = 1 038 240).
- Doubles as the evidence trail for §1.6 and for the PET vectorisation claim.

# New Indices and Adapters: DIRINT, PET, MRT observation adapters

> **Status:** all **accepted** for implementation (decided 2026-07-31), promoted
> out of `IDEAS.md`. This file carries the design detail; `TODO.md` remains the
> canonical list of accepted work and points here.

Three additions are planned:

| Item | What | Blocker status | Where the value is |
|---|---|---|---|
| **DIRINT** (Perez et al. 1992) | A correction on top of the existing `approximate_fdir_disc` | **None** — primary source is open access (verified) | Only with the stability index, on <= 1 h data (measured: −11.4 % RMSE) |
| **PET** (Höppe 1999) | Physiological Equivalent Temperature, the outdoor companion to UTCI | Citable equation set (see §2) | Vectorisation: pythermalcomfort is ~2700x slower than thermofeel's comparable solver (measured) |
| **MRT observation adapters** (§3) | Reconstruct the kernel's flux inputs from station GHI or cloud cover | **None** | MRT (and hence UTCI) at stations that lack the full flux set, through the *same* kernel as ERA5-HEAT |

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

# 2. PET — Höppe (1999) — IMPLEMENTED

**Status: shipped** (`thermofeel.calculate_pet`, this PR). The citation blocker
below was cleared by an external implementation brief that supplied the
reconciled equation set from Höppe (1999) and Walther & Goestchel (2018),
cross-checked against the authors' own released code. Validation lives in
`validation/pet/`; the model decisions (frozen clothing temperature, Woodcock
diffusion, constant alpha, whole-body activity watts, blood heat capacity) are
recorded in the function docstring and the CHANGELOG.

Delivered against the brief: identity to 3.6e-06 K, rms 1.00 K against Höppe's
published Table 1 (the two source papers differ by up to 1.31 K on those same
cases), rms 0.181 K against `pythermalcomfort.pet_steady`, and a 308x
throughput gain from the triangularised two-stage bisection — numpy-only, no
SciPy. The original Höppe skin-diffusion variant is deliberately **not**
implemented: the brief marked its own port UNVERIFIED, and shipping unverified
physics would breach `DESIGN.md` §5. The guide page carries the consequence,
that the Matzarakis/Mayer stress classes derive from that variant.

Remaining follow-ups: obtain **VDI 3787 Part 2** (still not held) to settle the
normative diffusion model, the blood heat capacity and the posture labelling;
and the §9 performance items beyond cache blocking (hoisting the stage-B
invariants, in-place residual arithmetic, an optional float32 path).

Original assessment, retained as the record of why it was deferred for 2.3.0:

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

# 3. MRT observation adapters (GHI and cloud-cover routes)

Source: an external implementation brief ("MRT observation adapters for
thermofeel", July 2026), reviewed and corrected against the code below. The
brief's physics, literature survey and uncertainty budget are adopted; three
factual corrections and one significant scope reduction are recorded here.

## 3.1 Decision

Add two **observation adapters** that reconstruct the radiative inputs
`calculate_mean_radiant_temperature` already consumes, and leave the MRT kernel
untouched:

| Adapter | Required inputs | Optional |
|---|---|---|
| **A — GHI** | global horizontal irradiance, 2 m temperature, humidity, `cossza`, day of year | albedo, skin temperature, cloud fraction |
| **B — cloud cover** | total cloud cover, 2 m temperature, humidity, `cossza`, day of year | albedo, skin temperature |

Adapter B is Adapter A prefixed by a cloud -> irradiance step, so they share the
whole downstream path. The architectural value is that station MRT stays
**definitionally comparable** with ERA5-HEAT and operational MRT: same kernel,
same angle factors, same `fp`, same `a_ir/eps_p`.

MENEX_2005 is implemented separately as a **legacy-compatibility function**
(§3.7), not as a default — the brief's §6 analysis of its defects is convincing
and is adopted.

**Wind is not an input.** Neither route needs it: it never appears in the
flux-reconstruction route, and in MENEX it appears only inside `Irc`, which
cancels in the MRT inversion. Wind remains a downstream UTCI input. The adapter
signatures must not demand it.

## 3.2 Correction 1 (blocking): the kernel signature is not what the brief assumes

The brief states that its `FluxSet` "mirrors the ERA5 names already in the
library (`ssrd`, `fdir`, `ssru`, `strd`, `stru`) so that
`calculate_mean_radiant_temperature` is reused with no change". **It does not.**
The shipped signature is seven arguments and takes **net** fluxes:

```python
calculate_mean_radiant_temperature(ssrd, ssr, dsrp, strd, fdir, strr, cossza)
```

and derives the upward components internally:

```python
dsw = ssrd - fdir      # diffuse down      rsw = ssrd - ssr    # reflected (up)
lur = strd - strr      # longwave up       Istar = dsrp        # direct normal
```

There is no `ssru`/`stru` in the API. An adapter emitting the brief's `FluxSet`
literally would not plug in. Required mapping:

| Brief | thermofeel | Conversion |
|---|---|---|
| `S↓` | `ssrd` | identity |
| `S↓_dir` (horizontal) | `fdir` | identity |
| `S↑` (upward SW) | `ssr` | `ssr = ssrd - S↑ = (1 - albedo)·ssrd` |
| `I*` (normal to beam) | `dsrp` | 7th kernel input; `approximate_dsrp(fdir, cossza)` |
| `L↓` | `strd` | identity |
| `L↑` (upward LW) | `strr` | `strr = strd - L↑` |
| `S↓_diffuse` | — | **not an adapter output**; the kernel computes `ssrd - fdir` itself |

The `FluxSet` container must therefore carry the seven kernel arguments (or be
dropped in favour of returning them directly). Suggest naming it after what it
is — the kernel's input tuple — rather than after ERA5 upward fluxes.

## 3.3 Correction 2: the shortwave path already exists (major scope reduction)

Adapter A's entire shortwave chain is already shipped and validated in 2.3.0:

| Brief proposes | Already in thermofeel |
|---|---|
| `kt = GHI/(E₀·f_ecc·sin γ)` | computed inside `approximations.approximate_fdir_erbs` (identical: `ssrd / (solar_constant · e0 · max(cossza, min_cossza))`) |
| `kd = Erbs(kt)`, `S↓_dir = (1−kd)·GHI` | **is** `approximate_fdir_erbs(ssrd, cossza, doy=doy)`, which returns `fdir` directly. Coefficients verified identical to the brief's transcription |
| `calculate_diffuse_fraction(kt, model="erbs")` | **not needed** — the kernel computes `dsw = ssrd - fdir` internally |
| `I* = S↓_dir / sin γ`, guarded | **is** `approximate_dsrp(fdir, cossza, threshold=0.1)`, shipped with the low-sun guard |

So **the only genuinely new physics is the longwave pair** (`strd`, `strr`) plus
Adapter B's cloud -> GHI step. That is a much smaller change than the brief
implies, and it means the shortwave side inherits the SURFRAD/pvlib validation
already done.

Two consequences:

- **`sin γ` is `cossza`.** Solar elevation and zenith are complementary
  (`sin(elevation) = cos(zenith)`), so the adapters take `cossza` like every
  other thermofeel function. Do **not** add a separate solar-elevation argument;
  the kernel already derives `gamma` internally for `fp`.
- **Default separation model should be DISC, not Erbs.** Our own SURFRAD
  validation (§1.2 harness) measured DISC at 43.1 W m⁻² hourly RMSE against
  Erbs's 54.0 — ~20 % better at every station. Better `fdir` means a better
  direct/diffuse split into the kernel. Expose
  `separation_model="disc"|"erbs"|"dirint"` (mapping onto shipped code, plus
  §1), default `"disc"`, with `"erbs"` retained for literature reproduction.
  The brief's suggestion to add Yang4 for sub-hourly use is noted but deferred:
  no implementation exists in-repo and it would need its own validation.

## 3.4 Correction 3: the kernel's Stefan-Boltzmann constant is not CODATA

The brief states the kernel uses "σ = CODATA". It does not — the shipped code
hard-codes `0.0000000567` (5.67e-8) against CODATA's 5.670374419e-8. The effect
is ~0.005 K at 300 K (MRT scales as σ^-0.25), i.e. negligible, but the brief's
§6.6 criticism of MENEX for using 5.667e-8 applies in weaker form to thermofeel
itself. Either state the constant honestly in the docs or promote it to a named
CODATA constant as a separate, clearly-flagged change — **not** silently inside
this work, since it would shift every existing MRT/UTCI regression value.

## 3.5 Adapter A — GHI available

```
fdir  = approximate_fdir_<separation>(ssrd, cossza, doy)      # existing
dsrp  = approximate_dsrp(fdir, cossza)                        # existing
ssr   = (1 - albedo) * ssrd                                   # trivial
strd  = eps_all * SIGMA * t2m**4                              # NEW
L_up  = eps_g * SIGMA * t_skin**4 + (1 - eps_g) * strd        # NEW
strr  = strd - L_up                                           # NEW
```

**Longwave down.** Clear-sky emissivity Brutsaert (1975),
`eps_c = 1.24·(e/T_a)^(1/7)` with `e` in hPa, `T_a` in K; cloud modifier
Unsworth & Monteith (1975), `eps_all = (1 − 0.84c)·eps_c + 0.84c`, clamped to
<= 1. `e` comes from the existing helpers — either
`calculate_saturation_vapour_pressure(td_k)` (exact: `e = e_s(T_d)`) or
`calculate_nonsaturation_vapour_pressure(t2_k, rh)`; pick one, document it, and
do not introduce a third vapour-pressure formula (see `IDEAS.md`, "vapour-pressure
formula audit").

Expose `clear_sky_emissivity="brutsaert"|"prata"|"dilley_obrien"` and
`cloud_emissivity="unsworth_monteith"`. When GHI is available, prefer
Crawford & Duchon (1999) for `c` — derived from the observed/clear-sky solar
ratio rather than an okta count (Staiger & Matzarakis 2010 evaluate exactly
this for human-biometeorology use).

**Skin temperature.** `t_skin = t2m` is the documented cheap default; accept an
optional `skin_temperature` so ERA5 `skt`, a station sensor or satellite LST can
be substituted. `eps_g` default 0.95, exposed.

## 3.6 Adapter B — cloud cover only

Prefix Adapter A with Kasten & Czeplak (1980):

```
G_cs = 910·sin γ − 30          W m-2, valid γ > 5 deg     ->  910·cossza − 30
K_c  = 1 − 0.75·(N/8)^3.4                                 ->  1 − 0.75·tcc^3.4
ssrd = K_c · G_cs
```

then Adapter A verbatim, with `c = tcc` for the longwave cloud modifier.

- **`G_cs` must be swappable.** The `910·sin γ − 30` fit is Hamburg-specific with
  no turbidity, water-vapour, altitude or aerosol dependence. Expose
  `clear_sky_model`; recommend CAMS **McClear** for production (already
  Copernicus infrastructure, keeps aerosol treatment consistent with the ECMWF
  chain) as an optional external source, never a hard dependency.
- **The 0.75 coefficient is a cloud-type average** (0.39 cirriform -> 0.84
  nimbostratus). Expose it; it is the dominant term in the uncertainty budget
  (§3.9).
- **Public boundary is a fraction 0–1**, with an okta helper. **WMO okta code 9**
  (sky obscured) is not a cloud amount: reject it or map to a documented
  sentinel; never coerce to 8.
- **`tcc_source` (`human`|`automatic`|`model`|`unknown`) is carried as metadata
  and must not branch the computation.** Adopted from the brief: automated
  systems report 0 oktas in 14–19 % of hours vs 1 % for human observers and
  8 oktas in 33–40 % vs 19–24 % (Smith, Bright & Crook 2017, 1.12 M
  station-hours); human/automated pairs agree exactly in only 39 % of hours
  (Wauben et al. 2006). The code path is identical; the error distribution is
  not.

## 3.7 MENEX_2005 — legacy only

**Decision: `legacy` is its own submodule, `thermofeel.legacy`**, namespace-only
and not re-exported at the top level — the same device already used for
`thermofeel.approximations`, and for the same reason: the call site should carry
the marker so the choice stays visible in the code that makes it.

The word "legacy" needs to be unambiguous, because it could be misread as
"deprecated thermofeel API scheduled for removal", which is *not* what it means
here. Proposed module docstring, to be used verbatim (**phrasing needs
maintainer sign-off**):

> **Historical models, retained for reproducibility.**
>
> `thermofeel.legacy` holds complete, faithful implementations of published
> human-biometeorology models that thermofeel does **not** recommend for new
> work, but which are needed to reproduce existing literature and to compare
> against results produced by older tooling.
>
> "Legacy" describes the **model**, not thermofeel's own API: nothing in this
> submodule is deprecated thermofeel code, and nothing here is scheduled for
> removal. These functions are supported, tested and documented to the same
> standard as the rest of the library — the difference is in what we recommend,
> not in what we maintain.
>
> Like `thermofeel.approximations`, this submodule is namespace-only and is not
> re-exported at the top level, so calling code always carries the `legacy.`
> marker. Each function's docstring states what supersedes it, and why.

Consequences for MENEX_2005: ship it there for reproducing BioKlima-based
literature, with the defects documented in the docstring and **no default use**:

- SolAlt bands cloudiness into four steps, so oktas 2/3/4 give *identical* MRT
  and the band edges are discontinuous — ~16 K step for a one-okta change at
  40 deg elevation.
- The nocturnal cloud sign is **inverted**: `L_a` is a clear-sky formula with no
  cloud term, so MENEX makes overcast nights *cooler* than clear ones. Nocturnal
  heat stress is where a heat index carries epidemiological weight; do not use
  MENEX at night.
- `T_g = 1.25·t` is dimensionally unprincipled (multiplicative on a Celsius
  value; pins `T_g = T_a` at 0 °C, warms ground below freezing), worth +4.6 K at
  40 °C.
- **Implement the 0.5-weighted MRT form**, and record the published-source
  erratum: the PDF prints the MRT equation without the 0.5 weights, which
  contradicts its own net-longwave equation and yields ~67 °C on a 15 °C day.
  The 0.5-weighted form is what the derived literature implements, and it makes
  MENEX two-hemisphere — convention-compatible with our kernel.
- Use MENEX's own constants verbatim (`sigma = 5.667e-8`, `0 °C = 273.0`) for
  bit-reproducibility, and say so; do not mix MENEX `K_t` with Kasten-Czeplak
  `G_cs` (they differ by 6–15 % and cross over near 35 deg, worth up to 10.4 K).

## 3.8 API sketch (corrected)

SI throughout: Kelvin in/out, irradiance W m⁻², cloud cover a 0–1 fraction at
the public boundary.

```python
# building blocks
calculate_clear_sky_ghi(cossza, model="kasten_czeplak")                  -> W m-2
calculate_ghi_from_cloud_cover(tcc, cossza, attenuation=0.75, ...)       -> W m-2
calculate_downward_longwave(t2m, vapour_pressure_hpa, cloud_fraction,
                            clear_sky="brutsaert",
                            cloud="unsworth_monteith")                   -> W m-2
calculate_upward_longwave(skin_temperature, lw_down, emissivity=0.95)    -> W m-2

# adapters -> the SEVEN kernel arguments (see 3.2)
reconstruct_mrt_inputs_from_ghi(ssrd, t2m, td2m, cossza, doy, *,
                                albedo=0.20, skin_temperature=None,
                                cloud_fraction=None,
                                separation_model="disc", ...)
reconstruct_mrt_inputs_from_cloud_cover(tcc, t2m, td2m, cossza, doy, *,
                                        albedo=0.20, skin_temperature=None,
                                        tcc_source="unknown", ...)

# convenience: adapter + existing kernel, unchanged
calculate_mean_radiant_temperature_from_ghi(...)                         -> K
calculate_mean_radiant_temperature_from_cloud_cover(...)                 -> K

# legacy, namespaced
legacy.menex2005_mrt_solalt(t2m, vapour_pressure_hpa, tcc, cossza)       -> K
legacy.menex2005_mrt_solglob(t2m, vapour_pressure_hpa, ghi, cossza)      -> K
```

Guards: reuse the existing `cossza` floor conventions; clamp `eps_all <= 1`;
clamp `kt` per standard hourly QC; reject okta code 9. NaN propagation and the
array contract follow `ROBUSTNESS.md` §5 like every other function.

## 3.9 Uncertainty — and what we must not claim

Adopted from the brief, at Ta = 25 °C / RH = 50 %:

| Source | MRT impact |
|---|---|
| Cloud type at fixed N | **12.6 K at γ=30°, 18.5 K at γ=60°** |
| ±2 okta observation uncertainty | up to 16 K (flux route) |
| Solar disc obscured or not | 20–30 K at intermediate N |
| Clear-sky model choice | ~10–15 % in GHI |
| Skin-temperature scheme | up to 4.6 K |
| Geometry convention (two-hemisphere vs six-directional) | ~3 K |

**Decision: the adapters return bare arrays, like every other thermofeel
function.** The source brief recommended emitting an uncertainty estimate
alongside MRT; we are not doing that. Every public function in the library
returns a plain `np.ndarray`, and returning a tuple or dataclass here would
fork the API contract for one family of functions. The uncertainty is handled
by *documentation* instead — a stated tolerance class per adapter, plus the two
statements below. If a machine-readable uncertainty is ever wanted, it belongs
in a separate function, not in the return value of these.

Two statements must appear in user-facing docs:

1. **Adapter B returns a conditional mean, not an instantaneous estimate.** An
   okta count cannot say whether the solar disc is obscured; at N = 4 the beam is
   either near-full or near-zero. Adequate for climatology and gap-filling, not
   for instantaneous values. If stations report **sunshine duration**, that is a
   far stronger beam constraint (WMO: direct normal >= 120 W m⁻²) and is worth a
   third adapter later.
2. **Adapter B cannot be represented as meeting the ISO 7726 ±2 °C comfort
   tolerance.** State the tolerance class per adapter.

The two-hemisphere convention carries a documented ~−3 K offset against
six-directional measurement (Holmer et al. 2018). That is a property of the
convention, not of the adapters, and any validation against six-directional or
globe-thermometer data must account for it.

## 3.10 Validation plan

**We already have the ideal dataset cached.** SURFRAD (used for the 2.3.0 fdir
campaign) measures every component the kernel needs:

| Kernel argument | SURFRAD column |
|---|---|
| `ssrd` | `dw_solar` |
| `ssr` | `netsolar` (measured) |
| `dsrp` | `direct_n` (measured, not approximated) |
| `strd` | `dw_ir` (measured) |
| `strr` | `netir` (measured) |
| `fdir` | `direct_n · cossza` |

So: compute a **reference MRT from fully measured fluxes** through the unmodified
kernel, then degrade the inputs to GHI-only (Adapter A) and record bias/RMSE per
station and per solar-elevation band — the brief's Tier 2, on data and a harness
that already exist (`validation/vlib/surfrad.py`, `validation/fdir/`).

- **Tier 1 (unit):** `fp` against VDI at γ = 0/45/90°; `K_c = 1` at N = 0 and
  `0.25` at N = 8; Erbs continuity at kt = 0.22 and 0.80; MRT -> `T_a` when all
  fluxes are blackbody emission at `T_a`; **seam test** — Adapter A fed with
  Adapter B's own GHI must reproduce Adapter B to machine precision.
- **Tier 2 (fixtures):** the SURFRAD degradation study above.
- **Adapter B cloud source — RESOLVED (feasibility verified 2026-07-31).**
  SURFRAD reports no cloud amount, so Adapter B must be validated against an
  independent observation. Deriving `c` from the measured clear-sky ratio is
  partly circular and is rejected. Use **METAR/ASOS sky cover from the Iowa
  State ASOS archive** (`mesonet.agron.iastate.edu`, public, no
  authentication, hourly, historical). Note this is *not* the same source as
  `examples/compute-obs.py`, which uses the NOAA Aviation Weather API — that is
  near-real-time only and cannot serve a 2023 validation period.

  Verified pairings (hourly sky cover successfully retrieved for all three):

  | SURFRAD site | ASOS | separation |
  |---|---|---:|
  | Desert Rock, NV | `DRA` | **0.8 km — effectively co-located** |
  | Bondville, IL | `CMI` (Champaign) | 8.2 km |
  | Goodwin Creek, MS | `UOX` (Oxford) | 34.0 km |

  Design consequences:

  - **Desert Rock is the primary validation site**, being co-located; the other
    two are supporting. Cloud cover is spatially variable, so the 34 km
    separation at Goodwin Creek is a genuine limitation and must be reported
    alongside its statistics, not buried.
  - METAR reports sky *condition codes*, not oktas. Use the standard mapping
    `CLR/SKC -> 0, FEW -> 2, SCT -> 4, BKN -> 6, OVC -> 8` oktas, and document
    that this quantisation is itself an error source on top of the ±2 okta
    observation uncertainty already in §3.9.
  - `VV` (vertical visibility, sky obscured) is the METAR analogue of WMO okta
    code 9: **reject it**, consistent with §3.6. Do not coerce it to 8.
  - METARs are reported near HH:53; aggregate SURFRAD to the same hour and
    match on the hour, reusing the existing hourly aggregation.
  - `tcc_source` for this data is `automatic` (ASOS ceilometer), which per §3.6
    is exactly the population that over-reports both 0 and 8 oktas — so the
    validation must state it, and the resulting error statistics are *not*
    transferable to human-observed SYNOP cloud cover.
- Worth stating in release notes: the ERA5-HEAT validation itself used exactly
  Adapter B's information content (Ta, Td, u10, TCC at 177 SYNOP stations, via
  RayMan Pro), i.e. the reference implementation being replaced here is a
  closed-source Windows GUI.

## 3.11 Risks and open questions

| | Risk | Mitigation |
|---|---|---|
| R1 | `FluxSet` as specified does not fit the kernel | §3.2 mapping is mandatory; seam test enforces it |
| R2 | Adapter B error is large (12–18 K cloud-type term) | Tolerance class in docs; conditional-mean framing; no ISO 7726 claim |
| R3 | Strategy-parameter sprawl (4+ pluggable schemes) | Keep the *default* path fixed and validated; strategies documented as opt-in |
| R4 | Duplicating shipped Erbs/dsrp code | §3.3 — reuse, do not reimplement |
| R5 | Silent σ change would shift every existing MRT/UTCI value | §3.4 — separate, flagged change if done at all |

Resolved (2026-07-31):

- **Uncertainty is not returned.** Bare arrays, as everywhere else in the
  library; tolerance class documented instead — see §3.9.
- **`legacy` is its own submodule**, `thermofeel.legacy`, namespace-only —
  see §3.7 for the agreed rationale and the proposed module docstring.

Still open:

- Adapter B validation source (see §3.10 open item).
- Sign-off on the exact `thermofeel.legacy` docstring wording in §3.7.

## 3.12 Definition of done

- [ ] Adapters emit the **seven** kernel arguments; kernel unchanged
- [ ] Every public function returns a bare `np.ndarray` (no uncertainty tuple)
- [ ] Shortwave path reuses `approximate_fdir_*` and `approximate_dsrp`
- [ ] Longwave helpers added with cited schemes and strategy parameters
- [ ] `thermofeel.legacy` created: namespace-only, not re-exported at top level,
      module docstring per §3.7
- [ ] MENEX in `thermofeel.legacy`, 0.5-weighted form, erratum documented,
      defects and superseding function named in its docstring
- [ ] Seam test (A fed B's GHI == B) passes to machine precision
- [ ] SURFRAD degradation study produces a published bias/RMSE table
- [ ] Tolerance class and conditional-mean caveat in the guide page
- [ ] `make all` green at 100 % coverage

# 4. Supporting: benchmark suite

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

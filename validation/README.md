# thermofeel 2.3.0 validation campaign

Scientific validation of the methods added in the 2.3.0 release (PR
[#55](https://github.com/ecmwf/thermofeel/pull/55)), against published
standards, independent peer-reviewed implementations, and ground
observations. Everything here is reproducible: each index directory contains
the validation code, the generated result tables (`results/*.csv`) and plots
(`plots/*.png`), and a README stating methodology, pre-registered acceptance
criteria, results and limitations.

**This directory is campaign material, not part of the shipped library.**

## Scope and verdicts

| Method | Directory | Primary evidence | Verdict |
|---|---|---|---|
| `calculate_pmv`, `calculate_ppd` | [`pmv_ppd/`](pmv_ppd/) | ISO 7730:2005 Annex D table (12/12 in tolerance); pythermalcomfort over 46 800-point ISO envelope (max ΔPMV 7.9e-3); 58 500 ASHRAE-DB-II field records vs real votes (bin r = 0.994) | **PASS** |
| `approximations.approximate_fdir_erbs` | [`fdir/`](fdir/) | NOAA SURFRAD obs, 3 stations × 4 months (hourly rRMSE 11–23 %, R² ≥ 0.957); pvlib machine-precision cross-check; global IFS open-data consistency | **PASS** |
| `approximations.approximate_fdir_disc` | [`fdir/`](fdir/) | as above (hourly rRMSE 9–18 %, R² ≥ 0.971 — the better default) | **PASS** |
| `approximations.earth_sun_distance_factor` | [`earth_sun_distance_factor/`](earth_sun_distance_factor/) | NREL SPA ephemeris 2015–2025 (max error 0.104 % ≡ 1.4 W m⁻² TOA); pvlib Spencer bit-level | **PASS** |
| `calculate_discomfort_index` | [`discomfort_index/`](discomfort_index/) | pythermalcomfort oracle to its 0.1 °C print resolution (residual = pure rounding noise); analytic identities exact; SURFRAD obs demonstration | **PASS** |
| `calculate_summer_simmer_index` | [`summer_simmer_index/`](summer_simmer_index/) | exact transcription + affine-THI identity; provenance spread across secondary sources quantified (±2–4 °F) | **PASS** (provenance: medium confidence, documented) |
| `calculate_relative_strain_index` | [`relative_strain_index/`](relative_strain_index/) | exact transcription; Magnus-variant sensitivity ≤ 0.0132 ≪ band width 0.10; behavioural gates; obs demonstration | **PASS** |
| `calculate_apparent_temperature_radiation` | [`apparent_temperature_radiation/`](apparent_temperature_radiation/) | exact fidelity to the BoM-published Steadman 1994 form; family coherence with the non-radiative AT (Q\* ≈ 0); obs demonstration | **PASS** |

## Methodology

Each method is validated on a tiered evidence ladder, using every rung
available for it:

1. **Standard conformance** — published validation tables with their own
   tolerances (ISO 7730 Annex D for PMV/PPD).
2. **Independent implementations** — element-wise agreement with
   peer-reviewed independent code on identical inputs (pythermalcomfort,
   pvlib), with rounding/limit switches disabled where the reference allows,
   and quantization accounted for where it does not.
3. **Independent transcriptions** — for indices *defined by* a published
   closed form (DI, SSI, RSI, AT-rad), a from-scratch transcription of the
   authority's formula, plus the formula's analytic identities and
   behavioural properties (monotonicity, boundary values, documented domain
   edges).
4. **Ground observations** — NOAA SURFRAD pyrheliometer/pyranometer data as
   independent truth for the radiation estimators; ASHRAE Global Thermal
   Comfort Database II field votes for PMV; SURFRAD met sensors for
   climatological face-validity of the heat indices.
5. **Pipeline consistency** — the fdir estimators run through the exact
   earthkit open-data pipeline of `examples/compute-thermal-indices.py` on a
   live global field, with physical-bounds gates.

**Pre-registration.** Every script hard-codes its acceptance criteria as
constants at the top (chosen from the standards' tolerances, reference print
resolutions, or published model-skill benchmarks — each README says which)
and exits non-zero on failure; the verdicts above are the scripts' own
output, not post-hoc judgement.

**Metrics** (`vlib/metrics.py`): error = predicted − observed;
MBE = mean error (bias); MAE; RMSE; rMBE/rRMSE = normalised by the observed
mean (%); R² = squared Pearson correlation. N is always reported.

## Data sources

| Dataset | Use | Access |
|---|---|---|
| NOAA SURFRAD (GML) | radiation ground truth; met observations | public HTTPS, no auth; credit NOAA GML |
| ECMWF IFS open data | global `ssrd` field (0.25°) | public, CC BY 4.0, attribute ECMWF |
| ISO 7730:2005 Annex D table | PMV/PPD conformance | via the MIT-licensed validation-data-comfort-models repository |
| ASHRAE Global Thermal Comfort Database II v2.1 | PMV field validation | public (Dryad doi:10.6078/D1F671 / GitHub) |
| NREL SPA (via pvlib) | Earth–Sun distance truth | computed locally |

Downloads are cached under `validation/_cache/` (git-ignored, ~200 MB fully
warmed). Generated `results/*.csv` and `plots/*.png` **are** versioned
(explicit `.gitignore` exceptions) so reviewers can read the evidence without
running anything.

## Running

```bash
pip install -e ".[validation]"          # earthkit-*, pvlib, pythermalcomfort, pandas, matplotlib
python validation/earth_sun_distance_factor/validate_esdf.py
python validation/pmv_ppd/validate_pmv_ppd.py
python validation/fdir/validate_fdir_surfrad.py
python validation/fdir/validate_fdir_global.py          # live open data: map differs run to run
python validation/discomfort_index/validate_di.py
python validation/summer_simmer_index/validate_ssi.py
python validation/relative_strain_index/validate_rsi.py
python validation/apparent_temperature_radiation/validate_at_radiation.py
```

Each script prints its evidence and ends with `RESULT: PASS`/`FAIL`
(exit code 0/1). All numbers quoted in the READMEs come from these runs; the
only run-to-run variation is the live open-data field in
`validate_fdir_global.py` (its physical gates must hold for any field).

## What this campaign does *not* claim

- It does not upgrade estimation-grade code: the fdir estimators remain
  clearly-labelled approximations with ~10–20 % errors (now *measured*, not
  guessed), and `.approximations.` namespacing keeps that visible at call
  sites.
- It does not validate Fanger's PMV *model* against votes — that is a
  75-year scientific literature; it validates thermofeel's *implementation*
  and reproduces the model's documented field behaviour on 58 500 records.
- It does not resolve the Summer Simmer Index's provenance; it bounds the
  ambiguity (±2–4 °F between secondary traditions) and documents it.
- The observational demonstrations for DI/RSI/AT-rad are face-validity
  studies (band occupancy, diurnal structure), not epidemiological
  validation — those indices have no instrument-measurable truth.

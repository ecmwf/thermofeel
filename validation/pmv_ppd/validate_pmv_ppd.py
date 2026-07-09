# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Validation of ``thermofeel.calculate_pmv`` and ``calculate_ppd``.

Three parts:

A. *Standard conformance*: the ISO 7730:2005 Annex D validation table
   (12 cases), as published in the standard and reproduced in the
   validation-data-comfort-models repository (Tartarini), with the source's
   stated tolerances |dPMV| <= 0.1 and |dPPD| <= 1.

B. *Independent implementation*: element-wise comparison against
   pythermalcomfort (CBE Berkeley, Tartarini & Schiavon 2020,
   https://doi.org/10.1016/j.softx.2020.100578) ``pmv_ppd_iso`` (ISO
   7730:2005 model, limits and rounding disabled) over a 46 800-point
   factorial grid spanning the ISO validity envelope, plus the PPD closed
   form fed with identical PMV.

C. *Field data*: PMV computed from ~50k measured indoor environments of the
   ASHRAE Global Thermal Comfort Database II (Foldvary Licina et al. 2018,
   https://doi.org/10.1016/j.buildenv.2018.06.022; data
   https://doi.org/10.6078/D1F671) against the occupants' recorded thermal
   sensation votes -- the classic Fanger scatter. Informational (no
   pass/fail): PMV-vs-vote discrepancies in field data are a documented
   property of the model itself, not of an implementation.

Run:  .venv-val/bin/python validation/pmv_ppd/validate_pmv_ppd.py
"""

import gzip
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[0]))

from vlib import metrics, style  # noqa: E402

import thermofeel as tf  # noqa: E402

CACHE = HERE.parents[0] / "_cache"
ISO_JSON_URL = (
    "https://raw.githubusercontent.com/FedericoTartarini/"
    "validation-data-comfort-models/main/ts_pmv_ppd.json"
)
DB2_URL = (
    "https://raw.githubusercontent.com/CenterForTheBuiltEnvironment/"
    "ashrae-db-II/master/v2.1.0/db_measurements_v2.1.0.csv.gz"
)

# Pre-registered acceptance criteria (see README.md)
TOL_ISO_PMV = 0.1  # A: tolerance stated by the validation table
TOL_ISO_PPD = 1.0
TOL_IMPL_PMV = 0.05  # B: max |dPMV| vs pythermalcomfort
TOL_IMPL_PPD = 1.0  # B: max |dPPD| vs pythermalcomfort
TOL_PPD_FORMULA = 1e-9  # B: PPD closed form on identical PMV


def _download(url: str, dest: Path) -> Path:
    import requests

    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"  downloading {url}")
    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(resp.content)
    return dest


def _tf_pmv(tdb, tr, vr, rh, met, clo):
    """thermofeel PMV from Celsius/relative-humidity inputs."""
    return tf.calculate_pmv(
        t2_k=tf.celsius_to_kelvin(np.asarray(tdb, dtype=float)),
        mrt_k=tf.celsius_to_kelvin(np.asarray(tr, dtype=float)),
        var=np.asarray(vr, dtype=float),
        rh=np.asarray(rh, dtype=float),
        met=np.asarray(met, dtype=float),
        clo=np.asarray(clo, dtype=float),
    )


def part_a_iso_table(plt):
    """ISO 7730:2005 Annex D conformance."""
    data = json.loads(_download(ISO_JSON_URL, CACHE / "ts_pmv_ppd.json").read_text())
    rows = [r for r in data["data"] if r.get("source") == "ISO 7730 2005"]
    assert len(rows) == 12, f"expected the 12 ISO rows, got {len(rows)}"

    recs = []
    for r in rows:
        i, o = r["inputs"], r["outputs"]
        pmv = float(_tf_pmv(i["tdb"], i["tr"], i["vr"], i["rh"], i["met"], i["clo"]))
        ppd = float(tf.calculate_ppd(pmv))
        recs.append(
            {
                **{k: i[k] for k in ("tdb", "tr", "rh", "vr", "met", "clo")},
                "pmv_iso": o["pmv"],
                "ppd_iso": o["ppd"],
                "pmv_tf": round(pmv, 3),
                "ppd_tf": round(ppd, 2),
                "dpmv": round(pmv - o["pmv"], 3),
                "dppd": round(ppd - o["ppd"], 2),
            }
        )
    df = pd.DataFrame(recs)
    metrics.write_results(df, HERE / "results" / "pmv_iso7730_annexd.csv")

    ok = (df.dpmv.abs() <= TOL_ISO_PMV) & (df.dppd.abs() <= TOL_ISO_PPD)
    print(
        f"part A: {int(ok.sum())}/12 ISO cases within tolerance "
        f"(max |dPMV| {df.dpmv.abs().max():.3f}, max |dPPD| {df.dppd.abs().max():.2f})"
    )

    fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.2))
    x = np.arange(len(df))
    axs[0].errorbar(
        x,
        df.pmv_iso,
        yerr=TOL_ISO_PMV,
        fmt="s",
        ms=4,
        color="k",
        label="ISO 7730 Annex D (tol. ±0.1)",
        capsize=2,
        lw=0.8,
    )
    axs[0].plot(x, df.pmv_tf, "o", ms=3.5, color="tab:red", label="thermofeel")
    axs[0].set_xlabel("Annex D case")
    axs[0].set_ylabel("PMV")
    axs[0].legend()
    axs[1].errorbar(
        x,
        df.ppd_iso,
        yerr=TOL_ISO_PPD,
        fmt="s",
        ms=4,
        color="k",
        label="ISO 7730 Annex D (tol. ±1)",
        capsize=2,
        lw=0.8,
    )
    axs[1].plot(x, df.ppd_tf, "o", ms=3.5, color="tab:red", label="thermofeel")
    axs[1].set_xlabel("Annex D case")
    axs[1].set_ylabel("PPD [%]")
    axs[1].legend()
    fig.suptitle("PMV / PPD vs the ISO 7730:2005 Annex D validation table")
    style.save_png(fig, HERE / "plots" / "pmv_iso7730_annexd.png")
    return bool(ok.all())


def part_b_pythermalcomfort(plt):
    """Element-wise agreement with pythermalcomfort over the ISO envelope."""
    from pythermalcomfort.models import pmv_ppd_iso

    tdb, dtr, vr, rh, met, clo = np.meshgrid(
        np.arange(10.0, 35.0, 2.0),
        np.array([-6.0, -3.0, 0.0, 3.0, 6.0]),
        np.array([0.05, 0.1, 0.2, 0.4, 0.8, 1.2]),
        np.array([20.0, 40.0, 60.0, 80.0]),
        np.array([0.8, 1.0, 1.2, 1.6, 2.0, 2.4]),
        np.array([0.3, 0.5, 0.8, 1.0, 1.5]),
        indexing="ij",
    )
    tdb, vr, rh, met, clo = (a.ravel() for a in (tdb, vr, rh, met, clo))
    tr = tdb + dtr.ravel()

    ours = _tf_pmv(tdb, tr, vr, rh, met, clo)
    ref = pmv_ppd_iso(
        tdb=tdb,
        tr=tr,
        vr=vr,
        rh=rh,
        met=met,
        clo=clo,
        model="7730-2005",
        limit_inputs=False,
        round_output=False,
    )
    dpmv = ours - ref.pmv
    ppd_ours = tf.calculate_ppd(ours)
    dppd = ppd_ours - ref.ppd
    # PPD closed form isolated: identical PMV into both implementations
    dppd_formula = tf.calculate_ppd(ref.pmv) - ref.ppd

    stats = [
        (
            {"comparison": "PMV grid vs pythermalcomfort (N=46800)"},
            metrics.error_stats(ours, ref.pmv),
        ),
        (
            {"comparison": "PPD grid vs pythermalcomfort"},
            metrics.error_stats(ppd_ours, ref.ppd),
        ),
    ]
    print(
        f"part B: grid max |dPMV| = {np.nanmax(np.abs(dpmv)):.2e}, "
        f"max |dPPD| = {np.nanmax(np.abs(dppd)):.2e}, "
        f"PPD-formula max |d| = {np.nanmax(np.abs(dppd_formula)):.2e}, "
        f"non-converged: {int(np.isnan(ours).sum())}"
    )

    fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.0))
    axs[0].hist(dpmv[np.isfinite(dpmv)], bins=80, color="tab:blue")
    axs[0].set_yscale("log")
    axs[0].set_xlabel("PMV(thermofeel) − PMV(pythermalcomfort)")
    axs[0].set_ylabel("count")
    axs[0].set_title(f"N = {dpmv.size}, max |Δ| = {np.nanmax(np.abs(dpmv)):.1e}")
    sc = axs[1].scatter(
        ref.pmv[::37], dpmv[::37], s=2, c=met[::37], cmap="plasma", alpha=0.6
    )
    axs[1].set_xlabel("PMV (pythermalcomfort)")
    axs[1].set_ylabel("ΔPMV")
    fig.colorbar(sc, ax=axs[1], label="met")
    fig.suptitle("Implementation agreement over the ISO envelope (factorial grid)")
    style.save_png(fig, HERE / "plots" / "pmv_vs_pythermalcomfort.png")

    ok = (
        np.nanmax(np.abs(dpmv)) <= TOL_IMPL_PMV
        and np.nanmax(np.abs(dppd)) <= TOL_IMPL_PPD
        and np.nanmax(np.abs(dppd_formula)) <= TOL_PPD_FORMULA
        and int(np.isnan(ours).sum()) == 0
    )
    return ok, stats


def part_c_ashrae_db2(plt):
    """Field data: ASHRAE Global Thermal Comfort Database II."""
    from pythermalcomfort.models import pmv_ppd_iso
    from pythermalcomfort.utilities import v_relative

    raw = _download(DB2_URL, CACHE / "db_measurements_v2.1.0.csv.gz").read_bytes()
    cols = ["ta", "tr", "vel", "rh", "met", "clo", "thermal_sensation"]
    df = pd.read_csv(io.BytesIO(gzip.decompress(raw)), usecols=cols, low_memory=False)
    n_total = len(df)
    df = df.dropna()
    # plausibility screen (instrument/entry errors); ISO-envelope-ish bounds
    df = df[
        df.ta.between(10, 40)
        & df.tr.between(5, 60)
        & df.vel.between(0, 4)
        & df.rh.between(5, 100)
        & df.met.between(0.6, 4)
        & df.clo.between(0, 3)
        & df.thermal_sensation.between(-3, 3)
    ]
    print(f"part C: DB II records: {n_total} total -> {len(df)} complete+plausible")

    # ISO 7730 relative air velocity (body-movement correction above 1 met)
    vr = v_relative(df.vel.to_numpy(), df.met.to_numpy())
    args = dict(
        tdb=df.ta.to_numpy(),
        tr=df.tr.to_numpy(),
        vr=vr,
        rh=df.rh.to_numpy(),
        met=df.met.to_numpy(),
        clo=df.clo.to_numpy(),
    )
    pmv_tf = _tf_pmv(
        args["tdb"], args["tr"], args["vr"], args["rh"], args["met"], args["clo"]
    )
    ref = pmv_ppd_iso(**args, model="7730-2005", limit_inputs=False, round_output=False)
    dpmv = pmv_tf - ref.pmv
    impl_stats = metrics.error_stats(pmv_tf, ref.pmv)
    print(
        f"part C: implementation agreement on field inputs: "
        f"max |dPMV| = {np.nanmax(np.abs(dpmv)):.2e} over N={len(df)}"
    )

    # classic Fanger diagnostics: bin computed PMV, compare with votes
    tsv = df.thermal_sensation.to_numpy()
    bins = np.arange(-3.25, 3.26, 0.5)
    ib = np.digitize(pmv_tf, bins) - 1
    rows = []
    for b in range(len(bins) - 1):
        sel = ib == b
        if sel.sum() < 100:  # need a stable bin mean
            continue
        centre = 0.5 * (bins[b] + bins[b + 1])
        dissat = np.mean(np.abs(tsv[sel]) >= 2) * 100.0
        rows.append(
            {
                "pmv_bin": centre,
                "N": int(sel.sum()),
                "mean_tsv": float(np.mean(tsv[sel])),
                "empirical_dissatisfied_pct": float(dissat),
                "ppd_fanger": float(tf.calculate_ppd(centre)),
            }
        )
    binned = pd.DataFrame(rows)
    r_bin = np.corrcoef(binned.pmv_bin, binned.mean_tsv)[0, 1]
    print(f"part C: bin-mean PMV vs bin-mean vote correlation r = {r_bin:.3f}")

    fig, axs = plt.subplots(1, 2, figsize=(7.8, 3.3))
    axs[0].plot([-3, 3], [-3, 3], "k--", lw=0.8, label="1:1")
    axs[0].scatter(
        binned.pmv_bin,
        binned.mean_tsv,
        s=np.sqrt(binned.N),
        color="tab:red",
        label="DB II bin mean (size ∝ √N)",
    )
    axs[0].set_xlabel("computed PMV (thermofeel), bin centre")
    axs[0].set_ylabel("mean observed thermal sensation vote")
    axs[0].set_title(f"PMV vs field votes (r = {r_bin:.3f})")
    axs[0].legend()
    pmv_c = np.linspace(-3, 3, 121)
    axs[1].plot(
        pmv_c,
        tf.calculate_ppd(pmv_c),
        color="k",
        label="PPD(PMV), ISO 7730 / thermofeel",
    )
    axs[1].scatter(
        binned.pmv_bin,
        binned.empirical_dissatisfied_pct,
        s=np.sqrt(binned.N),
        color="tab:red",
        label="DB II: % votes with |TSV| ≥ 2",
    )
    axs[1].set_xlabel("computed PMV (thermofeel), bin centre")
    axs[1].set_ylabel("dissatisfied [%]")
    axs[1].set_title("PPD curve vs field dissatisfaction")
    axs[1].legend()
    fig.suptitle(
        f"ASHRAE Global Thermal Comfort Database II (N = {len(df)} measurements)"
    )
    style.save_png(fig, HERE / "plots" / "pmv_field_ashrae_db2.png")

    metrics.write_results(binned, HERE / "results" / "pmv_field_ashrae_db2.csv")
    ok = np.nanmax(np.abs(dpmv)) <= TOL_IMPL_PMV
    return ok, [
        ({"comparison": f"PMV DB-II vs pythermalcomfort (N={len(df)})"}, impl_stats)
    ]


def main() -> int:
    plt = style.use_style()
    failures = []

    if not part_a_iso_table(plt):
        failures.append("A_iso7730_table")
    ok_b, stats_b = part_b_pythermalcomfort(plt)
    if not ok_b:
        failures.append("B_pythermalcomfort_grid")
    ok_c, stats_c = part_c_ashrae_db2(plt)
    if not ok_c:
        failures.append("C_db2_implementation")

    df = metrics.stats_frame(stats_b + stats_c)
    metrics.write_results(df, HERE / "results" / "pmv_agreement_stats.csv")

    print("\nRESULT:", "FAIL " + ",".join(failures) if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

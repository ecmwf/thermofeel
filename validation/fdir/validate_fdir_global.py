# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Global physical-consistency check of the ``fdir`` estimators (earthkit stack).

This leg exercises the estimators exactly as ``examples/compute-thermal-indices.py``
uses them operationally: a real IFS 0.25 deg open-data forecast field of
accumulated ``ssrd`` (public open data carries ``ssrd`` but *not* ``fdir``,
which is precisely the gap the approximations fill), converted to a step-mean
flux, with the step-mean cosine of the solar zenith angle from
``earthkit.meteo.solar.cos_solar_zenith_angle_integrated``.

Because open data has no ``fdir`` truth, this is a *consistency* check, not an
accuracy measurement (that is what the SURFRAD leg is for). Pre-registered
gates on the full global field:

- G1 bounds: estimates are finite and 0 <= fdir <= ssrd everywhere;
- G2 night: fdir == 0 wherever the step-mean cossza <= 0.065;
- G3 cross-model agreement (daytime, ssrd > 10 W m-2): Erbs and DISC estimates
  correlate with R^2 >= 0.95 and pooled |MBE| <= 15 W m-2.

Run:  .venv-val/bin/python validation/fdir/validate_fdir_global.py
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[0]))

from vlib import metrics, style  # noqa: E402

from thermofeel.approximations import (  # noqa: E402
    approximate_fdir_disc,
    approximate_fdir_erbs,
)

CACHE = HERE.parents[0] / "_cache" / "opendata"
STEP = 12  # hours; window [base, base+12] spans a full daylight hemisphere

# Pre-registered acceptance criteria (see README.md)
TOL_BOUNDS = 1e-9  # W m-2 slack on 0 <= fdir <= ssrd
TOL_XMODEL_R2 = 0.95
TOL_XMODEL_MBE = 15.0  # W m-2


def fetch():
    """Latest 00z IFS open-data ssrd+sp at STEP.

    Cached per calendar day: re-running on the same day reuses the GRIB (so
    the numbers are reproducible within a day); a later day fetches the then-
    latest run. Delete the cache file to force a refresh.
    """
    import datetime as dt

    from ecmwf.opendata import Client

    CACHE.mkdir(parents=True, exist_ok=True)
    target = CACHE / f"ssrd_sp_step{STEP}_{dt.date.today():%Y%m%d}.grib2"
    if not target.exists():
        Client("ecmwf").retrieve(
            type="fc",
            stream="oper",
            time=0,
            step=STEP,
            param=["ssrd", "sp"],
            target=str(target),
        )
    return target


def main() -> int:
    import earthkit.data as ekd
    from earthkit.meteo import solar

    plt = style.use_style()
    failures = []

    # index by shortName metadata (FieldList.sel is unreliable across readers;
    # same approach as examples/compute-thermal-indices.py)
    fields = {
        f.metadata("shortName"): f
        for f in ekd.from_source("file", str(fetch())).to_fieldlist()
    }
    ssrd_f = fields["ssrd"]
    sp_f = fields["sp"]
    base = ssrd_f.metadata("base_datetime")
    valid = ssrd_f.metadata("valid_datetime")
    print(f"IFS open data: base {base}, valid {valid} (step {STEP} h)")

    lat, lon = ssrd_f.geography.latlons()
    ssrd = ssrd_f.to_numpy(flatten=False) / (STEP * 3600.0)  # J m-2 -> W m-2
    p_hpa = sp_f.to_numpy(flatten=False) / 100.0

    cossza = solar.cos_solar_zenith_angle_integrated(base, valid, lat, lon)
    doy = float(valid.timetuple().tm_yday)

    fdir_erbs = approximate_fdir_erbs(ssrd, cossza, doy=doy)
    fdir_disc = approximate_fdir_disc(ssrd, cossza, doy, pressure_hpa=p_hpa)

    # ------------------------------------------------------------------ gates
    for name, arr in (("erbs", fdir_erbs), ("disc", fdir_disc)):
        finite = bool(np.isfinite(arr).all())
        lo = float(arr.min())
        over = float((arr - np.maximum(ssrd, 0.0)).max())
        night = (
            float(np.abs(arr[cossza <= 0.065]).max())
            if (cossza <= 0.065).any()
            else 0.0
        )
        print(
            f"G1/G2 {name}: finite={finite}, min={lo:.3g}, "
            f"max(fdir-ssrd)={over:.3g}, max|night|={night:.3g} W m-2"
        )
        if not finite or lo < -TOL_BOUNDS or over > TOL_BOUNDS:
            failures.append(f"G1_bounds_{name}")
        if night > TOL_BOUNDS:
            failures.append(f"G2_night_{name}")

    day = (cossza > 0.065) & (ssrd > 10.0)
    stats = metrics.error_stats(fdir_disc[day], fdir_erbs[day])
    print(
        f"G3 cross-model (N={stats['N']}): R2={stats['R2']:.4f}, "
        f"MBE={stats['MBE']:+.2f} W m-2, RMSE={stats['RMSE']:.2f} W m-2"
    )
    if stats["R2"] < TOL_XMODEL_R2 or abs(stats["MBE"]) > TOL_XMODEL_MBE:
        failures.append("G3_cross_model")
    table = metrics.stats_frame(
        [({"comparison": f"DISC vs Erbs, daytime, {valid}"}, stats)]
    )
    metrics.write_results(table, HERE / "results" / "fdir_global_stats.csv")

    # ------------------------------------------------------------------- maps
    import earthkit.plots as ekp

    for name, arr in (("erbs", fdir_erbs), ("disc", fdir_disc)):
        chart = ekp.Map()
        chart.contourf(
            x=lon,
            y=lat,
            z=arr,
            levels=np.linspace(0.0, 800.0, 17),
            extend="max",
        )
        chart.coastlines()
        chart.gridlines()
        chart.title(
            f"fdir estimated from IFS open-data ssrd, {name.upper()} "
            f"(step-mean {base:%Y-%m-%d %HZ} +{STEP}h)"
        )
        out = HERE / "plots" / f"fdir_global_{name}.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(out))
        print(f"  plot: {out} ({out.stat().st_size / 1024:.0f} KB)")

    fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.3))
    hb = axs[0].hexbin(
        fdir_erbs[day], fdir_disc[day], gridsize=60, bins="log", cmap="viridis"
    )
    axs[0].plot([0, 900], [0, 900], "r--", lw=0.9)
    axs[0].set_xlabel("fdir Erbs [W m$^{-2}$]")
    axs[0].set_ylabel("fdir DISC [W m$^{-2}$]")
    axs[0].set_title(f"daytime cells, R$^2$ = {stats['R2']:.3f}")
    fig.colorbar(hb, ax=axs[0], label="log$_{10}$ count")
    d = (fdir_disc - fdir_erbs)[day]
    axs[1].hist(d, bins=80, color="tab:blue")
    axs[1].set_yscale("log")
    axs[1].set_xlabel("DISC − Erbs [W m$^{-2}$]")
    axs[1].set_title(f"MBE {stats['MBE']:+.1f}, RMSE {stats['RMSE']:.1f} W m$^{{-2}}$")
    fig.suptitle(f"Cross-model agreement on the global IFS field, valid {valid}")
    style.save_png(fig, HERE / "plots" / "fdir_global_cross_model.png")

    print("\nRESULT:", "FAIL " + ",".join(failures) if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

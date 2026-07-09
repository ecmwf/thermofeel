# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Observational validation of ``approximate_fdir_erbs`` / ``approximate_fdir_disc``.

Ground truth: NOAA SURFRAD 1-minute observations (pyranometer global ``dw_solar``,
pyrheliometer direct-normal ``direct_n``, shaded diffuse), three climatically
contrasting stations x four seasonal months of 2023:

- Desert Rock, NV  (arid, high direct fraction)
- Bondville, IL    (continental, mixed skies)
- Goodwin Creek, MS (humid subtropical)

The models estimate the direct *horizontal* radiation (ECMWF ``fdir``) from the
global horizontal radiation alone; the reference is the independently measured
``direct_n * cos(zenith)``. Both models are *hourly* correlations (Erbs et al.
1982; Maxwell 1987), so the pass/fail gate is evaluated on hourly means; the
1-minute and 3-hourly (ECMWF-step-like) aggregations are reported alongside.

QC follows BSRN-style practice: NOAA QC flags == 0, zenith < 85 deg, GHI >= 5
W m-2, and a three-component closure check on GHI = diffuse + direct_n cos(z)
(8% for zenith < 75 deg, 15% for 75-85 deg, applied where GHI > 50 W m-2).

Pre-registered acceptance (per station, hourly, all-sky; thresholds calibrated
from the published skill of GHI-only separation models, see README):
rRMSE <= 40%, |rMBE| <= 15%, R^2 >= 0.85.

Run:  .venv-val/bin/python validation/fdir/validate_fdir_surfrad.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[0]))

from vlib import metrics, style, surfrad  # noqa: E402

from thermofeel.approximations import (  # noqa: E402
    SOLAR_CONSTANT,
    approximate_fdir_disc,
    approximate_fdir_erbs,
    earth_sun_distance_factor,
)

STATIONS = ("bon", "dra", "gwn")
MONTHS = (1, 4, 7, 10)
YEAR = 2023

# Pre-registered acceptance criteria (hourly, all-sky, per station)
TOL_RRMSE_PCT = 40.0
TOL_ABS_RMBE_PCT = 15.0
TOL_R2 = 0.85


def load_observations() -> pd.DataFrame:
    frames = []
    for sta in STATIONS:
        days = []
        for m in MONTHS:
            days += surfrad.month_days(YEAR, m)
        frames.append(surfrad.load(sta, days))
    df = pd.concat(frames, ignore_index=True)

    df["cossza"] = np.cos(np.radians(df.zen))
    df["fdir_obs"] = df.direct_n * df.cossza.clip(lower=0.0)
    df["doy"] = df.time.dt.dayofyear.astype(float)
    # station pressure [mb == hPa]; fall back to the standard atmosphere
    p_std = 1013.25 * (1.0 - 2.25577e-5 * df.elev) ** 5.25588
    df["p_hpa"] = df.pressure.fillna(p_std)

    # ------------------------------------------------------------------- QC
    ok = (
        np.isfinite(df.dw_solar)
        & np.isfinite(df.direct_n)
        & np.isfinite(df.diffuse)
        & (df.zen < 85.0)
        & (df.dw_solar >= 5.0)
    )
    closure = df.diffuse + df.fdir_obs - df.dw_solar
    rel = np.abs(closure) / df.dw_solar
    tight = (df.zen < 75.0) & (df.dw_solar > 50.0)
    loose = (df.zen >= 75.0) & (df.dw_solar > 50.0)
    ok &= np.where(tight, rel < 0.08, True)
    ok &= np.where(loose, rel < 0.15, True)
    n0 = len(df)
    df = df[ok].copy()
    print(
        f"QC: {n0} 1-min records -> {len(df)} valid daytime ({100 * len(df) / n0:.0f}%)"
    )

    e0 = earth_sun_distance_factor(df.doy.to_numpy())
    df["kt"] = np.clip(
        df.dw_solar / (SOLAR_CONSTANT * e0 * np.maximum(df.cossza, 0.065)), 0, 1
    )
    return df


def predict(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["fdir_erbs"] = approximate_fdir_erbs(
        df.dw_solar.to_numpy(), df.cossza.to_numpy(), doy=df.doy.to_numpy()
    )
    df["fdir_disc"] = approximate_fdir_disc(
        df.dw_solar.to_numpy(),
        df.cossza.to_numpy(),
        df.doy.to_numpy(),
        pressure_hpa=df.p_hpa.to_numpy(),
    )
    return df


def aggregate(df: pd.DataFrame, hours: int, min_minutes: int) -> pd.DataFrame:
    """Mean over `hours`-long windows of QC-valid minutes, then re-predict."""
    key = [df.station, df.time.dt.floor(f"{hours}h")]
    g = (
        df.groupby(key)
        .agg(
            n=("dw_solar", "size"),
            dw_solar=("dw_solar", "mean"),
            cossza=("cossza", "mean"),
            doy=("doy", "mean"),
            p_hpa=("p_hpa", "mean"),
            fdir_obs=("fdir_obs", "mean"),
            kt=("kt", "mean"),
        )
        .reset_index()
    )
    g = g[g.n >= min_minutes]
    return predict(g)


def cross_check_pvlib(df: pd.DataFrame) -> None:
    """Implementation fidelity vs pvlib on the same observations."""
    import pvlib

    zen = np.degrees(np.arccos(df.cossza.to_numpy()))
    ghi = df.dw_solar.to_numpy()
    doy = df.doy.to_numpy()

    # pvlib erbs returns dhi; direct horizontal = ghi - dhi. Match thermofeel's
    # solar constant by scaling: pvlib uses its own I0 via get_extra_radiation.
    out = pvlib.irradiance.erbs(ghi, zen, doy.astype(int))
    ref_erbs = np.clip(ghi - out["dhi"], 0, ghi)
    ours_erbs = approximate_fdir_erbs(
        ghi, df.cossza.to_numpy(), doy=doy, solar_constant=1366.1
    )
    d_erbs = np.nanmax(np.abs(ours_erbs - ref_erbs))

    disc = pvlib.irradiance.disc(ghi, zen, doy.astype(int), pressure=101325.0)
    ref_disc = np.clip(disc["dni"] * df.cossza.to_numpy(), 0, np.maximum(ghi, 0))
    ours_disc = approximate_fdir_disc(
        ghi,
        df.cossza.to_numpy(),
        doy,
        pressure_hpa=1013.25,
        solar_constant=1370.0,
    )
    d_disc = np.nanmax(np.abs(ours_disc - ref_disc))
    print(
        f"pvlib cross-check on SURFRAD inputs (matched constants): "
        f"max |d erbs| = {d_erbs:.2e} W m-2, max |d disc| = {d_disc:.2e} W m-2"
    )
    if d_erbs > 1e-8 or d_disc > 1e-8:
        raise AssertionError("pvlib cross-check failed")


def score(df: pd.DataFrame, label: str) -> list:
    rows = []
    for model in ("erbs", "disc"):
        for sta, g in df.groupby("station", observed=True):
            rows.append(
                (
                    {"resolution": label, "model": model, "station": sta},
                    metrics.error_stats(g[f"fdir_{model}"], g.fdir_obs),
                )
            )
        rows.append(
            (
                {"resolution": label, "model": model, "station": "ALL"},
                metrics.error_stats(df[f"fdir_{model}"], df.fdir_obs),
            )
        )
    return rows


def plot_scatter(plt, hourly: pd.DataFrame) -> None:
    fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.6), sharey=True)
    lim = 900.0
    for ax, model in zip(axs, ("erbs", "disc")):
        s = metrics.error_stats(hourly[f"fdir_{model}"], hourly.fdir_obs)
        hb = ax.hexbin(
            hourly.fdir_obs,
            hourly[f"fdir_{model}"],
            gridsize=60,
            bins="log",
            cmap="viridis",
            extent=(0, lim, 0, lim),
        )
        ax.plot([0, lim], [0, lim], "r--", lw=0.9)
        ax.set_xlabel("observed fdir [W m$^{-2}$]")
        ax.set_title(
            f"{model.upper()}  MBE {s['MBE']:+.0f}, RMSE {s['RMSE']:.0f} W m$^{{-2}}$"
            f" ({s['rRMSE_pct']:.0f}%), R$^2$ {s['R2']:.2f}"
        )
    axs[0].set_ylabel("estimated fdir [W m$^{-2}$]")
    fig.colorbar(hb, ax=axs, label="log$_{10}$ count")
    fig.suptitle(
        f"Hourly direct horizontal radiation vs SURFRAD, {YEAR} "
        f"(3 stations, 4 months, N = {len(hourly)})"
    )
    style.save_png(fig, HERE / "plots" / "fdir_surfrad_scatter_hourly.png")


def plot_kt_bias(plt, hourly: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    edges = np.arange(0.05, 0.90, 0.05)
    centres = 0.5 * (edges[:-1] + edges[1:])
    for model, color in (("erbs", "tab:blue"), ("disc", "tab:orange")):
        err = hourly[f"fdir_{model}"] - hourly.fdir_obs
        mbe, rmse = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = (hourly.kt >= lo) & (hourly.kt < hi)
            mbe.append(err[sel].mean() if sel.sum() >= 30 else np.nan)
            rmse.append(np.sqrt((err[sel] ** 2).mean()) if sel.sum() >= 30 else np.nan)
        ax.plot(centres, mbe, "-o", ms=3, color=color, label=f"{model.upper()} MBE")
        ax.plot(centres, rmse, "--", color=color, label=f"{model.upper()} RMSE")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xlabel("clearness index $k_t$ (hourly)")
    ax.set_ylabel("error [W m$^{-2}$]")
    ax.set_title("Error vs sky condition (hourly, all stations pooled)")
    ax.legend(ncols=2)
    style.save_png(fig, HERE / "plots" / "fdir_surfrad_kt_bias.png")


def plot_example_days(plt, df: pd.DataFrame) -> None:
    bon = df[(df.station == "bon") & (df.time.dt.month == 7)].copy()
    daily = bon.groupby(bon.time.dt.date).agg(kt=("kt", "mean"), n=("kt", "size"))
    daily = daily[daily.n >= 500]
    clear = daily.kt.idxmax()
    cloudy = (daily.kt - 0.45).abs().idxmin()

    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.2), sharey=True)
    for ax, day, name in ((axs[0], clear, "clear"), (axs[1], cloudy, "mixed/cloudy")):
        d = bon[bon.time.dt.date == day]
        t = d.time.dt.hour + d.time.dt.minute / 60.0
        ax.plot(t, d.dw_solar, color="0.6", lw=0.8, label="GHI (ssrd) obs")
        ax.plot(t, d.fdir_obs, color="k", lw=1.2, label="fdir obs")
        ax.plot(t, d.fdir_erbs, color="tab:blue", lw=1.0, label="Erbs")
        ax.plot(t, d.fdir_disc, color="tab:orange", lw=1.0, label="DISC")
        ax.set_title(f"Bondville {day} ({name}, $k_t$ = {daily.kt[day]:.2f})")
        ax.set_xlabel("hour (UTC)")
    axs[0].set_ylabel("irradiance [W m$^{-2}$]")
    axs[0].legend()
    fig.suptitle("Example days: 1-minute observations vs GHI-only estimates")
    style.save_png(fig, HERE / "plots" / "fdir_surfrad_example_days.png")


def main() -> int:
    plt = style.use_style()
    df = predict(load_observations())
    cross_check_pvlib(df)

    hourly = aggregate(df, 1, min_minutes=50)
    threeh = aggregate(df, 3, min_minutes=150)
    print(f"aggregations: {len(hourly)} hourly windows, {len(threeh)} 3-hourly windows")

    rows = score(df, "1-min") + score(hourly, "hourly") + score(threeh, "3-hourly")
    table = metrics.stats_frame(rows)
    metrics.write_results(table, HERE / "results" / "fdir_surfrad_stats.csv")

    # DISC solar-constant sensitivity (default 1361 vs Maxwell's 1370)
    alt = approximate_fdir_disc(
        hourly.dw_solar.to_numpy(),
        hourly.cossza.to_numpy(),
        hourly.doy.to_numpy(),
        pressure_hpa=hourly.p_hpa.to_numpy(),
        solar_constant=1370.0,
    )
    s_def = metrics.error_stats(hourly.fdir_disc, hourly.fdir_obs)
    s_alt = metrics.error_stats(alt, hourly.fdir_obs)
    print(
        f"DISC solar-constant sensitivity (hourly, pooled): "
        f"RMSE {s_def['RMSE']:.1f} (SC=1361) vs {s_alt['RMSE']:.1f} (SC=1370) W m-2"
    )

    plot_scatter(plt, hourly)
    plot_kt_bias(plt, hourly)
    plot_example_days(plt, df)

    failures = []
    gate = table[(table.resolution == "hourly") & (table.station != "ALL")]
    for _, r in gate.iterrows():
        checks = (
            r.rRMSE_pct <= TOL_RRMSE_PCT
            and abs(r.rMBE_pct) <= TOL_ABS_RMBE_PCT
            and r.R2 >= TOL_R2
        )
        flag = "ok" if checks else "FAIL"
        print(
            f"gate {r.model:4s} {r.station}: rRMSE {r.rRMSE_pct:5.1f}% "
            f"rMBE {r.rMBE_pct:+5.1f}% R2 {r.R2:.3f} -> {flag}"
        )
        if not checks:
            failures.append(f"{r.model}:{r.station}")

    print("\nRESULT:", "FAIL " + ",".join(failures) if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

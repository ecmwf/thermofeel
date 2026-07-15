# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Validation of ``thermofeel.calculate_discomfort_index`` (Thom / Giles).

Parts:

1. *Independent implementation*: element-wise comparison against
   pythermalcomfort ``discomfort_index`` (which implements the same Giles et
   al. 1990 Celsius/relative-humidity closed form) over a dense
   temperature-humidity grid AND over real observed (T, RH) pairs.

2. *Analytic identities* of the formula: DI(T, 100%) = T for all T, and
   DI(14.5 degC, RH) = 14.5 degC for all RH.

3. *Observational demonstration*: DI computed from NOAA SURFRAD 1-minute
   observations for July 2023 at a humid-subtropical station (Goodwin Creek,
   MS) and a desert station (Desert Rock, NV), summarised as occupancy of the
   Giles discomfort categories -- a face-validity check that the index
   behaves climatologically sensibly on real data (informational).

Run:  .venv-val/bin/python validation/discomfort_index/validate_di.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[0]))

from vlib import metrics, style, surfrad  # noqa: E402

import thermofeel as tf  # noqa: E402

# Pre-registered acceptance criteria (see README.md).
# pythermalcomfort's discomfort_index rounds its output to 1 decimal place and
# offers no round_output switch, so the sharpest available oracle statement is:
# the difference must be pure quantization noise -- bounded by half the quantum
# (0.05 degC) with an RMSE matching uniform rounding error (0.1/sqrt(12) degC).
TOL_ORACLE_MAX_C = 0.05 + 1e-9  # half of the 0.1 degC output quantum
ORACLE_RMSE_BAND_C = (0.024, 0.034)  # 0.1/sqrt(12) = 0.0289 degC
TOL_IDENTITY_C = 1e-12  # analytic identities [degC]

#: Giles et al. (1990) discomfort categories [degC]
BANDS = [
    (-np.inf, 21.0, "no discomfort"),
    (21.0, 24.0, "<50% feel discomfort"),
    (24.0, 27.0, ">50% feel discomfort"),
    (27.0, 29.0, "most feel discomfort"),
    (29.0, 32.0, "everyone feels severe stress"),
    (32.0, np.inf, "state of medical emergency"),
]


def di_celsius(t_c, rh):
    """thermofeel DI evaluated in Celsius."""
    return tf.kelvin_to_celsius(
        tf.calculate_discomfort_index(tf.celsius_to_kelvin(np.asarray(t_c)), rh)
    )


def band_occupancy(di_c: np.ndarray) -> pd.DataFrame:
    rows = []
    for lo, hi, name in BANDS:
        pct = 100.0 * np.mean((di_c >= lo) & (di_c < hi))
        rows.append({"band": f"[{lo}, {hi})", "category": name, "pct": pct})
    return pd.DataFrame(rows)


def main() -> int:
    from pythermalcomfort.models import discomfort_index

    plt = style.use_style()
    failures = []

    # -------------------------------------------------- 1. independent oracle
    t_c, rh = np.meshgrid(np.arange(-10.0, 50.1, 0.5), np.arange(0.0, 100.1, 1.0))
    ours = di_celsius(t_c.ravel(), rh.ravel())
    ref = discomfort_index(tdb=t_c.ravel(), rh=rh.ravel()).di
    d_grid = np.nanmax(np.abs(ours - ref))
    print(f"oracle grid  (N={ours.size}): max |dDI| = {d_grid:.2e} degC")

    obs = pd.concat(
        [
            surfrad.load(sta, surfrad.month_days(2023, 7))[
                ["time", "temp", "rh", "station"]
            ]
            for sta in ("gwn", "dra")
        ],
        ignore_index=True,
    ).dropna()
    ours_o = di_celsius(obs.temp.to_numpy(), obs.rh.to_numpy())
    ref_o = discomfort_index(tdb=obs.temp.to_numpy(), rh=obs.rh.to_numpy()).di
    d_obs = np.nanmax(np.abs(ours_o - ref_o))
    print(f"oracle obs   (N={len(obs)}): max |dDI| = {d_obs:.2e} degC")
    # the reference rounds to 0.1 degC: the residual must be pure rounding
    # noise (max <= 0.05 degC, RMSE = 0.1/sqrt(12) = 0.0289 degC), on both
    # the grid and the observation inputs
    rmse_g = float(np.sqrt(np.nanmean((ours - ref) ** 2)))
    rmse_o = float(np.sqrt(np.nanmean((ours_o - ref_o) ** 2)))
    print(
        f"residual vs reference rounding: RMSE grid = {rmse_g:.4f}, "
        f"obs = {rmse_o:.4f} degC "
        f"(pure 0.1-quantum rounding noise = 0.0289 degC)"
    )
    rmse_ok = all(
        ORACLE_RMSE_BAND_C[0] <= r <= ORACLE_RMSE_BAND_C[1] for r in (rmse_g, rmse_o)
    )
    if max(d_grid, d_obs) > TOL_ORACLE_MAX_C or not rmse_ok:
        failures.append("oracle")

    stats = [
        (
            {"comparison": "DI vs pythermalcomfort, grid"},
            metrics.error_stats(ours, ref),
        ),
        (
            {"comparison": "DI vs pythermalcomfort, SURFRAD obs Jul 2023"},
            metrics.error_stats(ours_o, ref_o),
        ),
    ]
    metrics.write_results(
        metrics.stats_frame(stats), HERE / "results" / "di_oracle_stats.csv"
    )

    # -------------------------------------------------- 2. analytic identities
    t = np.arange(-10.0, 50.1, 0.1)
    d1 = np.max(np.abs(di_celsius(t, np.full_like(t, 100.0)) - t))
    r = np.arange(0.0, 100.1, 0.5)
    d2 = np.max(np.abs(di_celsius(np.full_like(r, 14.5), r) - 14.5))
    print(f"identity DI(T,100)=T:      max |d| = {d1:.2e} degC")
    print(f"identity DI(14.5,RH)=14.5: max |d| = {d2:.2e} degC")
    if max(d1, d2) > TOL_IDENTITY_C:
        failures.append("identities")

    # -------------------------------------------------- 3. observational demo
    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.4))
    width = 0.38
    ticks = np.arange(len(BANDS))
    tables = {}
    for k, (sta, color) in enumerate((("gwn", "tab:red"), ("dra", "tab:orange"))):
        sel = obs[obs.station == sta]
        di = di_celsius(sel.temp.to_numpy(), sel.rh.to_numpy())
        occ = band_occupancy(di)
        tables[sta] = occ
        axs[0].bar(
            ticks + (k - 0.5) * width,
            occ.pct,
            width,
            color=color,
            label=f"{sta} ({surfrad.STATIONS[sta][0]})",
        )
    axs[0].set_xticks(ticks)
    axs[0].set_xticklabels([b[2] for b in BANDS], fontsize=6.5, rotation=25, ha="right")
    axs[0].set_ylabel("share of July 2023 minutes [%]")
    axs[0].set_title("Giles category occupancy (1-min obs)")
    axs[0].legend()

    gwn = obs[obs.station == "gwn"]
    di_gwn = di_celsius(gwn.temp.to_numpy(), gwn.rh.to_numpy())
    # diurnal cycle binned by the actual UTC hour of each observation
    hh = gwn.time.dt.hour.to_numpy()
    med = [np.nanmedian(di_gwn[hh == h]) for h in range(24)]
    q25 = [np.nanpercentile(di_gwn[hh == h], 25) for h in range(24)]
    q75 = [np.nanpercentile(di_gwn[hh == h], 75) for h in range(24)]
    axs[1].fill_between(range(24), q25, q75, alpha=0.3, color="tab:red")
    axs[1].plot(range(24), med, color="tab:red", label="DI median (IQR shaded)")
    axs[1].plot(
        range(24),
        [np.nanmedian(gwn.temp.to_numpy()[hh == h]) for h in range(24)],
        color="k",
        ls="--",
        lw=0.9,
        label="air temperature median",
    )
    for y, lbl in ((24, ">50% discomfort"), (27, "most"), (29, "severe")):
        axs[1].axhline(y, color="0.5", lw=0.6)
        axs[1].text(0.2, y + 0.1, lbl, fontsize=6.5, color="0.35")
    axs[1].set_xlabel("hour (UTC)")
    axs[1].set_ylabel("[degC]")
    axs[1].set_title("Goodwin Creek MS, July 2023 diurnal cycle")
    axs[1].legend(loc="lower right")
    fig.suptitle("Thom Discomfort Index on SURFRAD observations")
    style.save_png(fig, HERE / "plots" / "di_surfrad_demo.png")

    demo = pd.concat(
        [t.assign(station=s) for s, t in tables.items()], ignore_index=True
    )
    metrics.write_results(demo, HERE / "results" / "di_band_occupancy.csv")

    print("\nRESULT:", "FAIL " + ",".join(failures) if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

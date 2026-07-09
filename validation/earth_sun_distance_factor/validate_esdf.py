# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Validation of ``thermofeel.approximations.earth_sun_distance_factor``.

Two-tier validation:

1. *Implementation fidelity*: the function implements the Spencer (1971)
   Fourier series; it must agree with the independent pvlib implementation
   (``pvlib.irradiance.get_extra_radiation(..., method="spencer")`` with a
   unit solar constant) to machine precision.

2. *Physical accuracy*: the factor is compared against a high-precision
   ephemeris truth, (1 AU / R)^2 with R from the NREL Solar Position
   Algorithm (Reda & Andreas 2004) via
   ``pvlib.solarposition.nrel_earthsun_distance`` (stated uncertainty of the
   SPA R is ~1e-7 AU), daily at 12 UTC over 2015-2025 (three leap years).

Run:  .venv-val/bin/python validation/earth_sun_distance_factor/validate_esdf.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pvlib

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))  # repo root -> import thermofeel
sys.path.insert(0, str(HERE.parents[0]))  # validation/ -> import vlib

from vlib import metrics, style  # noqa: E402

from thermofeel.approximations import earth_sun_distance_factor  # noqa: E402

SOLAR_CONSTANT = 1361.0  # W m-2, matches thermofeel.approximations

# Pre-registered acceptance criteria (see README.md)
TOL_PVLIB_REL = 1e-12  # tier 1: same series, independent code
TOL_SPA_MAX_REL_PCT = 0.5  # tier 2: max |rel. error| vs SPA ephemeris


def main() -> int:
    plt = style.use_style()
    failures = []

    # ------------------------------------------------------------------ tier 1
    doy = np.arange(1.0, 367.0)
    ours = earth_sun_distance_factor(doy)
    pvlib_spencer = pvlib.irradiance.get_extra_radiation(
        doy, solar_constant=1.0, method="spencer"
    )
    rel = np.abs(ours - pvlib_spencer) / pvlib_spencer
    print(f"tier 1 (vs pvlib Spencer): max |rel diff| = {rel.max():.3e}")
    if rel.max() > TOL_PVLIB_REL:
        failures.append("tier1_pvlib_spencer")

    # ------------------------------------------------------------------ tier 2
    times = pd.date_range("2015-01-01 12:00", "2025-12-31 12:00", freq="D", tz="UTC")
    r_au = pvlib.solarposition.nrel_earthsun_distance(times).to_numpy()
    truth = 1.0 / r_au**2
    doy2 = times.dayofyear.to_numpy(dtype=float)
    pred = earth_sun_distance_factor(doy2)

    rel_err_pct = 100.0 * (pred - truth) / truth
    toa_err = SOLAR_CONSTANT * (pred - truth)  # impact on TOA normal irradiance

    stats = metrics.error_stats(pred, truth)
    stats_rows = [
        (
            {
                "comparison": "factor vs SPA (1/R^2), daily 2015-2025",
                "max_abs_rel_err_pct": float(np.max(np.abs(rel_err_pct))),
                "mean_abs_rel_err_pct": float(np.mean(np.abs(rel_err_pct))),
                "max_abs_TOA_err_Wm2": float(np.max(np.abs(toa_err))),
            },
            stats,
        )
    ]

    # leap-year edge: 29 Feb and 31 Dec (doy 366) of the leap years
    leap = times[times.dayofyear == 366]
    leap_err = rel_err_pct[times.dayofyear == 366]
    for t, e in zip(leap, leap_err):
        print(f"  doy 366 ({t.date()}): rel err {e:+.4f} %")

    max_rel = float(np.max(np.abs(rel_err_pct)))
    print(f"tier 2 (vs SPA):  max |rel err| = {max_rel:.4f} %")
    print(f"                  mean |rel err| = {np.mean(np.abs(rel_err_pct)):.4f} %")
    print(f"                  max |TOA err| = {np.max(np.abs(toa_err)):.3f} W m-2")
    if max_rel > TOL_SPA_MAX_REL_PCT:
        failures.append("tier2_spa")

    df = metrics.stats_frame(stats_rows)
    metrics.write_results(df, HERE / "results" / "esdf_stats.csv")

    # ------------------------------------------------------------------- plots
    fig, axs = plt.subplots(2, 1, figsize=(7, 5.4), sharex=False)
    yr = times.year == 2024  # leap year, shown as the annual cycle example
    axs[0].plot(times.dayofyear[yr], truth[yr], lw=2.2, color="k", label="SPA (truth)")
    axs[0].plot(
        times.dayofyear[yr],
        pred[yr],
        lw=1.0,
        color="tab:orange",
        label="thermofeel (Spencer 1971)",
    )
    axs[0].set_xlabel("day of year (2024)")
    axs[0].set_ylabel("Earth-Sun distance factor (1 AU / R)$^2$")
    axs[0].set_title("Annual cycle, leap year 2024")
    axs[0].legend()

    sc = axs[1].scatter(times.dayofyear, rel_err_pct, c=times.year, s=3, cmap="viridis")
    axs[1].axhline(0, color="k", lw=0.6)
    axs[1].set_xlabel("day of year (all days 2015-2025)")
    axs[1].set_ylabel("relative error [%]")
    axs[1].set_title(
        f"Error vs SPA ephemeris: max {max_rel:.3f} %, "
        f"max TOA impact {np.max(np.abs(toa_err)):.2f} W m$^{{-2}}$"
    )
    fig.colorbar(sc, ax=axs[1], label="year")
    style.save_png(fig, HERE / "plots" / "esdf_vs_spa.png")

    print("\nRESULT:", "FAIL " + ",".join(failures) if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

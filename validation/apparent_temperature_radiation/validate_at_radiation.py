# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Validation of ``thermofeel.calculate_apparent_temperature_radiation``.

The radiation-inclusive Steadman (1994) Apparent Temperature is *defined* by
the formula the Australian Bureau of Meteorology publishes operationally,
AT = Ta + 0.348 e - 0.70 ws + 0.70 Q/(ws + 10) - 4.25, with e the vapour
pressure from BoM's Magnus approximation. Validation therefore rests on:

1. *Independent transcription*: an element-wise comparison against a from-
   scratch transcription of the BoM-published formula (including BoM's
   vapour-pressure approximation) over a wide input grid, plus the exact
   reduction to AT = Ta + 0.348 e - 0.70 ws - 4.25 at Q = 0.

2. *Family coherence with the non-radiative form*: solving
   AT_rad(Q*) = AT_norad gives Q*(e, ws) = ((0.25 - 0.018 e)(ws + 10))/0.70;
   the identity is verified numerically and Q* is mapped over realistic
   (e, ws) -- showing the non-radiative form corresponds to a *small* implied
   radiation load (indoor/shade conditions), as Steadman designed.

3. *Observational demonstration*: diurnal AT with and without radiation at
   Desert Rock, NV (July 2023), using the SURFRAD *site* net radiation as a
   stand-in for Q. NOTE: Steadman's Q is the net radiation absorbed per unit
   body-surface area, not the surface energy-budget net radiation, so this is
   a qualitative demonstration of the radiation term's behaviour, not a
   quantitative validation of Q itself (informational).

Run:
  .venv-val/bin/python \
      validation/apparent_temperature_radiation/validate_at_radiation.py
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[0]))

from vlib import metrics, style, surfrad  # noqa: E402

import thermofeel as tf  # noqa: E402

# Pre-registered acceptance criteria (see README.md)
TOL_TRANSCRIPTION_C = 1e-9  # vs independent BoM transcription [degC]
TOL_IDENTITY_C = 1e-9  # AT_rad(Q*) = AT_norad [degC]


def at_rad_c(t_c, ws, rh, q):
    return tf.kelvin_to_celsius(
        tf.calculate_apparent_temperature_radiation(
            tf.celsius_to_kelvin(np.asarray(t_c, dtype=float)), ws, rh, q
        )
    )


def at_norad_c(t_c, ws, rh):
    return tf.kelvin_to_celsius(
        tf.calculate_apparent_temperature(
            tf.celsius_to_kelvin(np.asarray(t_c, dtype=float)), ws, rh
        )
    )


def bom_transcription(t_c, ws, rh, q):
    """BoM 'including radiation' AT, transcribed from
    http://www.bom.gov.au/info/thermal_stress/ (Steadman 1994)."""
    t_c = np.asarray(t_c, dtype=float)
    e = rh / 100.0 * 6.105 * np.exp(17.27 * t_c / (237.7 + t_c))
    return t_c + 0.348 * e - 0.70 * ws + 0.70 * q / (ws + 10.0) - 4.25


def main() -> int:
    plt = style.use_style()
    failures = []

    t_c, ws, rh, q = np.meshgrid(
        np.arange(-10.0, 45.1, 2.5),
        np.arange(0.0, 20.1, 1.0),
        np.arange(5.0, 100.1, 5.0),
        np.arange(-100.0, 1000.1, 50.0),
        indexing="ij",
    )
    t_c, ws, rh, q = (a.ravel() for a in (t_c, ws, rh, q))

    # ------------------------------------------------- 1. transcription check
    ours = at_rad_c(t_c, ws, rh, q)
    ref = bom_transcription(t_c, ws, rh, q)
    d1 = np.nanmax(np.abs(ours - ref))
    print(f"transcription check (N={ours.size}): max |dAT| = {d1:.2e} degC")
    if d1 > TOL_TRANSCRIPTION_C:
        failures.append("transcription")
    stats = [
        (
            {"comparison": "AT_rad vs independent BoM transcription"},
            metrics.error_stats(ours, ref),
        )
    ]

    # Q = 0 reduction
    sel = q == 0.0
    d_q0 = np.nanmax(
        np.abs(
            ours[sel]
            - (
                t_c[sel]
                + 0.348
                * (
                    rh[sel]
                    / 100.0
                    * 6.105
                    * np.exp(17.27 * t_c[sel] / (237.7 + t_c[sel]))
                )
                - 0.70 * ws[sel]
                - 4.25
            )
        )
    )
    print(f"Q=0 reduction: max |d| = {d_q0:.2e} degC")
    if d_q0 > TOL_TRANSCRIPTION_C:
        failures.append("q0_reduction")

    # --------------------------------------------------- 2. family coherence
    tt, wsx, rhx = np.meshgrid(
        np.arange(0.0, 45.1, 1.0),
        np.arange(0.0, 15.1, 0.5),
        np.arange(5.0, 100.1, 5.0),
        indexing="ij",
    )
    tt, wsx, rhx = tt.ravel(), wsx.ravel(), rhx.ravel()
    e = rhx / 100.0 * 6.105 * np.exp(17.27 * tt / (237.7 + tt))
    q_star = (0.25 - 0.018 * e) * (wsx + 10.0) / 0.70
    d2 = np.nanmax(np.abs(at_rad_c(tt, wsx, rhx, q_star) - at_norad_c(tt, wsx, rhx)))
    print(f"family coherence AT_rad(Q*) = AT_norad: max |d| = {d2:.2e} degC")
    print(
        f"  implied Q* range over the grid: {q_star.min():.0f} to "
        f"{q_star.max():.0f} W m-2 (median {np.median(q_star):.0f})"
    )
    if d2 > TOL_IDENTITY_C:
        failures.append("family_coherence")

    metrics.write_results(
        metrics.stats_frame(stats), HERE / "results" / "at_radiation_stats.csv"
    )

    # ------------------------------------------------------------------ plots
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.4))

    ee = np.arange(2.0, 40.1, 0.5)
    ww = np.arange(0.0, 15.05, 0.25)
    E, W = np.meshgrid(ee, ww)
    QS = (0.25 - 0.018 * E) * (W + 10.0) / 0.70
    cf = axs[0].contourf(E, W, QS, levels=14, cmap="RdYlBu_r")
    fig.colorbar(cf, ax=axs[0], label="Q* [W m$^{-2}$]")
    axs[0].set_xlabel("vapour pressure e [hPa]")
    axs[0].set_ylabel("10 m wind speed [m s$^{-1}$]")
    axs[0].set_title("Implied radiation load Q* where AT$_{rad}$ = AT$_{no-rad}$")

    # observational demonstration at Desert Rock, July 2023
    obs = surfrad.load("dra", surfrad.month_days(2023, 7))[
        ["time", "temp", "rh", "windspd", "totalnet"]
    ].dropna()
    at_r = at_rad_c(
        obs.temp.to_numpy(),
        obs.windspd.to_numpy(),
        obs.rh.to_numpy(),
        obs.totalnet.to_numpy(),
    )
    at_n = at_norad_c(obs.temp.to_numpy(), obs.windspd.to_numpy(), obs.rh.to_numpy())
    hh = obs.time.dt.hour.to_numpy()
    med_r = [np.nanmedian(at_r[hh == h]) for h in range(24)]
    med_n = [np.nanmedian(at_n[hh == h]) for h in range(24)]
    med_t = [np.nanmedian(obs.temp.to_numpy()[hh == h]) for h in range(24)]
    axs[1].plot(range(24), med_t, "k--", lw=0.9, label="air temperature")
    axs[1].plot(range(24), med_n, color="tab:blue", label="AT (no radiation)")
    axs[1].plot(
        range(24),
        med_r,
        color="tab:red",
        label="AT (radiation), Q = SITE net radiation:\nupper-bound stand-in, "
        "body-absorbed Q is far smaller",
    )
    axs[1].set_xlabel("hour (UTC; local = UTC-8)")
    axs[1].set_ylabel("[degC]")
    axs[1].set_title("Desert Rock NV, July 2023 medians (demonstration)")
    axs[1].legend()
    fig.suptitle("Apparent Temperature with radiation (Steadman 1994 / BoM)")
    style.save_png(fig, HERE / "plots" / "at_radiation_checks.png")

    print("\nRESULT:", "FAIL " + ",".join(failures) if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Validation of ``thermofeel.calculate_relative_strain_index`` (Lee & Henschel).

thermofeel implements the hectopascal closed form stated by Asghari et al.
(2020), RSI = (Ta - 21) / (58 - e), with e from thermofeel's BoM-constant
Magnus formula. Parts:

1. *Transcription check*: an independent reimplementation using thermofeel's
   own vapour-pressure convention must agree to machine precision.

2. *Reference-formula sensitivity*: Asghari et al. state the vapour pressure
   as e = (RH/100) 6.112 x 10^(7.5 Ta / (237.7 + Ta)) hPa, while thermofeel
   uses its library-wide BoM convention e = (RH/100) 6.105 exp(17.27 Ta /
   (237.7 + Ta)). Over the index's validity envelope (21-35 degC) the two
   Magnus variants differ by well under 1%, and the induced |dRSI| must stay
   below a pre-registered bound (0.02, i.e. one order of magnitude below the
   narrowest assessment band width... band widths are 0.10).

3. *Behavioural checks*: RSI strictly increases with temperature and humidity
   inside the validity envelope; the assessment-band crossing temperatures
   (Asghari's five levels) are tabulated per RH; the documented domain edge
   (denominator -> 0 as e -> 58 hPa near ~35-36 degC at saturation) is shown.

4. *Observational demonstration*: band occupancy at a humid vs an arid
   SURFRAD station in July 2023 (informational).

Run:  .venv-val/bin/python validation/relative_strain_index/validate_rsi.py
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

# Pre-registered acceptance criteria (see README.md)
TOL_TRANSCRIPTION = 1e-12
TOL_MAGNUS = 0.02  # |dRSI| from the Magnus-variant difference, in-domain
T_DOMAIN = (21.0, 35.0)  # validity envelope used for the gates [degC]

#: Asghari et al. (2020) assessment bands
BANDS = [
    (-np.inf, 0.15, "comfort"),
    (0.15, 0.25, "discomfort for sensitive"),
    (0.25, 0.35, "discomfort for all"),
    (0.35, 0.45, "heat-stroke risk >50%"),
    (0.45, np.inf, "hyperthermia risk for all"),
]


def rsi_thermofeel(t_c, rh):
    return tf.calculate_relative_strain_index(
        tf.celsius_to_kelvin(np.asarray(t_c, dtype=float)), rh
    )


def rsi_reimpl_bom(t_c, rh):
    """Independent transcription with thermofeel's BoM vapour pressure."""
    t_c = np.asarray(t_c, dtype=float)
    e = rh / 100.0 * 6.105 * np.exp(17.27 * t_c / (237.7 + t_c))
    return (t_c - 21.0) / (58.0 - e)


def rsi_asghari(t_c, rh):
    """The closed form with the vapour pressure exactly as printed in
    Asghari et al. (2020): e = (RH/100) 6.112 x 10^(7.5 T/(237.7+T))."""
    t_c = np.asarray(t_c, dtype=float)
    e = rh / 100.0 * 6.112 * 10.0 ** (7.5 * t_c / (237.7 + t_c))
    return (t_c - 21.0) / (58.0 - e)


def main() -> int:
    plt = style.use_style()
    failures = []

    t_c, rh = np.meshgrid(
        np.arange(T_DOMAIN[0], T_DOMAIN[1] + 0.01, 0.25),
        np.arange(10.0, 90.1, 1.0),
    )
    t_c, rh = t_c.ravel(), rh.ravel()

    # ------------------------------------------------- 1. transcription check
    d1 = np.nanmax(np.abs(rsi_thermofeel(t_c, rh) - rsi_reimpl_bom(t_c, rh)))
    print(f"transcription check (N={t_c.size}): max |dRSI| = {d1:.2e}")
    if d1 > TOL_TRANSCRIPTION:
        failures.append("transcription")

    # ------------------------------------------- 2. Magnus-variant sensitivity
    ours = rsi_thermofeel(t_c, rh)
    ref = rsi_asghari(t_c, rh)
    d2 = np.nanmax(np.abs(ours - ref))
    print(
        f"Magnus-variant sensitivity in-domain {T_DOMAIN} degC: "
        f"max |dRSI| = {d2:.4f} (bands are 0.10 wide)"
    )
    if d2 > TOL_MAGNUS:
        failures.append("magnus_sensitivity")
    stats = [
        (
            {"comparison": "RSI (BoM Magnus) vs Asghari-printed Magnus, in-domain"},
            metrics.error_stats(ours, ref),
        )
    ]
    metrics.write_results(
        metrics.stats_frame(stats), HERE / "results" / "rsi_formula_stats.csv"
    )

    # ------------------------------------------------- 3. behavioural checks
    # boundary identity: at Ta = 21 degC the numerator vanishes, so RSI = 0
    # for every humidity (and RH-monotonicity only holds strictly above it)
    rh_line = np.arange(0.0, 100.1, 1.0)
    d_bound = np.max(np.abs(rsi_thermofeel(np.full_like(rh_line, 21.0), rh_line)))
    print(f"boundary identity RSI(21 degC, RH) = 0: max |d| = {d_bound:.2e}")
    if d_bound > TOL_TRANSCRIPTION:
        failures.append("boundary_identity")

    grid_t = np.arange(T_DOMAIN[0] + 0.25, T_DOMAIN[1] + 0.01, 0.25)  # T > 21
    grid_rh = np.arange(10.0, 90.1, 1.0)
    tt, rr = np.meshgrid(grid_t, grid_rh)
    z = rsi_thermofeel(tt, rr)
    mono_t = np.all(np.diff(z, axis=1) > 0)
    mono_rh = np.all(np.diff(z, axis=0) > 0)
    print(f"strictly increasing in T: {mono_t}, in RH: {mono_rh} (21 < T <= 35)")
    if not (mono_t and mono_rh):
        failures.append("monotonicity")

    # band-crossing temperatures per RH (informational reference table)
    rows = []
    fine_t = np.arange(15.0, 45.001, 0.01)
    for rh_v in (20.0, 40.0, 60.0, 80.0, 100.0):
        vals = rsi_thermofeel(fine_t, np.full_like(fine_t, rh_v))
        rec = {"rh_pct": rh_v}
        for thr, name in (
            (0.15, "sensitive"),
            (0.25, "all"),
            (0.35, "heatstroke"),
            (0.45, "hyperthermia"),
        ):
            above = fine_t[(vals >= thr) & (fine_t > 21.0)]
            rec[f"T_at_RSI_{thr}"] = (
                round(float(above.min()), 1) if len(above) else np.nan
            )
        rows.append(rec)
    crossings = pd.DataFrame(rows)
    metrics.write_results(crossings, HERE / "results" / "rsi_band_crossings.csv")

    # ------------------------------------------------------------------ plots
    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.4))
    t_plot = np.arange(21.0, 39.001, 0.05)
    for rh_v, c in (
        (20, "tab:blue"),
        (50, "tab:green"),
        (80, "tab:orange"),
        (100, "tab:red"),
    ):
        axs[0].plot(
            t_plot,
            rsi_thermofeel(t_plot, np.full_like(t_plot, float(rh_v))),
            color=c,
            label=f"RH {rh_v}%",
        )
    for thr, _, name in [(b[0], b[1], b[2]) for b in BANDS[1:]]:
        axs[0].axhline(thr, color="0.6", lw=0.5)
        axs[0].text(21.2, thr + 0.005, name, fontsize=6, color="0.35")
    axs[0].set_ylim(-0.1, 1.2)
    axs[0].axvline(35, color="k", lw=0.6, ls=":")
    axs[0].text(35.2, 1.0, "validity edge\n(e -> 58 hPa at\nsaturation)", fontsize=6)
    axs[0].set_xlabel("air temperature [degC]")
    axs[0].set_ylabel("RSI [-]")
    axs[0].set_title("Response and assessment bands (Asghari et al. 2020)")
    axs[0].legend(loc="upper left")

    obs = pd.concat(
        [
            surfrad.load(s, surfrad.month_days(2023, 7))[["temp", "rh", "station"]]
            for s in ("gwn", "dra")
        ],
        ignore_index=True,
    ).dropna()
    # restrict to the index's warm-season domain (T > 21 degC, T <= 35 degC)
    obs = obs[(obs.temp > 21.0) & (obs.temp <= 35.0)]
    ticks = np.arange(len(BANDS))
    width = 0.38
    for k, (sta, color) in enumerate((("gwn", "tab:red"), ("dra", "tab:orange"))):
        sel = obs[obs.station == sta]
        v = rsi_thermofeel(sel.temp.to_numpy(), sel.rh.to_numpy())
        occ = [100.0 * np.mean((v >= lo) & (v < hi)) for lo, hi, _ in BANDS]
        axs[1].bar(
            ticks + (k - 0.5) * width,
            occ,
            width,
            color=color,
            label=f"{sta} ({surfrad.STATIONS[sta][0]})",
        )
    axs[1].set_xticks(ticks)
    axs[1].set_xticklabels([b[2].replace(" ", "\n", 1) for b in BANDS], fontsize=6.5)
    axs[1].set_ylabel("share of warm July minutes [%]")
    axs[1].set_title("Band occupancy on SURFRAD obs, T in (21, 35] degC")
    axs[1].legend()
    fig.suptitle("Relative Strain Index")
    style.save_png(fig, HERE / "plots" / "rsi_checks.png")

    print("\nRESULT:", "FAIL " + ",".join(failures) if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

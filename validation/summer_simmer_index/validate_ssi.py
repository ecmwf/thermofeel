# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Validation of ``thermofeel.calculate_summer_simmer_index`` (Pepi 1987).

The SSI has a known provenance problem: the original Weatherwise article is
not openly available, and public secondary sources are not fully consistent
with each other. thermofeel implements the *common 1987 closed form* and its
docstring/CHANGELOG say so explicitly. The validation therefore:

1. *Independent transcription check*: compares against the closed form
   SSI_F = 1.98 (T_F - (0.55 - 0.0055 RH)(T_F - 58)) - 56.83 transcribed
   independently from public sources (e.g. the vCalc equation page), over a
   dense grid -- guards against implementation/unit-conversion slips.

2. *Structural identity*: verifies SSI_F = 1.98 THI_F - 56.83 exactly, where
   THI_F is Thom's index evaluated in Fahrenheit with the 58 degF pivot, and
   quantifies the relationship to ``calculate_discomfort_index`` (14.5 degC
   pivot vs 58 degF = 14.44 degC: a small, bounded family offset).

3. *Provenance uncertainty quantification* (informational): evaluates the
   closed form at the anchor points quoted by public secondary sources
   (Pepi's narrative "90 degF / 50% RH feels like about 100 degF" vs the
   later "New SSI" table value 106 degF) and reports the discrepancies, so
   users see the documented medium confidence, not a false precision claim.

4. *Observational demonstration*: SSI vs DI on humid-subtropical July
   observations (SURFRAD Goodwin Creek) -- the affine family relation should
   be visible on real data (informational).

Run:  .venv-val/bin/python validation/summer_simmer_index/validate_ssi.py
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
TOL_TRANSCRIPTION_F = 1e-9  # closed form, independent transcription [degF]
TOL_IDENTITY_F = 1e-9  # SSI_F = 1.98*THI_F - 56.83 [degF]


def ssi_fahrenheit(t_c, rh):
    """thermofeel SSI evaluated in Fahrenheit."""
    return tf.kelvin_to_fahrenheit(
        tf.calculate_summer_simmer_index(
            tf.celsius_to_kelvin(np.asarray(t_c, dtype=float)), rh
        )
    )


def independent_ssi_f(t_f, rh):
    """The common 1987 closed form, transcribed independently (vCalc page)."""
    t_f = np.asarray(t_f, dtype=float)
    return 1.98 * (t_f - (0.55 - 0.0055 * rh) * (t_f - 58.0)) - 56.83


def main() -> int:
    plt = style.use_style()
    failures = []

    t_c, rh = np.meshgrid(np.arange(-10.0, 50.1, 0.5), np.arange(0.0, 100.1, 1.0))
    t_c, rh = t_c.ravel(), rh.ravel()
    t_f = t_c * 9.0 / 5.0 + 32.0

    # ------------------------------------------- 1. independent transcription
    ours_f = ssi_fahrenheit(t_c, rh)
    ref_f = independent_ssi_f(t_f, rh)
    d1 = np.nanmax(np.abs(ours_f - ref_f))
    print(f"transcription check (N={ours_f.size}): max |dSSI| = {d1:.2e} degF")
    if d1 > TOL_TRANSCRIPTION_F:
        failures.append("transcription")

    # ------------------------------------------------- 2. structural identity
    thi_f = t_f - (0.55 - 0.0055 * rh) * (t_f - 58.0)
    d2 = np.nanmax(np.abs(ours_f - (1.98 * thi_f - 56.83)))
    print(f"affine identity SSI_F = 1.98 THI_F - 56.83: max |d| = {d2:.2e} degF")
    if d2 > TOL_IDENTITY_F:
        failures.append("affine_identity")

    # family offset vs calculate_discomfort_index (14.5 degC vs 58 degF pivot)
    di_f = tf.kelvin_to_fahrenheit(
        tf.calculate_discomfort_index(tf.celsius_to_kelvin(t_c), rh)
    )
    fam = np.abs(di_f - thi_f)
    print(
        f"family offset |DI_F - THI_F(58F pivot)|: max = {fam.max():.3f} degF "
        f"(pivot 14.5 degC = 58.1 degF vs 58 degF)"
    )

    # ------------------------------------- 3. provenance anchors (informational)
    anchors = pd.DataFrame(
        [
            # T_F, RH, quoted SSI, source of the quote
            (90.0, 50.0, 100.0, "Pepi narrative (summersimmer.com)"),
            (90.0, 50.0, 106.0, "'New SSI' table (summersimmer.com)"),
            (100.0, 50.0, 120.0, "'New SSI' table (summersimmer.com)"),
        ],
        columns=["t_f", "rh", "ssi_quoted", "quote_source"],
    )
    anchors["ssi_closed_form"] = independent_ssi_f(
        anchors.t_f.to_numpy(), anchors.rh.to_numpy()
    ).round(1)
    anchors["difference"] = (anchors.ssi_closed_form - anchors.ssi_quoted).round(1)
    metrics.write_results(anchors, HERE / "results" / "ssi_provenance_anchors.csv")

    # ------------------------------------------------- 4. observational demo
    obs = surfrad.load("gwn", surfrad.month_days(2023, 7))[["temp", "rh"]].dropna()
    ssi_o = tf.kelvin_to_celsius(
        tf.calculate_summer_simmer_index(
            tf.celsius_to_kelvin(obs.temp.to_numpy()), obs.rh.to_numpy()
        )
    )
    di_o = tf.kelvin_to_celsius(
        tf.calculate_discomfort_index(
            tf.celsius_to_kelvin(obs.temp.to_numpy()), obs.rh.to_numpy()
        )
    )

    fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.3))
    hb = axs[0].hexbin(di_o, ssi_o, gridsize=50, bins="log", cmap="viridis")
    axs[0].set_xlabel("Discomfort Index [degC]")
    axs[0].set_ylabel("Summer Simmer Index [degC]")
    axs[0].set_title("Affine family relation on real obs (Goodwin Creek, Jul 2023)")
    fig.colorbar(hb, ax=axs[0], label="log$_{10}$ count")

    tt = np.arange(20.0, 46.0, 0.5)
    for rr, c in ((20, "tab:blue"), (50, "tab:green"), (80, "tab:red")):
        axs[1].plot(
            tt,
            tf.kelvin_to_celsius(
                tf.calculate_summer_simmer_index(
                    tf.celsius_to_kelvin(tt), np.full_like(tt, float(rr))
                )
            ),
            color=c,
            label=f"SSI, RH {rr}%",
        )
    axs[1].plot(tt, tt, "k--", lw=0.8, label="air temperature")
    axs[1].set_xlabel("air temperature [degC]")
    axs[1].set_ylabel("SSI [degC]")
    axs[1].set_title("Response curves")
    axs[1].legend()
    fig.suptitle("Summer Simmer Index (common 1987 closed form)")
    style.save_png(fig, HERE / "plots" / "ssi_checks.png")

    print("\nRESULT:", "FAIL " + ",".join(failures) if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

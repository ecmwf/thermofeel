# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Validation of ``thermofeel.calculate_pet``.

Four parts:

A. *Definitional identity*: PET is defined as the air temperature of the
   reference environment reproducing the body state, so feeding that
   environment to the function must return the air temperature exactly. This
   pins the reference environment, the clothing convention and the whole-body
   activity-unit convention simultaneously.

B. *Published reference values*: Hoeppe (1999) Table 1, together with the
   Walther & Goestchel (2018) recomputation of the same cases. The two papers
   disagree with each other by up to 1.3 K, which bounds what any
   implementation can claim.

C. *Independent implementation*: element-wise comparison against
   ``pythermalcomfort.pet_steady`` (v4) over a random sweep of the input
   space. This is a *fallible* oracle: it wraps two SciPy ``fsolve`` calls per
   point in a Python loop and emits non-convergence warnings on a
   non-negligible fraction of hot, humid inputs. It also makes several
   documented model choices that thermofeel deliberately does not follow, so
   exact agreement is neither expected nor desirable.

D. *Throughput*: the reason for the vectorised reformulation. Measured on this
   machine, against the same reference implementation.

Run:  .venv-val/bin/python validation/pet/validate_pet.py
"""

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[0]))

from vlib import metrics, style  # noqa: E402

import thermofeel as tf  # noqa: E402

# Pre-registered acceptance criteria (see README.md)
TOL_IDENTITY_K = 1e-4  # A: definitional identity
TOL_HOEPPE_RMS_K = 1.2  # B: rms vs Hoeppe's published table
TOL_PTC_RMS_K = 0.40  # C: rms vs pythermalcomfort
TOL_PTC_MAX_K = 2.0  # C: worst single point

#: Hoeppe (1999) Table 1, with the Walther & Goestchel (2018) recomputation.
#: (Ta [degC], Tmrt [degC], v [m/s], vpa [hPa], PET_hoeppe, PET_walther)
REFERENCE_CASES = [
    (21.0, 21.0, 0.1, 12.0, 21.0, 21.00),
    (-5.0, 40.0, 0.5, 2.0, 10.0, 10.09),
    (-5.0, -5.0, 5.0, 2.0, -13.0, -11.69),
    (30.0, 60.0, 1.0, 21.0, 43.0, 42.70),
    (30.0, 30.0, 1.0, 21.0, 29.0, 28.95),
]


def _pet_c(t_c, mrt_c, va, vpa=None, rh=None, **kw):
    """thermofeel PET in Celsius, from Celsius inputs."""
    out = tf.calculate_pet(
        tf.celsius_to_kelvin(np.asarray(t_c, dtype=float)),
        tf.celsius_to_kelvin(np.asarray(mrt_c, dtype=float)),
        np.asarray(va, dtype=float),
        rh=rh,
        vapour_pressure_hpa=vpa,
        **kw,
    )
    return tf.kelvin_to_celsius(out)


def part_a_identity(plt):
    t_c = np.arange(-20.0, 45.1, 0.5)
    pet = _pet_c(t_c, t_c, np.full_like(t_c, 0.1), vpa=np.full_like(t_c, 12.0))
    err = pet - t_c
    worst = float(np.max(np.abs(err)))
    print(f"part A: definitional identity, max |PET - Ta| = {worst:.2e} K")

    fig, ax = plt.subplots(figsize=(6.2, 2.6))
    ax.plot(t_c, err * 1e6, color="tab:blue")
    ax.set_xlabel("air temperature of the reference environment [degC]")
    ax.set_ylabel("PET − Ta [µK]")
    ax.set_title(
        f"Reference-environment identity: max |error| = {worst:.1e} K "
        f"(bisection tolerance)"
    )
    style.save_png(fig, HERE / "plots" / "pet_identity.png")
    return worst <= TOL_IDENTITY_K


def part_b_published(plt):
    rows = []
    for t_c, mrt_c, va, vpa, hoeppe, walther in REFERENCE_CASES:
        got = float(_pet_c([t_c], [mrt_c], [va], vpa=[vpa])[0])
        rows.append(
            {
                "Ta": t_c,
                "Tmrt": mrt_c,
                "v": va,
                "vpa": vpa,
                "pet_hoeppe_1999": hoeppe,
                "pet_walther_2018": walther,
                "pet_thermofeel": round(got, 2),
                "d_vs_hoeppe": round(got - hoeppe, 2),
                "d_vs_walther": round(got - walther, 2),
                "paper_spread": round(abs(hoeppe - walther), 2),
            }
        )
    df = pd.DataFrame(rows)
    metrics.write_results(df, HERE / "results" / "pet_published_cases.csv")
    rms = float(np.sqrt(np.mean(df.d_vs_hoeppe**2)))
    print(
        f"part B: rms vs Hoeppe (1999) = {rms:.2f} K; "
        f"the two papers differ by up to {df.paper_spread.max():.2f} K"
    )

    fig, ax = plt.subplots(figsize=(6.6, 3.0))
    x = np.arange(len(df))
    ax.plot(x, df.pet_hoeppe_1999, "s", ms=7, color="k", label="Höppe (1999)")
    ax.plot(
        x,
        df.pet_walther_2018,
        "d",
        ms=6,
        color="tab:green",
        label="Walther & Goestchel (2018)",
    )
    ax.plot(x, df.pet_thermofeel, "o", ms=5, color="tab:red", label="thermofeel")
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{r.Ta:.0f}/{r.Tmrt:.0f}\n{r.v}m/s" for r in df.itertuples()], fontsize=7
    )
    ax.set_ylabel("PET [degC]")
    ax.set_title("Published reference cases (Höppe 1999, Table 1)")
    ax.legend()
    style.save_png(fig, HERE / "plots" / "pet_published_cases.png")
    return rms <= TOL_HOEPPE_RMS_K


def part_c_pythermalcomfort(plt):
    from pythermalcomfort.models import pet_steady

    rng = np.random.default_rng(20260731)
    n = 500
    ta = rng.uniform(-20.0, 45.0, n)
    mrt = ta + rng.uniform(-10.0, 40.0, n)
    va = rng.uniform(0.1, 8.0, n)
    rh = rng.uniform(10.0, 95.0, n)
    clo = rng.uniform(0.3, 1.6, n)
    m_act = rng.uniform(60.0, 200.0, n)

    ours = _pet_c(ta, mrt, va, rh=rh, clo=clo, m_act_w=m_act)
    with warnings.catch_warnings():
        # pythermalcomfort's fsolve emits "not making good progress" on a
        # non-negligible fraction of hot/humid inputs; it is a fallible oracle.
        warnings.simplefilter("ignore")
        ref = np.array(
            [
                float(pet_steady(tdb=a, tr=m, v=v, rh=r, met=mm / 58.2, clo=c).pet)
                for a, m, v, r, c, mm in zip(ta, mrt, va, rh, clo, m_act)
            ]
        )

    stats = metrics.error_stats(ours, ref)
    d = ours - ref
    print(
        f"part C: vs pythermalcomfort (N={n}): rms = {stats['RMSE']:.3f} K, "
        f"max |d| = {np.max(np.abs(d)):.3f} K, median |d| = "
        f"{np.median(np.abs(d)):.3f} K"
    )

    fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.2))
    hb = axs[0].hexbin(ref, ours, gridsize=45, bins="log", cmap="viridis")
    lim = [min(ref.min(), ours.min()), max(ref.max(), ours.max())]
    axs[0].plot(lim, lim, "r--", lw=0.9)
    axs[0].set_xlabel("PET, pythermalcomfort [degC]")
    axs[0].set_ylabel("PET, thermofeel [degC]")
    axs[0].set_title(f"N = {n}, RMSE = {stats['RMSE']:.3f} K")
    fig.colorbar(hb, ax=axs[0], label="log$_{10}$ count")
    axs[1].hist(d, bins=60, color="tab:blue")
    axs[1].set_yscale("log")
    axs[1].set_xlabel("thermofeel − pythermalcomfort [K]")
    axs[1].set_title("residual (documented model-variant differences)")
    fig.suptitle("Independent-implementation comparison")
    style.save_png(fig, HERE / "plots" / "pet_vs_pythermalcomfort.png")

    metrics.write_results(
        metrics.stats_frame(
            [({"comparison": f"PET vs pythermalcomfort (N={n})"}, stats)]
        ),
        HERE / "results" / "pet_agreement_stats.csv",
    )
    return stats["RMSE"] <= TOL_PTC_RMS_K and np.max(np.abs(d)) <= TOL_PTC_MAX_K


def part_d_throughput():
    from pythermalcomfort.models import pet_steady

    rng = np.random.default_rng(0)
    n = 200_000
    ta = rng.uniform(-10.0, 40.0, n)
    mrt = ta + rng.uniform(-5.0, 30.0, n)
    va = rng.uniform(0.2, 8.0, n)
    rh = rng.uniform(20.0, 90.0, n)

    _pet_c(ta[:1000], mrt[:1000], va[:1000], rh=rh[:1000])  # warm up
    t0 = time.perf_counter()
    _pet_c(ta, mrt, va, rh=rh)
    tf_rate = n / (time.perf_counter() - t0)

    m = 300  # pythermalcomfort is far too slow to run the full array
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t0 = time.perf_counter()
        for a, mm, v, r in zip(ta[:m], mrt[:m], va[:m], rh[:m]):
            pet_steady(tdb=a, tr=mm, v=v, rh=r, met=80 / 58.2, clo=0.9)
        ptc_rate = m / (time.perf_counter() - t0)

    print(
        f"part D: thermofeel {tf_rate:,.0f} pt/s vs pythermalcomfort "
        f"{ptc_rate:,.0f} pt/s -> {tf_rate / ptc_rate:.0f}x"
    )
    o1280 = 6_599_680
    print(
        f"        one O1280 field ({o1280:,} pts): "
        f"{o1280 / tf_rate:,.0f} s vs {o1280 / ptc_rate / 3600:,.1f} core-hours"
    )
    df = pd.DataFrame(
        [
            {
                "implementation": "thermofeel.calculate_pet",
                "points_per_s": round(tf_rate),
                "o1280_field_s": round(o1280 / tf_rate, 1),
            },
            {
                "implementation": "pythermalcomfort.pet_steady",
                "points_per_s": round(ptc_rate),
                "o1280_field_s": round(o1280 / ptc_rate, 1),
            },
        ]
    )
    metrics.write_results(df, HERE / "results" / "pet_throughput.csv")


def main() -> int:
    plt = style.use_style()
    failures = []
    if not part_a_identity(plt):
        failures.append("A_identity")
    if not part_b_published(plt):
        failures.append("B_published")
    if not part_c_pythermalcomfort(plt):
        failures.append("C_pythermalcomfort")
    part_d_throughput()  # informational
    print("\nRESULT:", "FAIL " + ",".join(failures) if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

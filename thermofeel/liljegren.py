# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""
Liljegren Wet Bulb Globe Temperature - internal implementation.

Physically based WBGT after Liljegren et al. (2008), the "gold standard" for
deriving WBGT from standard meteorological variables, used operationally by KNMI
for its heat-force indicator. The globe and natural-wet-bulb sensors are each
solved from their steady-state energy balance by fixed-point iteration.

The constants and property functions are transcribed from Liljegren's reference
C code (MIT-licensed mirror github.com/mdljts/wbgt) and Kong & Huber (2022), and
validated bit-for-bit against that reference. This module is internal: the
public entry points are ``thermofeel.calculate_wbgt_liljegren``,
``thermofeel.calculate_wind_speed_2m_liljegren`` and
``thermofeel.calculate_heat_force``.

Reference: Liljegren et al. (2008) https://doi.org/10.1080/15459620802310770
See also: Kong and Huber (2022) https://doi.org/10.1029/2021EF002334
"""

import numpy as np
from numpy.typing import ArrayLike

from .helpers import celsius_to_kelvin, kelvin_to_celsius

# Physical constants (Liljegren et al. 2008; mdljts/wbgt header)
STEFANB = 5.6696e-8  # Stefan-Boltzmann constant [W m-2 K-4]
CP = 1003.5  # specific heat of dry air [J kg-1 K-1]
M_AIR = 28.97  # molar mass of dry air [g mol-1]
M_H2O = 18.015  # molar mass of water [g mol-1]
R_GAS = 8314.34  # universal gas constant [J kmol-1 K-1]
R_AIR = R_GAS / M_AIR  # gas constant for air [J kg-1 K-1]
PR = CP / (CP + 1.25 * R_AIR)  # Prandtl number
RATIO = CP * M_AIR / M_H2O  # psychrometric grouping
EMIS_WICK = 0.95
ALB_WICK = 0.4
D_WICK = 0.007  # wick diameter [m]
L_WICK = 0.0254  # wick length [m]
EMIS_GLOBE = 0.95
ALB_GLOBE = 0.05
D_GLOBE = 0.0508  # globe diameter [m]
EMIS_SFC = 0.999
ALB_SFC = 0.45
CZA_MIN = 0.00873  # cos(89.5 deg): below this the sun is treated as down
MIN_SPEED = 0.13  # internal floor on wind speed in the Reynolds number [m/s]
CONVERGENCE = 0.02  # iteration tolerance [K]
MAX_ITER = 500  # iteration cap
MIN_WIND_10M = 0.62  # KNMI minimum 10 m wind (~0.5 m/s at 2 m) [m/s]

# Pasquill-Gifford stability lookup for the 10 m -> 2 m wind profile
# (Liljegren et al. 2008; Kong & Huber 2022 PyWBGT). Rows are wind-speed bins,
# columns are solar-radiation / night bins; the value is the stability class.
LSRDT = np.array(
    [
        [1, 1, 2, 4, 0, 5, 6, 0],
        [1, 2, 3, 4, 0, 5, 6, 0],
        [2, 2, 3, 4, 0, 4, 4, 0],
        [3, 3, 4, 4, 0, 0, 0, 0],
        [3, 4, 4, 4, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
    ]
)
# Wind-profile power-law exponent per stability class (urban terrain).
URBAN_EXP = np.array([0.15, 0.15, 0.20, 0.25, 0.30, 0.30])


def esat(tk: ArrayLike) -> np.ndarray:
    """Saturation vapour pressure over liquid water [hPa], Buck (1981)."""
    y = (tk - 273.15) / (tk - 32.18)
    return 1.004 * 6.1121 * np.exp(17.502 * y)


def dew_point(e: ArrayLike) -> np.ndarray:
    """Dew-point temperature [K] from vapour pressure [hPa] (inverse of esat)."""
    z = np.log(e / (6.1121 * 1.004))
    return 273.15 + 240.97 * z / (17.502 - z)


def viscosity(tk: ArrayLike) -> np.ndarray:
    """Dynamic viscosity of air [kg m-1 s-1] (Bird, Stewart & Lightfoot)."""
    sigma = 3.617
    eps_kappa = 97.0
    tr = tk / eps_kappa
    omega = (tr - 2.9) / 0.4 * (-0.034) + 1.048
    return 2.6693e-6 * np.sqrt(M_AIR * tk) / (sigma * sigma * omega)


def thermal_cond(tk: ArrayLike) -> np.ndarray:
    """Thermal conductivity of air [W m-1 K-1] (Eucken relation)."""
    return (CP + 1.25 * R_AIR) * viscosity(tk)


def diffusivity(tk: ArrayLike, pair: ArrayLike) -> np.ndarray:
    """Diffusivity of water vapour in air [m2 s-1]; pair in hPa (BSL p.505)."""
    pcrit_air = 36.4
    pcrit_h2o = 218.0
    tcrit_air = 132.0
    tcrit_h2o = 647.3
    a = 3.640e-4
    b = 2.334
    pcrit13 = (pcrit_air * pcrit_h2o) ** (1.0 / 3.0)
    tcrit512 = (tcrit_air * tcrit_h2o) ** (5.0 / 12.0)
    tcrit12 = np.sqrt(tcrit_air * tcrit_h2o)
    mmix = np.sqrt(1.0 / M_AIR + 1.0 / M_H2O)
    patm = pair / 1013.25
    return a * (tk / tcrit12) ** b * pcrit13 * tcrit512 * mmix / patm * 1e-4


def evap(tk: ArrayLike) -> np.ndarray:
    """Latent heat of vaporisation [J kg-1], valid 283-313 K."""
    return (313.15 - tk) / 30.0 * (-71100.0) + 2.4073e6


def emis_atm(tk: ArrayLike, rh: ArrayLike) -> np.ndarray:
    """Clear-sky atmospheric emissivity; rh as fraction (Oke 2nd ed.)."""
    e = rh * esat(tk)
    return 0.575 * e**0.143


def h_sphere_in_air(tk: ArrayLike, pair: ArrayLike, speed: ArrayLike) -> np.ndarray:
    """Convective heat-transfer coefficient for the globe (sphere) [W m-2 K-1]."""
    density = pair * 100.0 / (R_AIR * tk)
    re = np.maximum(speed, MIN_SPEED) * density * D_GLOBE / viscosity(tk)
    nu = 2.0 + 0.6 * np.sqrt(re) * PR**0.3333
    return nu * thermal_cond(tk) / D_GLOBE


def h_cylinder_in_air(tk: ArrayLike, pair: ArrayLike, speed: ArrayLike) -> np.ndarray:
    """Convective heat-transfer coefficient for the wick (cylinder) [W m-2 K-1]."""
    a = 0.56
    b = 0.281
    c = 0.4
    density = pair * 100.0 / (R_AIR * tk)
    re = np.maximum(speed, MIN_SPEED) * density * D_WICK / viscosity(tk)
    nu = b * re ** (1.0 - c) * PR ** (1.0 - a)
    return nu * thermal_cond(tk) / D_WICK


def solve_globe(
    ta: ArrayLike,
    rh: ArrayLike,
    pair: ArrayLike,
    speed: ArrayLike,
    solar: ArrayLike,
    fdir: ArrayLike,
    cza: ArrayLike,
) -> np.ndarray:
    """Globe temperature [degC] by fixed-point iteration of the energy balance.

    ``rh`` is a fraction. NaN is returned where the iteration does not converge
    within the cap.
    """
    tsfc = ta
    emis = emis_atm(ta, rh)
    # Direct-beam geometry term; guarded so fdir == 0 contributes 0 even when the
    # sun is at/below the horizon (cza -> 0), avoiding 0 * inf.
    cza_safe = np.where(cza > CZA_MIN, cza, 1.0)
    beam = np.where(fdir > 0.0, fdir * (1.0 / (2.0 * cza_safe) - 1.0), 0.0)

    tg_prev = np.array(ta, dtype=float, copy=True)
    result = np.full(np.shape(ta), np.nan, dtype=float)
    converged = np.zeros(np.shape(ta), dtype=bool)
    for _ in range(MAX_ITER):
        tref = 0.5 * (tg_prev + ta)
        h = h_sphere_in_air(tref, pair, speed)
        tg_new = (
            0.5 * (emis * ta**4 + EMIS_SFC * tsfc**4)
            - h / (STEFANB * EMIS_GLOBE) * (tg_prev - ta)
            + solar
            / (2.0 * STEFANB * EMIS_GLOBE)
            * (1.0 - ALB_GLOBE)
            * (beam + 1.0 + ALB_SFC)
        ) ** 0.25
        now = (~converged) & (np.abs(tg_new - tg_prev) < CONVERGENCE)
        result = np.where(now, tg_new - 273.15, result)
        converged = converged | now
        tg_prev = np.where(converged, tg_prev, 0.9 * tg_prev + 0.1 * tg_new)
        if converged.all():
            break
    return result


def solve_wetbulb(
    ta: ArrayLike,
    rh: ArrayLike,
    pair: ArrayLike,
    speed: ArrayLike,
    solar: ArrayLike,
    fdir: ArrayLike,
    cza: ArrayLike,
    rad: float,
) -> np.ndarray:
    """Wet-bulb temperature [degC] by fixed-point iteration of the energy balance.

    ``rh`` is a fraction. With ``rad=1`` this is the natural wet-bulb temperature
    (the term entering WBGT); ``rad=0`` gives the psychrometric wet bulb. NaN
    where not converged.
    """
    tsfc = ta
    # Solar-zenith angle, guarded so tan(sza) stays finite when the sun is at or
    # below the horizon (fdir is already 0 there, so the term contributes 0).
    cza_safe = np.where(cza > CZA_MIN, cza, 1.0)
    sza = np.arccos(np.clip(cza_safe, -1.0, 1.0))
    emis = emis_atm(ta, rh)
    eair = rh * esat(ta)

    tw_prev = dew_point(eair)
    result = np.full(np.shape(ta), np.nan, dtype=float)
    converged = np.zeros(np.shape(ta), dtype=bool)
    for _ in range(MAX_ITER):
        tref = 0.5 * (tw_prev + ta)
        h = h_cylinder_in_air(tref, pair, speed)
        fatm = STEFANB * EMIS_WICK * (
            0.5 * (emis * ta**4 + EMIS_SFC * tsfc**4) - tw_prev**4
        ) + (1.0 - ALB_WICK) * solar * (
            (1.0 - fdir) * (1.0 + 0.25 * D_WICK / L_WICK)
            + fdir * ((np.tan(sza) / np.pi) + 0.25 * D_WICK / L_WICK)
            + ALB_SFC
        )
        ewick = esat(tw_prev)
        density = pair * 100.0 / (R_AIR * tref)
        sc = viscosity(tref) / (density * diffusivity(tref, pair))
        tw_new = (
            ta
            - evap(tref) / RATIO * (ewick - eair) / (pair - ewick) * (PR / sc) ** 0.56
            + (fatm / h) * rad
        )
        now = (~converged) & (np.abs(tw_new - tw_prev) < CONVERGENCE)
        result = np.where(now, tw_new - 273.15, result)
        converged = converged | now
        tw_prev = np.where(converged, tw_prev, 0.9 * tw_prev + 0.1 * tw_new)
        if converged.all():
            break
    return result


def wind_speed_2m(va: ArrayLike, cossza: ArrayLike, ssrd: ArrayLike) -> np.ndarray:
    """10 m -> 2 m wind speed via the Liljegren stability-dependent profile.

    ``va * (2/10)**p``, where the power-law exponent ``p`` comes from a
    Pasquill-Gifford stability class derived from solar elevation, incoming
    radiation and wind speed. The result is floored at MIN_SPEED (0.13 m/s).
    """
    va = np.asarray(va, dtype=float)
    cossza = np.asarray(cossza, dtype=float)
    ssrd = np.asarray(ssrd, dtype=float)

    daytime = cossza > 0.0

    # radiation / night column index of the stability table
    col = np.where(
        ssrd >= 925.0,
        0,
        np.where(ssrd >= 675.0, 1, np.where(ssrd >= 175.0, 2, 3)),
    )
    col = np.where(daytime, col, 5)

    # wind-speed row index of the stability table (day and night thresholds)
    row_day = np.where(
        va >= 6.0,
        4,
        np.where(va >= 5.0, 3, np.where(va >= 3.0, 2, np.where(va >= 2.0, 1, 0))),
    )
    row_night = np.where(va >= 2.5, 2, np.where(va >= 2.0, 1, 0))
    row = np.where(daytime, row_day, row_night)

    stability_class = LSRDT[row, col]
    exponent = URBAN_EXP[stability_class - 1]
    return np.maximum(va * (2.0 / 10.0) ** exponent, MIN_SPEED)


def wbgt(
    t2_k: ArrayLike,
    rh: ArrayLike,
    pressure: ArrayLike,
    speed_2m: ArrayLike,
    ssrd: ArrayLike,
    fdir: ArrayLike,
    cossza: ArrayLike,
) -> np.ndarray:
    """Liljegren WBGT [K] from the 2 m wind and the other meteorological inputs.

    ``rh`` is a percentage and ``speed_2m`` is the wind already scaled to 2 m.
    Applies the KNMI direct-beam guards (clamp fdir to [0, 0.9]; zero it below
    89.5 deg zenith), solves the globe and natural-wet-bulb energy balances, and
    combines them as ``WBGT = 0.7*Tnw + 0.2*Tg + 0.1*Ta``. NaN where the
    iteration does not converge.
    """
    t2_k = np.asarray(t2_k, dtype=float)
    rh = np.asarray(rh, dtype=float)
    pressure = np.asarray(pressure, dtype=float)
    ssrd = np.asarray(ssrd, dtype=float)
    fdir = np.asarray(fdir, dtype=float)
    cossza = np.asarray(cossza, dtype=float)

    rh_frac = rh / 100.0

    # KNMI direct-beam guards: clamp to [0, 0.9] and zero below the horizon.
    fdir = np.clip(fdir, 0.0, 0.9)
    fdir = np.where(cossza < CZA_MIN, 0.0, fdir)

    tg_c = solve_globe(t2_k, rh_frac, pressure, speed_2m, ssrd, fdir, cossza)
    tnwb_c = solve_wetbulb(t2_k, rh_frac, pressure, speed_2m, ssrd, fdir, cossza, 1.0)
    t2_c = kelvin_to_celsius(t2_k)

    wbgt_c = 0.1 * t2_c + 0.2 * tg_c + 0.7 * tnwb_c
    return celsius_to_kelvin(wbgt_c)

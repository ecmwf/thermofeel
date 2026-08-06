# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""
Physiological Equivalent Temperature (PET) -- the MEMI two-node solver.

PET is the air temperature of a standard indoor reference environment at which
the human heat balance closes with the *same core and skin temperature* as in
the outdoor environment being assessed (Hoeppe 1999). The underlying model is
MEMI (Munich Energy-balance Model for Individuals): a steady-state core/skin
two-node model plus a clothing node.

This module is internal and works in the units of the source papers (degC, hPa,
whole-body W). The public entry point is ``thermofeel.calculate_pet``, which
applies the library's Kelvin-in/Kelvin-out contract.

Vectorisation
-------------
The published implementations hand the 3x3 non-linear system in
(T_core, T_skin, T_clothing) to a per-point ``scipy.optimize.fsolve``, which
does not vectorise. The system is triangularisable, because the core-node
balance contains no clothing temperature and its forcing term contains no body
temperature at all:

    core     balance  ->  T_core     = f(T_skin)   closed form (quadratic)
    clothing balance  ->  T_clothing = g(T_skin)   1-D Newton (quartic)
    whole-body sum    ->  a scalar residual, monotone in T_skin

so the problem collapses to two nested 1-D *monotone* root finds, each solved by
bisection over the whole array at once: no per-element branching, no compaction,
and unconditional convergence given a valid bracket. The reformulation is
faithful, not an approximation -- it is the same second-order polynomial the
original Fortran solves (Walther & Goestchel 2018, Eqs. 47-53).

Model decisions
---------------
Documented in full in ``thermofeel.calculate_pet``; in brief, this follows the
VDI/Walther reference implementations: the clothing temperature is *frozen* at
its actual-environment value in the reference-environment stage, the Woodcock
clothing-aware evaporative resistance is used, the body-temperature weighting is
the constant alpha = 0.1, and activity metabolism is whole-body watts.

References
----------
- Hoeppe, P. (1999). The physiological equivalent temperature - a universal
  index for the biometeorological assessment of the thermal environment.
  Int J Biometeorol 43:71-75. https://doi.org/10.1007/s004840050118
- Walther, E. & Goestchel, Q. (2018). The P.E.T. comfort index: Questioning the
  model. Building and Environment 137:1-10.
  https://doi.org/10.1016/j.buildenv.2018.03.054
- Reference code: https://github.com/eddes/AREP
"""

import numpy as np

# --- physical constants -----------------------------------------------------
SBC = 5.67e-8  # Stefan-Boltzmann [W m-2 K-4] (as used by the reference codes)
EPS_SKIN = 0.99  # skin emissivity
EPS_CLO = 0.95  # clothing emissivity
L_VAP = 2.42e6  # latent heat of vaporisation [J kg-1]
C_AIR = 1010.0  # specific heat of air [J kg-1 K-1]
P0 = 1013.25  # reference pressure [hPa]

# --- physiological constants ------------------------------------------------
# RHO_CB is applied by every reference implementation as J L-1 K-1, although the
# source paper's nomenclature (rho_b = 1060 kg m-3, c_b = 4180 J kg-1 K-1) would
# give 4431. The reference codes define rho_b and then never use it. 3640 is
# kept here to reproduce them; see the docstring of calculate_pet.
RHO_CB = 3640.0  # effective volumetric heat capacity of blood [J L-1 K-1]
U_SK = 5.28  # skin tissue conductance [W m-2 K-1]
TC_SET = 36.6  # core set-point temperature [degC]
TSK_SET = 34.0  # skin set-point temperature [degC]
ALPHA = 0.1  # skin fraction of body mass (VDI: constant)
TBODY_SET = ALPHA * TSK_SET + (1.0 - ALPHA) * TC_SET  # 36.34 degC
Q_B_SET = 6.3  # set-point skin blood flow [L m-2 h-1]
Q_B_MAX = 90.0  # maximum skin blood flow [L m-2 h-1]
C_DIL = 75.0  # vasodilatation coefficient [L m-2 h-1 K-1]
C_STR = 0.5  # vasoconstriction coefficient [K-1]
C_SW = 304.94  # sweating coefficient [g m-2 h-1 K-1]
SW_MAX = 500.0  # sweat-rate cap [g m-2 h-1]
LR = 1.67  # Lewis ratio [K hPa-1]
I_M = 0.38  # Woodcock permeability index [-]

# --- the reference (indoor) environment that defines PET --------------------
REF_VPA = 12.0  # water vapour pressure [hPa]
REF_V = 0.1  # air velocity [m s-1]
REF_CLO = 0.9  # clothing insulation [clo]
REF_M_ACT = 80.0  # activity metabolism [W], whole body

# --- solver settings --------------------------------------------------------
# 26 bisections x 4 Newton gives ~5e-6 K against a fully converged reference,
# about four orders of magnitude tighter than any meaningful PET tolerance.
N_BISECT = 26
N_NEWTON = 4
# Brackets are deliberately wider than the reference implementations' [-40, 60]:
# at the cold end of the operational envelope (around -50 degC air with dry air
# and strong wind) the balance closes at a skin temperature just below -40 degC,
# and the narrow bracket would report NaN for a root that does exist. The
# solution there is physically meaningless - skin does not sit at -40 degC - but
# reporting it is consistent with the library's "out-of-range inputs return the
# raw value, the caller masks" contract, and keeps NaN meaning "no root at all".
TSK_BRACKET = (-70.0, 70.0)  # skin-temperature bracket, stage A [degC]
TX_BRACKET = (-90.0, 120.0)  # reference air-temperature bracket, stage B [degC]
# Cache-blocking: keeps the working set in L2. Measured optimum on the
# development machine; correctness is independent of this value.
CHUNK = 16384


def _svp_magnus_hpa(t_c):
    """Saturation vapour pressure over water [hPa], Magnus form.

    Same constants as ``calculate_nonsaturation_vapour_pressure`` (Bureau of
    Meteorology / VDI); evaluated in Celsius here because the solver works
    internally in the units of the source papers.
    """
    return 6.105 * np.exp(17.27 * t_c / (237.7 + t_c))


def _dubois_area(mass, height):
    """DuBois body-surface area [m2] (Walther & Goestchel Eq. 8)."""
    return 0.203 * mass**0.425 * height**0.725


def _basal_metabolism(mass, height, age, sex):
    """Basal metabolic rate [W], whole body."""
    if sex == "male":
        return (
            3.45
            * mass**0.75
            * (
                1.0
                + 0.004 * (30.0 - age)
                + 0.010 * (height * 100.0 / mass ** (1.0 / 3.0) - 43.4)
            )
        )
    return (
        3.19
        * mass**0.75
        * (
            1.0
            + 0.004 * (30.0 - age)
            + 0.018 * (height * 100.0 / mass ** (1.0 / 3.0) - 42.1)
        )
    )


def _geometry(clo, height, a_du):
    """Clothing geometry; depends only on (clo, height) so it is hoisted out of
    both root finds. Returns (f_cl, f_acl, a_clo, r_cl, h_cl)."""
    f_cl = 1.0 + 0.31 * clo  # Burton area increase factor
    f_acl = (173.51 * clo - 2.36 - 100.76 * clo * clo + 19.28 * clo**3) / 100.0
    a_clo = a_du * f_acl + a_du * (f_cl - 1.0)
    f_acl = np.minimum(f_acl, 1.0)
    r_cl = clo / 6.45  # clothing resistance [m2 K W-1]
    y = np.where(
        clo >= 2.0,
        1.0,
        np.where(clo > 0.6, (height - 0.2) / height, np.where(clo > 0.3, 0.5, 0.1)),
    )
    r1 = f_acl * a_du / (2.0 * np.pi * height * y)  # inner cylinder radius
    r2 = a_du * (f_cl - 1.0 + f_acl) / (2.0 * np.pi * height * y)  # outer radius
    h_cl = (
        2.0
        * np.pi
        * height
        * y
        * (r2 - r1)
        / (r_cl * np.log(r2 / r1) * a_clo)  # clothing conductance [W m-2 K-1]
    )
    return f_cl, f_acl, a_clo, r_cl, h_cl


def _core_from_skin(t_sk, q):
    """Closed-form root of the core-node balance ``K(Tc, Tsk) (Tc - Tsk) = q``.

    The conductance ``K`` is piecewise in ``Tc`` (set flow / vasodilated /
    saturated) and monotone increasing, so the left-hand side is monotone in
    ``Tc`` and exactly one branch is valid.
    """
    s = 1.0 / (1.0 + C_STR * np.maximum(TSK_SET - t_sk, 0.0))
    k_set = RHO_CB * Q_B_SET * s / 3600.0 + U_SK
    k_max = RHO_CB * Q_B_MAX / 3600.0 + U_SK
    tc_set_branch = t_sk + q / k_set
    tc_max_branch = t_sk + q / k_max
    # vasodilated branch: quadratic in u = Tc - TC_SET
    a1 = RHO_CB * C_DIL * s / 3600.0
    m = TC_SET - t_sk
    b = a1 * m + k_set
    disc = np.maximum(b * b - 4.0 * a1 * (k_set * m - q), 0.0)
    tc_dil_branch = TC_SET + (-b + np.sqrt(disc)) / (2.0 * a1)
    saturated = (Q_B_SET + C_DIL * np.maximum(tc_max_branch - TC_SET, 0.0)) * s >= (
        Q_B_MAX
    )
    return np.where(
        saturated,
        tc_max_branch,
        np.where(tc_set_branch <= TC_SET, tc_set_branch, tc_dil_branch),
    )


def _solve_tcl(t_sk, t_a, t_mrt, h_c, a_clo, a_du, f_eff, h_cl, t0):
    """Newton iteration on the clothing-node balance (a quartic in T_clothing).

    The balance is monotone decreasing in ``t_cl``, so Newton from any start in
    the physical range converges; N_NEWTON=4 reaches ~1e-6 K.
    """
    t_cl = t0
    k_c = h_c * a_clo / a_du
    k_r = f_eff * a_clo * EPS_CLO * SBC / a_du
    tmrt4 = (t_mrt + 273.15) ** 4
    for _ in range(N_NEWTON):
        tk = t_cl + 273.15
        f = k_c * (t_a - t_cl) + k_r * (tmrt4 - tk**4) + h_cl * (t_sk - t_cl)
        t_cl = t_cl - f / (-k_c - 4.0 * k_r * tk**3 - h_cl)
    return t_cl


def _h_conv(va, p_atm_hpa):
    """Convective heat transfer coefficient [W m-2 K-1] (standing posture)."""
    return (2.67 + 6.5 * va**0.67) * (p_atm_hpa / P0) ** 0.55


def _q_resp(h_met, t_a, vpa, p_atm_hpa):
    """Respiratory sensible + latent heat exchange [W m-2]."""
    t_exp = 0.47 * t_a + 21.0  # expired-air temperature [degC]
    rtv = h_met * 1.44e-6  # respiratory mass flow [kg m-2 s-1]
    vp_exp = 6.11 * 10.0 ** (7.45 * t_exp / (235.0 + t_exp))
    return (
        C_AIR * (t_a - t_exp) * rtv + 0.623 * L_VAP / p_atm_hpa * (vpa - vp_exp) * rtv
    )


def _balance(t_c, t_sk, t_cl, t_a, t_mrt, vpa, h_c, geom, a_du, f_eff, h_met, q_res):
    """Whole-body heat balance residual [W m-2]; zero when the body is in
    steady state. This is the sum of the core, skin and clothing node balances,
    in which the internal conduction terms cancel."""
    f_cl, f_acl, a_clo, r_cl, _ = geom

    # --- thermoregulation
    t_b = ALPHA * t_sk + (1.0 - ALPHA) * t_c
    m_sw = np.minimum(C_SW * np.maximum(t_b - TBODY_SET, 0.0), SW_MAX)
    e_sw = (L_VAP / 1000.0) * m_sw / 3600.0  # sweat evaporation [W m-2]

    # --- evaporation
    p_vsk = _svp_magnus_hpa(t_sk)
    f_pcl = 1.0 / (1.0 + 0.92 * h_c * r_cl)
    e_max = LR * h_c * f_pcl * (p_vsk - vpa)
    e_max = np.where(e_max == 0.0, 1e-3, e_max)  # reference-code guard
    w = np.minimum(e_sw / e_max, 1.0)  # skin wettedness
    r_ecl = (1.0 / (f_cl * h_c) + r_cl) / (LR * I_M)  # Woodcock resistance
    e_diff = (1.0 - w) * (p_vsk - vpa) / r_ecl
    e_vap = -(e_diff + np.maximum(e_sw, 0.0))

    # --- radiation and convection, bare skin and clothed fractions
    r_bare = (
        f_eff
        * (1.0 - f_acl)
        * EPS_SKIN
        * SBC
        * ((t_mrt + 273.15) ** 4 - (t_sk + 273.15) ** 4)
    )
    r_clo = (
        f_eff
        * (a_clo / a_du)
        * EPS_CLO
        * SBC
        * ((t_mrt + 273.15) ** 4 - (t_cl + 273.15) ** 4)
    )
    c_bare = h_c * (t_a - t_sk) * (1.0 - f_acl)
    c_clo = h_c * (t_a - t_cl) * a_clo / a_du

    return h_met + q_res + r_bare + r_clo + c_bare + c_clo + e_vap


def _pet_chunk(
    t_a, t_mrt, va, vpa, clo, m_act_w, a_du, m_bas, p_atm_hpa, f_eff, height
):
    """PET [degC] for one contiguous block of points."""
    # ---------------------------------------------------------------- stage A
    # Solve (T_core, T_skin, T_clothing) in the actual environment.
    h_c_a = _h_conv(va, p_atm_hpa)
    geom_a = _geometry(clo, height, a_du)
    h_met_a = (m_act_w + m_bas) / a_du
    q_res_a = _q_resp(h_met_a, t_a, vpa, p_atm_hpa)
    forcing = h_met_a + q_res_a  # independent of every body temperature

    def residual_a(t_sk):
        t_c = _core_from_skin(t_sk, forcing)
        t_cl = _solve_tcl(
            t_sk,
            t_a,
            t_mrt,
            h_c_a,
            geom_a[2],
            a_du,
            f_eff,
            geom_a[4],
            0.5 * (t_sk + t_a),
        )
        return _balance(
            t_c,
            t_sk,
            t_cl,
            t_a,
            t_mrt,
            vpa,
            h_c_a,
            geom_a,
            a_du,
            f_eff,
            h_met_a,
            q_res_a,
        )

    lo = np.full(t_a.shape, TSK_BRACKET[0])
    hi = np.full(t_a.shape, TSK_BRACKET[1])
    # residual_a is monotone decreasing in t_sk: a root exists iff r(lo) >= 0
    # and r(hi) <= 0. Elements failing this are reported as NaN rather than
    # silently returning a bracket edge.
    bracketed_a = (residual_a(lo) >= 0.0) & (residual_a(hi) <= 0.0)
    for _ in range(N_BISECT):
        mid = 0.5 * (lo + hi)
        above = residual_a(mid) > 0.0
        lo = np.where(above, mid, lo)
        hi = np.where(above, hi, mid)
    t_sk = 0.5 * (lo + hi)
    t_c = _core_from_skin(t_sk, forcing)
    t_cl = _solve_tcl(
        t_sk,
        t_a,
        t_mrt,
        h_c_a,
        geom_a[2],
        a_du,
        f_eff,
        geom_a[4],
        0.5 * (t_sk + t_a),
    )

    # ---------------------------------------------------------------- stage B
    # Find the reference-environment air temperature reproducing (T_core,
    # T_skin). The clothing temperature is deliberately FROZEN at its
    # actual-environment value: see calculate_pet's docstring.
    h_c_r = _h_conv(np.asarray(REF_V), p_atm_hpa)
    geom_r = _geometry(np.asarray(REF_CLO), height, a_du)
    h_met_r = (REF_M_ACT + m_bas) / a_du

    def residual_b(tx):
        return _balance(
            t_c,
            t_sk,
            t_cl,
            tx,
            tx,
            REF_VPA,
            h_c_r,
            geom_r,
            a_du,
            f_eff,
            h_met_r,
            _q_resp(h_met_r, tx, REF_VPA, p_atm_hpa),
        )

    lo = np.full(t_a.shape, TX_BRACKET[0])
    hi = np.full(t_a.shape, TX_BRACKET[1])
    # residual_b is monotone increasing in tx
    bracketed_b = (residual_b(lo) <= 0.0) & (residual_b(hi) >= 0.0)
    for _ in range(N_BISECT):
        mid = 0.5 * (lo + hi)
        below = residual_b(mid) < 0.0
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)

    pet_c = 0.5 * (lo + hi)
    return np.where(bracketed_a & bracketed_b, pet_c, np.nan)


def pet(
    t_a,
    t_mrt,
    va,
    vpa,
    m_act_w,
    clo,
    height,
    mass,
    age,
    sex,
    p_atm_hpa,
    f_eff,
):
    """PET [degC] from inputs in the units of the source papers.

    Internal entry point; ``thermofeel.calculate_pet`` is the public one.
    """
    t_a, t_mrt, va, vpa, clo, m_act_w, height, mass, age, p_atm_hpa, f_eff = (
        np.broadcast_arrays(
            *[
                np.asarray(x, dtype=float)
                for x in (
                    t_a,
                    t_mrt,
                    va,
                    vpa,
                    clo,
                    m_act_w,
                    height,
                    mass,
                    age,
                    p_atm_hpa,
                    f_eff,
                )
            ]
        )
    )

    a_du = _dubois_area(mass, height)
    m_bas = _basal_metabolism(mass, height, age, sex)

    # Elements that cannot be evaluated: non-finite inputs, or a clothing
    # insulation of zero or less, at which the cylinder geometry is singular
    # (the reference implementations divide by zero there).
    valid = (
        np.isfinite(t_a)
        & np.isfinite(t_mrt)
        & np.isfinite(va)
        & np.isfinite(vpa)
        & np.isfinite(clo)
        & np.isfinite(m_act_w)
        & np.isfinite(a_du)
        & np.isfinite(m_bas)
        & (clo > 0.0)
        & (va >= 0.0)
    )
    # substitute a benign value on invalid elements so the arithmetic below
    # cannot raise; those elements are overwritten with NaN at the end
    safe_clo = np.where(valid, clo, REF_CLO)
    args = [np.where(valid, a, 0.0) for a in (t_a, t_mrt, va, vpa, m_act_w)]

    out = np.empty(t_a.shape, dtype=float)
    flat = out.reshape(-1)
    n = flat.size
    fa, fm, fv, fp, fw = (a.reshape(-1) for a in args)
    fc = safe_clo.reshape(-1)
    f_adu, f_mbas, f_p, f_feff, f_h = (
        np.broadcast_to(x, t_a.shape).reshape(-1)
        for x in (a_du, m_bas, p_atm_hpa, f_eff, height)
    )
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        for start in range(0, n, CHUNK):
            s = slice(start, min(start + CHUNK, n))
            flat[s] = _pet_chunk(
                fa[s],
                fm[s],
                fv[s],
                fp[s],
                fc[s],
                fw[s],
                f_adu[s],
                f_mbas[s],
                f_p[s],
                f_feff[s],
                f_h[s],
            )
    return np.where(valid, out, np.nan)

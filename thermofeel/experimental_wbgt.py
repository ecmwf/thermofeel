# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import numpy as np

import thermofeel as thermofeel

# Lijigren WBGT
global Pr, rair, diamglobe, emisglobe, diamwick, emiswick, ratio, cp

mair = 28.97  # molecular weight of dry air (grams per mole)
mh2o = 18.015  # molecular weight of water vapor (grams per mole)
rgas = 8314.34  # ideal gas constant (J/kg mol · K)
cp = 1003.5  # Specific heat capacity of air at constant pressure (J·kg-1·K-1)
stefanb = 0.000000056696  # stefan-boltzmann constant
ratio = cp * mair * (mh2o ** (-1))
rair = rgas * (mair ** (-1))
Pr = cp * ((cp + 1.25 * rair) ** (-1))  # Prandtl number

# globe constants
diamglobe = 0.0508  # diameter of globe (m)
emisglobe = 0.95  # emissivity of globe
albglobe = 0.05  # albedo of globe

# wick constants
emiswick = 0.95  # emissivity of the wick
albwick = 0.4  # albedo of the wick
diamwick = 0.007  # diameter of the wick
lenwick = 0.0254  # length of the wick

# surface constant
albsfc = 0.45


def viscosity(t2m):
    """
    Calculate viscosity
    :param t2m: (float array) 2m temperature [K]
    returns  air viscosity [kg/(m/s)]
    https://www.tandfonline.com/doi/abs/10.1080/15459620802310770

    """
    omega = 1.2945 - t2m / 1141.176470588
    visc = 0.0000026693 * (np.sqrt(28.97 * t2m)) * ((13.082689 * omega) ** (-1))
    return visc


def thermcond(t2m, cp, rair):
    """
    Calculate thermal conductivity
    :param t2m: (float array) 2m temperature [K]
    :param cp: (int) Specific heat capacity of air at constant pressure [J^kg-1^K-1]
    :param rair: (int) physical constant
    returns  thermal conductivity
    https://www.tandfonline.com/doi/abs/10.1080/15459620802310770

    """
    tc = (cp + 1.25 * rair) * viscosity(t2m)
    return tc


def emisatm(t2m, rh, ps):
    """

    :param t2m: (float array) 2m temperature [K]
    :param rh:
    :param ps: (float array) Surface Pressure [Pa]
    returns atmospheric emissivity

    https://www.tandfonline.com/doi/abs/10.1080/15459620802310770
    """
    esat = thermofeel.calculate_saturation_vapour_pressure(t2m)

    e = rh * 0.01 * (esat * 0.01)
    emis_atm = 0.575 * (e**0.143)
    return emis_atm


def diffusivity(t2m, ps):
    """

    :param t2m: (float array) 2m temperature [K]
    :param ps: (float array) Surface Pressure [Pa]
    returns diffusivity of water vapour in air [m2/s]
    https://www.tandfonline.com/doi/abs/10.1080/15459620802310770
    """
    diff = (
        2.471773765165648e-05
        * ((t2m * 0.0034210563748421257) ** 2.334)
        * ((ps / 101325) ** (-1))
    )
    return diff


def h_evap(t2m):
    """

    :param t2m:  (float array) 2m temperature [K]
    returns heat of evaporation (J/(kg K))
    """
    hevap = 1665134.5 + 2370.0 * t2m

    return hevap


def h_sphere_and_cylinder_in_air(t2m, ps, va, diamglobe, diamwick, Pr, cp, rair):
    """

    :param t2m: (float array) 2m temperature [K]
    :param ps: (float array) Surface Pressure [Pa]
    :param va: (float array) 10 wind speed [m/s]
    :param diamglobe: constant
    :param diamwick: constant
    :param Pr: constant
    :param cp: constant
    :param rair: constant

    returns convective heat tranfer coefficient for flow around a sphere [W/m2 K]
    https://www.tandfonline.com/doi/abs/10.1080/15459620802310770
    """
    thermcon = thermcond(t2m, cp, rair)
    density = ps * ((rair * t2m) ** (-1))
    Re_globe = va * density * diamglobe * ((viscosity(t2m)) ** (-1))
    Nu_globe = 2 + 0.6 * np.sqrt(Re_globe) * np.power(Pr, 0.3333)
    h_globe = Nu_globe * thermcon * (diamglobe ** (-1))

    Re_wick = va * density * diamwick * ((viscosity(t2m)) ** (-1))
    Nu_wick = 0.281 * (Re_wick**0.6) * (Pr**0.44)
    h_wick = Nu_wick * thermcon * (diamwick ** (-1))

    return h_globe, h_wick


def bgt_lijigren(
    t2m, rh, ps, ssrd, fdir, cossza, eatm, diamglobe, diamwick, Pr, va, cp, rair
):
    """

    :param t2m: (float array) 2m temperature [K]
    :param rh: Humidity Ratio by Vapor Partial Pressure
    :param ps: (float array) Surface Pressure [Pa]
    :param ssrd: is surface solar radiation downwards [J/m^-2]
    :param fdir: is total sky radiation [J/m^-2]
    :param cossza: is cosine of the solar zenith angle []
    :param eatm: constant
    :param diamglobe: constant
    :param diamwick: constant
    :param Pr: constant
    :param va: 10 meter wind speed [m/s]
    :param cp: constant
    :param rair: constant

    returns globe temperature from lijigren
    https://www.tandfonline.com/doi/abs/10.1080/15459620802310770
    """
    # equation of Tg that needs to be solved by iteration
    i = 0
    diamglobe = diamglobe
    diamwick = diamwick
    Pr = Pr
    va = va
    tsfc = t2m
    albsfc = 0.45
    albglobe = 0.05
    eglobe = 0.0508
    tglobe_prev = t2m  # first guess is 2m temperature
    tref = 0.5 * (tglobe_prev + t2m)
    esfc = 0.999
    h = h_sphere_and_cylinder_in_air(tref, ps, va, diamglobe, diamwick, Pr, cp, rair)
    h = h[0]
    while i < 50:
        tg = np.power(
            0.5 * (eatm * np.power(t2m, 4))
            + esfc * np.power(tsfc, 4)
            - h / (stefanb * eglobe) * (tglobe_prev - t2m)
            + ssrd
            / (2.0 * stefanb * eglobe)
            * (1.0 - albglobe)
            * (fdir * (1.0 / (2.0 * cossza) - 1.0) + 1.0 + albsfc),
            0.25,
        )
        tg_filter = np.where(np.abs(tg - tglobe_prev))
        tglobe_prev[tg_filter] = 0.9 * tglobe_prev[tg_filter] + 0.1 * tg[tg_filter]
        i = i + 1
    return tg


def wbt_lijigren(
    t2m, td, rh, ps, va, ssrd, fdir, cossza, rair, ratio, diamglobe, diamwick, Pr, cp
):
    """

    :param t2m: (float array) 2m temperature [K]
    :param td: (float array) 2m dew point temperature [K]
    :param rh: Humidity Ratio by Vapor Partial Pressure
    :param ps: (float array) Surface Pressure [Pa]
    :param va: 10 meter wind speed [m/s]
    :param ssrd: is surface solar radiation downwards [J/m^-2]
    :param fdir: is total sky radiation [J/m^-2]
    :param cossza: is cosine of the solar zenith angle
    :param rair: constant
    :param ratio: constant
    :param diamglobe: constant
    :param diamwick: constant
    :param Pr: constant
    :param cp: constant

    returns wet bulb temperature lijigren

    https://www.tandfonline.com/doi/abs/10.1080/15459620802310770
    """
    a = 0.56
    i = 0

    tref = t2m
    tsfc = t2m
    sza = np.arccos(cossza)
    eair = rh * thermofeel.calculate_saturation_vapour_pressure(t2m)
    twb_prev = td
    h = h_sphere_and_cylinder_in_air(tref, ps, va, diamglobe, diamwick, Pr, cp, rair)
    h = h[0]
    tref = 0.5 * (twb_prev + t2m)
    emis_atm = emisatm(t2m, rh, ps)
    emiswick = 0.95
    albwick = 0.4
    dwick = 0.007
    lwick = 0.0254
    emisfc = 0.999

    while i < 50:
        Fatm = stefanb * emiswick * (
            0.5 * (emis_atm * np.power(t2m, 4) + emisfc * np.power(tsfc, 4))
            - np.power(twb_prev, 4)
        ) + (1.0 - albwick) * ssrd * (
            (1.0 - fdir) * (1.0 + 0.25 * dwick / lwick)
            + fdir * ((np.tan(sza) / np.pi) + 0.25 * dwick / lwick)
            + albsfc
        )
        ewick = thermofeel.calculate_saturation_vapour_pressure(twb_prev)
        density = ps * 100 / (rair * tref)
        Sc = viscosity(tref) / (density * diffusivity(t2m, ps))
        twb = (
            t2m
            - h_evap(t2m) / ratio * (ewick - eair) / (ps - ewick) * np.power(Pr / Sc, a)
            + (Fatm / h * -1)
        )
        twb_filter = np.where(np.abs(twb - twb_prev))
        twb_prev[twb_filter] = 0.9 * twb_prev[twb_filter] + 0.1 * twb[twb_filter]

        i = i + 1

    return twb

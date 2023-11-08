import math

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from netCDF4 import Dataset, date2num, num2date

from thermofeel import (
    calculate_cos_solar_zenith_angle,
    calculate_relative_humidity_percent,
    calculate_saturation_vapour_pressure_multiphase,
)

# from datetime import datetime


to_radians = math.pi / 180.0

# physical constants
global M_AIR, M_H2O, R_GAS, Cp, STEFANB, RATIO, R_AIR, Pr
M_AIR = 28.97  # molecular weight of dry air (grams per mole)
M_H2O = 18.015  # molecular weight of water vapor (grams per mole)
R_GAS = 8314.34  # ideal gas constant (J/kg mol · K)
Cp = 1003.5  # Specific heat capacity of air at constant pressure (J·kg-1·K-1)
STEFANB = 0.000000056696  # stefan-boltzmann constant
RATIO = Cp * M_AIR / M_H2O
R_AIR = R_GAS * (M_AIR ** (-1))
Pr = Cp * ((Cp + 1.25 * R_AIR) ** (-1))  # Prandtl number

# globe constants
global D_GLOBE, EMIS_GLOBE, ALB_GLOBE
D_GLOBE = 0.0508  # diameter of globe (m)
EMIS_GLOBE = 0.95  # emissivity of globe
ALB_GLOBE = 0.05  # albedo of globe

# wick constants
global EMIS_WICK, ALB_WICK, D_WICK, L_WICK
EMIS_WICK = 0.95  # emissivity of the wick
ALB_WICK = 0.4  # albedo of the wick
D_WICK = 0.007  # diameter of the wick
L_WICK = 0.0254  # length of the wick

# surface constant
global EMIS_SFC, ALB_SFC
EMIS_SFC = 0.999
ALB_SFC = 0.45

global MAX_ITER, CONVERGENCE
MAX_ITER = 50
CONVERGENCE = 0.02

global MISSING_VALUE
MISSING_VALUE = -9999.0
MIN_SPEED = 0.13


def viscosity(t2m):
    """
    Calculate viscosity
    :param t2m: (float array) 2m temperature [K]
    returns  air viscosity [kg/(m/s)]

    """
    omega = 1.2945 - t2m / 1141.176470588
    visc = 0.0000026693 * (np.sqrt(28.97 * t2m)) * ((13.082689 * omega) ** (-1))
    return visc


def thermal_cond(tk):
    """
    Calculate thermal conductivity
    :param tk: (float array) temperature [K]
    :param cp: (int) Specific heat capacity of air at constant pressure (J·kg-1·K-1)
    :param R_AIR: (int) physical constant
    returns  air viscosity [kg/(m/s)]
    """
    tc = (Cp + 1.25 * R_AIR) * viscosity(tk)
    return tc


def emis_atm(t2m, rhpc):
    """
    Purpose: calculate the atmospheric emissivity.
    Reference: Oke (2nd edition), page 373.
    :param t2m: 2m temperature [K]
    :param rh: relative humidity must be in % [0-100]
    :return:
    """
    liquid = 0
    esat = calculate_saturation_vapour_pressure_multiphase(t2m, liquid)

    rh = rhpc * 0.01
    e = rh * esat
    emis_atm = 0.575 * (e**0.143)
    return emis_atm


def diffusivity(tk, ps):
    # tk: air temperature (K)
    # ps: surface pressure (Pa)
    # return diffusivity of water vapor in air (m2/s)

    diff = (
        2.471773765165648e-05
        * ((tk * 0.0034210563748421257) ** 2.334)
        * ((ps / 101325) ** (-1))
    )
    return diff


def evap(tk):
    # tk: air temperature (K)
    # return heat of evaporation (J/(kg K))
    hevap = 1665134.5 + 2370.0 * tk

    return hevap


def h_sphere_and_cylinder_in_air(t2m, ps, va, diamglobe, diamwick, Pr, cp, R_AIR):
    # tas: air temperature (K)
    # ps: surface pressure (Pa)
    # sfcwind: 2 meter wind (m/s)
    # return convective heat tranfer coefficient for flow around a sphere (W/(m2 K))
    # h = h_sphere_and_cylinder_in_air(tref,ps,va,diamglobe,diamwick,Pr,cp,R_AIR)
    thermcon = thermal_cond(t2m)
    density = ps * ((R_AIR * t2m) ** (-1))
    Re_globe = va * density * diamglobe * ((viscosity(t2m)) ** (-1))
    Nu_globe = 2 + 0.6 * np.sqrt(Re_globe) * np.power(Pr, 0.3333)
    h_globe = Nu_globe * thermcon * (diamglobe ** (-1))

    Re_wick = va * density * diamwick * ((viscosity(t2m)) ** (-1))
    Nu_wick = 0.281 * (Re_wick**0.6) * (Pr**0.44)
    h_wick = Nu_wick * thermcon * (diamwick ** (-1))

    return h_globe, h_wick


def h_sphere_in_air(diameter, Tair, Pair, speed):
    """
    Calculate the convective heat transfer coefficient, W/(m2 K), for flow around a sphere.
    diameter, sphere diameter, m
        Tair, air temperature, K
        Pair, barometric pressure, mb
        speed, fluid (air) speed, m/s
    Reference: Bird, Stewart, and Lightfoot (BSL), page 409.
    """
    density = Pair * 100.0 / (R_AIR * Tair)
    Re = max(speed, MIN_SPEED) * density * diameter / viscosity(Tair)
    Nu = 2.0 + 0.6 * math.sqrt(Re) * pow(Pr, 0.3333)
    return Nu * thermal_cond(Tair) / diameter


def h_cylinder_in_air(diameter, length, Tair, Pair, speed):
    """
    diameter, cylinder diameter, [m]
        length,   cylinder length, [m] - does not seem to be used
        Tair, air temperature, [K]
        Pair, barometric pressure, [hPa]
        speed,	fluid (wind) speed, [m/s]
    """
    # parameters from Bedingfield and Drew
    a = 0.56
    b = 0.281
    c = 0.4
    density = (Pair * 100.0) / (R_AIR * Tair)
    Re = max(speed, MIN_SPEED) * density * diameter / viscosity(Tair)
    Nu = b * pow(Re, (1.0 - c)) * pow(Pr, (1.0 - a))
    return Nu * thermal_cond(Tair) / diameter


def bgt_lijigren(t2m, rh, ps, ssrd, fdir, cossza, va):
    """
    Purpose: to calculate the globe temperature.
    Method by J. Liljegren, Argonne National Laboratory
    """
    Tair = t2m
    speed = va
    Pair = ps
    cza = cossza
    solar = ssrd

    Tsfc = Tair
    Tglobe_prev = Tair
    # first guess is the air temperature

    converged = False
    iter = 0
    while not converged and iter < MAX_ITER:
        iter = iter + 1
        Tref = 0.5 * (
            Tglobe_prev + Tair
        )  # evaluate properties at the average temperature
        h = h_sphere_in_air(D_GLOBE, Tref, Pair, speed)
        Tglobe_new = pow(
            0.5 * (emis_atm(Tair, rh) * pow(Tair, 4.0) + EMIS_SFC * pow(Tsfc, 4.0))
            - h / (STEFANB * EMIS_GLOBE) * (Tglobe_prev - Tair)
            + solar
            / (2.0 * STEFANB * EMIS_GLOBE)
            * (1.0 - ALB_GLOBE)
            * (fdir * (1.0 / (2.0 * cza) - 1.0) + 1.0 + ALB_SFC),
            0.25,
        )

        if abs(Tglobe_new - Tglobe_prev) < CONVERGENCE:
            converged = True

        Tglobe_prev = 0.9 * Tglobe_prev + 0.1 * Tglobe_new

    if converged:
        return Tglobe_new - 273.15
    else:
        return MISSING_VALUE

    # bgt = bgt_lijigren(t2m, rh, ps, ssrd/3600, fdir_frac, cossza, eatm, diamglobe, diamwick, Pr, va, cp, R_AIR)
    #
    # equation of Tg that needs to be solved by iteration
    # i = 0
    # diamglobe = diamglobe
    # diamwick= diamwick
    # Pr = Pr
    # va = va
    # tsfc = t2m
    # albsfc = 0.45
    # albglobe = 0.05
    # eglobe = 0.95
    # tglobe_prev = np.copy(t2m) #first guess is 2m temperature
    # tref = 0.5 * (tglobe_prev + t2m)
    # esfc = 0.999
    # h = h_sphere_and_cylinder_in_air(tref,ps,va,diamglobe,diamwick,Pr,cp,R_AIR)
    # h = h[0]
    # # while i < 50:
    # #     tg = np.power(0.5*( eatm * np.power(t2m,4)) + esfc * np.power(tsfc,4) \
    # #                 - h/(STEFANB*eglobe)*(tglobe_prev - t2m) \
    # #                     + ssrd/(2.*STEFANB*eglobe)*(1.-albglobe)\
    # #                         *(fdir*(1./(2.*cossza)-1.)+1.+ albsfc), 0.25)
    # #     tg_filter = np.where(np.abs(tg - tglobe_prev) < tg)
    # #     tglobe_prev[tg_filter] = 0.9 * tglobe_prev[tg_filter] + 0.1* tg[tg_filter]
    # #     print(np.average(tg))
    # #     i = i + 1
    # # tg = np.power(0.5*(eatm * np.power(t2m, 4)) + esfc * np.power(tsfc, 4) - h/(STEFANB*eglobe)*(tglobe_prev - t2m)
    #      + ssrd/(2.*STEFANB*eglobe)*(1.-albglobe) * (fdir*(1./(2.*cossza)-1.)+1.+ albsfc), 0.25)
    # tg = np.power(0.5 * ((1+eatm) * np.power(t2m, 4)) - h/(STEFANB*eglobe)*(tglobe_prev - t2m)
    #    + ssrd/(2.*STEFANB*eglobe)*(1.-albglobe) * (fdir*(1./(2.*cossza)-1.)+1.+ albsfc), 0.25)
    # tg_filter = np.where(np.abs(tg - tglobe_prev) < tg)
    # tglobe_prev[tg_filter] = 0.9 * tglobe_prev[tg_filter] + 0.1* tg[tg_filter]
    # print(np.average(tg))
    # return tg


def wbt_lijigren(t2m, td, rh, ps, va, ssrd, fdir, cossza):
    """
    Purpose: to calculate the natural wet bulb temperature.
    Method by J. Liljegren, Argonne National Laboratory
    """

    Tair = t2m
    Pair = ps
    speed = va
    solar = ssrd

    radiative = 1
    liquid = 0

    a = 0.56  # from Bedingfield and Drew

    Tsfc = Tair
    sza = math.acos(cossza)
    eair = rh * calculate_saturation_vapour_pressure_multiphase(t2m, liquid)
    Tdew = td
    Twb_prev = Tdew  # first guess is the dew point temperature

    converged = False
    iter = 0
    while not converged and iter < MAX_ITER:
        iter = iter + 1
        Tref = 0.5 * (Twb_prev + Tair)  # evaluate properties at the average temperature
        h = h_cylinder_in_air(D_WICK, L_WICK, Tref, Pair, speed)
        Fatm = STEFANB * EMIS_WICK * (
            0.5 * (emis_atm(Tair, rh) * pow(Tair, 4.0) + EMIS_SFC * pow(Tsfc, 4.0))
            - pow(Twb_prev, 4.0)
        ) + (1.0 - ALB_WICK) * solar * (
            (1.0 - fdir) * (1.0 + 0.25 * D_WICK / L_WICK)
            + fdir * ((math.tan(sza) / math.pi) + 0.25 * D_WICK / L_WICK)
            + ALB_SFC
        )
        ewick = calculate_saturation_vapour_pressure_multiphase(Twb_prev, liquid)
        density = Pair * 100.0 / (R_AIR * Tref)

        Sc = viscosity(Tref) / (density * diffusivity(Tref, Pair))  # Schmidt number
        print(f"Sc {Sc}")

        Twb_new = (
            Tair
            - evap(Tref) / RATIO * (ewick - eair) / (Pair - ewick) * pow(Pr / Sc, a)
            + (Fatm / h * radiative)
        )

        print(f"Twb_prev {Twb_prev} Twb_new {Twb_new} Increment {Twb_new-Twb_prev}")

        if abs(Twb_new - Twb_prev) < CONVERGENCE:
            converged = True
        Twb_prev = 0.9 * Twb_prev + 0.1 * Twb_new

    if converged:
        return Twb_new - 273.15
    else:
        return MISSING_VALUE

    # h = h_sphere_and_cylinder_in_air(tref,ps,va,diamglobe,diamwick,Pr,cp,R_AIR)
    # h = h[1]
    # tref = 0.5 * (twb_prev + t2m)
    # emis_atm = emis_atm(t2m,rh,ps)
    # emiswick =	0.95
    # albwick =	0.4
    # dwick = 0.007
    # lwick =	0.0254
    # emisfc = 0.999

    # converged = False

    # while i < 50:
    #     Fatm = STEFANB * emiswick *  ( 0.5*( emis_atm *np.power(t2m,4) + emisfc*np.power(tsfc,4) )
    #                                   - np.power(twb_prev,4) ) + (1.-albwick) * ssrd * \


# 	       ( (1.-fdir)*(1.+0.25*dwick/lwick) + fdir*((np.tan(sza)/np.pi)+0.25*dwick/lwick) + albsfc )
#     ewick = calculate_saturation_vapour_pressure(twb_prev)
#     density = ps * 100 / (R_AIR * tref)
#     Sc = viscosity(tref)/(density* diffusivity(t2m, ps))
#     twb = t2m - evap(t2m) / ratio * (ewick - eair)/(ps - ewick) * np.power(Pr/Sc,a) \
#         + (Fatm/h * -1)
#    # twb_filter = np.where(np.abs(twb - twb_prev) < twb)
#     #twb_prev[twb_filter] = 0.9 * twb_prev[twb_filter] + 0.1* twb[twb_filter]

#     i= i + 1

# return twb


def main():
    filename = "March2013LijiCake.nc"
    print("Loading NetCDF file {filename} --  Go and have a coffee ....")

    fh2 = Dataset(filename, mode="r")
    # print(fh2.variables)
    lons = fh2.variables["longitude"][:]  # degrees
    lats = fh2.variables["latitude"][:]  # degrees
    nctime = fh2.variables["time"]  # get values
    dtime = num2date(nctime[1], nctime.units)
    # print(dtime)
    ssrd = fh2.variables["ssrd"][1]  # Jm-2
    fdir = fh2.variables["fdir"][1]  # Jm-2

    lon_mg, lat_mg = np.meshgrid(lons, lats)
    va = np.sqrt(fh2.variables["u10"][1] ** 2 + fh2.variables["v10"][1] ** 2)  # ms-1
    t2m = fh2.variables["t2m"][1]  # K
    td = fh2.variables["d2m"][1]  # K
    ps = fh2.variables["sp"][1] / 100  # hPa
    # rh = 0.62198 * calculate_saturation_vapour_pressure(t2m) / (ps - calculate_saturation_vapour_pressure(t2m))
    rh = calculate_relative_humidity_percent(t2m, td)  # %
    # print(rh)
    # eatm = emis_atm(t2m, rh)  # unitless
    # print(eatm)
    cossza = calculate_cos_solar_zenith_angle(
        lat=lat_mg, lon=lon_mg, y=dtime.year, m=dtime.month, d=dtime.day, h=dtime.hour
    )  # unitless

    #
    # NOTE: CHECK THIS MAKES SENSE OR HAS IMPACT IN RESULTS
    #
    cossza = np.where(
        cossza == 0, 0.25, cossza
    )  # TODO: Check it tolerance checking is needed instead!

    #
    # NOTE: CHECK THIS MAKES SENSE OR HAS IMPACT IN RESULTS
    #
    fdir_frac = np.where(fdir > 0, fdir / ssrd, 0)  # unitless

    ssrd = ssrd / 3600.0
    bgt = np.zeros(len(cossza))
    wbt = np.zeros(len(cossza))
    for i in range(len(cossza)):
        t2m_ = float(t2m[i][1])
        td_ = float(td[i][1])
        rh_ = float(rh[i][1])
        ps_ = float(ps[i][1])
        ssrd_ = float(ssrd[i][1])
        fdir_frac_ = float(fdir_frac[i][1])
        cossza_ = float(cossza[i][1])
        va_ = float(va[i][1])
        bgt[i] = bgt_lijigren(t2m_, rh_, ps_, ssrd_, fdir_frac_, cossza_, va_)
        wbt[i] = wbt_lijigren(t2m_, td_, rh_, ps_, va_, ssrd_, fdir_frac_, cossza_)
        print(f"i {i} -> bgt {bgt[i]} wbt {wbt[i]}")

    WBGTvar = 0.7 * (wbt - 273.15) + 0.2 * (bgt - 273.15) + 0.1 * (t2m - 273.15)

    # print (WBGTvar)
    # h = h_sphere_and_cylinder_in_air(t2m, ps, va,diamglobe,diamwick,Pr,cp,R_AIR)

    ########################

    # fig = plt.figure()
    # ax = plt.axes(projection=ccrs.AlbersEqualArea())
    # ax.coastlines()
    # filled_c = plt.pcolormesh(lon_mg, lat_mg, (h[0]),
    #                           transform=ccrs.PlateCarree(),
    #                           cmap="RdBu_r")
    # fig.colorbar(filled_c, orientation="horizontal",extend="both")
    # plt.savefig("h0.png")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax = plt.axes(projection=ccrs.AlbersEqualArea())
    # ax.coastlines()
    # filled_c = plt.pcolormesh(lon_mg, lat_mg, (h[1]),
    #                           transform=ccrs.PlateCarree(),
    #                           cmap="RdBu_r")
    # fig.colorbar(filled_c, orientation="horizontal",extend="both")
    # plt.savefig("h1.png")
    # plt.clf()

    fig = plt.figure()
    ax = plt.axes(projection=ccrs.AlbersEqualArea())
    ax.coastlines()
    filled_c = plt.pcolormesh(
        lon_mg,
        lat_mg,
        bgt,
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
        vmin=273,
        vmax=350,
    )
    fig.colorbar(filled_c, orientation="horizontal", extend="both")
    plt.savefig("bgt.png")
    plt.clf()

    # fig = plt.figure()
    # ax = plt.axes(projection=ccrs.PlateCarree())
    # ax.coastlines()
    # filled_c = plt.pcolormesh(lon_mg, lat_mg, (wbt),
    #                           transform=ccrs.PlateCarree(),
    #                           cmap="RdBu_r")
    # fig.colorbar(filled_c, orientation="horizontal",extend="both")
    # plt.savefig("wbt.png")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax = plt.axes(projection=ccrs.PlateCarree())
    # ax.coastlines()
    # filled_c = plt.pcolormesh(lon_mg, lat_mg, (WBGTvar),
    #                           transform=ccrs.PlateCarree(),
    #                           cmap="RdBu_r")
    # fig.colorbar(filled_c, orientation="horizontal",extend="both")
    # plt.savefig("wbgtvar.png")
    # plt.clf()
    #
    # fig = plt.figure()
    # ax = plt.axes(projection=ccrs.PlateCarree())
    # ax.coastlines()
    # filled_c = plt.pcolormesh(lon_mg, lat_mg, (t2m),
    #                           transform=ccrs.PlateCarree(),
    #                           cmap="RdBu_r")
    # fig.colorbar(filled_c, orientation="horizontal",extend="both")
    # plt.savefig("t2m.png")
    # plt.clf()

    time = pd.date_range(start="2013-03-01", end="2013-03-31").to_pydatetime().tolist()

    ncout = Dataset("wbgt2013liji.nc", "w", format="NETCDF4")
    ncout.createDimension("latitude", len(lat_mg))
    ncout.createDimension("longitude", len(lon_mg))
    ncout.createDimension("time", len(time))

    ncout_latitude = ncout.createVariable("latitude", "f4", ("latitude",))
    ncout_latitude[:] = lat_mg
    ncout_longitude = ncout.createVariable("longitude", "f4", ("longitude",))
    ncout_longitude[:] = lon_mg

    units = "hours in march 2013"
    ncout_time = ncout.createVariable("time", "i4", ("time",))
    ncout_time[:] = date2num(time, units)
    ncout_time.units = units

    ncout_dist = ncout.createVariable("WBGT", "f4", ("time", "latitude", "longitude"))
    ncout_dist.long_name = "WBGT"
    ncout_dist.units = "C"
    ncout_dist[:] = WBGTvar

    ncout.close()


if __name__ == "__main__":
    main()

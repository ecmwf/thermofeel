# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""
  thermofeel is a library to calculate human thermal comfort indexes.

    Currently calculates the thermal indexes:
    * Universal Thermal Climate Index
    * Mean Radiant Temperature
    * Mean Radiant Temperature from Wet Bulb Globe Temperature
    * Heat Index Simplified
    * Heat Index Adjusted
    * Humidex
    * Apparent Temperature
    * Wind Chill
    * Net Effective Temperature

    In support of the above indexes, it also calculates:
    * Solar Declination Angle
    * Solar Zenith Angle
    * Relative Humidity Percentage
    * Saturation vapour pressure
    * Wet Bulb Globe Temperature Simple
    * Wet Bulb Globe Temperature

  """

import math

import numpy as np

from .helpers import (
    __wrap,
    farenheit_to_celcius,
    kelvin_to_farenheit,
    kPa_to_hPa,
    to_julian_date,
    to_radians,
)


# solar declination angle [degrees] + time correction for solar angle
def solar_declination_angle(jd, h):
    g = (360 / 365.25) * (jd + (h / 24))  # fractional year g in degrees
    while g > 360:
        g = g - 360
    grad = g * to_radians
    # declination in [degrees]
    d = (
        0.396372
        - 22.91327 * math.cos(grad)
        + 4.025430 * math.sin(grad)
        - 0.387205 * math.cos(2 * grad)
        + 0.051967 * math.sin(2 * grad)
        - 0.154527 * math.cos(3 * grad)
        + 0.084798 * math.sin(3 * grad)
    )
    # time correction in [ h.degrees ]
    tc = (
        0.004297
        + 0.107029 * math.cos(grad)
        - 1.837877 * math.sin(grad)
        - 0.837378 * math.cos(2 * grad)
        - 2.340475 * math.sin(2 * grad)
    )
    return d, tc


def calculate_relative_humidity_percent(t2m, td):
    """
    Calculate relative humidity in percent
    :param t2m: (float array) 2m temperature [K]
    :param td: (float array) dew point temperature [K]

    returns relative humidity [%]
    """
    t2m = __wrap(t2m)
    td = __wrap(td)

    t2m = kelvin_to_celcius(t2m)
    td = kelvin_to_celcius(td)

    # saturated vapour pressure
    es = 6.11 * 10.0 ** (7.5 * t2m / (237.3 + t2m))
    # vapour pressure
    e = 6.11 * 10.0 ** (7.5 * td / (237.3 + td))
    rh = (e / es) * 100
    return rh


def calculate_saturation_vapour_pressure(t2m):
    """
    Calculate saturation vapour pressure over water
    :param t2m: (float array) 2m temperature [K]
     returns relative humidity [hPa]
    http://www.thunderscientific.com/tech_info/reflibrary/its90formulas.pdf
    """

    tk = __wrap(t2m)

    g = [
        -2.8365744e3,
        -6.028076559e3,
        1.954263612e1,
        -2.737830188e-2,
        1.6261698e-5,
        7.0229056e-10,
        -1.8680009e-13,
        2.7150305,
    ]
    ess = g[7] * np.log(tk)
    for i in range(7):
        ess += g[i] * np.power(tk, (i - 2))

    ess = np.exp(ess) * 0.01  # hPa

    return ess


def calculate_cos_solar_zenith_angle(h, lat, lon, y, m, d):
    """
    calculate solar zenith angle
    :param lat: (float array) latitude [degrees]
    :param lon: (float array) longitude [degrees]
    :param y: year [int]
    :param m: month [int]
    :param d: day [int]
    :param h: hour [int]

    https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1002/2015GL066868

    see also:
    http://answers.google.com/answers/threadview/id/782886.html

    returns cosine of the solar zenith angle
    """

    # convert to julian days counting from the beginning of the year
    jd_ = to_julian_date(d, m, y)  # julian date of data
    jd11_ = to_julian_date(1, 1, y)  # julian date 1st Jan
    jd = jd_ - jd11_ + 1  # days since start of year

    # declination angle + time correction for solar angle
    d, tc = solar_declination_angle(jd, h)
    drad = d * to_radians

    latrad = lat * to_radians

    sindec_sinlat = np.sin(drad) * np.sin(latrad)
    cosdec_coslat = np.cos(drad) * np.cos(latrad)

    # solar hour angle [h.deg]
    sharad = ((h - 12) * 15 + lon + tc) * to_radians
    csza = sindec_sinlat + cosdec_coslat * np.cos(sharad)

    return np.clip(csza, 0, None)


def calculate_cos_solar_zenith_angle_integrated(lat, lon, y, m, d, h, base, step):
    """
    calculate solar zenith angle
    :param lat: (int array) latitude [degrees]
    :param lon: (int array) longitude [degrees]
    :param y: year [int]
    :param m: month [int]
    :param d: day [int]
    :param h: hour [int]
    :param base: base time of forecast enum [0,6,12,18]
    :param step: step interval of forecast enum [0,1,3,6,..]

    https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1002/2015GL066868

    returns cosine of the solar zenith angle [degrees]
    """

    # OLD CODE
    # maxx = 0
    # base_offset = 0
    # h_offset = 0
    # # step_offset = 0
    # if base == 0:
    #     base_offset = 0
    # if base == 6:
    #     base_offset = 594
    # if base == 18:
    #     base_offset = 1806
    # if base == 12:
    #     base_offset = 1212  # 1200+12

    # if step == 1:
    #     h_offset = 0.5
    #     maxx = 32.0
    # if step == 3:
    #     h_offset = 1.5
    #     maxx = 14.0
    # if step == 6:
    #     h_offset = 3
    #     maxx = 100.0

    # hh = h - base_offset - h_offset

    # POTENTIAL BUG in OLD CODE ABOVE -- POSSIBLE FIX BELOW -- Needs checking
    # This ensures we can integrate but doesn't use hours offsets
    hh = h
    maxx = 100.0

    # convert to julian days counting from the beginning of the year
    jd_ = to_julian_date(d, m, y)  # julian date of data
    jd11_ = to_julian_date(1, 1, y)  # julian date 1st Jan
    jd = jd_ - jd11_ + 1  # days since start of year

    # declination angle + time correction for solar angle
    d, tc = solar_declination_angle(jd, hh)
    drad = d * to_radians

    # solar hour angle
    sha = (hh - 12) * 15 + (lon + 180) / 2 + tc  # solar hour angle
    sharad = sha * to_radians

    latrad = lat * to_radians
    lonrad = (lon + 180) / 2 * to_radians
    accumulationperiod = step  # hours

    # 2*3.1415=day in rad; /24=day hour in rad; *accumulationperiod/2;
    # = half of the (radiation) accumulation period - SUN IN THE MIDDLE!
    zhalftimestep = (2 * np.pi) / 24 * accumulationperiod / 2
    zsolartimestart = sharad - zhalftimestep
    zsolartimeend = sharad + zhalftimestep
    ztandec = math.sin(drad) / max((math.cos(drad), 1.0e-12))
    zcoshouranglesunset = (
        -ztandec * np.sin(latrad) / np.clip(np.cos(latrad), 1.0e-12, None)
    )
    zsindecsinlat = math.sin(drad) * np.sin(latrad)
    zcosdeccoslat = math.cos(drad) * np.cos(latrad)

    def solar_zenith_angle_average(
        lonrad,
        zcoshouranglesunset,
        zsindecsinlat,
        zcosdeccoslat,
        zsolartimestart,
        zsolartimeend,
        sha,
    ):
        # start and end hour
        def horizon(val, rrange, lonrad, zsolartimestart, zsolartimeend, maxx):
            if val < rrange * math.pi:
                zhouranglestart = zsolartimestart + lonrad - ((rrange - 1.0) * math.pi)
                zhourangleend = zsolartimeend + lonrad - ((rrange - 1.0) * math.pi)
            else:
                if rrange < maxx:
                    zhouranglestart, zhourangleend = horizon(
                        val, rrange + 2.0, lonrad, zsolartimestart, zsolartimeend, maxx
                    )
                else:
                    zhouranglestart = zsolartimestart + lonrad - (maxx + 1.0) * math.pi
                    zhourangleend = zsolartimeend + lonrad - (maxx + 1.0) * math.pi

            return (zhouranglestart, zhourangleend)

        # calculating the solar zenith angle

        PMU0 = 0
        if zcoshouranglesunset <= 1:
            zhouranglestart, zhourangleend = horizon(
                sha * math.pi / 180 + lonrad,
                2,
                lonrad,
                zsolartimestart,
                zsolartimeend,
                maxx,
            )

            if zcoshouranglesunset >= -1:
                zhouranglesunset = math.acos(zcoshouranglesunset)
                if (
                    zhourangleend <= -zhouranglesunset
                    or zhouranglestart >= zhouranglesunset
                ):
                    PMU0 = 0
                zhouranglestart = max(
                    -zhouranglesunset, min(zhouranglestart, zhouranglesunset)
                )
                zhourangleend = max(
                    -zhouranglesunset, min(zhourangleend, zhouranglesunset)
                )

            if (zhourangleend - zhouranglestart) > 1.0e-8:
                PMU0 = max(
                    0.0,
                    zsindecsinlat
                    + (
                        zcosdeccoslat
                        * (math.sin(zhourangleend) - math.sin(zhouranglestart))
                    )
                    / (zhourangleend - zhouranglestart),
                )
            else:
                PMU0 = 0.0

        return PMU0

    return np.vectorize(solar_zenith_angle_average)(
        lonrad,
        zcoshouranglesunset,
        zsindecsinlat,
        zcosdeccoslat,
        zsolartimestart,
        zsolartimeend,
        sha,
    )


def calculate_mean_radiant_temperature(ssrd, ssr, fdir, strd, strr, cossza):
    """
    mrt - Mean Radiant Temperature
    :param ssrd: is surface solar radiation downwards [J/m^-2]
    :param ssr: is surface net solar radiation [J/m^-2]
    :param fdir: is Total sky direct solar radiation at surface [J/m^-2]
    :param strd: is Surface thermal radiation downwards [J/m^-2]
    :param strr: is Surface net thermal radiation [J/m^-2]
    :param cossza: is cosine of solar zenith angle [degrees]

    returns Mean Radiant Temperature [K]
    """
    dsw = ssrd - fdir
    rsw = ssrd - ssr
    lur = strd - strr

    # calculate fp projected factor area

    gamma = np.arcsin(cossza) * 180 / np.pi
    fp = 0.308 * np.cos(to_radians * gamma * 0.998 - (gamma * gamma / 50000))

    # filter statement for solar zenith angle
    csza_filter1 = np.where((cossza > 0.01))
    # print(csza_filter1)
    fdir[csza_filter1] = fdir[csza_filter1] / cossza[csza_filter1]

    # calculate mean radiant temperature
    mrt = np.power(
        (
            (1 / 0.0000000567)
            * (
                0.5 * strd
                + 0.5 * lur
                + (0.7 / 0.97) * (0.5 * dsw + 0.5 * rsw + fp * fdir)
            )
        ),
        0.25,
    )

    return mrt


def calculate_utci(t2_k, va_ms, mrt_k, e_hPa):
    """
    UTCI
    :param t2m: (float array) is 2m temperature [K]
    :param va: (float array) is wind speed at 10 meters [m/s]
    :param mrt:(float array) is mean radiant temperature [K]
    :param ehPa: (float array) is water vapour pressure [hPa]

    Calculate UTCI with a 6th order polynomial approximation according to:
    Brode, P. et al. Deriving the operational procedure for the
    Universal Thermal Climate Index (UTCI). Int J Biometeorol (2012) 56: 48.1

    returns UTCI [°C]

    """
    t2 = __wrap(t2_k)
    va = __wrap(va_ms)
    mrt_kw = __wrap(mrt_k)
    ehPa = __wrap(e_hPa)

    rh = ehPa / 10.0  # rh in kPa

    t2m = kelvin_to_celcius(t2)  # polynomial approx. is in Celsius
    mrt = kelvin_to_celcius(mrt_kw)  # polynomial approx. is in Celsius

    e_mrt = np.subtract(mrt, t2m)

    t2m2 = t2m * t2m
    t2m3 = t2m2 * t2m
    t2m4 = t2m3 * t2m
    t2m5 = t2m4 * t2m
    t2m6 = t2m5 * t2m

    va2 = va * va
    va3 = va2 * va
    va4 = va3 * va
    va5 = va4 * va
    va6 = va5 * va

    e_mrt2 = e_mrt * e_mrt
    e_mrt3 = e_mrt2 * e_mrt
    e_mrt4 = e_mrt3 * e_mrt
    e_mrt5 = e_mrt4 * e_mrt
    e_mrt6 = e_mrt5 * e_mrt

    rh2 = rh * rh
    rh3 = rh2 * rh
    rh4 = rh3 * rh
    rh5 = rh4 * rh
    rh6 = rh5 * rh

    utci = (
        t2m
        + 6.07562052e-01
        + -2.27712343e-02 * t2m
        + 8.06470249e-04 * t2m2
        + -1.54271372e-04 * t2m3
        + -3.24651735e-06 * t2m4
        + 7.32602852e-08 * t2m5
        + 1.35959073e-09 * t2m6
        + -2.25836520e00 * va
        + 8.80326035e-02 * t2m * va
        + 2.16844454e-03 * t2m2 * va
        + -1.53347087e-05 * t2m3 * va
        + -5.72983704e-07 * t2m4 * va
        + -2.55090145e-09 * t2m5 * va
        + -7.51269505e-01 * va2
        + -4.08350271e-03 * t2m * va2
        + -5.21670675e-05 * t2m2 * va2
        + 1.94544667e-06 * t2m3 * va2
        + 1.14099531e-08 * t2m4 * va2
        + 1.58137256e-01 * va3
        + -6.57263143e-05 * t2m * va3
        + 2.22697524e-07 * t2m2 * va3
        + -4.16117031e-08 * t2m3 * va3
        + -1.27762753e-02 * va4
        + 9.66891875e-06 * t2m * va4
        + 2.52785852e-09 * t2m2 * va4
        + 4.56306672e-04 * va5
        + -1.74202546e-07 * t2m * va5
        + -5.91491269e-06 * va6
        + 3.98374029e-01 * e_mrt
        + 1.83945314e-04 * t2m * e_mrt
        + -1.73754510e-04 * t2m2 * e_mrt
        + -7.60781159e-07 * t2m3 * e_mrt
        + 3.77830287e-08 * t2m4 * e_mrt
        + 5.43079673e-10 * t2m5 * e_mrt
        + -2.00518269e-02 * va * e_mrt
        + 8.92859837e-04 * t2m * va * e_mrt
        + 3.45433048e-06 * t2m2 * va * e_mrt
        + -3.77925774e-07 * t2m3 * va * e_mrt
        + -1.69699377e-09 * t2m4 * va * e_mrt
        + 1.69992415e-04 * va2 * e_mrt
        + -4.99204314e-05 * t2m * va2 * e_mrt
        + 2.47417178e-07 * t2m2 * va2 * e_mrt
        + 1.07596466e-08 * t2m3 * va2 * e_mrt
        + 8.49242932e-05 * va3 * e_mrt
        + 1.35191328e-06 * t2m * va3 * e_mrt
        + -6.21531254e-09 * t2m2 * va3 * e_mrt
        + -4.99410301e-06 * va4 * e_mrt
        + -1.89489258e-08 * t2m * va4 * e_mrt
        + 8.15300114e-08 * va5 * e_mrt
        + 7.55043090e-04 * e_mrt2
        + -5.65095215e-05 * t2m * e_mrt2
        + -4.52166564e-07 * t2m * e_mrt2
        + 2.46688878e-08 * t2m3 * e_mrt2
        + 2.42674348e-10 * t2m4 * e_mrt2
        + 1.54547250e-04 * va * e_mrt2
        + 5.24110970e-06 * t2m * va * e_mrt2
        + -8.75874982e-08 * t2m2 * va * e_mrt2
        + -1.50743064e-09 * t2m3 * va * e_mrt2
        + -1.56236307e-05 * va2 * e_mrt2
        + -1.33895614e-07 * t2m * va2 * e_mrt2
        + 2.49709824e-09 * t2m2 * va2 * e_mrt2
        + 6.51711721e-07 * va3 * e_mrt2
        + 1.94960053e-09 * t2m * va3 * e_mrt2
        + -1.00361113e-08 * va4 * e_mrt2
        + -1.21206673e-05 * e_mrt3
        + -2.18203660e-07 * t2m * e_mrt3
        + 7.51269482e-09 * t2m2 * e_mrt3
        + 9.79063848e-11 * t2m3 * e_mrt3
        + 1.25006734e-06 * va * e_mrt3
        + -1.81584736e-09 * t2m * va * e_mrt3
        + -3.52197671e-10 * t2m2 * va * e_mrt3
        + -3.36514630e-08 * va2 * e_mrt3
        + 1.35908359e-10 * t2m * va2 * e_mrt3
        + 4.17032620e-10 * va3 * e_mrt3
        + -1.30369025e-09 * e_mrt4
        + 4.13908461e-10 * t2m * e_mrt4
        + 9.22652254e-12 * t2m2 * e_mrt4
        + -5.08220384e-09 * va * e_mrt4
        + -2.24730961e-11 * t2m * va * e_mrt4
        + 1.17139133e-10 * va2 * e_mrt4
        + 6.62154879e-10 * e_mrt5
        + 4.03863260e-13 * t2m * e_mrt5
        + 1.95087203e-12 * va * e_mrt5
        + -4.73602469e-12 * e_mrt6
        + 5.12733497e00 * rh
        + -3.12788561e-01 * t2m * rh
        + -1.96701861e-02 * t2m2 * rh
        + 9.99690870e-04 * t2m3 * rh
        + 9.51738512e-06 * t2m4 * rh
        + -4.66426341e-07 * t2m5 * rh
        + 5.48050612e-01 * va * rh
        + -3.30552823e-03 * t2m * va * rh
        + -1.64119440e-03 * t2m2 * va * rh
        + -5.16670694e-06 * t2m3 * va * rh
        + 9.52692432e-07 * t2m4 * va * rh
        + -4.29223622e-02 * va2 * rh
        + 5.00845667e-03 * t2m * va2 * rh
        + 1.00601257e-06 * t2m2 * va2 * rh
        + -1.81748644e-06 * t2m3 * va2 * rh
        + -1.25813502e-03 * va3 * rh
        + -1.79330391e-04 * t2m * va3 * rh
        + 2.34994441e-06 * t2m2 * va3 * rh
        + 1.29735808e-04 * va4 * rh
        + 1.29064870e-06 * t2m * va4 * rh
        + -2.28558686e-06 * va5 * rh
        + -3.69476348e-02 * e_mrt * rh
        + 1.62325322e-03 * t2m * e_mrt * rh
        + -3.14279680e-05 * t2m2 * e_mrt * rh
        + 2.59835559e-06 * t2m3 * e_mrt * rh
        + -4.77136523e-08 * t2m4 * e_mrt * rh
        + 8.64203390e-03 * va * e_mrt * rh
        + -6.87405181e-04 * t2m * va * e_mrt * rh
        + -9.13863872e-06 * t2m2 * va * e_mrt * rh
        + 5.15916806e-07 * t2m3 * va * e_mrt * rh
        + -3.59217476e-05 * va2 * e_mrt * rh
        + 3.28696511e-05 * t2m * va2 * e_mrt * rh
        + -7.10542454e-07 * t2m2 * va2 * e_mrt * rh
        + -1.24382300e-05 * va3 * e_mrt * rh
        + -7.38584400e-09 * t2m * va3 * e_mrt * rh
        + 2.20609296e-07 * va4 * e_mrt * rh
        + -7.32469180e-04 * e_mrt2 * rh
        + -1.87381964e-05 * t2m * e_mrt2 * rh
        + 4.80925239e-06 * t2m2 * e_mrt2 * rh
        + -8.75492040e-08 * t2m3 * e_mrt2 * rh
        + 2.77862930e-05 * va * e_mrt2 * rh
        + -5.06004592e-06 * t2m * va * e_mrt2 * rh
        + 1.14325367e-07 * t2m2 * va * e_mrt2 * rh
        + 2.53016723e-06 * va2 * e_mrt2 * rh
        + -1.72857035e-08 * t2m * va2 * e_mrt2 * rh
        + -3.95079398e-08 * va3 * e_mrt2 * rh
        + -3.59413173e-07 * e_mrt3 * rh
        + 7.04388046e-07 * t2m * e_mrt3 * rh
        + -1.89309167e-08 * t2m2 * e_mrt3 * rh
        + -4.79768731e-07 * va * e_mrt3 * rh
        + 7.96079978e-09 * t2m * va * e_mrt3 * rh
        + 1.62897058e-09 * va2 * e_mrt3 * rh
        + 3.94367674e-08 * e_mrt4 * rh
        + -1.18566247e-09 * t2m * e_mrt4 * rh
        + 3.34678041e-10 * va * e_mrt4 * rh
        + -1.15606447e-10 * e_mrt5 * rh
        + -2.80626406e00 * rh2
        + 5.48712484e-01 * t2m * rh2
        + -3.99428410e-03 * t2m2 * rh2
        + -9.54009191e-04 * t2m3 * rh2
        + 1.93090978e-05 * t2m4 * rh2
        + -3.08806365e-01 * va * rh2
        + 1.16952364e-02 * t2m * va * rh2
        + 4.95271903e-04 * t2m2 * va * rh2
        + -1.90710882e-05 * t2m3 * va * rh2
        + 2.10787756e-03 * va2 * rh2
        + -6.98445738e-04 * t2m * va2 * rh2
        + 2.30109073e-05 * t2m2 * va2 * rh2
        + 4.17856590e-04 * va3 * rh2
        + -1.27043871e-05 * t2m * va3 * rh2
        + -3.04620472e-06 * va4 * rh2
        + 5.14507424e-02 * e_mrt * rh2
        + -4.32510997e-03 * t2m * e_mrt * rh2
        + 8.99281156e-05 * t2m2 * e_mrt * rh2
        + -7.14663943e-07 * t2m3 * e_mrt * rh2
        + -2.66016305e-04 * va * e_mrt * rh2
        + 2.63789586e-04 * t2m * va * e_mrt * rh2
        + -7.01199003e-06 * t2m2 * va * e_mrt * rh2
        + -1.06823306e-04 * va2 * e_mrt * rh2
        + 3.61341136e-06 * t2m * va2 * e_mrt * rh2
        + 2.29748967e-07 * va3 * e_mrt * rh2
        + 3.04788893e-04 * e_mrt2 * rh2
        + -6.42070836e-05 * t2m * e_mrt2 * rh2
        + 1.16257971e-06 * t2m2 * e_mrt2 * rh2
        + 7.68023384e-06 * va * e_mrt2 * rh2
        + -5.47446896e-07 * t2m * va * e_mrt2 * rh2
        + -3.59937910e-08 * va2 * e_mrt2 * rh2
        + -4.36497725e-06 * e_mrt3 * rh2
        + 1.68737969e-07 * t2m * e_mrt3 * rh2
        + 2.67489271e-08 * va * e_mrt3 * rh2
        + 3.23926897e-09 * e_mrt4 * rh2
        + -3.53874123e-02 * rh3
        + -2.21201190e-01 * t2m * rh3
        + 1.55126038e-02 * t2m2 * rh3
        + -2.63917279e-04 * t2m3 * rh3
        + 4.53433455e-02 * va * rh3
        + -4.32943862e-03 * t2m * va * rh3
        + 1.45389826e-04 * t2m2 * va * rh3
        + 2.17508610e-04 * va2 * rh3
        + -6.66724702e-05 * t2m * va2 * rh3
        + 3.33217140e-05 * va3 * rh3
        + -2.26921615e-03 * e_mrt * rh3
        + 3.80261982e-04 * t2m * e_mrt * rh3
        + -5.45314314e-09 * t2m2 * e_mrt * rh3
        + -7.96355448e-04 * va * e_mrt * rh3
        + 2.53458034e-05 * t2m * va * e_mrt * rh3
        + -6.31223658e-06 * va2 * e_mrt * rh3
        + 3.02122035e-04 * e_mrt2 * rh3
        + -4.77403547e-06 * t2m * e_mrt2 * rh3
        + 1.73825715e-06 * va * e_mrt2 * rh3
        + -4.09087898e-07 * e_mrt3 * rh3
        + 6.14155345e-01 * rh4
        + -6.16755931e-02 * t2m * rh4
        + 1.33374846e-03 * t2m2 * rh4
        + 3.55375387e-03 * va * rh4
        + -5.13027851e-04 * t2m * va * rh4
        + 1.02449757e-04 * va2 * rh4
        + -1.48526421e-03 * e_mrt * rh4
        + -4.11469183e-05 * t2m * e_mrt * rh4
        + -6.80434415e-06 * va * e_mrt * rh4
        + -9.77675906e-06 * e_mrt2 * rh4
        + 8.82773108e-02 * rh5
        + -3.01859306e-03 * t2m * rh5
        + 1.04452989e-03 * va * rh5
        + 2.47090539e-04 * e_mrt * rh5
        + 1.48348065e-03 * rh6
    )

    # print(f"utci {utci}")

    utci_filtert2m = np.where(t2m >= 70)
    utci_filtert2m2 = np.where(t2m <= -70)
    utci_filterva = np.where(17 <= va)
    utci_filterva2 = np.where(0 >= va)
    utci_filterrh = np.where(5 < rh)
    utci_filtere_mrt = np.where(e_mrt >= 100.0)
    utci_filtere_mrt2 = np.where(e_mrt <= -30)

    utci[utci_filtert2m] = -9999
    # print(f"utci f1 {utci}")
    utci[utci_filtert2m2] = -9999
    # print(f"utci f2 {utci}")
    utci[utci_filterva] = -9999
    # print(f"utci f3 {utci}")
    utci[utci_filterva2] = -9999
    # print(f"utci f4 {utci}")
    utci[utci_filterrh] = -9999
    # print(f"utci f5 {utci}")
    utci[utci_filtere_mrt] = -9999
    # print(f"utci f6 {utci}")
    utci[utci_filtere_mrt2] = -9999
    # print(f"utci f7 {utci}")

    # print(f"utci {utci}")

    return utci


def calculate_wbgts(t2m):
    """
    wgbts - Wet Bulb Globe Temperature Simple
    :param t2m: 2m temperature [K]
    :param rh: relative humidity [pa]

    https://link.springer.com/article/10.1007/s00484-011-0453-2
    http://www.bom.gov.au/info/thermal_stress/#approximation
    https://www.jstage.jst.go.jp/article/indhealth/50/4/50_MS1352/_pdf

    returns Wet Bulb Globe Temperature [°C]
    """
    t2m = __wrap(t2m)
    rh = calculate_saturation_vapour_pressure(t2m)
    rh = kPa_to_hPa(rh)
    t2m = kelvin_to_celcius(t2m)
    wbgts = 0.567 * t2m + 0.393 * rh + 3.38
    return wbgts


def calculate_wbgt(t2m, mrt, va):
    """
    calculate wet bulb globe temperature
    :param t2m: 2m temperature [K]
    :param mrt: mean radiant temperature [K]
    :param va: wind speed at 10 meters [m/s]

    returns wet bulb globe temperature [°C]
    """
    t2m = __wrap(t2m)
    mrt = __wrap(mrt)
    va = __wrap(va)

    f = (1.1e8 * va ** 0.6) / (0.98 * 0.15 ** 0.4)
    a = f / 2
    b = -f * t2m - mrt ** 4
    rt1 = 3 ** (1 / 3)
    rt2 = np.sqrt(3) * np.sqrt(27 * a ** 4 - 16 * b ** 3) + 9 * a ** 2
    rt3 = 2 * 2 ** (2 / 3) * b
    a = a.clip(min=0)
    wbgt_quartic = -1 / 2 * np.sqrt(
        rt3 / (rt1 * rt2 ** (1 / 3)) + (2 ** (1 / 3) * rt2 ** (1 / 3)) / 3 ** (2 / 3)
    ) + 1 / 2 * np.sqrt(
        (4 * a)
        / np.sqrt(
            rt3 / (rt1 * rt2 ** (1 / 3))
            + (2 ** (1 / 3) * rt2 ** (1 / 3)) / 3 ** (2 / 3)
        )
        - (2 ** (1 / 3) * rt2 ** (1 / 3)) / 3 ** (2 / 3)
        - rt3 / (rt1 * rt2 ** (1 / 3))
    )
    wbgt_quartic = kelvin_to_celcius(wbgt_quartic)
    return wbgt_quartic


def calculate_mrt_from_wbgt(t2m, wbgt, va):
    """
    calculate mean radiant temperature from wet bulb globe temperature
    :param t2m: 2m temperature [K]
    :param wbgt: wet bulb globe temperature in kelvin [K]
    :param va: wind speed at 10 meters [m/s]
    returns mean radiant temperature [K]
    """
    t2m = __wrap(t2m)
    wbgt = __wrap(wbgt)
    va = __wrap(va)

    f = (1.1e8 * va ** 0.6) / (0.98 * 0.15 ** 0.4)
    wbgt4 = wbgt ** 4
    mrtc = wbgt4 + f * (wbgt - t2m)
    mrtc2 = np.sqrt(np.sqrt(mrtc))
    return kelvin_to_celcius(mrtc2)


def calculate_humidex(t2m, td):
    """
    humidex - heat index used by the Canadian Meteorological Service
    :param t2m: 2m temperature [K]
    :param td: dew point temperature [K]

    returns humidex [°C]
    """
    t2m = __wrap(t2m)
    td = __wrap(td)
    e = 6.11 * np.exp(5417.7530 * ((1 / t2m) - (1 / td)))
    h = 0.5555 * (e - 10.0)
    humidex = (t2m + h) - 273.15
    return humidex


def calculate_net_effective_temperature(t2m, va, td):
    """
    Net Effective Temperature used in Hong Kong, Poland and Germany
    :param t2m: 2m temperature [K]
    :param td: 2m dew point temperature [K]
    :param rh: Relative Humidity [pa]
    :param va: Wind speed at 10 meters [m/s]

    returns net effective temperature [°C]
    """
    rh = calculate_relative_humidity_percent(t2m, td)
    t2m = __wrap(t2m)
    va = __wrap(va)
    rh = __wrap(rh)
    t2m = kelvin_to_celcius(t2m)
    rh = kPa_to_hPa(rh)
    ditermeq = 1 / 1.76 + 1.4 * va ** 0.75
    net = 37 - (37 - t2m / 0.68 - 0.0014 * rh + ditermeq) - 0.29 * t2m * (1 - 0.01 * rh)
    return net


def calculate_apparent_temperature(t2m, va, rh=None):
    """
    Apparent Temperature version without radiation
    :param t2m: 2m Temperature [K]
    :param rh: Relative Humidity [pa]

    returns apparent temperature [K]
    """
    if rh is None:
        rh = calculate_saturation_vapour_pressure(t2m)
    t2m = __wrap(t2m)
    va = __wrap(va)
    rh = __wrap(rh)
    va = va * 4.87 / np.log10(67.8 * 10 - 5.42)  # converting to 2m, ~1.2m wind speed
    at = t2m + 0.33 * rh - 0.7 * va - 4
    at = kelvin_to_celcius(at)
    return at


def calculate_wind_chill(t2m, va):
    """
    Wind Chill
    :param t2m: 2m Temperature [K]
    :param va: wind speed at 10 meters [m/s]

    returns wind chill [°C]
    """
    t2m = __wrap(t2m)
    va = __wrap(va)
    t2m = kelvin_to_celcius(t2m)
    va = va * 2.23694  # convert to miles per hour
    windchill = 13.12 + 0.6215 * t2m - 11.37 * va ** 0.16 + 0.3965 + t2m + va ** 0.16
    return windchill


def calculate_heat_index_simplified(t2m, rh=None):
    """
    Heat Index
       :param t2m: np.array 2m temperature [K]
       :param rh: Relative Humidity [pa]
       returns heat index [°C]
    """
    t2m = __wrap(t2m)
    if rh is None:
        rh = calculate_saturation_vapour_pressure(t2m)
    t2m = kelvin_to_celcius(t2m)
    rh = kPa_to_hPa(rh)

    hiarray = [
        -8.784695,
        1.61139411,
        2.338549,
        0.14611605,
        1.2308094e-2,
        2.211732e-3,
        7.2546e-4,
        3.58e-6,
    ]

    hi = (
        -hiarray[0]
        + hiarray[1] * t2m
        + hiarray[2] * rh
        - hiarray[3] * t2m * rh
        - hiarray[4] * rh ** 2
        + hiarray[5] * t2m ** 2 * rh
        + hiarray[6] * t2m * rh ** 2
        - hiarray[7] * t2m ** 2 * rh ** 2
    )

    return hi


def calculate_heat_index_adjusted(t2m, td):
    """
    Heat Index adjusted
       :param t2m: np.array 2m temperature [K]
       :param td: np.array 2m dewpoint temperature  [K]
       returns heat index [°C]
    """
    t2m = __wrap(t2m)
    td = __wrap(td)

    rh = calculate_relative_humidity_percent(t2m, td)
    t2m = kelvin_to_farenheit(t2m)

    hiarray = [
        42.379,
        2.04901523,
        10.1433312,
        0.22475541,
        0.00683783,
        0.05481717,
        0.00122874,
        0.00085282,
        0.00000199,
    ]

    hi = (
        -hiarray[0]
        + hiarray[1] * t2m
        + hiarray[2] * rh
        - hiarray[3] * t2m * rh
        - hiarray[4] * t2m ** 2
        - hiarray[5] * rh ** 2
        + hiarray[6] * t2m ** 2 * rh
        + hiarray[7] * t2m * rh ** 2
        - hiarray[8] * t2m ** 2 * rh ** 2
    )

    hi_filter1 = np.where(t2m > 80)
    hi_filter2 = np.where(t2m < 112)
    hi_filter3 = np.where(rh <= 13)
    hi_filter4 = np.where(t2m < 87)
    hi_filter5 = np.where(rh > 85)
    hi_filter6 = np.where(t2m < 80)

    adjustment1 = (
        (13 - rh[hi_filter1 and hi_filter2 and hi_filter3])
        / 4
        * np.sqrt(
            17 - np.abs(t2m[hi_filter1 and hi_filter2 and hi_filter3] - 0.95) / 17
        )
    )
    adjustment2 = (rh[hi_filter1 and hi_filter4 and hi_filter5] - 85) * (
        (87 - t2m[hi_filter1 and hi_filter4 and hi_filter5]) / 5
    )
    adjustment3 = 0.5 * (
        t2m[hi_filter6]
        + 61.0
        + ((t2m[hi_filter6] - 68.0) * 1.2)
        + (rh[hi_filter6] * 0.094)
    )

    hi[hi_filter1 and hi_filter2 and hi_filter3] = (
        hi[hi_filter1 and hi_filter2 and hi_filter3] - adjustment1
    )
    hi[hi_filter1 and hi_filter4 and hi_filter5] = (
        hi[hi_filter1 and hi_filter4 and hi_filter5] - adjustment2
    )
    hi[hi_filter6] = adjustment3
    hi = farenheit_to_celcius(hi)
    return hi


# Helpers

# convert celcius to kelvin
def celcius_to_kelvin(tc):
    tk = tc + 273.15
    return tk


# convert kelvin to celcius
def kelvin_to_celcius(tk):
    tc = tk - 273.15
    return tc

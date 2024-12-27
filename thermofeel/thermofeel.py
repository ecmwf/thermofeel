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
    * Apparent Temperature
    * Heat Index Adjusted
    * Heat Index Simplified
    * Humidex
    * Normal Effective Temperature
    * Wet Bulb Globe Temperature
    * Wet Bulb Globe Temperature Simple
    * Wind Chill

    In support of the above indexes, it also calculates:
    * Globe Temperature
    * Mean Radiant Temperature
    * Mean Radiant Temperature from Globe Temperature
    * Relative Humidity Percentage
    * Saturation vapour pressure
    * Wet Bulb Temperature

    To calculate the cos of the solar zenith angle, we suggest to use the
    earthkit-meteo library (github.com:ecmwf/earthkit-meteo.git)
  """

import math

import numpy as np

from .helpers import (
    celsius_to_kelvin,
    fahrenheit_to_kelvin,
    kelvin_to_celsius,
    kelvin_to_fahrenheit,
)

to_radians = math.pi / 180.0


def calculate_relative_humidity_percent(t2_k, td_k):
    """
    Relative Humidity in percent
        :param t2_k: (float array) 2m temperature [K]
        :param td_k: (float array) dew point temperature [K]
        returns relative humidity [%]
    """

    t2_c = kelvin_to_celsius(t2_k)
    td_c = kelvin_to_celsius(td_k)
    # saturated vapour pressure
    es = 6.11 * 10.0 ** (7.5 * t2_c / (237.3 + t2_c))
    # vapour pressure
    e = 6.11 * 10.0 ** (7.5 * td_c / (237.3 + td_c))
    rh = (e / es) * 100
    return rh


def calculate_saturation_vapour_pressure(t2_k):
    """
    Saturation vapour pressure over water
        :param t2_k: (float array) 2m temperature [K]
        returns saturation vapor pressure over water in the pure phase [hPa] == [mBar]
    Reference: Hardy (1998)
    https://www.decatur.de/javascript/dew/resources/its90formulas.pdf
    """

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
    ess = g[7] * np.log(t2_k)
    for i in range(7):
        ess = ess + g[i] * np.power(t2_k, (i - 2))

    ess = np.exp(ess) * 0.01  # hPa

    return ess


def calculate_saturation_vapour_pressure_multiphase(t2_k, phase):
    """
    Saturation vapour pressure over liquid water and ice
        :param t2_k: (float array) 2m temperature [K]
        :param phase: 0 over liquid water and 1 over ice
        returns pressure of water vapor over a surface of liquid water or ice [hPa] == [mBar]
    Reference: ECMWF IFS Documentation CY45R1 - Part IV : Physical processes (2018) pp. 116
    https://doi.org/10.21957/4whwo8jw0
    https://metview.readthedocs.io/en/latest/api/functions/saturation_vapour_pressure.html
    """
    T0 = 273.16  # triple point of water 273.16 K (0.01 °C) at 611.73 Pa
    es = np.zeros_like(t2_k)
    y = (t2_k - T0) / (t2_k - 32.19)  # over liquid water
    es[phase == 0] = 6.1121 * np.exp(17.502 * y[phase == 0])
    y = (t2_k - T0) / (t2_k + 0.7)  # over ice
    es[phase == 1] = 6.1121 * np.exp(22.587 * y[phase == 1])

    return es


def calculate_nonsaturation_vapour_pressure(t2_k, rh):
    """
    Non saturated vapour pressure
        :param t2_k: (float array) 2m temperature [K]
        :param rh: (float array) relative humidity percentage [%]
        returns non saturated vapor pressure [hPa] == [mBar]
    Reference: Bureau of Meteorology (2010)
    http://www.bom.gov.au/info/thermal_stress/#approximation
    """

    t2_c = kelvin_to_celsius(t2_k)
    ens = rh / 100 * 6.105 * np.exp(17.27 * t2_c / (237.7 + t2_c))

    return ens


def scale_windspeed(va, h):
    """
    Scaling wind speed from 10 metres to height h
        :param va: (float array) 10m wind speed [m/s]
        :param h: (float array) height at which wind speed needs to be scaled [m]
        returns wind speed at height h
    Reference: Bröde et al. (2012)
    https://doi.org/10.1007/s00484-011-0454-1
    """
    c = 1 / np.log10(10 / 0.01)  #
    c = 0.333333333333
    vh = va * np.log10(h / 0.01) * c

    return vh


def approximate_dsrp(fdir, cossza, threshold=0.1):
    """
    Helper function to approximate dsrp from fdir and cossza
    Note that the function introduces large errors as cossza approaches zero.
    Only use if dsrp is not available in your dataset.
        :param fdir: (float array) total sky direct solar radiation at surface [W m-2]
        :param cossza: (float array) cosine of solar zenith angle [dimentionless]
        returns direct radiation from the Sun [W m-2]
    """
    # filter statement for solar zenith angle to avoid division by zero.
    csza_filter1 = np.where((cossza > threshold))
    dsrp = np.copy(fdir)  # dsrp = fdir for cossza <= 0.01, equals to fdir
    dsrp[csza_filter1] = dsrp[csza_filter1] / cossza[csza_filter1]
    return dsrp


def calculate_dew_point_from_relative_humidity(rh, t2_k):
    """
    Dew point temperature at 2m from relative humidity in percent
        :param rh: (float array) relative humidity [%]
        :param t2_k: (float array) 2m temperature [K]
        returns dew point temperature [K]
    Reference: Alduchov and Eskridge (1996)
    https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2
    """
    t2_c = kelvin_to_celsius(t2_k)
    td_c = (
        243.04
        * (np.log(rh / 100) + ((17.625 * t2_c) / (243.04 + t2_c)))
        / (17.625 - np.log(rh / 100) - ((17.625 * t2_c) / (243.04 + t2_c)))
    )
    td_k = celsius_to_kelvin(td_c)
    return td_k


def calculate_mean_radiant_temperature(ssrd, ssr, dsrp, strd, fdir, strr, cossza):
    """
    MRT - Mean Radiant Temperature
        :param ssrd: (float array) surface solar radiation downwards [W m-2]
        :param ssr: (float array) surface net solar radiation [W m-2]
        :param dsrp: (float array) direct solar radiation [W m-2]
        :param strd: (float array) surface thermal radiation downwards [W m-2]
        :param fdir: (float array) total sky direct solar radiation at surface [W m-2]
        :param strr: (float array) surface net thermal radiation [W m-2]
        :param cossza: (float array) cosine of solar zenith angle [dimentionless]
        returns mean radiant temperature [K]
    Reference: Di Napoli et al. (2020)
    https://link.springer.com/article/10.1007/s00484-020-01900-5
    """

    dsw = ssrd - fdir
    rsw = ssrd - ssr
    lur = strd - strr
    # Istar = dsrp

    # calculate fp projected factor area

    gamma = np.arcsin(cossza) * 180 / np.pi
    fp = 0.308 * np.cos(to_radians * gamma * (0.998 - gamma * gamma / 50000))

    # calculate mean radiant temperature
    mrt = np.power(
        (
            (1 / 0.0000000567)
            * (
                0.5 * strd
                + 0.5 * lur
                + (0.7 / 0.97) * (0.5 * dsw + 0.5 * rsw + fp * dsrp)
            )
        ),
        0.25,
    )

    return mrt


def calculate_utci_polynomial(t2m, mrt, va, wvp):
    """
    Helper function to calculate the UTCI polynomial approximation
        :param t2_k: (float array) is 2m temperature [K]
        :param mrt: (float array) is mean radiant temperature [K]
        :param va: (float array) is wind speed at 10 meters [m/s]
        :param wvp: (float array) is water vapour pressure [kPa]
    returns UTCI [K]
    Reference: Brode et al. (2012)
    https://doi.org/10.1007/s00484-011-0454-1
    """
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

    wvp2 = wvp * wvp
    wvp3 = wvp2 * wvp
    wvp4 = wvp3 * wvp
    wvp5 = wvp4 * wvp
    wvp6 = wvp5 * wvp

    varh2 = va * wvp2
    va2_rh = va2 * wvp
    va2_e_mrt = va2 * e_mrt
    e_mrt_rh = e_mrt * wvp
    e_mrt_rh2 = e_mrt * wvp2
    e_mrt2_rh = e_mrt2 * wvp
    e_mrt2_rh2 = e_mrt2 * wvp2
    e_mrt_rh3 = e_mrt * wvp3
    va_e_mrt = va * e_mrt
    va_e_mrt2 = va * e_mrt2
    va_rh = va * wvp
    t2m_va = t2m * va
    e_mrt3_rh = e_mrt3 * wvp
    e_mrt4_rh = e_mrt4 * wvp

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
        + -2.00518269e-02 * va_e_mrt
        + 8.92859837e-04 * t2m * va_e_mrt
        + 3.45433048e-06 * t2m2 * va_e_mrt
        + -3.77925774e-07 * t2m3 * va_e_mrt
        + -1.69699377e-09 * t2m4 * va_e_mrt
        + 1.69992415e-04 * va2_e_mrt
        + -4.99204314e-05 * t2m * va2_e_mrt
        + 2.47417178e-07 * t2m2 * va2_e_mrt
        + 1.07596466e-08 * t2m3 * va2_e_mrt
        + 8.49242932e-05 * va3 * e_mrt
        + 1.35191328e-06 * t2m * va3 * e_mrt
        + -6.21531254e-09 * t2m2 * va3 * e_mrt
        + -4.99410301e-06 * va4 * e_mrt
        + -1.89489258e-08 * t2m * va4 * e_mrt
        + 8.15300114e-08 * va5 * e_mrt
        + 7.55043090e-04 * e_mrt2
        + -5.65095215e-05 * t2m * e_mrt2
        + -4.52166564e-07 * t2m2 * e_mrt2
        + 2.46688878e-08 * t2m3 * e_mrt2
        + 2.42674348e-10 * t2m4 * e_mrt2
        + 1.54547250e-04 * va_e_mrt2
        + 5.24110970e-06 * t2m * va_e_mrt2
        + -8.75874982e-08 * t2m2 * va_e_mrt2
        + -1.50743064e-09 * t2m3 * va_e_mrt2
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
        + -1.81584736e-09 * t2m_va * e_mrt3
        + -3.52197671e-10 * t2m2 * va * e_mrt3
        + -3.36514630e-08 * va2 * e_mrt3
        + 1.35908359e-10 * t2m * va2 * e_mrt3
        + 4.17032620e-10 * va3 * e_mrt3
        + -1.30369025e-09 * e_mrt4
        + 4.13908461e-10 * t2m * e_mrt4
        + 9.22652254e-12 * t2m2 * e_mrt4
        + -5.08220384e-09 * va * e_mrt4
        + -2.24730961e-11 * t2m_va * e_mrt4
        + 1.17139133e-10 * va2 * e_mrt4
        + 6.62154879e-10 * e_mrt5
        + 4.03863260e-13 * t2m * e_mrt5
        + 1.95087203e-12 * va * e_mrt5
        + -4.73602469e-12 * e_mrt6
        + 5.12733497e00 * wvp
        + -3.12788561e-01 * t2m * wvp
        + -1.96701861e-02 * t2m2 * wvp
        + 9.99690870e-04 * t2m3 * wvp
        + 9.51738512e-06 * t2m4 * wvp
        + -4.66426341e-07 * t2m5 * wvp
        + 5.48050612e-01 * va_rh
        + -3.30552823e-03 * t2m * va_rh
        + -1.64119440e-03 * t2m2 * va_rh
        + -5.16670694e-06 * t2m3 * va_rh
        + 9.52692432e-07 * t2m4 * va_rh
        + -4.29223622e-02 * va2_rh
        + 5.00845667e-03 * t2m * va2_rh
        + 1.00601257e-06 * t2m2 * va2_rh
        + -1.81748644e-06 * t2m3 * va2_rh
        + -1.25813502e-03 * va3 * wvp
        + -1.79330391e-04 * t2m * va3 * wvp
        + 2.34994441e-06 * t2m2 * va3 * wvp
        + 1.29735808e-04 * va4 * wvp
        + 1.29064870e-06 * t2m * va4 * wvp
        + -2.28558686e-06 * va5 * wvp
        + -3.69476348e-02 * e_mrt_rh
        + 1.62325322e-03 * t2m * e_mrt_rh
        + -3.14279680e-05 * t2m2 * e_mrt_rh
        + 2.59835559e-06 * t2m3 * e_mrt_rh
        + -4.77136523e-08 * t2m4 * e_mrt_rh
        + 8.64203390e-03 * va * e_mrt_rh
        + -6.87405181e-04 * t2m_va * e_mrt_rh
        + -9.13863872e-06 * t2m2 * va * e_mrt_rh
        + 5.15916806e-07 * t2m3 * va * e_mrt_rh
        + -3.59217476e-05 * va2 * e_mrt_rh
        + 3.28696511e-05 * t2m * va2 * e_mrt_rh
        + -7.10542454e-07 * t2m2 * va2 * e_mrt_rh
        + -1.24382300e-05 * va3 * e_mrt_rh
        + -7.38584400e-09 * t2m * va3 * e_mrt_rh
        + 2.20609296e-07 * va4 * e_mrt_rh
        + -7.32469180e-04 * e_mrt2_rh
        + -1.87381964e-05 * t2m * e_mrt2_rh
        + 4.80925239e-06 * t2m2 * e_mrt2_rh
        + -8.75492040e-08 * t2m3 * e_mrt2_rh
        + 2.77862930e-05 * va * e_mrt2_rh
        + -5.06004592e-06 * t2m_va * e_mrt2_rh
        + 1.14325367e-07 * t2m2 * va * e_mrt2_rh
        + 2.53016723e-06 * va2 * e_mrt2_rh
        + -1.72857035e-08 * t2m * va2 * e_mrt2_rh
        + -3.95079398e-08 * va3 * e_mrt2_rh
        + -3.59413173e-07 * e_mrt3_rh
        + 7.04388046e-07 * t2m * e_mrt3_rh
        + -1.89309167e-08 * t2m2 * e_mrt3_rh
        + -4.79768731e-07 * va * e_mrt3_rh
        + 7.96079978e-09 * t2m_va * e_mrt3_rh
        + 1.62897058e-09 * va2 * e_mrt3_rh
        + 3.94367674e-08 * e_mrt4_rh
        + -1.18566247e-09 * t2m * e_mrt4_rh
        + 3.34678041e-10 * va * e_mrt4_rh
        + -1.15606447e-10 * e_mrt5 * wvp
        + -2.80626406e00 * wvp2
        + 5.48712484e-01 * t2m * wvp2
        + -3.99428410e-03 * t2m2 * wvp2
        + -9.54009191e-04 * t2m3 * wvp2
        + 1.93090978e-05 * t2m4 * wvp2
        + -3.08806365e-01 * varh2
        + 1.16952364e-02 * t2m * varh2
        + 4.95271903e-04 * t2m2 * varh2
        + -1.90710882e-05 * t2m3 * varh2
        + 2.10787756e-03 * va2 * wvp2
        + -6.98445738e-04 * t2m * va2 * wvp2
        + 2.30109073e-05 * t2m2 * va2 * wvp2
        + 4.17856590e-04 * va3 * wvp2
        + -1.27043871e-05 * t2m * va3 * wvp2
        + -3.04620472e-06 * va4 * wvp2
        + 5.14507424e-02 * e_mrt_rh2
        + -4.32510997e-03 * t2m * e_mrt_rh2
        + 8.99281156e-05 * t2m2 * e_mrt_rh2
        + -7.14663943e-07 * t2m3 * e_mrt_rh2
        + -2.66016305e-04 * va * e_mrt_rh2
        + 2.63789586e-04 * t2m_va * e_mrt_rh2
        + -7.01199003e-06 * t2m2 * va * e_mrt_rh2
        + -1.06823306e-04 * va2 * e_mrt_rh2
        + 3.61341136e-06 * t2m * va2 * e_mrt_rh2
        + 2.29748967e-07 * va3 * e_mrt_rh2
        + 3.04788893e-04 * e_mrt2_rh2
        + -6.42070836e-05 * t2m * e_mrt2_rh2
        + 1.16257971e-06 * t2m2 * e_mrt2_rh2
        + 7.68023384e-06 * va * e_mrt2_rh2
        + -5.47446896e-07 * t2m_va * e_mrt2_rh2
        + -3.59937910e-08 * va2 * e_mrt2_rh2
        + -4.36497725e-06 * e_mrt3 * wvp2
        + 1.68737969e-07 * t2m * e_mrt3 * wvp2
        + 2.67489271e-08 * va * e_mrt3 * wvp2
        + 3.23926897e-09 * e_mrt4 * wvp2
        + -3.53874123e-02 * wvp3
        + -2.21201190e-01 * t2m * wvp3
        + 1.55126038e-02 * t2m2 * wvp3
        + -2.63917279e-04 * t2m3 * wvp3
        + 4.53433455e-02 * va * wvp3
        + -4.32943862e-03 * t2m_va * wvp3
        + 1.45389826e-04 * t2m2 * va * wvp3
        + 2.17508610e-04 * va2 * wvp3
        + -6.66724702e-05 * t2m * va2 * wvp3
        + 3.33217140e-05 * va3 * wvp3
        + -2.26921615e-03 * e_mrt_rh3
        + 3.80261982e-04 * t2m * e_mrt_rh3
        + -5.45314314e-09 * t2m2 * e_mrt_rh3
        + -7.96355448e-04 * va * e_mrt_rh3
        + 2.53458034e-05 * t2m_va * e_mrt_rh3
        + -6.31223658e-06 * va2 * e_mrt_rh3
        + 3.02122035e-04 * e_mrt2 * wvp3
        + -4.77403547e-06 * t2m * e_mrt2 * wvp3
        + 1.73825715e-06 * va * e_mrt2 * wvp3
        + -4.09087898e-07 * e_mrt3 * wvp3
        + 6.14155345e-01 * wvp4
        + -6.16755931e-02 * t2m * wvp4
        + 1.33374846e-03 * t2m2 * wvp4
        + 3.55375387e-03 * va * wvp4
        + -5.13027851e-04 * t2m_va * wvp4
        + 1.02449757e-04 * va2 * wvp4
        + -1.48526421e-03 * e_mrt * wvp4
        + -4.11469183e-05 * t2m * e_mrt * wvp4
        + -6.80434415e-06 * va * e_mrt * wvp4
        + -9.77675906e-06 * e_mrt2 * wvp4
        + 8.82773108e-02 * wvp5
        + -3.01859306e-03 * t2m * wvp5
        + 1.04452989e-03 * va * wvp5
        + 2.47090539e-04 * e_mrt * wvp5
        + 1.48348065e-03 * wvp6
    )

    return utci


def calculate_utci(t2_k, va, mrt, td_k=None, ehPa=None):
    """
    UTCI - Universal Thermal Climate Index
        :param t2_k: (float array) is 2m temperature [K]
        :param va: (float array) is wind speed at 10 meters [m/s]
        :param mrt: (float array) is mean radiant temperature [K]
        :param td_k: (float array) is 2m dew point temperature [K]
        :param ehPa: (float array) is water vapour pressure [hPa]
    returns UTCI [K]
    Reference: Brode et al. (2012)
    https://doi.org/10.1007/s00484-011-0454-1
    """

    if ehPa is not None:
        wvp = ehPa / 10.0  # water vapour pressure in kPa
    else:
        if td_k is not None:
            rh_pc = calculate_relative_humidity_percent(t2_k, td_k)
            ehPa = calculate_saturation_vapour_pressure(t2_k) * rh_pc / 100.0
            wvp = ehPa / 10.0  # water vapour pressure in kPa
        else:
            raise ValueError("Missing input ehPa or td_k")

    t2_c = kelvin_to_celsius(t2_k)  # polynomial approx. is in Celsius
    mrt_c = kelvin_to_celsius(mrt)  # polynomial approx. is in Celsius

    utci = calculate_utci_polynomial(t2_c, mrt_c, va, wvp)
    utci_k = celsius_to_kelvin(utci)

    return utci_k


def calculate_wbgt_simple(t2_k, rh):
    """
    WBGT - Wet Bulb Globe Temperature computed by a the simpler algorithm
        :param t2_k: (float array) 2m temperature [K]
        :param rh: (float array) relative humidity percentage [%]
        returns Wet Bulb Globe Temperature [K]
    Reference: ACSM (1984)
    https://doi.org/10.1080/00913847.1984.11701899
    See also: http://www.bom.gov.au/info/thermal_stress/#approximation
    https://www.jstage.jst.go.jp/article/indhealth/50/4/50_MS1352/_pdf
    """
    t2_c = kelvin_to_celsius(t2_k)
    e = calculate_nonsaturation_vapour_pressure(t2_k, rh)
    wbgt = 0.567 * t2_c + 0.393 * e + 3.94
    wbgt_k = celsius_to_kelvin(wbgt)

    return wbgt_k


def calculate_wbt(t2_k, rh):
    """
    Wet Bulb Temperature
        :param t2_k: (float array) 2m temperature [K]
        :param rh: (float array) relative humidity percentage [%]
        returns wet bulb temperature [K]
    Reference: Stull (2011)
    https://doi.org/10.1175/JAMC-D-11-0143.1
    """
    t2_c = kelvin_to_celsius(t2_k)
    tw = (
        t2_c * np.arctan(0.151977 * np.sqrt(rh + 8.313659))
        + np.arctan(t2_c + rh)
        - np.arctan(rh - 1.676331)
        + 0.00391838 * (rh) ** (3 / 2) * np.arctan(0.023101 * rh)
        - 4.686035
    )
    tw_k = celsius_to_kelvin(tw)

    return tw_k


def calculate_bgt(t2_k, mrt, va):
    """
    Globe temperature
        :param t2_k: (float array) 2m temperature [K]
        :param mrt: (float array) mean radiant temperature [K]
        :param va: (float array) wind speed at 10 meters [m/s]
        returns globe temperature [K]
    Reference: Guo et al. 2018
    https://doi.org/10.1016/j.enbuild.2018.08.029
    """
    v = scale_windspeed(
        va, 1.1
    )  # formula requires wind speed at 1.1m (i.e., at the level of the globe)

    # a = 1
    d = (1.1e8 * v**0.6) / (0.95 * 0.15**0.4)
    e = -(mrt**4) - d * t2_k

    q = 12 * e
    s = 27 * (d**2)
    delta = ((s + np.sqrt(s**2 - 4 * (q**3))) / 2) ** (1 / 3)
    Q = 0.5 * np.sqrt((1 / 3) * (delta + q / delta))

    bgt = -Q + 0.5 * np.sqrt(-4 * (Q**2) + d / Q)

    # f = (1.1e8 * va**0.6) / (0.95 * 0.15**0.4)
    # a = f / 2
    # b = -f * t2_k - mrt**4
    # rt1 = 3 ** (1 / 3)
    # rt2 = np.sqrt(3) * np.sqrt(27 * a**4 - 16 * b**3) + 9 * a**2
    # rt3 = 2 * 2 ** (2 / 3) * b
    # a = a.clip(min=0)
    # bgt = -1 / 2 * np.sqrt(
    #     rt3 / (rt1 * rt2 ** (1 / 3)) + (2 ** (1 / 3) * rt2 ** (1 / 3)) / 3 ** (2 / 3)
    # ) + 1 / 2 * np.sqrt(
    #     (4 * a)
    #     / np.sqrt(
    #         rt3 / (rt1 * rt2 ** (1 / 3))
    #         + (2 ** (1 / 3) * rt2 ** (1 / 3)) / 3 ** (2 / 3)
    #     )
    #     - (2 ** (1 / 3) * rt2 ** (1 / 3)) / 3 ** (2 / 3)
    #     - rt3 / (rt1 * rt2 ** (1 / 3))
    # )

    return bgt


def calculate_wbgt(t2_k, mrt, va, td_k):
    """
    WBGT - Wet Bulb Globe Temperature
        :param t2_k: (float array) 2m temperature [K]
        :param mrt: (float array) mean radiant temperature [K]
        :param va: (float array) wind speed at 10 meters [m/s]
        :param td_k: (float array) dew point temperature [K]
        returns wet bulb globe temperature [K]
    Reference: Stull (2011)
    https://doi.org/10.1175/JAMC-D-11-0143.1
    See also: http://www.bom.gov.au/info/thermal_stress/
    """

    bgt_k = calculate_bgt(t2_k, mrt, va)
    bgt_c = kelvin_to_celsius(bgt_k)

    rh = calculate_relative_humidity_percent(t2_k, td_k)
    t2_c = kelvin_to_celsius(t2_k)
    tw_k = calculate_wbt(t2_k, rh)
    tw_c = kelvin_to_celsius(tw_k)

    wbgt = 0.7 * tw_c + 0.2 * bgt_c + 0.1 * t2_c
    wbgt_k = celsius_to_kelvin(wbgt)

    return wbgt_k


def calculate_mrt_from_bgt(t2_k, bgt_k, va):
    """
    Mean radiant temperature from globe temperature
        :param t2_k: (float array) 2m temperature [K]
        :param bgt_k: (float array) globe temperature [K]
        :param va: (float array) wind speed at 10 meters [m/s]
        returns mean radiant temperature [K]
    Reference: Brimicombe et al. (2023)
    https://doi.org/10.1029/2022GH000701
    """
    v = scale_windspeed(
        va, 1.1
    )  # formula requires wind speed at 1.1m (i.e., at the level of the globe)
    f = (1.1e8 * v**0.6) / (0.95 * 0.15**0.4)
    bgt4 = bgt_k**4
    mrtc = bgt4 + f * (bgt_k - t2_k)
    mrtc2 = np.sqrt(np.sqrt(mrtc))

    return mrtc2


def calculate_humidex(t2_k, td_k):
    """
        Humidex
        :param t2_k: (float array) 2m temperature [K]
        :param td_k: (float array) dew point temperature [K]
        returns humidex [K]
    Reference: Blazejczyk et al. (2012)
    https://doi.org/10.1007/s00484-011-0453-2
    """
    vp = 6.11 * np.exp(5417.7530 * ((1 / 273.16) - (1 / td_k)))  # vapour pressure [hPa]
    h = 0.5555 * (vp - 10.0)
    humidex = t2_k + h

    return humidex


def calculate_normal_effective_temperature(t2_k, va, rh):
    """
    NET - Normal Effective Temperature
        :param t2_k: (float array) 2m temperature [K]
        :param va: (float array) wind speed at 10 meters [m/s]
        :param rh: (float array) relative humidity percentage [%]
        returns normal effective temperature [K]
    Reference: Li and Chan (2006)
    https://doi.org/10.1017/S1350482700001602
    """
    t2_k = kelvin_to_celsius(t2_k)
    v = scale_windspeed(va, 1.2)  # formula requires wind speed at 1.2m
    ditermeq = 1 / (1.76 + 1.4 * v**0.75)
    net = (
        37
        - ((37 - t2_k) / (0.68 - 0.0014 * rh + ditermeq))
        - 0.29 * t2_k * (1 - 0.01 * rh)
    )
    net_k = celsius_to_kelvin(net)

    return net_k


def calculate_apparent_temperature(t2_k, va, rh):
    """
    Apparent Temperature - version without radiation
        :param t2_k: (float array) 2m temperature [K]
        :param va: (float array) wind speed at 10 meters [m/s]
        :param rh: (float array) relative humidity percentage [%]
        returns apparent temperature [K]
    Reference: Steadman (1984)
    https://doi.org/10.1175/1520-0450(1984)023%3C1674:AUSOAT%3E2.0.CO;2
    See also: http://www.bom.gov.au/info/thermal_stress/#atapproximation
    """
    t2_c = kelvin_to_celsius(t2_k)
    e = calculate_nonsaturation_vapour_pressure(t2_k, rh)
    at = t2_c + 0.33 * e - 0.7 * va - 4
    at_k = celsius_to_kelvin(at)

    return at_k


def calculate_wind_chill(t2_k, va):
    """
    Wind Chill
        :param t2_k: (float array) 2m Temperature [K]
        :param va: (float array) wind speed at 10 meters [m/s]
        returns wind chill [K]
        Computation is only valid for temperatures between -50°C and 5°C and wind speeds between 5km/h and 80km/h.
        For input values outside those ranges, computed results not be considered valid.
    Reference: Blazejczyk et al. (2012)
    https://doi.org/10.1007/s00484-011-0453-2
    See also: https://web.archive.org/web/20130627223738/http://climate.weatheroffice.gc.ca/prods_servs/normals_documentation_e.html  # noqa
    """
    t2_c = kelvin_to_celsius(t2_k)  # kelvin_to_celsius(tk)
    v = va * 3.6  # convert to kilometers per hour
    windchill = 13.12 + 0.6215 * t2_c - 11.37 * v**0.16 + 0.3965 * t2_c * v**0.16
    windchill_k = celsius_to_kelvin(windchill)

    return windchill_k


def calculate_heat_index_simplified(t2_k, rh):
    """
    Heat Index
        :param t2m: (float array) 2m temperature [K]
        :param rh: (float array) relative humidity [%]
        returns heat index [K]
    Reference: Blazejczyk et al. (2012)
    https://doi.org/10.1007/s00484-011-0453-2
    """

    t2_c = kelvin_to_celsius(t2_k)

    hiarray = [
        8.784695,
        1.61139411,
        2.338549,
        0.14611605,
        1.2308094e-2,
        1.6424828e-2,
        2.211732e-3,
        7.2546e-4,
        3.582e-6,
    ]
    hi = np.copy(t2_c)

    hi_filter1 = np.where(t2_c > 20)

    hi[hi_filter1] = (
        -hiarray[0]
        + hiarray[1] * t2_c[hi_filter1]
        + hiarray[2] * rh[hi_filter1]
        - hiarray[3] * t2_c[hi_filter1] * rh[hi_filter1]
        - hiarray[4] * t2_c[hi_filter1] ** 2
        - hiarray[5] * rh[hi_filter1] ** 2
        + hiarray[6] * t2_c[hi_filter1] ** 2 * rh[hi_filter1]
        + hiarray[7] * t2_c[hi_filter1] * rh[hi_filter1] ** 2
        - hiarray[8] * t2_c[hi_filter1] ** 2 * rh[hi_filter1] ** 2
    )

    hi_k = celsius_to_kelvin(hi)

    return hi_k


def calculate_heat_index_adjusted(t2_k, td_k):
    """
    Heat Index adjusted
       :param t2_k: (float array) 2m temperature [K]
       :param td_k: (float array) 2m dewpoint temperature  [K]
       returns heat index [K]
    Reference: https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml
    """

    rh = calculate_relative_humidity_percent(t2_k, td_k)
    t2_f = kelvin_to_fahrenheit(t2_k)

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

    hi_initial = 0.5 * (t2_f + 61 + ((t2_f - 68) * 1.2) + (rh * 0.094))

    hi = (
        -hiarray[0]
        + hiarray[1] * t2_f
        + hiarray[2] * rh
        - hiarray[3] * t2_f * rh
        - hiarray[4] * t2_f**2
        - hiarray[5] * rh**2
        + hiarray[6] * t2_f**2 * rh
        + hiarray[7] * t2_f * rh**2
        - hiarray[8] * t2_f**2 * rh**2
    )

    hi_filter1 = np.where(t2_f > 80)
    hi_filter2 = np.where(t2_f < 112)
    hi_filter3 = np.where(rh <= 13)
    hi_filter4 = np.where(t2_f < 87)
    hi_filter5 = np.where(rh > 85)
    hi_filter6 = np.where(t2_f < 80)
    hi_filter7 = np.where((hi_initial + t2_f) / 2 < 80)

    f_adjust1 = hi_filter1 and hi_filter2 and hi_filter3
    f_adjust2 = hi_filter1 and hi_filter4 and hi_filter5

    adjustment1 = (
        (13 - rh[f_adjust1]) / 4 * np.sqrt(17 - np.abs(t2_f[f_adjust1] - 95) / 17)
    )

    adjustment2 = (rh[f_adjust2] - 85) / 10 * ((87 - t2_f[f_adjust2]) / 5)

    adjustment3 = 0.5 * (
        t2_f[hi_filter6]
        + 61.0
        + ((t2_f[hi_filter6] - 68.0) * 1.2)
        + (rh[hi_filter6] * 0.094)
    )

    hi[f_adjust1] = hi[f_adjust1] - adjustment1

    hi[f_adjust2] = hi[f_adjust2] + adjustment2

    hi[hi_filter6] = adjustment3

    hi[hi_filter7] = hi_initial[hi_filter7]

    hi_k = fahrenheit_to_kelvin(hi)

    return hi_k

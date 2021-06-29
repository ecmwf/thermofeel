# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from math import cos, pi, sin

import numpy as np

to_radians = pi / 180


def __julian_date(d, m, y):
    return (
        d
        - 32075
        + 1461 * (y + 4800 + (m - 14) / 12) / 4
        + 367 * (m - 2 - (m - 14) / 12 * 12) / 12
        - 3 * ((y + 4900 + (m - 14) / 12) / 100) / 4
    )


# solar declination angle [degrees] + time correction for solar angle
def __declination_angle(jd, h):
    g = (360 / 365.25) * (jd + (h / 24))  # fractional year g in degrees
    while g > 360:
        g = g - 360
    grad = g * to_radians
    # declination in [degrees]
    d = (
        0.396372
        - 22.91327 * cos(grad)
        + 4.02543 * sin(grad)
        - 0.387205 * cos(2 * grad)
        + 0.051967 * sin(2 * grad)
        - 0.154527 * cos(3 * grad)
        + 0.084798 * sin(3 * grad)
    )
    # time correction in [ h.degrees ]
    tc = (
        0.004297
        + 0.107029 * cos(grad)
        - 1.837877 * sin(grad)
        - 0.837378 * cos(2 * grad)
        - 2.340475 * sin(2 * grad)
    )
    return d, tc


# validate data convert int float to numpy array
def __wrap(variable):
    if isinstance(variable, int) or isinstance(variable, float):
        variable = np.array(variable)
        return variable
    if variable is None:
        raise ValueError
    else:
        return variable


# convert farenheit to kelvin
def __farenheit_to_kelvin(t2m):
    t2m = 5 * (t2m - 273) / 9 + 32
    return t2m


def __kelvin_to_farenheit(t2m):
    t2m = (t2m - 273.15) * 9 / 5 + 32
    return t2m


def __farenheit_to_celcius(t2m):
    t2m = (t2m - 32) * 5 / 9
    return t2m


# convert kelvin to celcius
def __kelvin_to_celcius(t2m):
    t2m = t2m - 273.15
    return t2m


# convert celcius to kelvin
def __celcius_to_kelvin(t2m):
    t2m = t2m + 273.15
    return t2m


# convert from pa to hpa for e (relative humidity)
def __pa_to_hpa(rh):
    rh = rh / 10
    return rh

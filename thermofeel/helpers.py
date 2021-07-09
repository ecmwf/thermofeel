# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import math

import numpy as np

to_radians = math.pi / 180


def to_julian_date(d, m, y):
    return (
        d
        - 32075
        + 1461 * (y + 4800 + (m - 14) / 12) / 4
        + 367 * (m - 2 - (m - 14) / 12 * 12) / 12
        - 3 * ((y + 4900 + (m - 14) / 12) / 100) / 4
    )


# validate data convert int float to numpy array
def __wrap(variable):
    if isinstance(variable, int) or isinstance(variable, float):
        variable = np.array(variable)
        return variable
    if variable is None:
        raise ValueError
    else:
        return variable


# convert Fahrenheit to Kelvin
# def __fahrenheit_to_kelvin(tf):
#     tk = 5 * (tf - 273) / 9 + 32
#     return tk


# convert Kelvin to Fahrenheit
def kelvin_to_fahrenheit(tk):
    tf = (tk - 273.15) * 9 / 5 + 32
    return tf


# convert Fahrenheit to Celsius
def fahrenheit_to_celsius(tf):
    tc = (tf - 32) * 5 / 9
    return tc


# convert from kPa to hPa for e (saturation water vapour pressure)
def kPa_to_hPa(rh_kpa):
    rh_hpa = rh_kpa / 10
    return rh_hpa

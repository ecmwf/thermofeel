# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


# Converters


# convert Celsius to Kelvin
def celsius_to_kelvin(tc):
    tk = tc + 273.15
    return tk


# convert Kelvin to Celsius
def kelvin_to_celsius(tk):
    tc = tk - 273.15
    return tc


# convert Kelvin to Fahrenheit
def kelvin_to_fahrenheit(tk):
    tf = (tk - 273.15) * 9 / 5 + 32
    return tf


# convert Fahrenheit to Celsius
def fahrenheit_to_celsius(tf):
    tc = (tf - 32) * 5 / 9
    return tc


# convert Fahrenheit to Kelvin
def fahrenheit_to_kelvin(tf):
    tk = (tf + 459.67) * 5 / 9
    return tk

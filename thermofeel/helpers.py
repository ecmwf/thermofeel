# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import functools
import math
import time
import os

import numpy as np

to_radians = math.pi / 180


def timer(func):
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        tic = time.perf_counter()
        value = func(*args, **kwargs)
        toc = time.perf_counter()
        elapsed_time = toc - tic
        print(f"{func} elapsed time: {elapsed_time:0.6f} s")
        return value

    return wrapper_timer


optnumba_jit_functions = {}


def optnumba_jit(_func=None, *, nopython=True, nogil=True, parallel=True):
    def decorator_optnumba(func):
        @functools.wraps(func)
        def jited_function(*args, **kwargs):

            global optnumba_jit_functions

            if func in optnumba_jit_functions:
                return optnumba_jit_functions[func](*args, **kwargs)

            if os.environ.get('THERMOFEEL_NO_NUMBA'):
                optnumba_jit_functions[func] = func
            else:
                try:
                    import numba

                    print(
                        f"Numba trying to compile {func}, args: nopython {nopython} nogil {nogil} parallel {parallel}"
                    )
                    optnumba_jit_functions[func] = numba.jit(
                        nopython=nopython, nogil=nogil, parallel=parallel
                    )(func)
                except Exception as e:
                    print(
                        f"Numba compilation failed for {func}, reverting to pure python code -- Exception caught: {e}"
                    )
                    optnumba_jit_functions[func] = func

            assert (
                func in optnumba_jit_functions
                and optnumba_jit_functions[func] is not None
            )

            return optnumba_jit_functions[func](*args, **kwargs)

        return jited_function

    if _func is None:
        return decorator_optnumba
    else:
        return decorator_optnumba(_func)


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

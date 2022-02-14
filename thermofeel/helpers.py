# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import functools
import math
import os
import time

to_radians = math.pi / 180

func_timers = {}


def timer(func):
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        entry = func_timers.get(func, {})

        tic = time.perf_counter()
        value = func(*args, **kwargs)
        toc = time.perf_counter()
        elapsed = toc - tic

        total = elapsed + entry.get("elapsed", 0)
        entry["elapsed"] = total
        entry["calls"] = entry.get("calls", 0) + 1
        func_timers[func] = entry

        # print(f"{func} elapsed time: {elapsed:0.6f} s")

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

            if os.environ.get("THERMOFEEL_NO_NUMBA"):
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


# convert from kPa to hPa for e (saturation water vapour pressure)
def kPa_to_hPa(rh_kpa):
    rh_hpa = rh_kpa / 10
    return rh_hpa

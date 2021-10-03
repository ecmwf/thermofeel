# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import sys

import eccodes
import numpy as np
from grib import decode_grib

from thermofeel.thermofeel import (
    calculate_cos_solar_zenith_angle,
    calculate_cos_solar_zenith_angle_integrated,
)


# Compute based on Numerical integration
def calc_cossza(message, begin, end):

    lats = message["lats"]
    lons = message["lons"]

    # print(lats.size)
    # print(lons.size)

    assert lats.size == lons.size

    dt = message["forecast_datetime"]

    # print(dt.year, dt.month, dt.day, dt.hour)

    # in hours
    h_begin = begin
    h_end = end
    nsplits = 2 * (end - begin)  # 2 x (3-0) = 6   |---|---|---|---|---|---|
    #                                                 0                       3
    assert nsplits > 0

    time_steps = np.linspace(h_begin, h_end, num=nsplits)

    integral = np.zeros(lats.size)
    cossza = np.zeros(lats.size)
    last_t = -1
    for s in range(len(time_steps) - 1):
        # print(f"interval {s+1}")
        # simpsons rule
        ti = time_steps[s]
        tf = time_steps[s + 1]
        t = [ti, (tf + ti) / 2, tf]
        print(t)
        w = ((tf - ti) / 6) * np.array([1, 4, 1])

        for n in range(len(w)):
            time_h = t[n]
            if time_h != last_t:  # don't recompute if last evaluation is same point
                cossza = calculate_cos_solar_zenith_angle(
                    lat=lats, lon=lons, y=dt.year, m=dt.month, d=dt.day, h=time_h
                )
            last_t = time_h
            integral += w[n] * cossza

    integral /= end - begin

    return integral


# uses the calculate_cos_solar_zenith_angle_integrated
# this is NOT the numerical integration
def calc_cossza_int(message, begin, end):

    lats = message["lats"]
    lons = message["lons"]
    assert lats.size == lons.size

    dt = message["forecast_datetime"]
    # print(dt.year, dt.month, dt.day, dt.hour)

    integral = calculate_cos_solar_zenith_angle_integrated(
        lat=lats,
        lon=lons,
        y=dt.year,
        m=dt.month,
        d=dt.day,
        h=dt.hour,
        tbegin=begin,
        tend=end,
    )

    return integral


def main():

    msgs = decode_grib(sys.argv[1], True)

    output = open(sys.argv[2], "wb")

    print(msgs)
    step_begin = 0
    for m in msgs:
        step_end = int(m["step"])
        print(f"Interval [{step_begin},{step_end}]")

        # integrated_cossza = calc_cossza(m, step_begin, step_end)
        integrated_cossza = calc_cossza_int(m, step_begin, step_end)

        handle = eccodes.codes_clone(m["grib"])
        eccodes.codes_set_values(handle, integrated_cossza)
        eccodes.codes_write(handle, output)
        eccodes.codes_release(handle)

        step_begin = step_end


if __name__ == "__main__":
    sys.exit(main())

# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import sys

import eccodes
from grib import decode_grib

from thermofeel.thermofeel import calculate_cos_solar_zenith_angle_integrated


# this is specific to ECMWF's forecast
def forecast_step_interval(step):
    assert step % 3 == 0  # steps are multiples of 3
    if step == 0:
        return (0, 0)
    if step <= 144:
        return (step - 3, step)
    if step > 144:
        return (step - 6, step)


def main():
    msgs = decode_grib(sys.argv[1], True)

    output = open(sys.argv[2], "wb")

    for m in msgs:
        lats = m["lats"]
        lons = m["lons"]
        assert lats.size == lons.size

        dt = m["base_datetime"]

        ftime = int(m["time"] / 100)

        step_begin, step_end = forecast_step_interval(int(m["step"]))

        print(f"Date {dt} -- Time {ftime} -- Interval [{step_begin},{step_end}]")

        cossza = calculate_cos_solar_zenith_angle_integrated(
            lat=lats,
            lon=lons,
            y=dt.year,
            m=dt.month,
            d=dt.day,
            h=ftime,
            tbegin=step_begin,
            tend=step_end,
        )

        # encode results in GRIB
        handle = eccodes.codes_clone(m["grib"])
        eccodes.codes_set_values(handle, cossza)
        eccodes.codes_write(handle, output)
        eccodes.codes_release(handle)


if __name__ == "__main__":
    sys.exit(main())

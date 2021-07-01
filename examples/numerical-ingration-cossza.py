# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import sys

import numpy as np
from grib import decode_grib

from thermofeel.thermofeel import calculate_cos_solar_zenith_angle


def calc_cossza(message):

    lats = message["lats"]
    lons = message["lons"]
    assert lats.size == lons.size

    print(lats.size)

    dt = message["forecast_datetime"]

    print(dt.year, dt.month, dt.day, dt.hour)

    shape = (message["Nj"], message["Ni"])

    latsmat = np.reshape(lats, shape)
    lonsmat = np.reshape(lons, shape)

    h_begin = 0
    h_end = 6
    nsplits = 10

    time_steps = np.linspace(h_begin, h_end, num=nsplits * h_end - 1)

    integral = np.zeros(shape)

    for s in range(len(time_steps) - 1):
        print(s)
        # simpsons rule
        ti = time_steps[s]
        tf = time_steps[s + 1]
        t = [ti, (tf + ti) / 2, tf]
        print(t)
        w = ((tf - ti) / 6) * np.array([1, 4, 1])

        for n in range(len(w)):
            cossza = calculate_cos_solar_zenith_angle(
                lat=lats, lon=lons, y=dt.year, m=dt.month, d=dt.day, h=t[n]
            )
            integral += w[n] * np.reshape(cossza, shape)

    print(np.max(integral))

    fname = sys.argv[2]
    print("=> ", fname)
    np.savez(fname, lats=latsmat, lons=lonsmat, values=integral)


def main():
    msgs = decode_grib(sys.argv[1])
    print(msgs)
    for m in msgs:
        calc_cossza(m)


if __name__ == "__main__":
    sys.exit(main())

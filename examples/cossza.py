# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))

import thermofeel

from thermofeel.thermofeel import *
from grib import *

def calc_cossza(message):
    
    lats = message["lats"]
    lons = message["lons"]
    assert lats.size == lons.size

    print(lats.size)

    dt = message["datetime"]
    date = message["date"]
    time = message["time"]
    step = message["step"]

    print(dt.year, dt.month, dt.day, dt.hour)

    shape = (message["Nj"], message["Ni"])

    latsmat = np.reshape(lats, shape)
    lonsmat = np.reshape(lons, shape)

    # vectorised computation
    cossza = calculate_cos_solar_zenith_angle(lat=lats, lon=lons, y=dt.year, m=dt.month, d=dt.day, h=dt.hour)
    valsmat = np.reshape(cossza, shape)
    
    print(valsmat)

    new_msg = message.copy()
    new_msg['values'] = valsmat
    new_msg['paramId'] = "214001" # cossza from GRIB database

    return new_msg

def main():
    fout = open(sys.argv[2], 'wb')
    try:
        msgs = decode_grib(sys.argv[1], keep=True)
        print(msgs)
        for m in msgs:
            n = calc_cossza(m)
            encode_grib(n, fout)

    except eccodes.CodesInternalError as err:
        if eccodes.VERBOSE:
            eccodes.traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write(err.msg + '\n')

        return 1


if __name__ == "__main__":
    sys.exit(main())

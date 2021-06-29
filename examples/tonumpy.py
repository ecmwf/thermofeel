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


def save(message):

    lats = message["lats"]
    lons = message["lons"]
    vals = message["values"]

    assert lats.size == lons.size
    assert lats.size == vals.size

    print(lats.size)

    shape = (message["Nj"], message["Ni"])

    latsmat = np.reshape(lats, shape)
    lonsmat = np.reshape(lons, shape)
    valsmat = np.reshape(vals, shape)

    np.savez(sys.argv[2], lats=latsmat, lons=lonsmat, values=valsmat)


def main():
    msgs = decode_grib(sys.argv[1])
    for m in msgs:
        save(m)


if __name__ == "__main__":
    sys.exit(main())

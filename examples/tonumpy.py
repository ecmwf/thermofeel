import sys
import numpy as np

from thermofeel.grib import *
from datetime import datetime, timedelta, timezone

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
    try:
        msgs = decode_grib(sys.argv[1])
        for m in msgs:
            save(m)

    except eccodes.CodesInternalError as err:
        if eccodes.VERBOSE:
            eccodes.traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write(err.msg + '\n')
        return 1    


if __name__ == "__main__":
    sys.exit(main())

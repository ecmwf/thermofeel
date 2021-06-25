import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))

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

    # scalar computation
    # cossza = []
    # for i in range(len(lats)):
    #     v = calculate_cos_solar_zenith_angle_integrated(lat=lats[i], lon=lons[i], y=dt.year, m=dt.month, d=dt.day, h=dt.hour, base=time / 100, step=3)
    #     # v = calculate_cos_solar_zenith_angle(lat=lats[i], lon=lons[i], y=dt.year, m=dt.month, d=dt.day, h=dt.hour)
    #     cossza.append(v)

    shape = (message["Nj"], message["Ni"])

    latsmat = np.reshape(lats, shape)
    lonsmat = np.reshape(lons, shape)

    h_begin = 0
    h_end = 6
    nsplits = 10
     
    time_steps = np.linspace(h_begin, h_end, num=nsplits*h_end-1)

    integral = np.zeros(shape)

    for s in range(len(time_steps) - 1):
        print(s)
        # simpsons rule
        ti = time_steps[s]
        tf = time_steps[s+1]
        t = [ti, (tf + ti) / 2, tf]
        print(t)
        w = ((tf - ti) / 6) * np.array([1, 4, 1])

        for n in range(len(w)):
            cossza = calculate_cos_solar_zenith_angle(lat=lats, lon=lons, y=dt.year, m=dt.month, d=dt.day, h=t[n])
            integral += w[n] * np.reshape(cossza, shape)

    print(np.max(integral))

    fname = sys.argv[2]
    print('=> ', fname)
    np.savez(fname, lats=latsmat, lons=lonsmat, values=integral)



def main():
    try:
        msgs = decode_grib(sys.argv[1])
        print(msgs)
        for m in msgs:
            calc_cossza(m)

    except eccodes.CodesInternalError as err:
        if eccodes.VERBOSE:
            eccodes.traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write(err.msg + '\n')

        return 1


if __name__ == "__main__":
    sys.exit(main())

# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import sys

import numpy as np
from grib import decode_grib, encode_grib

from thermofeel.thermofeel import (
    calculate_cos_solar_zenith_angle,
    calculate_mean_radiant_temperature,
    calculate_utci,
    calculate_wbgts,
)


def calc_10m_wind_speed(messages):

    m10u = messages["10u"]
    m10v = messages["10v"]

    u = m10u["values"]
    v = m10v["values"]

    speed = np.sqrt(u ** 2 + v ** 2)

    # print(f"speed --> {speed}")

    m = m10u.copy()
    m["values"] = speed
    m["paramId"] = "207"
    m["shortName"] = "10si"

    messages["10si"] = m


def calc_mrt(messages):

    ssrd = messages["ssrd"]["values"]
    ssr = messages["ssr"]["values"]
    fdir = messages["fdir"]["values"]
    strd = messages["strd"]["values"]
    strr = messages["str"]["values"]
    cossza = messages["uvcossza"]["values"]

    mrt = calculate_mean_radiant_temperature(ssrd, ssr, fdir, strd, strr, cossza)

    m = messages["ssrd"].copy()
    m["values"] = mrt
    m["paramId"] = "261002"
    m["shortName"] = "mrt"

    messages["mrt"] = m


def calc_utci(messages):

    t2m = messages["2t"]["values"]
    speed = messages["10si"]["values"]
    mrt = messages["mrt"]["values"]

    utci = calculate_utci(t2m, speed, mrt)

    m = messages["2t"].copy()
    m["values"] = utci
    m["paramId"] = "261001"  # from GRIB database
    m["shortName"] = "utci"
    m["edition"] = 2

    messages["utci"] = m


def calc_wbgt(messages):

    t2m = messages["2t"]["values"]

    wbgt = calculate_wbgts(t2m)

    m = messages["2t"].copy()
    m["values"] = wbgt
    m["paramId"] = "000000"  # no currently on GRIB database
    m["shortName"] = "wbgt"

    messages["wbgt"] = m


def calc_cossza(messages):

    message = messages["2t"]  # reference for all common data

    lats = message["lats"]
    lons = message["lons"]
    assert lats.size == lons.size

    dt = message["datetime"]

    cossza = calculate_cos_solar_zenith_angle(
        lat=lats, lon=lons, y=dt.year, m=dt.month, d=dt.day, h=dt.hour
    )

    mcossza = message.copy()
    mcossza["values"] = cossza
    mcossza["paramId"] = "214001"  # cossza from GRIB database
    mcossza["shortName"] = "uvcossza"

    messages["uvcossza"] = mcossza


def main():
    inpath = sys.argv[1]
    outpath = sys.argv[2]

    msgs = {}

    for m in decode_grib(inpath, keep=True):
        name = m["shortName"]
        msgs[name] = m

    calc_10m_wind_speed(msgs)
    calc_cossza(msgs)
    calc_mrt(msgs)
    calc_utci(msgs)
    calc_wbgt(msgs)

    # for m in msgs:
    #     print(f"{m} --> {msgs[m]}")

    fout = open(outpath, "wb")
    encode_grib(msgs["utci"], fout)

    return 0


if __name__ == "__main__":
    sys.exit(main())

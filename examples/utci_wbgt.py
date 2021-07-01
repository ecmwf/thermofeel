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
    calculate_cos_solar_zenith_angle_integrated,
    calculate_mean_radiant_temperature,
    calculate_relative_humidity_percent,
    calculate_saturation_vapour_pressure,
    calculate_utci,
    calculate_wbgt,
    celcius_to_kelvin,
)


def calc_10m_wind_speed(messages):

    m10u = messages["10u"]
    m10v = messages["10v"]

    u = m10u["values"]
    v = m10v["values"]

    va = np.sqrt(u ** 2 + v ** 2)

    print(f"va --> {va}")

    m = m10u.copy()
    m["values"] = va
    m["paramId"] = "207"
    m["shortName"] = "10si"

    messages["10si"] = m


def calc_mrt(messages):

    step = messages["ssrd"]["step"]

    factor = 1.0 / (step * 3600.0)  # deaccumulate forecast values

    ssrd = messages["ssrd"]["values"] * factor
    ssr = messages["ssr"]["values"] * factor
    fdir = messages["fdir"]["values"] * factor
    strd = messages["strd"]["values"] * factor
    strr = messages["str"]["values"] * factor
    cossza = messages["uvcossza"]["values"] * factor

    mrt = calculate_mean_radiant_temperature(ssrd, ssr, fdir, strd, strr, cossza)

    print(f"mrt --> {mrt}")

    m = messages["ssrd"].copy()
    m["values"] = mrt
    m["paramId"] = "261002"
    m["shortName"] = "mrt"

    messages["mrt"] = m


def calc_utci(messages):

    t2m = messages["2t"]["values"]  # Kelvin
    t2d = messages["2d"]["values"]  # Kelvin
    va = messages["10si"]["values"]  # m/s
    mrt = messages["mrt"]["values"]  # Kelvin

    # calculate ehPa
    rh_pc = calculate_relative_humidity_percent(t2m, t2d)
    print(f"rh_pc --> {rh_pc}")
    ehPa = calculate_saturation_vapour_pressure(t2m) * rh_pc / 100.0
    print(f"ehPa --> {ehPa}")

    utci = calculate_utci(t2m, va, mrt, ehPa)  # in Celsius

    # convert to Kelvin
    missingValueFilter = np.where(utci != -9999)
    utci[missingValueFilter] = utci[missingValueFilter] + 273.15

    # TODO: REMOVE THIS HACK
    missingValueFilter = np.where(utci == -9999)
    utci[missingValueFilter] = 9999
    print(f"utci --> {utci}")

    m = messages["2t"].copy()
    m["values"] = utci
    # TODO: REMOVE THIS HACK -- writes UTCI as 2t
    # m["paramId"] = "261001"  # from GRIB database
    # m["shortName"] = "utci"
    # m["edition"] = 2

    messages["utci"] = m


def calc_wbgt(messages):

    t2m = messages["2t"]["values"]  # Kelvin
    mrt = messages["mrt"]["values"]  # Kelvin
    va = messages["10si"]["values"]  # m/s

    wbgt = calculate_wbgt(t2m, mrt, va)
    wbgt = celcius_to_kelvin(wbgt)

    print(f"wbgt --> {wbgt}")

    m = messages["2t"].copy()
    m["values"] = wbgt
    # no currently on GRIB database -- so we write this as 2t
    # m["paramId"] = "000000"
    # m["shortName"] = "wbgt"

    messages["wbgt"] = m


def calc_cossza(messages):

    message = messages["2t"]  # reference for all common data
    step = message["step"]
    base = message["time"] / 100

    print(f"base {base} step {step}")

    lats = message["lats"]
    lons = message["lons"]
    assert lats.size == lons.size

    dt = message["base_datetime"]
    print(f"date {dt}")

    cossza = calculate_cos_solar_zenith_angle_integrated(
        lat=lats,
        lon=lons,
        y=dt.year,
        m=dt.month,
        d=dt.day,
        h=dt.hour,
        base=base,
        step=step,
    )

    print(f"cossza --> {cossza}")

    mcossza = message.copy()
    mcossza["values"] = cossza
    mcossza["paramId"] = "214001"  # cossza from GRIB database
    mcossza["shortName"] = "uvcossza"

    messages["uvcossza"] = mcossza


def main():
    inpath = sys.argv[1]
    # outpath = sys.argv[2]

    msgs = {}

    for m in decode_grib(inpath, keep=True):
        name = m["shortName"]
        msgs[name] = m

    calc_10m_wind_speed(msgs)
    calc_cossza(msgs)
    calc_mrt(msgs)
    calc_utci(msgs)
    calc_wbgt(msgs)

    utci_grib = open("utci.grib", "wb")
    encode_grib(msgs["utci"], utci_grib)

    wbgt_grib = open("wbgt.grib", "wb")
    encode_grib(msgs["wbgt"], wbgt_grib)

    return 0


if __name__ == "__main__":
    sys.exit(main())

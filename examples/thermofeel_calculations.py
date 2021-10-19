# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
import metview as mv
import thermofeel
import numpy as np
from datetime import datetime, timedelta, timezone
from math import floor
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import eccodes
import sys
import math


def decode_grib(fpath, keep=False):
    messages = []
    i = 0
    with open(fpath, "rb") as f:
        while True:
            msg = eccodes.codes_any_new_from_file(f)
            if msg is None:
                break

            i += 1
            # print("message ", i)

            md = dict()

            # decode metadata

            # loop metadata key-values
            # it = eccodes.codes_keys_iterator_new(msg, 'mars')
            # while eccodes.codes_keys_iterator_next(it):
            #     k = eccodes.codes_keys_iterator_get_name(it)
            #     v = eccodes.codes_get_string(msg, k)
            #     print("%s = %s" % (k, v))
            # eccodes.codes_keys_iterator_delete(it)

            md["paramId"] = eccodes.codes_get_string(msg, "paramId")
            md["shortName"] = eccodes.codes_get_string(msg, "shortName")

            md["Ni"] = eccodes.codes_get_long(msg, "Ni")
            md["Nj"] = eccodes.codes_get_long(msg, "Nj")

            md["time"] = eccodes.codes_get_long(msg, "time")
            md["date"] = eccodes.codes_get_string(msg, "date")
            md["step"] = eccodes.codes_get_double(msg, "step")

            sname = md["shortName"]
            step = md["step"]

            print(f"reading grib {sname} step {step}")

            ldate = eccodes.codes_get_long(msg, "date")
            yyyy = floor(ldate / 10000)
            mm = floor((ldate - (yyyy * 10000)) / 100)
            dd = ldate - (yyyy * 10000) - mm * 100

            md["base_datetime"] = datetime(yyyy, mm, dd, tzinfo=timezone.utc)

            forecast_datetime = (
                    datetime(yyyy, mm, dd, tzinfo=timezone.utc)
                    + timedelta(minutes=60 * md["time"] / 100)
                    + timedelta(minutes=60 * md["step"])
            )

            md["forecast_datetime"] = forecast_datetime

            # decode data
            # get the lats, lons, values
            md["lats"] = eccodes.codes_get_double_array(msg, "latitudes")
            # print(lats)
            md["lons"] = eccodes.codes_get_double_array(msg, "longitudes")
            # print(lons)
            md["values"] = eccodes.codes_get_double_array(msg, "values")
            # print(values)

            if keep:
                md["grib"] = msg  # dont close message and keep it for writing
            else:
                eccodes.codes_release(msg)

            messages.append(md)

    f.close()
    return messages

def calc_cossza_int(message, begin, end):

    lats = message["lats"]
    lons = message["lons"]
    assert lats.size == lons.size

    dt = message["base_datetime"]
    # print(dt.year, dt.month, dt.day, dt.hour)

    integral = thermofeel.calculate_cos_solar_zenith_angle_integrated(
        lat=lats,
        lon=lons,
        y=dt.year,
        m=dt.month,
        d=dt.day,
        h=dt.hour,
        tbegin=begin,
        tend=end
    )

    return integral

def calc_apparent_temp(message):

    variables = message["shortName"]
    t2m = variables["2t"]["values"]
    u10 = variables["10u"]["values"]
    v10 = variables["10v"]["values"]
    va = np.sqrt(u10 ** 2 + v10 ** 2 )

    at = thermofeel.calculate_apparent_temperature(t2m=t2m,va=va)

    return at

def calc_mean_radiant_temp(message,begin,end):

    variables = message["shortName"]
    step = message["step"]

    factor = 1.0 / (step * 3600.0)
    ssrd = variables["ssrd"]["values"]
    ssr = variables["ssr"]["values"]
    fdir = variables["fdir"]["values"]
    strd = variables["strd"]["values"]
    strr = variables["str"]["values"]
    cossza = calc_cossza_int(message=message,begin=begin,end=end)

    mrt = thermofeel.calculate_mean_radiant_temperature(ssrd = ssrd * factor,
                                                        ssr = ssr * factor,
                                                        fdir = fdir * factor,
                                                        strd = strd * factor,
                                                        strr = strr * factor,
                                                        cossza = cossza)
    return mrt

def calculate_utci(message,begin,end):

    variables= message["shortName"]
    t2m = variables["2t"]["values"]
    u10 = variables["10u"]["values"]
    v10 = variables["10v"]["values"]
    va = np.sqrt(u10 ** 2 + v10 ** 2)
    mrt = calc_mean_radiant_temp(message=message,begin=begin,end=end)
    t2d = variables["2d"]["values"]
    rh_pc = thermofeel.calculate_relative_humidity_percent(t2m, t2d)
    ehPa = thermofeel.calculate_saturation_vapour_pressure(t2m) * rh_pc / 100.0
    utci = calculate_utci(t2_k = t2m,
                          va_ms = va,
                          mrt_k = mrt,
                          e_hPa= ehPa)
    return utci

def main():

    msgs = decode_grib("agrib", True)

    output = open("anothergrib", "wb")

    for m in msgs:

        lats = m["lats"]
        lons = m["lons"]
        assert lats.size == lons.size

        dt = m["base_datetime"]

        ftime = int(m["time"] / 100)

        step_begin = ftime
        step_end = ftime + int(m["step"])
        utci = calculate_utci(message=m,begin=step_begin,end=step_end)
        print(f"Date {dt} -- Time {ftime} -- Interval [{step_begin},{step_end}]")



        # encode results in GRIB
        handle = eccodes.codes_clone(m["grib"])
        eccodes.codes_set_values(handle, utci)
        eccodes.codes_write(handle, output)
        eccodes.codes_release(handle)

if __name__ == "__main__":
    sys.exit(main())


# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import math
import sys
from datetime import datetime, timedelta, timezone

import eccodes
import numpy as np

import thermofeel as thermofeel


def mydebug(name, values):
    print(
        f"{name} avg {np.average(values)} max {np.max(values)} "
        f"min {np.min(values)} stddev {np.std(values, dtype=np.float64)}"
    )


def decode_grib(fpath):

    # print(f"decoding file {fpath}")

    prev_step = None
    prev_number = None

    msgcount = 0
    messages = {}

    with open(fpath, "rb") as f:

        while True:

            msg = eccodes.codes_any_new_from_file(f)

            if msg is None:  # end of file, stop iterating
                # print(f"yielding {len(messages)} messages")
                yield messages
                for k, m in messages.items():
                    grib = m["grib"]
                    eccodes.codes_release(grib)
                messages = {}
                break

            md = dict()
            msgcount += 1

            # loop metadata key-values
            it = eccodes.codes_keys_iterator_new(msg, "mars")
            while eccodes.codes_keys_iterator_next(it):
                k = eccodes.codes_keys_iterator_get_name(it)
                v = eccodes.codes_get_string(msg, k)
                md[k] = v
            eccodes.codes_keys_iterator_delete(it)

            # change types
            step = int(md["step"])
            number = md.get("number", None)

            # on new step or number, return/yield group of messages accumulated so far
            # and ensure proper cleanup of memory

            stop = (prev_step is not None and step != prev_step) or (
                prev_number is not None and number != prev_number
            )

            if stop:
                # print(f"yielding {len(messages)} messages")
                yield messages
                for k, m in messages.items():
                    grib = m["grib"]
                    eccodes.codes_release(grib)
                messages = {}

            prev_number = number
            prev_step = step

            # print(f"message {msgcount} mars metadata: {md}")

            # aggregate messages on step, number, assuming they are contiguous

            md["paramId"] = eccodes.codes_get_string(msg, "paramId")
            md["shortName"] = eccodes.codes_get_string(msg, "shortName")

            md["Ni"] = eccodes.codes_get_long(msg, "Ni")
            md["Nj"] = eccodes.codes_get_long(msg, "Nj")

            md["time"] = eccodes.codes_get_long(msg, "time")
            md["date"] = eccodes.codes_get_string(msg, "date")
            md["step"] = step

            sname = md["shortName"]

            # print(f"message {msgcount} step {step} number {number} param {sname}")

            ldate = eccodes.codes_get_long(msg, "date")
            yyyy = math.floor(ldate / 10000)
            mm = math.floor((ldate - (yyyy * 10000)) / 100)
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

            md["grib"] = msg  # keep grib open

            assert sname not in messages

            messages[sname] = md

    f.close()


def calc_cossza_int(messages, begin, end):

    dt = messages["2t"]["base_datetime"]
    lats = messages["2t"]["lats"]
    lons = messages["2t"]["lons"]
    assert lats.size == lons.size

    # print(dt.year, dt.month, dt.day, dt.hour)

    integral = thermofeel.calculate_cos_solar_zenith_angle_integrated(
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


def calc_apparent_temp(messages):

    t2m = messages["2t"]["values"]
    u10 = messages["10u"]["values"]
    v10 = messages["10v"]["values"]
    va = np.sqrt(u10 ** 2 + v10 ** 2)

    at = thermofeel.calculate_apparent_temperature(t2m=t2m, va=va)

    return at


def calc_humidex(messages):
    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    humidex = thermofeel.calculate_humidex(t2m=t2m, td=td)

    return humidex


def calc_rela_humid_perc(messages):
    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    rhp = thermofeel.calculate_relative_humidity_percent(t2m=t2m, td=td)

    return rhp


def calc_mrt(messages, cossza):

    step = messages["2t"]["step"]

    factor = 1.0 / (step * 3600.0)

    ssrd = messages["ssrd"]["values"]
    ssr = messages["ssr"]["values"]
    fdir = messages["fdir"]["values"]
    strd = messages["strd"]["values"]
    strr = messages["str"]["values"]

    mrt = thermofeel.calculate_mean_radiant_temperature(
        ssrd=ssrd * factor,
        ssr=ssr * factor,
        fdir=fdir * factor,
        strd=strd * factor,
        strr=strr * factor,
        cossza=cossza * factor,
    )

    return mrt


def calc_utci(messages, mrt):

    t2m = messages["2t"]["values"]
    u10 = messages["10u"]["values"]
    v10 = messages["10v"]["values"]
    t2d = messages["2d"]["values"]

    va = np.sqrt(u10 ** 2 + v10 ** 2)
    mydebug("va", va)

    rh_pc = thermofeel.calculate_relative_humidity_percent(t2m, t2d)
    mydebug("rh_pc", rh_pc)

    ehPa = thermofeel.calculate_saturation_vapour_pressure(t2m) * rh_pc / 100.0
    mydebug("ehPa", ehPa)

    utci = thermofeel.calculate_utci(t2_k=t2m, va_ms=va, mrt_k=mrt, e_hPa=ehPa)

    return utci


def check_messages(msgs):
    assert "2t" in msgs
    assert "2d" in msgs
    assert "10u" in msgs
    assert "10v" in msgs
    assert "ssrd" in msgs
    assert "ssr" in msgs
    assert "fdir" in msgs
    assert "str" in msgs
    assert "strd" in msgs

    # check grids all compatible
    lats = msgs["2t"]["lats"]
    lons = msgs["2t"]["lons"]

    assert lats.size == lons.size

    ftime = msgs["2t"]["forecast_datetime"]

    for k, m in msgs.items():
        nlats = m["lats"].size
        nlons = m["lons"].size
        assert nlats == lats.size
        assert nlons == lons.size
        assert ftime == m["forecast_datetime"]


def output_grib(output, msg, paramid, values, missing=None):
    # encode results in GRIB
    grib = msg["grib"]
    handle = eccodes.codes_clone(grib)
    eccodes.codes_set_long(handle, "edition", 2)
    eccodes.codes_set_string(handle, "paramId", paramid)
    eccodes.codes_set_values(handle, values)
    if missing is not None:
        eccodes.codes_set_double(handle, "missingValue", missing)
    eccodes.codes_write(handle, output)
    eccodes.codes_release(handle)


def main():

    output = open(sys.argv[2], "wb")

    for msgs in decode_grib(sys.argv[1]):

        check_messages(msgs)

        # print(f"loaded {len(msgs)} parameters: {list(msgs.keys())}")

        msg = msgs["2t"]

        dt = msg["base_datetime"]

        ftime = int(msg["time"] / 100)

        step_begin = ftime
        step_end = ftime + msg["step"]

        # print(f"Date {dt} -- Time {ftime} -- Interval [{step_begin},{step_end}]")

        t2 = msgs["2t"]["values"]  # Kelvin
        mydebug("2t", t2)

        # humidex = calc_humidex(messages=msgs)
        # apparenttemp = calc_apparent_temp(messages=msgs)
        # rhp = calc_rela_humid_perc(messages=msgs)

        cossza = calc_cossza_int(messages=msgs, begin=step_begin, end=step_end)
        # mydebug("cossza", cossza)

        mrt = calc_mrt(messages=msgs, cossza=cossza)
        # mydebug("mrt", mrt)

        utci = calc_utci(messages=msgs, mrt=mrt)
        utci = thermofeel.celsius_to_kelvin(utci)
        # mydebug("utci", utci)

        # output_grib(output, msg, "167", t2)
        # output_grib(output,msg,"157",rhp)
        # output_grib(output,msg,"260255",apparenttemp)
        # output_grib(output,msg,"260004",humidex) #heat index parameter ID
        output_grib(output, msg, "214001", cossza)
        output_grib(output, msg, "261001", utci, float(-9999))
        output_grib(output, msg, "261002", mrt)


if __name__ == "__main__":
    sys.exit(main())

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

UTCI_MIN_VALUE = thermofeel.celsius_to_kelvin(-80)
UTCI_MAX_VALUE = thermofeel.celsius_to_kelvin(90)
MISSING_VALUE = -9999.0

############################################################################################################


def field_stats(name, values):
    print(
        f"{name} avg {np.nanmean(values)} max {np.nanmax(values)} "
        f"min {np.nanmin(values)} stddev {np.nanstd(values, dtype=np.float64)}"
    )


############################################################################################################

lats = None
lons = None


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
            # md["lats"] = eccodes.codes_get_double_array(msg, "latitudes")
            # print(lats)
            # md["lons"] = eccodes.codes_get_double_array(msg, "longitudes")
            # print(lons)
            global lats
            if lats is None:
                lats = eccodes.codes_get_double_array(msg, "latitudes")
            global lons
            if lons is None:
                lons = eccodes.codes_get_double_array(msg, "longitudes")

            md["values"] = eccodes.codes_get_double_array(msg, "values")
            # print(values)

            md["grib"] = msg  # keep grib open

            assert sname not in messages

            messages[sname] = md

    f.close()


@thermofeel.timer
def calc_cossza_int(dt, begin, end):

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
        integration_order=2,
    )

    return integral

<<<<<<< HEAD:examples/thermofeel_calculations.py

@thermofeel.timer
def calc_heat_index_ad(messages):
    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    hia = thermofeel.calculate.heat_index_adjusted(t2m=t2m, td=td)


@thermofeel.timer
def calc_apparent_temp(messages, va):
    t2m = messages["2t"]["values"]

    at = thermofeel.calculate_apparent_temperature(t2m=t2m, va=va)

    return at


@thermofeel.timer
def calc_humidex(messages):
    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    humidex = thermofeel.calculate_humidex(t2m=t2m, td=td)

    return humidex


@thermofeel.timer
def calc_rela_humid_perc(messages):
    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    rhp = thermofeel.calculate_relative_humidity_percent(t2m=t2m, td=td)

    return rhp

=======
>>>>>>> 5fbe388a33497eb423c603c897c8483f5ad24571:examples/compute_utci.py

@thermofeel.timer
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


@thermofeel.timer
def calc_va(messages):
    u10 = messages["10u"]["values"]
    v10 = messages["10v"]["values"]

    return np.sqrt(u10 ** 2 + v10 ** 2)


@thermofeel.optnumba_jit
def calc_ehPa_(rh_pc, svp):
    return svp * rh_pc * 0.01  # / 100.0


@thermofeel.timer
def calc_ehPa(t2m, t2d):
    rh_pc = thermofeel.calculate_relative_humidity_percent(t2m, t2d)
    svp = thermofeel.calculate_saturation_vapour_pressure(t2m)
    ehPa = calc_ehPa_(rh_pc, svp)
    return ehPa


@thermofeel.timer
def calc_utci_in_kelvin(t2m, va, mrt, ehPa):
    utci = thermofeel.calculate_utci(t2_k=t2m, va_ms=va, mrt_k=mrt, e_hPa=ehPa)
    return thermofeel.celsius_to_kelvin(utci)


@thermofeel.timer
def filter_utci(t2m, va, mrt, ehPa, utci):
    e_mrt = np.subtract(mrt, t2m)

    misses = np.where(t2m >= thermofeel.celsius_to_kelvin(70))
    t = np.where(t2m <= thermofeel.celsius_to_kelvin(-70))
    misses = np.union1d(t, misses)

    t = np.where(va >= 25.0)  # 90kph
    misses = np.union1d(t, misses)

    t = np.where(ehPa > 50.0)
    misses = np.union1d(t, misses)

    t = np.where(e_mrt >= 100.0)
    misses = np.union1d(t, misses)

    t = np.where(e_mrt <= -30)
    misses = np.union1d(t, misses)

    print(f"utci missing values: {len(misses)}")

    utci[misses] = MISSING_VALUE

    return misses


@thermofeel.timer
def validate_utci(utci, misses):

    utci[misses] = np.nan

    bads = 0
    for i in range(len(utci)):
        v = utci[i]
        if not np.isnan(v) and (v < UTCI_MIN_VALUE or v > UTCI_MAX_VALUE):
            bads += 1
            print("UTCI [", i, "] = ", utci[i], " : lat/lon ", lats[i], lons[i])

    nmisses = len(misses)
    if nmisses > 0 or bads > 0:
        print(f"UTCI => MISS {nmisses} BADS {bads}")

    utci[misses] = MISSING_VALUE


@thermofeel.timer
def calc_utci(messages, mrt, va):

    t2m = messages["2t"]["values"]
    t2d = messages["2d"]["values"]

    ehPa = calc_ehPa(t2m, t2d)
    utci = calc_utci_in_kelvin(t2m, va, mrt, ehPa)

    filter_utci(t2m, va, mrt, ehPa, utci)

    # validate_utci(utci, filter_utci(t2m, va, mrt, ehPa, utci))

    return utci


@thermofeel.timer
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

    assert lats.size == lons.size

    ftime = msgs["2t"]["forecast_datetime"]

    for k, m in msgs.items():
        assert lats.size == m["values"].size
        assert ftime == m["forecast_datetime"]


# @thermofeel.timer
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


@thermofeel.timer
<<<<<<< HEAD:examples/thermofeel_calculations.py
def output_gribs(output, msg, cossza, mrt, utci, apparenttemp, humidex, hia):

    # output_grib(output, msg, "167", t2)
    # output_grib(output,msg,"157",rhp)
    output_grib(output, msg, "260255", apparenttemp)
    output_grib(output, msg, "260005", humidex)  # wind chill parameter ID
    output_grib(output, msg, "26004", hia)
    # output_grib(output, msg, "260005", windchill)
=======
def output_gribs(output, msg, cossza, mrt, utci):
>>>>>>> 5fbe388a33497eb423c603c897c8483f5ad24571:examples/compute_utci.py
    output_grib(output, msg, "214001", cossza)
    output_grib(output, msg, "261001", utci, missing=MISSING_VALUE)
    output_grib(output, msg, "261002", mrt)


cossza = None
last_step_end = 0


@thermofeel.timer
def process_step(msgs, output):

    check_messages(msgs)

    # print(f"loaded {len(msgs)} parameters: {list(msgs.keys())}")

    msg = msgs["2t"]

    step = msg["step"]
    time = msg["time"]
    dt = msgs["2t"]["base_datetime"]

    ftime = int(time / 100)

    step_begin = ftime
    step_end = ftime + step

    print(
        f"dt {dt.date().isoformat()} time {time} step {step} - [{step_begin},{step_end}]"
    )

    global cossza
    global last_step_end
    if cossza is None:
        print(f"[{step_begin},{step_end}]")
        cossza = calc_cossza_int(dt=dt, begin=step_begin, end=step_end)
    else:
        print(f"[{last_step_end},{step_end}]")
        cossza += calc_cossza_int(dt=dt, begin=last_step_end, end=step_end)

    last_step_end = step_end

    mrt = calc_mrt(messages=msgs, cossza=cossza)
    va = calc_va(messages=msgs)
    utci = calc_utci(messages=msgs, mrt=mrt, va=va)

<<<<<<< HEAD:examples/thermofeel_calculations.py
    humidex = calc_humidex(messages=msgs)
    # windchill = calc_windchill(messages=msgs, va=va)
    apparenttemp = calc_apparent_temp(messages=msgs, va=va)
    hia = calc_heat_index_ad(messages=msgs)
    output_gribs(
        output=output,
        msg=msg,
        cossza=cossza,
        mrt=mrt,
        utci=utci,
        apparenttemp=apparenttemp,
        hia=hia,
        humidex=humidex,
    )
=======
    output_gribs(output=output, msg=msg, cossza=cossza, mrt=mrt, utci=utci)
>>>>>>> 5fbe388a33497eb423c603c897c8483f5ad24571:examples/compute_utci.py


def main():

    print(f"Thermofeel version: {thermofeel.__version__}")
    print(f"Python version: {sys.version}")
    print(f"Numpy version: {np.version.version}")
    np.show_config()

    output = open(sys.argv[2], "wb")

    print("----------------------------------------")
    for msgs in decode_grib(sys.argv[1]):
        process_step(msgs, output)
        print("----------------------------------------")


if __name__ == "__main__":
    sys.exit(main())

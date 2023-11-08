# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Note: This script is intended only as example usage of thermofeel library.
#       It is designed to be used with ECMWF forecast data.
#       The function ifs_step_intervals() is used to calculate the time interval based on the forecast step.
#       This is particular to the IFS model and ECMWF's NWP operational system.

import argparse
import math
import sys
from datetime import datetime, timedelta, timezone

import eccodes
import numpy as np

import thermofeel as thermofeel

UTCI_MIN_VALUE = thermofeel.celsius_to_kelvin(-80)
UTCI_MAX_VALUE = thermofeel.celsius_to_kelvin(90)
MISSING_VALUE = -9999.0

###########################################################################################################

lats = None
lons = None

results = {}
misses = {}

###########################################################################################################


def field_stats(name, values):
    if name in misses:
        values[misses[name]] = np.nan

    print(
        f"{name} min {np.nanmin(values)} max {np.nanmax(values)} "
        f"avg {np.nanmean(values)} stddev {np.nanstd(values, dtype=np.float64)} "
        f"missing {np.count_nonzero(np.isnan(values))}"
    )

    if name in misses:
        values[misses[name]] = MISSING_VALUE


###########################################################################################################


def decode_grib(fpath):
    print(f"decoding file {fpath}")

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

            # assert sname not in messages

            messages[sname] = md

    f.close()


###########################################################################################################


# @thermofeel.timer
def calc_cossza_int(messages):
    dt = messages["2t"]["base_datetime"]
    time, step, begin, end = timestep_interval(messages)

    cossza = thermofeel.calculate_cos_solar_zenith_angle_integrated(
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

    return cossza


# @thermofeel.timer
def approximate_dsrp(messages):
    """
    In the absence of dsrp, approximate it with fdir and cossza.
    Note this introduces some amoutn of error as cossza approaches zero
    """
    fdir = messages["fdir"]["values"]
    cossza = calc_field("cossza", calc_cossza_int, messages)

    dsrp = thermofeel.approximate_dsrp(fdir, cossza)

    return dsrp


# @thermofeel.timer
def calc_heatx(messages):
    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    heatx = thermofeel.calculate.heat_index_adjusted(t2m=t2m, td=td)

    return heatx


# @thermofeel.timer
def calc_aptmp(messages):
    t2m = messages["2t"]["values"]

    ws = calc_field("ws", calc_ws, messages)

    aptmp = thermofeel.calculate_apparent_temperature(t2m=t2m, va=ws)

    return aptmp


# @thermofeel.timer
def calc_humidex(messages):
    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    humidex = thermofeel.calculate_humidex(t2m=t2m, td=td)

    return humidex


# @thermofeel.timer
def calc_rhp(messages):
    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    rhp = thermofeel.calculate_relative_humidity_percent(t2m=t2m, td=td)

    return rhp


# @thermofeel.timer
def calc_mrt(messages):
    time, step, begin, end = timestep_interval(messages)

    assert begin < end

    cossza = calc_field("cossza", calc_cossza_int, messages)

    # will use dsrp if available, otherwise approximate it
    dsrp = calc_field("dsrp", approximate_dsrp, messages)

    seconds_in_time_step = (end - begin) * 3600  # steps are in hours

    f = 1.0 / float(seconds_in_time_step)

    ssrd = messages["ssrd"]["values"]
    ssr = messages["ssr"]["values"]
    fdir = messages["fdir"]["values"]
    dsrp = messages["dsrp"]["values"]
    strd = messages["strd"]["values"]
    strr = messages["str"]["values"]

    mrt = thermofeel.calculate_mean_radiant_temperature(
        ssrd * f, ssr * f, dsrp * f, strd * f, fdir * f, strr * f, cossza * f
    )

    return mrt


def calc_field(name, func, messages):
    if name in results:
        return results[name]

    values = func(messages)

    field_stats(name, values)
    results[name] = values

    return values


# @thermofeel.timer
def calc_ws(messages):
    u10 = messages["10u"]["values"]
    v10 = messages["10v"]["values"]

    ws = np.sqrt(u10**2 + v10**2)

    return ws


def compute_ehPa_(rh_pc, svp):
    return svp * rh_pc * 0.01  # / 100.0


# @thermofeel.timer
def compute_ehPa(t2m, t2d):
    rh_pc = thermofeel.calculate_relative_humidity_percent(t2m, t2d)
    svp = thermofeel.calculate_saturation_vapour_pressure(t2m)
    ehPa = compute_ehPa_(rh_pc, svp)
    return ehPa


# @thermofeel.timer
def compute_utci_in_kelvin(t2m, ws, mrt, ehPa):
    utci = thermofeel.calculate_utci(t2_k=t2m, va_ms=ws, mrt_k=mrt, e_hPa=ehPa)
    return thermofeel.celsius_to_kelvin(utci)


# @thermofeel.timer
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

    return misses


# @thermofeel.timer
def validate_utci(utci, misses):
    utci[misses] = np.nan

    field_stats("utci", utci)

    out_of_bounds = 0
    for i in range(len(utci)):
        v = utci[i]
        if not np.isnan(v) and (v < UTCI_MIN_VALUE or v > UTCI_MAX_VALUE):
            out_of_bounds += 1
            print("UTCI [", i, "] = ", utci[i], " : lat/lon ", lats[i], lons[i])

    nmisses = len(misses)
    if nmisses > 0 or out_of_bounds > 0:
        print(f"UTCI => MISS {nmisses} out_of_bounds {out_of_bounds}")

    utci[misses] = MISSING_VALUE


# @thermofeel.timer
def calc_utci(messages):
    t2m = messages["2t"]["values"]
    t2d = messages["2d"]["values"]

    ws = calc_field("ws", calc_ws, messages)
    mrt = calc_field("mrt", calc_mrt, messages)

    ehPa = compute_ehPa(t2m, t2d)
    utci = compute_utci_in_kelvin(t2m, ws, mrt, ehPa)

    missing = filter_utci(t2m, ws, mrt, ehPa, utci)
    misses["utci"] = missing

    # validate_utci(utci, missing)

    return utci


# @thermofeel.timer
def calc_wbgt(messages):
    t2m = messages["2t"]["values"]  # Kelvin
    t2d = messages["2d"]["values"]

    ws = calc_field("ws", calc_ws, messages)
    mrt = calc_field("mrt", calc_mrt, messages)

    wbgt = thermofeel.calculate_wbgt(t2m, mrt, ws, t2d)
    wbgt = thermofeel.celsius_to_kelvin(wbgt)

    return wbgt


# @thermofeel.timer
def calc_bgt(messages):
    t2m = messages["2t"]["values"]  # Kelvin

    ws = calc_field("ws", calc_ws, messages)
    mrt = calc_field("mrt", calc_mrt, messages)

    bgt = thermofeel.calculate_bgt(t2m, mrt, ws)
    bgt = thermofeel.celsius_to_kelvin(bgt)

    return bgt


# @thermofeel.timer
def calc_wbt(messages):
    t2m = messages["2t"]["values"]
    t2m = thermofeel.kelvin_to_celcius(t2m)

    rhp = calc_field("rhp", calc_rhp, messages)

    wbt = thermofeel.calculate_relative_humidity_percent(t2m=t2m, rh=rhp)

    return wbt


# @thermofeel.timer
def calc_net(messages):
    t2m = messages["2t"]["values"]  # Kelvin
    t2d = messages["2d"]["values"]

    ws = calc_field("ws", calc_ws, messages)

    net = thermofeel.calculate_net_effective_temperature(t2m, ws, t2d)
    net = thermofeel.celsius_to_kelvin(net)

    return net


# @thermofeel.timer
def calc_windchill(messages):
    t2m = messages["2t"]["values"]

    ws = calc_field("ws", calc_ws, messages)

    windchill = thermofeel.calculate_wind_chill(t2m, ws)

    return windchill


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
    """Encode field in GRIB2"""
    grib = msg["grib"]
    handle = eccodes.codes_clone(grib)
    eccodes.codes_set_long(handle, "edition", 2)
    eccodes.codes_set_string(handle, "paramId", paramid)
    eccodes.codes_set_values(handle, values)
    if missing is not None:
        eccodes.codes_set_double(handle, "missingValue", missing)
    eccodes.codes_write(handle, output)
    eccodes.codes_release(handle)


def ifs_step_intervals(step):
    """Computes the time integration interval for the IFS forecasting system given a forecast output step"""
    # assert step != 0 and step is not None

    # assert step > 0
    assert step <= 360
    if step > 0:
        return step - 3
    else:
        return step
    # if step <= 144:
    #    assert step % 3 == 0
    #   return step - 3
    # else:
    #    if step <= 360:
    #        assert step % 6 == 0
    #        return step - 6


def timestep_interval(messages):
    msg = messages["2t"]
    step = msg["step"]  # end of the forecast integration
    time = msg["time"]
    ftime = int(time / 100)  # forecast time in hours
    integration_start = ifs_step_intervals(step)  # start of forecast integration step
    step_begin = ftime + integration_start
    step_end = ftime + step
    return time, step, step_begin, step_end


# @thermofeel.timer
def process_step(args, msgs, output):
    check_messages(msgs)

    # print(f"loaded {len(msgs)} parameters: {list(msgs.keys())}")
    template = msgs["2t"]

    dt = msgs["2t"]["base_datetime"]
    time, step, step_begin, step_end = timestep_interval(msgs)
    print(
        f"dt {dt.date().isoformat()} time {time} step {step} - [{step_begin},{step_end}]"
    )

    global results
    global misses

    results = {}
    misses = {}

    # Windspeed - shortName ws
    if args.ws:
        ws = calc_field("ws", calc_ws, msgs)
        output_grib(output, template, "10", ws)

    # Cosine of Solar Zenith Angle - shortName uvcossza - ECMWF product
    # TODO: 214001 only exists for GRIB1 -- but here we use it for GRIB2 (waiting for WMO)
    if args.cossza:
        cossza = calc_field("cossza", calc_cossza_int, msgs)
        output_grib(output, template, "214001", cossza)

    # Mean Radiant Temperature - shortName mrt - ECMWF product
    if args.mrt:
        mrt = calc_field("mrt", calc_mrt, msgs)
        output_grib(output, template, "261002", mrt)

    # Univeral Thermal Climate Index - shortName utci - ECMWF product
    if args.utci:
        utci = calc_field("utci", calc_utci, msgs)
        output_grib(output, template, "261001", utci, missing=MISSING_VALUE)

    # Heat Index (adjusted) - shortName heatx - ECMWF product
    if args.heatx:
        heatx = calc_field("heatx", calc_heatx, msgs)
        output_grib(output, template, "260004", heatx)

    # Wind Chill factor - shortName wcf - ECMWF product
    if args.windchill:
        windchill = calc_field("windchill", calc_windchill, msgs)
        output_grib(output, template, "260005", windchill)

    # Apparent Temperature - shortName aptmp - ECMWF product
    if args.aptmp:
        aptmp = calc_field("aptmp", calc_aptmp, msgs)
        output_grib(output, template, "260255", aptmp)

    # Relative humidity percent at 2m - shortName 2r - ECMWF product
    if args.rhp:
        rhp = calc_field("rhp", calc_rhp, msgs)
        output_grib(output, template, "260242", rhp)

    # Humidex - shortName hx - TO BE RELEASED as ECMWF product
    # TODO: 261016 is experimental GRIB code, update once WMO publishes
    if args.humidex:
        humidex = calc_field("hmdx", calc_humidex, msgs)
        output_grib(output, template, "261016", humidex)

    # Normal Effective Temperature - shortName nefft - TO BE RELEASED as ECMWF product
    # TODO: 212002 is experimental GRIB code, update once WMO publishes
    if args.net:
        net = calc_field("net", calc_net, msgs)
        output_grib(output, template, "261018", net)

    # Globe Temperature - shortName gt
    # TODO: 212003 is experimental GRIB code, update once WMO publishes
    if args.bgt:
        bgt = calc_field("bgt", calc_bgt, msgs)
        output_grib(output, template, "261015", bgt)

    # Wet-bulb potential temperature - shortName wbt - TO BE RELEASED as ECMWF product
    # TODO: 212004 is experimental GRIB code, update once WMO publishes
    if args.wbt:
        wbt = calc_field("wbt", calc_wbt, msgs)
        output_grib(output, template, "261022", wbt)

    # Wet Bulb Globe Temperature - shortName wbgt - TO BE RELEASED as ECMWF product
    if args.wbgt:  #
        wbgt = calc_field("wbgt", calc_wbgt, msgs)
        output_grib(output, template, "261014", wbgt)

    # effective temperature 261017
    # standard effective temperature 261019

    return step


def command_line_options():
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="input file with GRIB messages")
    parser.add_argument("output", help="output file with GRIB messages")

    parser.add_argument(
        "--ws", help="compute wind speed from components", action="store_true"
    )
    parser.add_argument(
        "--cossza",
        help="compute Cosine of Solar Zenith Angle (cossza)",
        action="store_true",
    )
    parser.add_argument("--mrt", help="compute mrt", action="store_true")
    parser.add_argument(
        "--utci",
        help="compute UTCI Universal Thermal Climate Index",
        action="store_true",
    )
    parser.add_argument(
        "--heatx", help="compute Heat Index (adjusted)", action="store_true"
    )
    parser.add_argument(
        "--windchill", help="compute Windchill factor", action="store_true"
    )
    parser.add_argument(
        "--aptmp", help="compute Apparent Temperature", action="store_true"
    )
    parser.add_argument(
        "--rhp", help="compute relative humidity percent", action="store_true"
    )

    parser.add_argument("--humidex", help="compute humidex", action="store_true")
    parser.add_argument(
        "--net", help="compute net effective temperature", action="store_true"
    )

    # TODO: these outputs are not yet in WMO GRIB2 recognised parameters
    parser.add_argument(
        "--wbgt", help="compute Wet Bulb Globe Temperature", action="store_true"
    )
    parser.add_argument("--bgt", help="compute  Globe Temperature", action="store_true")
    parser.add_argument(
        "--wbt", help="compute Wet Bulb Temperature", action="store_true"
    )

    args = parser.parse_args()

    return args


def main():
    args = command_line_options()

    print(f"Thermofeel version: {thermofeel.__version__}")
    print(f"Python version: {sys.version}")
    print(f"Numpy version: {np.version.version}")
    # np.show_config()

    output = open(args.output, "wb")

    steps = []

    print("----------------------------------------")
    for msgs in decode_grib(args.input):
        step = process_step(args, msgs, output)
        steps.append(step)
        print("----------------------------------------")

    print(f"\nProcessed steps: {steps}\n")

    print("Performance summary:")
    print("--------------------")
    for func, stats in thermofeel.func_timers.items():
        assert stats["calls"] > 0
        average = stats["elapsed"] / stats["calls"]
        print(func, "->", stats, f"average {average:0.6f} s")


if __name__ == "__main__":
    sys.exit(main())

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

############################################################################################################


def field_stats(name, values):

    print(
        f"{name} avg {np.nanmean(values)} max {np.nanmax(values)} "
        f"min {np.nanmin(values)} stddev {np.nanstd(values, dtype=np.float64)} "
        f"missing {np.count_nonzero(np.isnan(values))}"
    )


############################################################################################################

lats = None
lons = None


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

            assert sname not in messages

            messages[sname] = md

    f.close()


@thermofeel.timer
def calc_cossza_int(messages, results):

    if "cossza" in results:
        return results["cossza"]

    dt = messages["2t"]["base_datetime"]
    time, step, begin, end = timestep_interval(messages)

    # print(dt.year, dt.month, dt.day, dt.hour)
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

    field_stats("cossza", cossza)
    results["cossza"] = cossza

    return cossza


@thermofeel.timer
def calc_heatx(messages, results):
    if "heatx" in results:
        return results["heatx"]

    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    heatx = thermofeel.calculate.heat_index_adjusted(t2m=t2m, td=td)

    field_stats("heatx", heatx)
    results["heatx"] = heatx

    return heatx


@thermofeel.timer
def calc_aptmp(messages, results):
    if "aptmp" in results:
        return results["aptmp"]

    t2m = messages["2t"]["values"]

    ws = calc_ws(messages, results)

    aptmp = thermofeel.calculate_apparent_temperature(t2m=t2m, va=ws)

    field_stats("aptmp", aptmp)
    results["aptmp"] = aptmp

    return aptmp


@thermofeel.timer
def calc_humidex(messages, results):
    if "humidex" in results:
        return results["humidex"]

    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    humidex = thermofeel.calculate_humidex(t2m=t2m, td=td)

    field_stats("humidex", humidex)
    results["humidex"] = humidex

    return humidex


@thermofeel.timer
def calc_rhp(messages, results):
    if "rhp" in results:
        return results["rhp"]

    t2m = messages["2t"]["values"]
    td = messages["2d"]["values"]

    rhp = thermofeel.calculate_relative_humidity_percent(t2m=t2m, td=td)

    field_stats("rhp", rhp)
    results["rhp"] = rhp

    return rhp


@thermofeel.timer
def calc_mrt(messages, results):
    if "mrt" in results:
        return results["mrt"]

    time, step, begin, end = timestep_interval(messages)

    assert begin < end

    cossza = calc_cossza_int(messages, results)

    seconds_since_start_forecast = step * 3600
    seconds_in_time_step = (end - begin) * 3600

    f1 = 1.0 / float(seconds_since_start_forecast)
    f2 = 1.0 / float(seconds_in_time_step)

    ssrd = messages["ssrd"]["values"]
    ssr = messages["ssr"]["values"]
    fdir = messages["fdir"]["values"]
    strd = messages["strd"]["values"]
    strr = messages["str"]["values"]

    mrt = thermofeel.calculate_mean_radiant_temperature(
        ssrd=ssrd * f1,  # de-accumulate since forecast start
        ssr=ssr * f1,
        fdir=fdir * f1,
        strd=strd * f1,
        strr=strr * f1,
        cossza=cossza * f2,  # de-accumulate time step integration
    )

    field_stats("mrt", mrt)
    results["mrt"] = mrt

    return mrt


@thermofeel.timer
def calc_ws(messages, results):
    if "ws" in results:
        return results["ws"]

    u10 = messages["10u"]["values"]
    v10 = messages["10v"]["values"]
    ws = np.sqrt(u10 ** 2 + v10 ** 2)
    results["ws"] = ws
    return ws


@thermofeel.optnumba_jit
def compute_ehPa_(rh_pc, svp):
    return svp * rh_pc * 0.01  # / 100.0


@thermofeel.timer
def compute_ehPa(t2m, t2d):
    rh_pc = thermofeel.calculate_relative_humidity_percent(t2m, t2d)
    svp = thermofeel.calculate_saturation_vapour_pressure(t2m)
    ehPa = compute_ehPa_(rh_pc, svp)
    return ehPa


# @thermofeel.timer
def compute_utci_in_kelvin(t2m, ws, mrt, ehPa):
    utci = thermofeel.calculate_utci(t2_k=t2m, va_ms=ws, mrt_k=mrt, ehPa=ehPa)
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

    utci[misses] = np.nan
    field_stats("utci", utci)
    utci[misses] = MISSING_VALUE

    return misses


@thermofeel.timer
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


@thermofeel.timer
def calc_utci(messages, results):

    if "utci" in results:
        return results["utci"]

    t2m = messages["2t"]["values"]
    t2d = messages["2d"]["values"]

    ws = calc_ws(messages, results)
    mrt = calc_mrt(messages, results)

    ehPa = compute_ehPa(t2m, t2d)
    utci = compute_utci_in_kelvin(t2m, ws, mrt, ehPa)

    filter_utci(t2m, ws, mrt, ehPa, utci)
    # validate_utci(utci, filter_utci(t2m, va, mrt, ehPa, utci))

    results["utci"] = utci

    return utci


@thermofeel.timer
def calc_wbgt(messages, results):

    if "wbgt" in results:
        return results["wbgt"]

    t2m = messages["2t"]["values"]  # Kelvin
    t2d = messages["2d"]["values"]

    ws = calc_ws(messages, results)
    mrt = calc_mrt(messages, results)

    wbgt = thermofeel.calculate_wbgt(t2m, mrt, ws, t2d)
    wbgt = thermofeel.celsius_to_kelvin(wbgt)

    field_stats("wbgt", wbgt)
    results["wbgt"] = wbgt

    return wbgt


@thermofeel.timer
def calc_bgt(messages, results):

    if "bgt" in results:
        return results["bgt"]

    t2m = messages["2t"]["values"]  # Kelvin

    ws = calc_ws(messages, results)
    mrt = calc_mrt(messages, results)

    bgt = thermofeel.calculate_bgt(t2m, mrt, ws)
    bgt = thermofeel.celsius_to_kelvin(bgt)

    field_stats("bgt", bgt)
    results["bgt"] = bgt

    return bgt


@thermofeel.timer
def calc_wbt(messages, results):
    if "wbt" in results:
        return results["wbt"]

    t2m = messages["2t"]["values"]
    t2m = thermofeel.kelvin_to_celcius(t2m)

    rh = calc_rhp(messages, results)

    wbt = thermofeel.calculate_relative_humidity_percent(t2m=t2m, rh=rh)

    field_stats("wbt", wbt)
    results["wbt"] = wbt

    return wbt


@thermofeel.timer
def calc_net(messages, results):

    if "net" in results:
        return results["net"]

    t2m = messages["2t"]["values"]  # Kelvin
    t2d = messages["2d"]["values"]

    ws = calc_ws(messages, results)

    net = thermofeel.calculate_net_effective_temperature(t2m, ws, t2d)
    net = thermofeel.celsius_to_kelvin(net)

    field_stats("net", net)
    results["net"] = net

    return net


@thermofeel.timer
def calc_windchill(messages, results):
    if "windchill" in results:
        return results["windchill"]

    t2m = messages["2t"]["values"]

    ws = calc_ws(messages, results)

    windchill = thermofeel.calculate_wind_chill(t2m, ws)

    field_stats("windchill", windchill)
    results["windchill"] = windchill

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


@thermofeel.timer
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


def ifs_step_intervals(step):
    """Computes the time integration interval for the IFS forecasting system given a forecast output step"""
    assert step != 0 and step is not None

    assert step > 0
    assert step <= 360

    if step <= 144:
        assert step % 3 == 0
        return step - 3
    else:
        if step <= 360:
            assert step % 6 == 0
            return step - 6


def timestep_interval(messages):
    msg = messages["2t"]
    step = msg["step"]  # end of the forecast integration
    time = msg["time"]
    ftime = int(time / 100)  # forecast time in hours
    integration_start = ifs_step_intervals(step)  # start of forecast integration step
    step_begin = ftime + integration_start
    step_end = ftime + step
    return time, step, step_begin, step_end


@thermofeel.timer
def process_step(args, msgs, output):

    check_messages(msgs)

    # print(f"loaded {len(msgs)} parameters: {list(msgs.keys())}")
    template = msgs["2t"]

    dt = msgs["2t"]["base_datetime"]
    time, step, step_begin, step_end = timestep_interval(msgs)
    print(
        f"dt {dt.date().isoformat()} time {time} step {step} - [{step_begin},{step_end}]"
    )

    results = {}

    # Windspeed - shortName ws
    if args.ws:
        ws = calc_ws(msgs, results)
        output_grib(output, template, "10", ws)

    # Cosine of Solar Zenith Angle - shortName uvcossza
    if args.cossza:
        cossza = calc_cossza_int(msgs, results)
        output_grib(output, template, "214001", cossza)

    # Mean Radiant Temperature - shortName mrt
    if args.mrt:
        mrt = calc_mrt(msgs, results)
        output_grib(output, template, "261002", mrt)

    # Univeral Thermal Climate Index - shortName utci
    if args.utci:
        utci = calc_utci(msgs, results)
        output_grib(output, template, "261001", utci, missing=MISSING_VALUE)

    # Heat Index (adjusted) - shortName heatx
    if args.heatx:
        heatx = calc_heatx(msgs, results)
        output_grib(output, template, "260004", heatx)

    # Wet Bulb Globe Temperature - shortName wbgt
    if args.wbgt:  #
        wbgt = calc_wbgt(msgs, results)
        output_grib(output, template, "260004", wbgt)

    # Wind Chill factor - shortName wcf
    if args.windchill:
        windchill = calc_windchill(msgs, results)
        output_grib(output, template, "260005", windchill)

    # Apparent Temperature - shortName aptmp
    if args.aptmp:
        aptmp = calc_aptmp(msgs, results)
        output_grib(output, template, "260255", aptmp)

    # Relative humidity percent at 2m - shortName 2r
    if args.rhp:
        rhp = calc_rhp(msgs, results)
        output_grib(output, template, "260242", rhp)

    # Humidex - shortName hx
    # TODO: 212001 is experimental GRIB code, update once WMO publishes
    if args.humidex:
        humidex = calc_humidex(msgs, results)
        output_grib(output, template, "212001", humidex)

    # Net Effective Temperature - shortName net
    # TODO: 212002 is experimental GRIB code, update once WMO publishes
    if args.net:
        net = calc_net(msgs, results)
        output_grib(output, template, "212002", net)

    # Globe Temperature - shortName gt
    # TODO: 212003 is experimental GRIB code, update once WMO publishes
    if args.bgt:
        bgt = calc_bgt(msgs, results)
        output_grib(output, template, "212003", bgt)

    # Wet Bulb Temperature - shortName wbt
    # TODO: 212004 is experimental GRIB code, update once WMO publishes
    if args.wbt:
        wbt = calc_wbt(msgs, results)
        output_grib(output, template, "212004", wbt)

    return step


def main():

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
        "--wbgt", help="compute Wet Bulb Globe Temperature", action="store_true"
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

    # TODO: these outputs are not yet in WMO GRIB2 recognised parameters
    parser.add_argument("--humidex", help="compute humidex", action="store_true")
    parser.add_argument(
        "--net", help="compute net effective temperature", action="store_true"
    )
    parser.add_argument("--bgt", help="compute  Globe Temperature", action="store_true")
    parser.add_argument(
        "--wbt", help="compute Wet Bulb Temperature", action="store_true"
    )

    args = parser.parse_args()

    print(f"Thermofeel version: {thermofeel.__version__}")
    print(f"Python version: {sys.version}")
    print(f"Numpy version: {np.version.version}")
    np.show_config()

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

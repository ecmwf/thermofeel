import sys
import eccodes
import datetime
import numpy as np
import pandas as pd
from ThermoFeel import *

def decode_grib(fpath):
    messages = []
    i = 0
    with open(fpath, "rb") as f:
        while True:
            msg = eccodes.codes_any_new_from_file(f)
            if msg is None:
                break

            i += 1
            print("message %d ----------------------------------", i)

            md = dict()

            # decode metadata

            it = eccodes.codes_keys_iterator_new(msg, 'ls')
            while eccodes.codes_keys_iterator_next(it):
                datedata = eccodes.codes_get_string(msg, "dataDate")
                stepRange = eccodes.codes_get_string(msg, "stepRange")
            eccodes.codes_keys_iterator_delete(it)

            dateobj = pd.to_datetime(datedata, format='%Y%m%d')
            forecast_datetime = dateobj + \
                pd.to_timedelta(int(stepRange), unit='h')

            md["date"] = datedata
            md["step"] = stepRange
            md["datetime"] = forecast_datetime

            # decode data
            # get the lats, lons, values
            md["lats"] = eccodes.codes_get_double_array(msg, "latitudes")
            # print(lats)
            md["lons"] = eccodes.codes_get_double_array(msg, "longitudes")
            # print(lons)
            md["values"] = eccodes.codes_get_double_array(msg, "values")
            # print(values)
            eccodes.codes_release(msg)

            messages.append(md)

    f.close()
    return messages

def run_cossza(values):
    lat = values[1]
    lon = values[2]
    h = pd.to_datetime(values[0], format='%Y%m%d %H:%m:%s')
    print(h)
    hours = h.hour
    print(hours)
    years = h.year
    months = h.month
    days = h.day
    base = hours[0]
    step = hours[1]
    c1 = []
    c2 = []
    for x in range(len(hours)):
        for i in range(len(lat)):
            c1.append(
                calculate_solar_zenith_angle(lat=lat[i], lon=lon[i], y=years[x], m=months[x], d=days[x], h=hours[x]))
            # c2.append(
            #     calculate_solar_zenith_angle_f(lat=lat[i], lon=lon[i], y=years[x], m=months[x], d=days[x], h=hours[x], base=base, step=step))
    # print(c1, c2)
    # return(c1,c2)
    return c1


def run_heatindex(values):
    output = calculate_heat_index(values[0])
    return(output)


def main():
    try:
        v = decode_grib(sys.argv[1])
        print(v)
        # run_cossza(v)
        # run_heatindex(decode_singleparam(sys.argv[1]))
    except eccodes.CodesInternalError as err:
        if eccodes.VERBOSE:
            eccodes.traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write(err.msg + '\n')

        return 1


if __name__ == "__main__":
    sys.exit(main())

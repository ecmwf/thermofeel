import sys
import eccodes
import datetime
import numpy as np
import pandas as pd
from ThermoFeel import *

def decode_singleparam(fpath):
    values_list = []
    lat_list = []
    lon_list = []
    i = 0
    with open(fpath, "rb") as f:
        while True:
            msg = eccodes.codes_any_new_from_file(f)
            if msg is None:
                break

            i += 1
            print("message %d ----------------------------------", i)

            # print metadata key-values
            it = eccodes.codes_keys_iterator_new(msg, 'ls')
            while eccodes.codes_keys_iterator_next(it):
                k = eccodes.codes_keys_iterator_get_name(it)
                v = eccodes.codes_get_string(msg, k)
                print("%s = %s" % (k, v))
            eccodes.codes_keys_iterator_delete(it)

            # get the lats, lons, values
            lats = eccodes.codes_get_double_array(msg, "latitudes")
            # print(lats)
            lons = eccodes.codes_get_double_array(msg, "longitudes")
            # print(lons)
            values = eccodes.codes_get_double_array(msg, "values")
            # print(values)
            lat_list.append(lats)
            lon_list.append(lons)
            values_list.append(values)
            eccodes.codes_release(msg)

    f.close()
    return values_list,lat_list,lon_list


def decode_csza(fpath):
    values_list = []
    lat_list = []
    long_list = []
    i = 0
    with open(fpath, "rb") as f:
        while True:
            msg = eccodes.codes_any_new_from_file(f)
            if msg is None:
                break

            i += 1
            print("message %d ----------------------------------", i)

            # print metadata key-values
            it = eccodes.codes_keys_iterator_new(msg, 'ls')
            while eccodes.codes_keys_iterator_next(it):
                datedata = eccodes.codes_get_string(msg, "dataDate")
                stepRange = eccodes.codes_get_string(msg, "stepRange")

            eccodes.codes_keys_iterator_delete(it)

            # get the lats, lons, values
            lats = eccodes.codes_get_double_array(msg, "latitudes")
            # print(lats)
            lons = eccodes.codes_get_double_array(msg, "longitudes")
            # print(lons)
            #values = eccodes.codes_get_double_array(msg, "values")
            # print(values)
            #datedata = eccodes.codes_get_double_array(msg, "dataDate")
            #stepRange = eccodes.codes_get_double_array(msg, "stepRange")
            dateobj = pd.to_datetime(datedata, format='%Y%m%d')
            forecast_datetime = dateobj + pd.to_timedelta(int(stepRange), unit='h')
            values_list.append(forecast_datetime)
            lat_list.append(lats)
            long_list.append(lons)
            eccodes.codes_release(msg)

    f.close()
    return values_list, lat_list, long_list
def run_cossza(values):
    h = pd.to_datetime(values[0], format='%Y%m%d %H:%m:%s')
    hours = h.hour
    years = h.year
    months = h.month
    days = h.day
    lon = values[2]
    lat = values[1]
    base = hours[0]
    step = hours[1]
    c1 = []
    c2 = []
    for x in range(len(hours)):
        for i in range(len(lat)):
            c1.append(
                calculate_solar_zenith_angle(lat=lat[i], lon=lon[i], y=years[x], m=months[x], d=days[x], h=hours[x]))
            c2.append(
                calculate_solar_zenith_angle_f(lat=lat[i], lon=lon[i], y=years[x], m=months[x], d=days[x], h=hours[x], base=base, step=step))
    print(c1, c2)
    return(c1,c2)
def run_heatindex(values):
    output = calculate_heat_index(values[0])
    return(output)
def main():
    try:
        #run_cossza(decode_csza(sys.argv[1]))
        run_heatindex(decode_singleparam(sys.argv[1]))
    except eccodes.CodesInternalError as err:
        if eccodes.VERBOSE:
            eccodes.traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write(err.msg + '\n')

        return 1


if __name__ == "__main__":
    sys.exit(main())

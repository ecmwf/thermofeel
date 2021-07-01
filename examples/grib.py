# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from datetime import datetime, timedelta, timezone
from math import floor

import eccodes


def encode_grib(msg, fout):
    # print(f"{msg}")
    handle = eccodes.codes_clone(msg["grib"])
    if "edition" in msg:
        eccodes.codes_set_long(handle, "edition", msg["edition"])

    paramId = msg["paramId"]
    print(f"grib {paramId}")
    # eccodes.codes_set_long(handle, "missingValue", -9999)
    eccodes.codes_set_string(handle, "paramId", paramId)
    eccodes.codes_set_values(handle, msg["values"])

    eccodes.codes_write(handle, fout)
    eccodes.codes_release(handle)


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

            sname = md["shortName"]

            print(f"reading grib {sname}")

            md["Ni"] = eccodes.codes_get_long(msg, "Ni")
            md["Nj"] = eccodes.codes_get_long(msg, "Nj")

            md["time"] = eccodes.codes_get_long(msg, "time")
            md["date"] = eccodes.codes_get_string(msg, "date")
            md["step"] = eccodes.codes_get_double(msg, "step")

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

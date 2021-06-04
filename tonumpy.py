import sys
import eccodes
import numpy as np

from datetime import datetime, timedelta, timezone

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

            # loop metadata key-values
            # it = eccodes.codes_keys_iterator_new(msg, 'mars')
            # while eccodes.codes_keys_iterator_next(it):
            #     k = eccodes.codes_keys_iterator_get_name(it)
            #     v = eccodes.codes_get_string(msg, k)
            #     print("%s = %s" % (k, v))
            # eccodes.codes_keys_iterator_delete(it)

            # time_l = eccodes.codes_get_long(msg, "time")
            # time_f = eccodes.codes_get_double(msg, "time")

            md['Ni'] = eccodes.codes_get_long(msg, "Ni")
            md['Nj'] = eccodes.codes_get_long(msg, "Nj")

            md["time"] = eccodes.codes_get_long(msg, "time")
            md["date"] = eccodes.codes_get_string(msg, "date")
            md["step"] = eccodes.codes_get_double(msg, "step")

            ldate = eccodes.codes_get_long(msg, "date")
            yyyy = floor(ldate/10000)
            mm = floor((ldate-(yyyy*10000))/100)
            dd = ldate-(yyyy*10000)-mm*100

            forecast_datetime = \
                datetime(yyyy, mm, dd, tzinfo=timezone.utc) \
                + timedelta(minutes=60 * md["time"] / 100) \
                + timedelta(minutes=60 * md["step"])

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


def save(message):
    
    lats = message["lats"]
    lons = message["lons"]
    vals = message["values"]

    assert lats.size == lons.size
    assert lats.size == vals.size

    print(lats.size)

    shape = (message["Nj"], message["Ni"])

    latsmat = np.reshape(lats, shape)
    lonsmat = np.reshape(lons, shape)
    valsmat = np.reshape(vals, shape)

    np.savez(sys.argv[2], lats=latsmat, lons=lonsmat, values=valsmat)

def main():
    try:
        msgs = decode_grib(sys.argv[1])
        for m in msgs:
            save(m)

    except eccodes.CodesInternalError as err:
        if eccodes.VERBOSE:
            eccodes.traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write(err.msg + '\n')
        return 1    


if __name__ == "__main__":
    sys.exit(main())

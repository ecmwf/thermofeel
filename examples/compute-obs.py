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
import csv
import sys

import dateutil.parser
import numpy as np

import thermofeel as thermofeel

#######################################################################################################################

# INPUTS
# BAR_dif.tab   DIF [W/m**2]  -- ssrd = ssrdir + ssrdiff --> ssrdiff = ssrd - fdir
# BAR_dir.tab   DIR [W/m**2]  --
# BAR_lwd.tab   LWD [W/m**2]  -- strd   (thermal radiadtion downwards)
# BAR_lwu.tab   LWU [W/m**2]  -- stru = strd - strr   thermal radiadtion upwards
# BAR_rh.tab    RH  [%]       -- use to compute td with calculate_dew_point_from_relative_humidity(rh, t2m)
# BAR_swd.tab   SWD [W/m**2]  -- ssrd
# BAR_swu.tab   SWU [W/m**2]  -- ssr = swd - swu --> swu = swd - ssr --> swu = ssrd - ssr
# BAR_t2m.tab   T2 [°C]       -- t2m

# ssrd = swd
# ssr = ssrd - swu
# dsrp =
# strd = lwd
# fdir =
# strr = lwd - lwu
# cossza


# Five stations with wind data, IDS: BAR  IZA  NYA  PAY   TAT

# files with names <STATIONID>_{dif,dir,lwd,lwu,rh,swd,swu,t2m,wind}.tab

# OUTPUTS, CSV file with each station
# cossza
# mrt
# wbgt

#######################################################################################################################


def command_line_options():
    parser = argparse.ArgumentParser()

    parser.add_argument("output", help="output filename to dump results in CSV format")
    parser.add_argument("--basepath", help="basename for input files")
    parser.add_argument(
        "--stationid", help="station 3 char identifier (e.g. TAT for Tateno) "
    )

    parser.add_argument("--max", help="max number of items", type=int, default=0)

    args = parser.parse_args()

    return args


def read_part_dataset(obsdata, varname, filepath, nmax=0):
    with open(filepath, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # skip first line

        nitems = 0
        for record in reader:
            x = record[0]
            y = record[1]
            if y.startswith("ø"):
                y = float(y[1:])

            # print(x, y)
            if x in obsdata:
                obsdata[x].update([(varname, y)])
            else:
                obsdata[x] = dict([(varname, y)])

            nitems = nitems + 1

            if nmax and nitems >= nmax:
                break

        print(f"{filepath}: {nitems} record(s)")


def station_location(id):
    # IZA () * Latitude:  * Longitude:
    # NYA () * Latitude: 78.922700 * Longitude:
    #
    return {
        "PAY": {
            "location": "Payerne",
            "lat": 46.812300,
            "lon": 6.942200,
        },  # 'url': 'https://doi.pangaea.de/10.1594/PANGAEA.671565' },
        "TAT": {
            "location": "Tateno",
            "lat": 36.058100,
            "lon": 140.125800,
        },  # 'url': 'https://doi.pangaea.de/10.1594/PANGAEA.673090' },
        "BAR": {
            "location": "Barrow",
            "lat": 71.323000,
            "lon": -156.607000,
        },  # 'url': 'https://doi.pangaea.de/10.1594/PANGAEA.789103' },
        "IZA": {
            "location": "   ",
            "lat": 28.309350,
            "lon": -16.499260,
        },  # 'url': 'https://doi.pangaea.de/10.1594/PANGAEA.871032' },
        "NYA": {
            "location": "   ",
            "lat": 78.922700,
            "lon": 11.927300,
        },  # 'url': 'https://doi.pangaea.de/10.1594/PANGAEA.670762' },
    }[id]


def main():
    args = command_line_options()

    print(f"Thermofeel version: {thermofeel.__version__}")
    # print(f"Python version: {sys.version}")
    # print(f"Numpy version: {np.version.version}")
    # np.show_config()

    station = station_location(args.stationid)
    location = station["location"]
    lat = station["lat"]
    lon = station["lon"]

    print(f"stationid {args.stationid} in {location} @ lat {lat} lon {lon}")

    obsdata = {}

    read_part_dataset(
        obsdata, "rh", args.basepath + args.stationid + "_rh.tab", args.max
    )
    read_part_dataset(
        obsdata, "t2m", args.basepath + args.stationid + "_t2m.tab", args.max
    )
    read_part_dataset(
        obsdata, "wind", args.basepath + args.stationid + "_wind.tab", args.max
    )
    read_part_dataset(
        obsdata, "dif", args.basepath + args.stationid + "_dif.tab", args.max
    )
    read_part_dataset(
        obsdata, "dir", args.basepath + args.stationid + "_dir.tab", args.max
    )
    read_part_dataset(
        obsdata, "lwd", args.basepath + args.stationid + "_lwd.tab", args.max
    )
    read_part_dataset(
        obsdata, "lwu", args.basepath + args.stationid + "_lwu.tab", args.max
    )
    read_part_dataset(
        obsdata, "swd", args.basepath + args.stationid + "_swd.tab", args.max
    )
    read_part_dataset(
        obsdata, "swu", args.basepath + args.stationid + "_swu.tab", args.max
    )

    print(f"total entries {len(obsdata.keys())}")

    completes = []
    for k in obsdata.keys():
        entry = obsdata[k]

        if (
            "rh" in entry
            and "t2m" in entry
            and "wind" in entry
            and "dif" in entry
            and "dir" in entry
            and "lwd" in entry
            and "lwu" in entry
            and "swd" in entry
            and "swu" in entry
        ):
            entry.update(station)
            entry.update({"stationid": args.stationid})
            entry.update({"datetime": k})

            dt = dateutil.parser.isoparse(k)
            # print(k,  entry)
            # print(dt)

            cossza = thermofeel.calculate_cos_solar_zenith_angle(
                lat=lat, lon=lon, y=dt.year, m=dt.month, d=dt.day, h=dt.hour
            )

            rh = float(entry["rh"])  # percent
            t2m = thermofeel.celsius_to_kelvin(float(entry["t2m"]))  # K
            va = float(entry["wind"])  # m/s
            dir = float(entry["dir"])  # m/s
            swd = float(entry["swd"])  # W/m**2
            swu = float(entry["swu"])  # W/m**2
            lwd = float(entry["lwd"])  # W/m**2
            lwu = float(entry["lwu"])  # W/m**2

            rh = np.array([rh])
            t2m = np.array([t2m])
            va = np.array([va])
            dir = np.array([dir])
            swd = np.array([swd])
            swu = np.array([swu])
            lwd = np.array([lwd])
            lwu = np.array([lwu])

            ssrd = swd
            ssr = ssrd - swu
            dsrp = dir
            strd = lwd
            strr = lwd - lwu
            fdir = dsrp * cossza

            td = thermofeel.calculate_dew_point_from_relative_humidity(rh, t2m)
            mrt = thermofeel.calculate_mean_radiant_temperature(
                ssrd, ssr, dsrp, strd, fdir, strr, cossza
            )
            wbgt = thermofeel.calculate_wbgt(t2m, mrt, va, td)

            # print(f"cossza={cossza} mrt {mrt} mrt_C {thermofeel.kelvin_to_celsius(mrt)} wbgt {wbgt}")

            entry.update(
                {"cossza": cossza, "td": td[0], "mrt": mrt[0], "wbgt": wbgt[0]}
            )

            completes.append(entry)

    print(f"complete entries {len(completes)}")

    if len(completes):
        try:
            outfile = args.output
            columns = [
                "stationid",
                "location",
                "lat",
                "lon",
                "datetime",
                "rh",
                "t2m",
                "wind",
                "dif",
                "dir",
                "lwd",
                "lwu",
                "swd",
                "swu",
                "cossza",
                "td",
                "mrt",
                "wbgt",
            ]
            print(f"creating output CSV {outfile} with columns {columns}")
            with open(outfile, "w") as f:
                writer = csv.DictWriter(f, fieldnames=columns)
                writer.writeheader()
                for e in completes:
                    writer.writerow(e)
        except IOError:
            print("I/O error")

    else:
        print("No complete set of observations to output")


if __name__ == "__main__":
    sys.exit(main())

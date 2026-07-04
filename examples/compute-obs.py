# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Compute thermofeel indices from near-real-time surface station observations.

This example fetches live weather-station reports from the internet and computes
the thermofeel thermal-comfort indices that need only temperature, humidity and
wind:

    Heat Index (adjusted), Humidex, Apparent Temperature, Wind Chill,
    Discomfort Index (Thom), Summer Simmer Index, Relative Strain Index.

The radiation-based indices (mean radiant temperature, UTCI, WBGT) are
intentionally NOT computed here: they require downwelling/upwelling short- and
long-wave radiation, which the near-real-time station sources below do not
provide.

Sources
-------
metar (default)
    NOAA Aviation Weather METAR API - global airport observations, no
    authentication required. https://aviationweather.gov/data/api/
meteostat (fallback)
    The Meteostat Python package (optional; ``pip install meteostat``). Imported
    lazily so the METAR path works even when Meteostat is not installed.

This script is intended only as example usage of the thermofeel library.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone

import numpy as np

import thermofeel

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

# Curated global stations. Values are (city / airport, latitude, longitude).
# Used for --list and to supply coordinates to the Meteostat source.
DEFAULT_STATIONS = {
    "EGLL": ("London Heathrow", 51.477, -0.461),
    "KJFK": ("New York JFK", 40.639, -73.764),
    "LFPG": ("Paris Charles de Gaulle", 49.015, 2.534),
    "RJTT": ("Tokyo Haneda", 35.553, 139.781),
    "YSSY": ("Sydney Kingsford Smith", -33.946, 151.177),
    "OMDB": ("Dubai Intl", 25.253, 55.365),
    "FACT": ("Cape Town Intl", -33.970, 18.602),
    "SBGR": ("Sao Paulo Guarulhos", -23.432, -46.470),
}

# A single default keeps a bare invocation cheap; pass --station to widen it and
# --list to see the curated set.
DEFAULT_STATION = "EGLL"
DEFAULT_HOURS = 6

METAR_URL = "https://aviationweather.gov/api/data/metar"
STATIONINFO_URL = "https://aviationweather.gov/api/data/stationinfo"
USER_AGENT = f"thermofeel/{thermofeel.__version__} (examples/compute-obs.py)"

KNOTS_TO_MS = 0.514444  # 1 knot in m/s
KMH_TO_MS = 1.0 / 3.6  # 1 km/h in m/s

SOURCE_LABELS = {
    "metar": "NOAA Aviation Weather METAR API",
    "meteostat": "Meteostat",
}

# Short reason shown when an index cannot be computed for lack of input.
MISSING_REASON = {
    "Heat Index (adjusted)": "needs dew point or humidity",
    "Humidex": "needs dew point or humidity",
    "Apparent Temperature": "needs wind and humidity",
    "Wind Chill": "needs wind",
    "Discomfort Index": "needs humidity",
    "Summer Simmer Index": "needs humidity",
    "Relative Strain Index": "needs humidity",
}

# Mapping from printed index label to CSV column name.
INDEX_CSV_KEYS = {
    "Heat Index (adjusted)": "heat_index_c",
    "Humidex": "humidex_c",
    "Apparent Temperature": "apparent_temp_c",
    "Wind Chill": "wind_chill_c",
    "Discomfort Index": "discomfort_index_c",
    "Summer Simmer Index": "summer_simmer_c",
    "Relative Strain Index": "relative_strain",
}

CSV_FIELDS = [
    "source",
    "icao",
    "name",
    "time_utc",
    "lat",
    "lon",
    "t2_c",
    "td_c",
    "rh",
    "wind_ms",
    "pressure_hpa",
    *INDEX_CSV_KEYS.values(),
]

NOTES = (
    "Notes:\n"
    "  - Indices use temperature, humidity and wind only; the radiation-based\n"
    "    indices (MRT, UTCI, WBGT) need radiation observations that these\n"
    "    near-real-time sources do not provide.\n"
    "  - Each index is meaningful only within its own validity range (e.g. Wind\n"
    "    Chill for cold, windy air; Humidex / Heat Index for warm, humid air).\n"
    "    Values outside those ranges are shown for information only.\n"
    "  - 'derived' humidity is converted between dew point and relative humidity\n"
    "    with thermofeel when the source reports only one of the two."
)

EPILOG = """\
examples:
  compute-obs.py --list                          # curated stations (offline)
  compute-obs.py --station EGLL,KJFK,RJTT        # latest METARs + indices
  compute-obs.py --station-info EGLL             # NOAA station metadata
  compute-obs.py --source meteostat --station EGLL --output out.csv
"""


# --------------------------------------------------------------------------- #
# Observation container
# --------------------------------------------------------------------------- #


@dataclass
class Observation:
    """One parsed near-real-time observation, normalised to SI-friendly units."""

    icao: str
    name: str
    time_utc: datetime | None
    t2_c: float | None  # air temperature [degC]
    td_c: float | None  # dew point [degC]
    rh: float | None  # relative humidity [%]
    va: float | None  # wind speed [m/s]
    pressure_hpa: float | None  # pressure [hPa]
    source: str
    lat: float | None
    lon: float | None


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #


def _num(value) -> float | None:
    """Coerce a value to float, mapping None/NaN/unparseable to None."""
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(result):
        return None
    return result


def _scalar(array) -> float:
    """Return the first element of a (1-element) array as a Python float."""
    return float(np.asarray(array, dtype=float).reshape(-1)[0])


def _fmt(value: float | None, nd: int = 1) -> str:
    """Format an optional float for a fixed-width table."""
    if value is None:
        return "n/a"
    return f"{value:.{nd}f}"


def _round(value: float | None, nd: int = 2) -> float | None:
    return None if value is None else round(value, nd)


def split_ids(text: str) -> list[str]:
    """Split a comma-separated list of ICAO ids into upper-case tokens."""
    return [part.strip().upper() for part in text.split(",") if part.strip()]


def fmt_time(dt: datetime | None) -> str:
    if dt is None:
        return "unknown"
    return dt.strftime("%Y-%m-%d %H:%M UTC")


# --------------------------------------------------------------------------- #
# Index computation
# --------------------------------------------------------------------------- #


def evaluate_observation(obs: Observation) -> dict:
    """Compute the non-radiation indices for one observation.

    Scalars are wrapped in ``np.array([...])`` because the thermofeel functions
    are vectorised over NumPy. Returns a dict with the (possibly derived) inputs
    and an ordered ``label -> (value, unit)`` mapping; a value of ``None`` marks
    an index whose required input was missing.
    """
    result = {
        "t2_c": obs.t2_c,
        "td_c": None,
        "td_derived": False,
        "rh": None,
        "rh_derived": False,
        "va": obs.va,
        "pressure_hpa": obs.pressure_hpa,
        "indices": OrderedDict(),
    }
    if obs.t2_c is None:
        return result

    t2_k = thermofeel.celsius_to_kelvin(np.array([obs.t2_c], dtype=float))

    td_k = None
    if obs.td_c is not None:
        td_k = thermofeel.celsius_to_kelvin(np.array([obs.td_c], dtype=float))

    rh = None
    if obs.rh is not None:
        rh = np.array([obs.rh], dtype=float)

    # Bridge dew point <-> relative humidity so every humidity-based index can be
    # computed whichever the source reports (METAR gives dew point; Meteostat
    # usually gives relative humidity).
    if rh is None and td_k is not None:
        rh = thermofeel.calculate_relative_humidity_percent(t2_k, td_k)
        result["rh_derived"] = True
    if td_k is None and rh is not None:
        td_k = thermofeel.calculate_dew_point_from_relative_humidity(rh, t2_k)
        result["td_derived"] = True

    va = None
    if obs.va is not None:
        va = np.array([obs.va], dtype=float)

    if td_k is not None:
        result["td_c"] = _scalar(thermofeel.kelvin_to_celsius(td_k))
    if rh is not None:
        result["rh"] = _scalar(rh)

    def to_c(kelvin_array) -> float:
        return _scalar(thermofeel.kelvin_to_celsius(kelvin_array))

    idx = result["indices"]

    # Heat Index (adjusted) and Humidex need temperature + dew point.
    if td_k is not None:
        idx["Heat Index (adjusted)"] = (
            to_c(thermofeel.calculate_heat_index_adjusted(t2_k, td_k)),
            "degC",
        )
        idx["Humidex"] = (to_c(thermofeel.calculate_humidex(t2_k, td_k)), "degC")
    else:
        idx["Heat Index (adjusted)"] = (None, "degC")
        idx["Humidex"] = (None, "degC")

    # Apparent Temperature needs temperature + wind + humidity.
    if va is not None and rh is not None:
        idx["Apparent Temperature"] = (
            to_c(thermofeel.calculate_apparent_temperature(t2_k, va, rh)),
            "degC",
        )
    else:
        idx["Apparent Temperature"] = (None, "degC")

    # Wind Chill needs temperature + wind.
    if va is not None:
        idx["Wind Chill"] = (to_c(thermofeel.calculate_wind_chill(t2_k, va)), "degC")
    else:
        idx["Wind Chill"] = (None, "degC")

    # Discomfort, Summer Simmer and Relative Strain need temperature + humidity.
    if rh is not None:
        idx["Discomfort Index"] = (
            to_c(thermofeel.calculate_discomfort_index(t2_k, rh)),
            "degC",
        )
        idx["Summer Simmer Index"] = (
            to_c(thermofeel.calculate_summer_simmer_index(t2_k, rh)),
            "degC",
        )
        idx["Relative Strain Index"] = (
            _scalar(thermofeel.calculate_relative_strain_index(t2_k, rh)),
            "",
        )
    else:
        idx["Discomfort Index"] = (None, "degC")
        idx["Summer Simmer Index"] = (None, "degC")
        idx["Relative Strain Index"] = (None, "")

    return result


# --------------------------------------------------------------------------- #
# METAR source (NOAA Aviation Weather, no auth)
# --------------------------------------------------------------------------- #


def fetch_metar_reports(icaos: list[str], hours: int) -> dict[str, dict]:
    """Fetch METAR JSON for the given ICAO ids and keep the newest per station.

    A single request covers every requested station (kind to the ~100/min rate
    limit). Returns a dict ``ICAO -> newest raw report`` (max ``obsTime``).
    """
    import requests

    params = {"ids": ",".join(icaos), "format": "json", "hours": hours}
    response = requests.get(
        METAR_URL, params=params, headers={"User-Agent": USER_AGENT}, timeout=30
    )
    response.raise_for_status()
    if response.status_code == 204 or not response.text.strip():
        return {}

    newest: dict[str, dict] = {}
    for report in response.json():
        icao = report.get("icaoId")
        if icao is None:
            continue
        obs_time = report.get("obsTime") or 0
        current = newest.get(icao)
        if current is None or obs_time > (current.get("obsTime") or 0):
            newest[icao] = report
    return newest


def observation_from_metar(icao: str, report: dict) -> Observation:
    """Convert a raw METAR report into a normalised Observation."""
    obs_time = report.get("obsTime")
    time_utc = (
        datetime.fromtimestamp(obs_time, tz=timezone.utc)
        if obs_time is not None
        else None
    )

    wspd_knots = _num(report.get("wspd"))
    va = wspd_knots * KNOTS_TO_MS if wspd_knots is not None else None

    # Prefer sea-level pressure when reported, otherwise the altimeter setting.
    pressure = _num(report.get("slp"))
    if pressure is None:
        pressure = _num(report.get("altim"))

    lat = _num(report.get("lat"))
    lon = _num(report.get("lon"))
    if (lat is None or lon is None) and icao in DEFAULT_STATIONS:
        _, lat, lon = DEFAULT_STATIONS[icao]

    return Observation(
        icao=icao,
        name=report.get("name") or icao,
        time_utc=time_utc,
        t2_c=_num(report.get("temp")),
        td_c=_num(report.get("dewp")),
        rh=None,  # computed from temperature + dew point by thermofeel
        va=va,
        pressure_hpa=pressure,
        source="metar",
        lat=lat,
        lon=lon,
    )


def collect_metar(icaos: list[str], hours: int) -> list[Observation]:
    try:
        newest = fetch_metar_reports(icaos, hours)
    except ImportError:
        print("  requests is not installed - install it with 'pip install requests'")
        return []
    except Exception as exc:  # network / HTTP / JSON errors
        print(f"  METAR request failed: {exc}")
        return []

    observations = []
    for icao in icaos:
        report = newest.get(icao)
        if report is None:
            print(f"  {icao}: no METAR in the last {hours}h - skipping")
            continue
        observations.append(observation_from_metar(icao, report))
    return observations


def print_station_info(icaos: list[str]) -> None:
    """Print NOAA station metadata for the given ICAO ids."""
    try:
        import requests
    except ImportError:
        print("requests is not installed - install it with 'pip install requests'")
        return

    params = {"ids": ",".join(icaos), "format": "json"}
    try:
        response = requests.get(
            STATIONINFO_URL,
            params=params,
            headers={"User-Agent": USER_AGENT},
            timeout=30,
        )
        response.raise_for_status()
    except Exception as exc:
        print(f"station info request failed: {exc}")
        return

    if response.status_code == 204 or not response.text.strip():
        print("no station info returned")
        return

    stations = response.json()
    if not stations:
        print("no station info returned")
        return

    header = (
        f"  {'ICAO':<6} {'Site':<30} {'Ctry':<5} {'lat':>9} {'lon':>10} {'elev':>6}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))
    for station in stations:
        print(
            f"  {str(station.get('icaoId', '')):<6} "
            f"{str(station.get('site', '')):<30.30} "
            f"{str(station.get('country', '')):<5} "
            f"{_fmt(_num(station.get('lat')), 3):>9} "
            f"{_fmt(_num(station.get('lon')), 3):>10} "
            f"{_fmt(_num(station.get('elev')), 0):>6}"
        )


# --------------------------------------------------------------------------- #
# Meteostat source (optional fallback)
# --------------------------------------------------------------------------- #


def _find_column(frame, name: str):
    """Return the column label matching ``name`` (columns may be enum or str)."""
    for column in frame.columns:
        if getattr(column, "value", column) == name:
            return column
    return None


def _column_value(row, name: str):
    """Return a row value by parameter name (labels may be enum or str)."""
    for label in row.index:
        if getattr(label, "value", label) == name:
            return row[label]
    return None


def _timestamp_to_utc(stamp) -> datetime | None:
    """Convert a pandas Timestamp (tz-naive UTC) to a tz-aware datetime."""
    try:
        pydt = stamp.to_pydatetime()
    except AttributeError:
        return None
    if pydt.tzinfo is None:
        return pydt.replace(tzinfo=timezone.utc)
    return pydt.astimezone(timezone.utc)


def fetch_meteostat_observation(
    icao: str, name: str, lat: float, lon: float
) -> tuple[Observation | None, str | None]:
    """Fetch the latest hourly observation from Meteostat for a location.

    Returns ``(Observation, None)`` on success or ``(None, message)`` otherwise.
    Meteostat is imported lazily so the METAR path works without it installed.
    """
    try:
        import meteostat
    except ImportError:
        return None, (
            "meteostat is not installed - run 'pip install meteostat' to use "
            "--source meteostat"
        )

    # Meteostat's bulk providers key off numeric station ids and want
    # timezone-naive UTC datetimes. Resolve the nearest station to the curated
    # coordinates (avoids maintaining an ICAO -> Meteostat id table), then pull a
    # short recent window and keep the newest row that has a temperature.
    now = datetime.now(timezone.utc).replace(tzinfo=None)
    start = now - timedelta(days=2)
    point = meteostat.Point(lat, lon)

    try:
        nearby = meteostat.stations.nearby(point, radius=100000, limit=1)
    except Exception as exc:
        return None, f"meteostat station lookup failed: {exc}"
    if nearby is None or nearby.empty:
        return None, "meteostat found no station near the coordinates"

    station_id = nearby.index[0]
    station_name = nearby.iloc[0].get("name") or name

    parameters = [
        meteostat.Parameter.TEMP,
        meteostat.Parameter.DWPT,
        meteostat.Parameter.RHUM,
        meteostat.Parameter.WSPD,
        meteostat.Parameter.PRES,
    ]
    try:
        frame = meteostat.hourly(station_id, start, now, parameters=parameters).fetch()
    except Exception as exc:
        return None, f"meteostat data request failed: {exc}"
    if frame is None or frame.empty:
        return None, "meteostat returned no recent hourly data"

    temp_col = _find_column(frame, "temp")
    if temp_col is None:
        return None, "meteostat data has no temperature column"
    valid = frame.dropna(subset=[temp_col])
    if valid.empty:
        return None, "meteostat data has no valid temperature"

    row = valid.iloc[-1]
    wspd_kmh = _num(_column_value(row, "wspd"))
    va = wspd_kmh * KMH_TO_MS if wspd_kmh is not None else None

    obs = Observation(
        icao=icao,
        name=f"{station_name} (via Meteostat)",
        time_utc=_timestamp_to_utc(valid.index[-1]),
        t2_c=_num(_column_value(row, "temp")),
        td_c=_num(_column_value(row, "dwpt")),
        rh=_num(_column_value(row, "rhum")),
        va=va,
        pressure_hpa=_num(_column_value(row, "pres")),
        source="meteostat",
        lat=lat,
        lon=lon,
    )
    return obs, None


def collect_meteostat(icaos: list[str]) -> list[Observation]:
    observations = []
    for icao in icaos:
        curated = DEFAULT_STATIONS.get(icao)
        if curated is None:
            print(
                f"  {icao}: not in the curated list - --source meteostat needs "
                "coordinates; use --source metar or a curated station (--list)"
            )
            continue
        name, lat, lon = curated
        obs, error = fetch_meteostat_observation(icao, name, lat, lon)
        if error is not None:
            print(f"  {icao}: {error}")
            continue
        observations.append(obs)
    return observations


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #


def list_stations() -> None:
    print("Curated default stations:")
    print(f"  {'ICAO':<6} {'City / airport':<26} {'lat':>9} {'lon':>10}")
    print(f"  {'-' * 6:<6} {'-' * 26:<26} {'-' * 9:>9} {'-' * 10:>10}")
    for icao, (city, lat, lon) in DEFAULT_STATIONS.items():
        print(f"  {icao:<6} {city:<26} {lat:>9.3f} {lon:>10.3f}")


def print_observation(obs: Observation, ev: dict) -> None:
    width = 68
    print("=" * width)
    print(f"{obs.icao}  {obs.name}")
    print("-" * width)
    print(f"  source          : {SOURCE_LABELS.get(obs.source, obs.source)}")
    print(f"  observation time: {fmt_time(obs.time_utc)}")
    if obs.lat is not None and obs.lon is not None:
        print(f"  location        : lat {obs.lat:.3f}, lon {obs.lon:.3f}")
    print("-" * width)

    td_note = "  (derived from RH)" if ev["td_derived"] else ""
    rh_note = "  (derived from Td)" if ev["rh_derived"] else ""
    print("  observed inputs")
    print(f"    air temperature   T   : {_fmt(ev['t2_c']):>7} degC")
    print(f"    dew point         Td  : {_fmt(ev['td_c']):>7} degC{td_note}")
    print(f"    relative humidity RH  : {_fmt(ev['rh']):>7} %{rh_note}")
    print(f"    wind speed        va  : {_fmt(ev['va']):>7} m/s")
    print(f"    pressure          p   : {_fmt(ev['pressure_hpa']):>7} hPa")
    print("-" * width)

    print("  thermal comfort indices")
    for label, (value, unit) in ev["indices"].items():
        if value is None:
            print(f"    {label:<24}: {'n/a':>8}  ({MISSING_REASON[label]})")
        elif unit:
            print(f"    {label:<24}: {value:8.1f} {unit}")
        else:
            print(f"    {label:<24}: {value:8.2f}  (dimensionless)")
    print("=" * width)


def csv_row(obs: Observation, ev: dict) -> dict:
    row = {
        "source": obs.source,
        "icao": obs.icao,
        "name": obs.name,
        "time_utc": fmt_time(obs.time_utc),
        "lat": obs.lat,
        "lon": obs.lon,
        "t2_c": _round(ev["t2_c"]),
        "td_c": _round(ev["td_c"]),
        "rh": _round(ev["rh"]),
        "wind_ms": _round(ev["va"]),
        "pressure_hpa": _round(ev["pressure_hpa"]),
    }
    for label, (value, _unit) in ev["indices"].items():
        row[INDEX_CSV_KEYS[label]] = _round(value)
    return row


def write_csv(path: str, rows: list[dict]) -> None:
    try:
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
    except OSError as exc:
        print(f"could not write CSV {path}: {exc}")
        return
    print(f"\nWrote {len(rows)} row(s) to {path}")


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def collect_observations(
    source: str, icaos: list[str], hours: int
) -> list[Observation]:
    if source == "metar":
        return collect_metar(icaos, hours)
    return collect_meteostat(icaos)


def command_line_options(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute thermofeel non-radiation thermal-comfort indices from "
            "near-real-time surface station observations."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EPILOG,
    )
    parser.add_argument(
        "--source",
        choices=["metar", "meteostat"],
        default="metar",
        help="observation source (default: metar - NOAA Aviation Weather, no auth)",
    )
    parser.add_argument(
        "--station",
        default=DEFAULT_STATION,
        help=(
            f"comma-separated ICAO id(s) to query (default: {DEFAULT_STATION}). "
            "Any ICAO works with --source metar; --source meteostat needs a "
            "curated station (see --list)"
        ),
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="print the curated default stations and exit (works offline)",
    )
    parser.add_argument(
        "--station-info",
        metavar="ICAO[,ICAO...]",
        help="fetch NOAA station metadata for the given ICAO id(s) and exit",
    )
    parser.add_argument(
        "--output",
        metavar="PATH.csv",
        help="also write the computed results to a CSV file",
    )
    parser.add_argument(
        "--hours",
        type=int,
        default=DEFAULT_HOURS,
        help=f"METAR look-back window in hours (default: {DEFAULT_HOURS})",
    )
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = command_line_options(argv)

    print(f"thermofeel version: {thermofeel.__version__}")

    if args.list:
        list_stations()
        return 0

    if args.station_info:
        print_station_info(split_ids(args.station_info))
        return 0

    icaos = split_ids(args.station)
    if not icaos:
        print("no station ids given (use --station ICAO[,ICAO...] or --list)")
        return 1

    print(f"source: {SOURCE_LABELS[args.source]}")
    print(f"stations: {', '.join(icaos)}\n")

    observations = collect_observations(args.source, icaos, args.hours)
    if not observations:
        print("\nNo observations could be retrieved.")
        return 1

    rows = []
    print()
    for obs in observations:
        ev = evaluate_observation(obs)
        print_observation(obs, ev)
        rows.append(csv_row(obs, ev))

    print()
    print(NOTES)

    if args.output:
        write_csv(args.output, rows)

    return 0


if __name__ == "__main__":
    sys.exit(main())

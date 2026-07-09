# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Access to NOAA SURFRAD surface-radiation observations.

SURFRAD (Surface Radiation Budget Network, NOAA GML,
https://gml.noaa.gov/grad/surfrad/) provides 1-minute (3-minute before 2009)
surface irradiance and meteorological observations at seven US stations, with
BSRN-style quality control. Daily ASCII files are public, no authentication:

    https://gml.noaa.gov/aftp/data/radiation/surfrad/{Station_Dir}/{YYYY}/
        {sta}{yy}{jjj}.dat

Each data row (after two header lines) holds date/time fields, the solar
zenith angle, and 20 (value, qc) pairs; qc == 0 means the value passed all
checks, and missing values are -9999.9. Format reference:
https://gml.noaa.gov/aftp/data/radiation/surfrad/README (and the realtime
README). Radiation is in W m-2, temperature degC, RH %, wind m/s, pressure mb.
"""

import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd

BASE_URL = "https://gml.noaa.gov/aftp/data/radiation/surfrad"
CACHE_DIR = Path(__file__).resolve().parents[1] / "_cache" / "surfrad"

#: station code -> (server directory, latitude, longitude, elevation [m])
STATIONS = {
    "bon": ("Bondville_IL", 40.05192, -88.37309, 230.0),
    "tbl": ("Boulder_CO", 40.12498, -105.23680, 1689.0),
    "dra": ("Desert_Rock_NV", 36.62373, -116.01947, 1007.0),
    "fpk": ("Fort_Peck_MT", 48.30783, -105.10170, 634.0),
    "gwn": ("Goodwin_Creek_MS", 34.2547, -89.8729, 98.0),
    "psu": ("Penn_State_PA", 40.72012, -77.93085, 376.0),
    "sxf": ("Sioux_Falls_SD", 43.73403, -96.62328, 473.0),
}

# The 20 (value, qc) pairs that follow the 8 date/time+zenith columns, in file
# order (SURFRAD daily-file format).
_PAIR_FIELDS = [
    "dw_solar",  # global downwelling solar (GHI) [W m-2]
    "uw_solar",
    "direct_n",  # direct-normal solar (pyrheliometer) [W m-2]
    "diffuse",  # downwelling diffuse solar [W m-2]
    "dw_ir",
    "dw_casetemp",
    "dw_dometemp",
    "uw_ir",
    "uw_casetemp",
    "uw_dometemp",
    "uvb",
    "par",
    "netsolar",  # net solar (dw_solar - uw_solar) [W m-2]
    "netir",  # net infrared [W m-2]
    "totalnet",  # net radiation (netsolar + netir) [W m-2]
    "temp",  # 10 m air temperature [degC]
    "rh",  # relative humidity [%]
    "windspd",  # 10 m wind speed [m s-1]
    "winddir",
    "pressure",  # station pressure [mb]
]
_MISSING = -9999.9


def _fetch(station: str, day: dt.date, timeout: float = 60.0) -> Path | None:
    """Download (or reuse from cache) one SURFRAD daily file.

    Returns the local path, or None if the file is unavailable (HTTP error).
    """
    import requests

    directory, *_ = STATIONS[station]
    name = f"{station}{day:%y}{day.timetuple().tm_yday:03d}.dat"
    local = CACHE_DIR / directory / str(day.year) / name
    if local.exists() and local.stat().st_size > 0:
        return local
    url = f"{BASE_URL}/{directory}/{day.year}/{name}"
    resp = requests.get(url, timeout=timeout)
    if resp.status_code != 200 or len(resp.content) < 1000:
        print(f"  !! unavailable: {url} (HTTP {resp.status_code})")
        return None
    local.parent.mkdir(parents=True, exist_ok=True)
    local.write_bytes(resp.content)
    return local


def _parse(path: Path) -> pd.DataFrame:
    """Parse one SURFRAD daily file into a QC-aware DataFrame."""
    raw = np.loadtxt(path, skiprows=2)
    if raw.ndim != 2 or raw.shape[1] != 8 + 2 * len(_PAIR_FIELDS):
        raise ValueError(f"{path}: unexpected column count {raw.shape}")
    out = {
        "time": pd.to_datetime(
            {
                "year": raw[:, 0].astype(int),
                "month": raw[:, 2].astype(int),
                "day": raw[:, 3].astype(int),
                "hour": raw[:, 4].astype(int),
                "minute": raw[:, 5].astype(int),
            },
            utc=True,
        ),
        "zen": raw[:, 7],
    }
    for k, field in enumerate(_PAIR_FIELDS):
        value = raw[:, 8 + 2 * k]
        qc = raw[:, 9 + 2 * k]
        value = np.where((value <= _MISSING) | (qc != 0), np.nan, value)
        out[field] = value
    return pd.DataFrame(out)


def load(station: str, days: list[dt.date]) -> pd.DataFrame:
    """Load SURFRAD observations for one station over the given days.

    Files are downloaded to the campaign cache on first use. Only QC-passing
    values are kept (failed QC -> NaN). Adds station metadata columns.
    """
    frames = []
    misses = 0
    for day in days:
        path = _fetch(station, day)
        if path is None:
            misses += 1
            continue
        frames.append(_parse(path))
    if not frames:
        raise RuntimeError(f"no SURFRAD data for {station}")
    df = pd.concat(frames, ignore_index=True)
    directory, lat, lon, elev = STATIONS[station]
    df["station"] = station
    df["lat"], df["lon"], df["elev"] = lat, lon, elev
    print(
        f"  {station} ({directory}): {len(frames)} days, {len(df)} records"
        + (f", {misses} missing days" if misses else "")
    )
    return df


def month_days(year: int, month: int) -> list[dt.date]:
    """All calendar days of one month."""
    first = dt.date(year, month, 1)
    nxt = dt.date(year + month // 12, month % 12 + 1, 1)
    return [first + dt.timedelta(days=i) for i in range((nxt - first).days)]

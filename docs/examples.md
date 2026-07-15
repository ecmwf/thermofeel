# Examples

Runnable examples live in the
[`examples/`](https://github.com/ecmwf/thermofeel/tree/main/examples) directory
of the repository:

- `compute-heat-force.py` — Liljegren WBGT and the KNMI 0–10 heat-force scale
  from hard-coded arrays (**self-contained**, no extra dependencies).
- `compute-thermal-indices.py` — compute indices on **gridded ECMWF forecast
  fields** fetched online with [earthkit-data](https://github.com/ecmwf/earthkit-data),
  and write the results to GRIB and/or NetCDF.
- `compute-obs.py` — compute indices from **near-real-time station observations**
  (global airport METARs, no account needed).
- `thermofeel-earthkit-demo.ipynb` — a notebook that fetches open data, computes
  indices and draws maps with
  [earthkit-plots](https://github.com/ecmwf/earthkit-plots).

Except for `compute-heat-force.py`, the examples use the ECMWF **earthkit** stack
and a few other libraries. Install them with the `examples` extra:

```bash
pip install "thermofeel[examples]"
```

The two authenticated online sources of `compute-thermal-indices.py`
(`--source polytope` and `--source mars`) additionally need
`pip install polytope-client ecmwf-api-client` and ECMWF credentials.

## Heat force from standard variables (self-contained)

The functions are vectorised over NumPy, so the same call works on a single
point or on a full gridded field.

```python
import numpy as np

import thermofeel as thermofeel

t2_k = np.array([298.15, 303.15, 308.15])  # 2 m air temperature [K]
rh = np.array([60.0, 70.0, 40.0])  # relative humidity [%]
pressure = np.array([1013.0, 1010.0, 1013.0])  # surface pressure [hPa]
va = np.array([3.0, 2.0, 1.0])  # 10 m wind speed [m/s]
ssrd = np.array([400.0, 600.0, 800.0])  # instantaneous SSRD [W/m2]
fdir = np.array([0.5, 0.6, 0.7])  # direct-beam fraction [0-1]
cossza = np.array([0.6, 0.8, 0.9])  # cosine of solar zenith angle

wbgt_k = thermofeel.calculate_wbgt_liljegren(
    t2_k, rh, pressure, va, ssrd, fdir, cossza
)
heat_force = thermofeel.calculate_heat_force(wbgt_k)
```

The cosine of the solar zenith angle can be obtained from
[earthkit-meteo](https://github.com/ecmwf/earthkit-meteo):

```python
from earthkit.meteo import solar

cossza = solar.cos_solar_zenith_angle(...)
```

## Gridded indices from ECMWF forecasts

`compute-thermal-indices.py` reads ECMWF forecast fields with earthkit-data,
computes a selection of indices with thermofeel, and writes them out. It offers
four data **sources**:

```console
$ python examples/compute-thermal-indices.py --source opendata --step 6 \
      --output indices.nc --output-format netcdf
```

- `--source opendata` (default, **no account**) — the free
  [ECMWF open data](https://www.ecmwf.int/en/forecasts/datasets/open-data).
  Open data does **not** include direct solar radiation (`fdir`), so the
  radiation-based indices (mean radiant temperature, UTCI, WBGT, PMV) are
  skipped with a warning; the temperature/humidity/wind indices are computed in
  full. Pass `--approximate-fdir[=erbs|disc]` to estimate `fdir` from global
  radiation — via the Erbs (1982) decomposition (default) or the DISC (Maxwell
  1987) model, from
  [`thermofeel.approximations`](guide/approximations.md) — and obtain **approximate**
  MRT/UTCI/WBGT (a demonstration only, not validation-grade).
- `--source file --input FORECAST.grib` — a local GRIB that already contains the
  full field set (including `fdir`); computes every index exactly.
- `--source polytope` and `--source mars` — near-real-time / archived ECMWF
  operational forecasts (with `fdir`) via
  [Polytope](https://github.com/ecmwf/polytope-client) or the MARS Web API.
  These require ECMWF credentials.

Output goes to GRIB and/or NetCDF (`--output-format grib|netcdf|both`). Indices
that have official ECMWF GRIB parameters are encoded with them; the newer
indices (discomfort, summer simmer, relative strain, radiation apparent
temperature, PMV/PPD) have no registered code and are written with
**experimental local** GRIB parameters, accompanied by a warning.

## Live station observations

`compute-obs.py` pulls the latest observations for airport weather stations and
computes the temperature/humidity/wind indices (station feeds carry no
radiation, so MRT/UTCI/WBGT are out of scope here).

```console
$ python examples/compute-obs.py --list                 # curated global stations
$ python examples/compute-obs.py --station EGLL,KJFK,RJTT
$ python examples/compute-obs.py --source meteostat --station EGLL
```

The default source is the free, no-account
[NOAA Aviation Weather METAR API](https://aviationweather.gov/data/api/); any
ICAO station identifier works, and [Meteostat](https://meteostat.net/) is
available as a fallback (`--source meteostat`).

## Maps notebook

`thermofeel-earthkit-demo.ipynb` ties it together: fetch open data with
earthkit-data, compute indices with thermofeel, and plot them on maps with
earthkit-plots (plus an optional station overlay from the METAR feed).

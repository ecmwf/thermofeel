# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Compute thermofeel thermal-comfort indices from ECMWF forecast fields.

This example modernises the original raw-``eccodes`` script to the **earthkit
1.0** stack (``earthkit-data`` + ``earthkit-meteo``). It fetches ECMWF surface
forecast fields from one of four selectable sources, computes thermofeel
indices with vectorised NumPy, and writes the results to NetCDF and/or GRIB via
``earthkit-data``.

Sources (``--source``)
----------------------
opendata  ECMWF open data (no authentication, ``ecmwf-open-data``). The latest
          run is auto-selected. Open data does **not** include direct solar
          radiation (``fdir``), so the radiation indices (MRT / UTCI / WBGT /
          BGT / PMV) are skipped -- unless ``--approximate-fdir[=erbs|disc]`` is
          given, which estimates ``fdir`` from ``ssrd`` via the Erbs (1982) or
          DISC (Maxwell 1987) model from ``thermofeel.approximations`` (a
          demonstration approximation, not validation-grade).
file      A local GRIB file (``--input PATH``); the full index set is computed
          when ``fdir`` is present.
polytope  ECMWF Polytope (needs ``polytope-client`` + a Polytope key).
mars      ECMWF MARS (needs ``ecmwf-api-client`` + ``~/.ecmwfapirc``).

Unit contract (thermofeel is SI in / SI out)
--------------------------------------------
Temperatures are in kelvin, wind speed in m/s, relative humidity in %, radiation
fluxes in W/m^2. thermofeel returns kelvin for every temperature-like index.

Radiation accumulation simplification
-------------------------------------
ECMWF surface radiation fields (``ssrd``, ``ssr``, ``strd``, ``str``, ``fdir``)
are **accumulated** from the start of the forecast (J/m^2). This example
converts them to a step-mean flux (W/m^2) by dividing by ``step * 3600`` seconds
and, consistently, uses the step-mean cosine of the solar zenith angle
(``earthkit.meteo.solar.cos_solar_zenith_angle_integrated`` returns the *mean*
over the interval in earthkit-meteo 1.0). This ignores the finer IFS
accumulation-interval structure and is a deliberate simplification for the
example. Radiation indices are skipped at step 0 (the window would be zero).
"""

import argparse
import importlib.util
import os
import sys
from dataclasses import dataclass
from typing import Callable, Optional

import earthkit.data as ekd
import numpy as np
from earthkit.meteo import solar, wind

import thermofeel as tf
from thermofeel.approximations import approximate_fdir_disc, approximate_fdir_erbs

###########################################################################################################
# Input fields requested from each source.
#
# Open data exposes fields by short name and does NOT include ``fdir``.
OPENDATA_PARAMS = ["2t", "2d", "10u", "10v", "sp", "ssrd", "ssr", "strd", "str"]
# MARS / Polytope are requested by paramId; 228021 == fdir (direct solar
# radiation), which is what unlocks the radiation indices.
#   167=2t 168=2d 165=10u 166=10v 134=sp 169=ssrd 176=ssr 175=strd 177=str
MARS_PARAMS = "167/168/165/166/134/169/176/175/177/228021"

# GRIB2 local-use ranges (192-254 per WMO) used to encode the experimental
# indices that have no WMO paramId. They round-trip as an ECMWF-local parameter
# (paramId reads back as 0 / "unknown"); the nominal private paramId documented
# alongside each index (261101...) is for human reference only.
LOCAL_DISCIPLINE = 192
LOCAL_CATEGORY = 192
LOCAL_NUMBER_BASE = 191  # parameterNumber = LOCAL_NUMBER_BASE + local_number
NOMINAL_PRIVATE_PARAMID_BASE = 261100  # nominal paramId = base + local_number

###########################################################################################################


def field_stats(name, values):
    """Print min/max/mean/std/NaN-count for a field (kept from the original)."""
    values = np.asarray(values, dtype=float)
    n_nan = int(np.count_nonzero(np.isnan(values)))
    if values.size == n_nan:
        print(f"  {name:34s} all NaN ({values.size} points)")
        return
    print(
        f"  {name:34s} "
        f"min {np.nanmin(values):+10.3f}  max {np.nanmax(values):+10.3f}  "
        f"mean {np.nanmean(values):+10.3f}  "
        f"std {np.nanstd(values, dtype=np.float64):8.3f}  nan {n_nan}"
    )


def banner(lines):
    """Print a prominent boxed warning/notice block."""
    width = max(len(line) for line in lines)
    bar = "!" * (width + 6)
    print(bar)
    for line in lines:
        print(f"!! {line.ljust(width)} !!")
    print(bar)


###########################################################################################################
# Computation context: raw fields + cached derived building blocks.


class Env:
    """Holds the decoded input fields and derives (and caches) the shared
    building blocks (wind speed, relative humidity, cossza, MRT, ...)."""

    def __init__(
        self, fields, template, lat, lon, base_time, valid_time, step, approximate_fdir
    ):
        self._f = fields  # shortName -> 2D ndarray of values
        self.template = template  # a source Field, reused as GRIB clone template
        self.lat = lat
        self.lon = lon
        self.base_time = base_time
        self.valid_time = valid_time
        self.step = int(step)
        self.approximate_fdir = approximate_fdir
        self._cache = {}
        self.t2 = fields["2t"]
        self.td = fields["2d"]
        self.u10 = fields.get("10u")
        self.v10 = fields.get("10v")

    def has(self, name):
        return name in self._f

    def get(self, name):
        return self._f[name]

    def _cached(self, key, fn):
        if key not in self._cache:
            self._cache[key] = fn()
        return self._cache[key]

    # --- derived building blocks -------------------------------------------
    def ws(self):
        return self._cached("ws", lambda: wind.speed(self.u10, self.v10))

    def rh(self):
        return self._cached(
            "rh", lambda: tf.calculate_relative_humidity_percent(self.t2, self.td)
        )

    def cossza(self):
        # Step-mean cosine of the solar zenith angle over [base_time, valid_time]
        # (earthkit-meteo 1.0 returns the interval MEAN in [0, 1]).
        return self._cached(
            "cossza",
            lambda: solar.cos_solar_zenith_angle_integrated(
                self.base_time, self.valid_time, self.lat, self.lon
            ),
        )

    def mean_flux(self, name):
        """Accumulated J/m^2 -> step-mean flux W/m^2 (see module docstring)."""
        return self.get(name) / (self.step * 3600.0)

    def fdir_flux(self):
        def _f():
            if self.has("fdir"):
                return self.mean_flux("fdir")
            # Estimate fdir from ssrd via thermofeel.approximations. Guarded by
            # the caller so this only runs when --approximate-fdir was requested;
            # self.approximate_fdir holds the chosen method ("erbs" | "disc").
            ssrd = self.mean_flux("ssrd")
            doy = self.valid_time.timetuple().tm_yday
            if self.approximate_fdir == "disc":
                pressure_hpa = self.get("sp") / 100.0 if self.has("sp") else 1013.25
                return approximate_fdir_disc(
                    ssrd, self.cossza(), doy, pressure_hpa=pressure_hpa
                )
            return approximate_fdir_erbs(ssrd, self.cossza(), doy=doy)

        return self._cached("fdir_flux", _f)

    def dsrp(self):
        # Direct solar radiation perpendicular to the beam, approximated from
        # fdir and cossza (thermofeel helper); no dsrp field is requested.
        return self._cached(
            "dsrp", lambda: tf.approximate_dsrp(self.fdir_flux(), self.cossza())
        )

    def mrt(self):
        def _f():
            return tf.calculate_mean_radiant_temperature(
                self.mean_flux("ssrd"),
                self.mean_flux("ssr"),
                self.dsrp(),
                self.mean_flux("strd"),
                self.fdir_flux(),
                self.mean_flux("str"),
                self.cossza(),
            )

        return self._cached("mrt", _f)

    def pmv(self):
        return self._cached(
            "pmv",
            lambda: tf.calculate_pmv(
                self.t2, self.mrt(), self.pmv_air_velocity(), rh=self.rh()
            ),
        )

    def pmv_air_velocity(self):
        # PMV wants the relative air velocity felt at the body, NOT the 10 m
        # meteorological wind. There is no such NWP field, so the 10 m wind speed
        # is used here as a coarse outdoor demonstration proxy.
        return self.ws()

    def body_net_radiation(self):
        # Steadman's radiation apparent-temperature wants the net radiation
        # absorbed per unit body-surface area (a caller-supplied quantity). Here
        # we demonstrate with the surface net all-wave flux (ssr + str), which is
        # only a rough proxy for the body-scale term.
        return self.mean_flux("ssr") + self.mean_flux("str")


###########################################################################################################
# Index compute functions. Each takes an Env and returns a NumPy array.


def _relative_humidity(env):
    return tf.calculate_relative_humidity_percent(env.t2, env.td)


def _heat_index(env):
    return tf.calculate_heat_index_adjusted(env.t2, env.td)


def _heat_index_simplified(env):
    return tf.calculate_heat_index_simplified(env.t2, env.rh())


def _humidex(env):
    return tf.calculate_humidex(env.t2, env.td)


def _apparent_temperature(env):
    return tf.calculate_apparent_temperature(env.t2, env.ws(), env.rh())


def _wind_chill(env):
    return tf.calculate_wind_chill(env.t2, env.ws())


def _normal_effective_temperature(env):
    return tf.calculate_normal_effective_temperature(env.t2, env.ws(), env.rh())


def _discomfort_index(env):
    return tf.calculate_discomfort_index(env.t2, env.rh())


def _summer_simmer_index(env):
    return tf.calculate_summer_simmer_index(env.t2, env.rh())


def _relative_strain_index(env):
    return tf.calculate_relative_strain_index(env.t2, env.rh())


def _cossza(env):
    return env.cossza()


def _mean_radiant_temperature(env):
    return env.mrt()


def _utci(env):
    return tf.calculate_utci(t2_k=env.t2, va=env.ws(), mrt=env.mrt(), td_k=env.td)


def _wet_bulb_temperature(env):
    return tf.calculate_wbt(env.t2, env.rh())


def _globe_temperature(env):
    return tf.calculate_bgt(env.t2, env.mrt(), env.ws())


def _wbgt(env):
    return tf.calculate_wbgt(env.t2, env.mrt(), env.ws(), env.td)


def _apparent_temperature_radiation(env):
    return tf.calculate_apparent_temperature_radiation(
        env.t2, env.ws(), env.rh(), env.body_net_radiation()
    )


def _pmv(env):
    return env.pmv()


def _ppd(env):
    return tf.calculate_ppd(env.pmv())


###########################################################################################################
# Index registry.
#
# category:  "basic"      -> computable for every source (no fdir needed): the
#                            2t/2d/10u/10v indices, plus cossza (pure solar
#                            geometry) and wet-bulb temperature (2t + rh only).
#            "radiation"  -> needs the surface radiation fluxes incl. fdir (exact
#                            or Erbs/DISC-approximated); skipped without fdir / at
#                            step 0: MRT, UTCI, WBGT, BGT, PMV and the radiation
#                            apparent temperature.
# GRIB code: ``paramid`` set  -> WMO/ECMWF paramId, clean encoding.
#            ``local_number`` -> experimental local GRIB2 code (see constants).


@dataclass(frozen=True)
class IndexSpec:
    key: str  # internal id + NetCDF variable name
    long_name: str
    units: str
    category: str
    compute: Callable
    paramid: Optional[int] = None
    local_number: Optional[int] = None

    @property
    def cli(self):
        return "--" + self.key.replace("_", "-")

    @property
    def dest(self):
        return self.key

    @property
    def experimental(self):
        return self.paramid is None

    @property
    def parameter_number(self):
        return (
            None if self.local_number is None else LOCAL_NUMBER_BASE + self.local_number
        )

    @property
    def nominal_paramid(self):
        if self.local_number is None:
            return None
        return NOMINAL_PRIVATE_PARAMID_BASE + self.local_number


INDEX_SPECS = [
    # --- basic (non-radiation) indices -------------------------------------
    IndexSpec(
        "relative_humidity",
        "Relative humidity",
        "%",
        "basic",
        _relative_humidity,
        paramid=260242,
    ),
    IndexSpec(
        "heat_index", "Heat index (adjusted)", "K", "basic", _heat_index, paramid=260004
    ),
    IndexSpec(
        "heat_index_simplified",
        "Heat index (simplified)",
        "K",
        "basic",
        _heat_index_simplified,
        local_number=1,
    ),
    IndexSpec("humidex", "Humidex", "K", "basic", _humidex, paramid=261016),
    IndexSpec(
        "apparent_temperature",
        "Apparent temperature",
        "K",
        "basic",
        _apparent_temperature,
        paramid=260255,
    ),
    IndexSpec("wind_chill", "Wind chill", "K", "basic", _wind_chill, paramid=260005),
    IndexSpec(
        "normal_effective_temperature",
        "Normal effective temperature",
        "K",
        "basic",
        _normal_effective_temperature,
        paramid=261018,
    ),
    IndexSpec(
        "discomfort_index",
        "Discomfort index (Thom)",
        "K",
        "basic",
        _discomfort_index,
        local_number=2,
    ),
    IndexSpec(
        "summer_simmer_index",
        "Summer simmer index",
        "K",
        "basic",
        _summer_simmer_index,
        local_number=3,
    ),
    IndexSpec(
        "relative_strain_index",
        "Relative strain index",
        "1",
        "basic",
        _relative_strain_index,
        local_number=4,
    ),
    # cossza is pure solar geometry (lat/lon/time) and wet-bulb temperature needs
    # only 2t + relative humidity, so neither depends on fdir -- they are
    # computable for every source, hence "basic".
    IndexSpec(
        "cossza",
        "Cosine of solar zenith angle (step mean)",
        "1",
        "basic",
        _cossza,
        paramid=214001,
    ),
    IndexSpec(
        "wet_bulb_temperature",
        "Wet-bulb temperature",
        "K",
        "basic",
        _wet_bulb_temperature,
        paramid=261022,
    ),
    # --- radiation-dependent indices ---------------------------------------
    IndexSpec(
        "mean_radiant_temperature",
        "Mean radiant temperature",
        "K",
        "radiation",
        _mean_radiant_temperature,
        paramid=261002,
    ),
    IndexSpec(
        "utci",
        "Universal Thermal Climate Index",
        "K",
        "radiation",
        _utci,
        paramid=261001,
    ),
    IndexSpec(
        "globe_temperature",
        "Globe temperature",
        "K",
        "radiation",
        _globe_temperature,
        paramid=261015,
    ),
    IndexSpec(
        "wbgt", "Wet-bulb globe temperature", "K", "radiation", _wbgt, paramid=261014
    ),
    IndexSpec(
        "apparent_temperature_radiation",
        "Apparent temperature (with radiation)",
        "K",
        "radiation",
        _apparent_temperature_radiation,
        local_number=5,
    ),
    IndexSpec("pmv", "Predicted mean vote", "1", "radiation", _pmv, local_number=6),
    IndexSpec(
        "ppd",
        "Predicted percentage of dissatisfied",
        "%",
        "radiation",
        _ppd,
        local_number=7,
    ),
]

SPECS_BY_KEY = {spec.key: spec for spec in INDEX_SPECS}


###########################################################################################################
# Field decoding / selection.


def load_fields(fl, step):
    """Index a fieldlist by GRIB shortName.

    earthkit-data 1.0's ``FieldList.sel(param=...)`` / ``sel(shortName=...)`` is
    unreliable across readers, so we index robustly on the ``shortName``
    metadata. When a file carries several steps, the field matching ``step`` is
    preferred. Returns ``(values_by_name, field_objs_by_name)``.
    """
    values = {}
    field_objs = {}
    for f in fl:
        name = f.metadata("shortName")
        try:
            f_step = int(f.metadata("step"))
        except Exception:
            f_step = None
        prefer = f_step == step
        if name in field_objs and not prefer:
            continue
        values[name] = f.to_numpy()
        field_objs[name] = f
    return values, field_objs


def build_env(fl, args):
    values, field_objs = load_fields(fl, args.step)

    missing = [p for p in ("2t", "2d") if p not in values]
    if missing:
        raise SystemExit(
            f"ERROR: required field(s) {missing} not found in the source; "
            f"available: {sorted(values)}"
        )

    template = field_objs.get("2t", next(iter(field_objs.values())))
    lat, lon = template.geography.latlons()
    base_time = template.metadata("base_datetime")
    valid_time = template.metadata("valid_datetime")
    step = int(template.metadata("step"))

    print(f"\nDecoded {len(values)} field(s): {', '.join(sorted(values))}")
    print(
        f"base time {base_time}  step {step} h  valid {valid_time}  "
        f"grid {lat.shape[0]}x{lat.shape[1]}"
    )
    print("\nInput field stats:")
    for name in sorted(values):
        field_stats(name, values[name])

    return Env(
        values, template, lat, lon, base_time, valid_time, step, args.approximate_fdir
    )


###########################################################################################################
# Source fetchers.


def _require_backend(module, pip_name, source):
    if importlib.util.find_spec(module) is None:
        print(
            f"ERROR: --source {source} needs the optional '{pip_name}' backend, "
            f"which is not installed.\n"
            f"       Install it with:  pip install {pip_name}"
        )
        return False
    return True


def parse_grid(grid):
    parts = [p.strip() for p in grid.replace(",", "/").split("/")]
    try:
        if len(parts) != 2:
            raise ValueError
        float(parts[0])  # validate numeric
        float(parts[1])
    except ValueError:
        raise SystemExit(
            f"ERROR: --grid expects 'dx/dy' (e.g. 0.25/0.25), got {grid!r}"
        ) from None
    # Return the "dx/dy" STRING form (not a [dx, dy] list): both MARS and
    # Polytope accept it, whereas Polytope rejects a two-element list with equal
    # values ("Duplicate values found in list for key 'grid'"), which would
    # break the common square-grid case (e.g. the 0.25/0.25 default).
    return f"{parts[0]}/{parts[1]}"


def mars_request(args):
    return {
        "class": "od",
        "stream": "oper",
        "type": "fc",
        "levtype": "sfc",
        "date": args.date if args.date else -1,
        "time": args.time if args.time is not None else "00",
        "step": str(args.step),
        "param": MARS_PARAMS,
        "grid": parse_grid(args.grid),
    }


def fetch_opendata(args):
    request = dict(
        type="fc",
        stream="oper",
        levtype="sfc",
        step=args.step,
        param=OPENDATA_PARAMS,
    )
    if args.date:
        request["date"] = args.date
    if args.time is not None:
        request["time"] = args.time
    print(f"Fetching ECMWF open data (no authentication): {request}")
    return ekd.from_source("ecmwf-open-data", request=request).to_fieldlist()


def fetch_file(args):
    if not args.input:
        raise SystemExit("ERROR: --source file requires --input PATH")
    if not os.path.exists(args.input):
        raise SystemExit(f"ERROR: input file not found: {args.input}")
    print(f"Reading GRIB file: {args.input}")
    return ekd.from_source("file", args.input).to_fieldlist()


def fetch_polytope(args):
    if not _require_backend("polytope", "polytope-client", "polytope"):
        raise SystemExit(2)
    request = mars_request(args)
    print(f"Requesting via ECMWF Polytope (ecmwf-mars): {request}")
    try:
        src = ekd.from_source(
            "polytope",
            "ecmwf-mars",
            request=request,
            stream=False,
            address="polytope.ecmwf.int",
        )
        return src.to_fieldlist()
    except Exception as exc:
        print(
            "ERROR: Polytope request failed. This usually means a missing or "
            "invalid Polytope key.\n"
            "       Obtain a key at https://polytope.ecmwf.int and store it in "
            "~/.polytopeapirc\n"
            f"       Underlying error: {exc}"
        )
        raise SystemExit(2) from exc


def fetch_mars(args):
    if not _require_backend("ecmwfapi", "ecmwf-api-client", "mars"):
        raise SystemExit(2)
    request = mars_request(args)
    print(f"Requesting via ECMWF MARS: {request}")
    try:
        return ekd.from_source("mars", request=request).to_fieldlist()
    except Exception as exc:
        print(
            "ERROR: MARS request failed. This usually means missing credentials.\n"
            "       Create ~/.ecmwfapirc with your API key "
            "(https://api.ecmwf.int/v1/key/)\n"
            f"       Underlying error: {exc}"
        )
        raise SystemExit(2) from exc


SOURCES = {
    "opendata": fetch_opendata,
    "file": fetch_file,
    "polytope": fetch_polytope,
    "mars": fetch_mars,
}


###########################################################################################################
# Radiation availability + index selection.


def radiation_status(env):
    """Return (ok, mode, reason) describing whether the radiation indices can be
    computed for this Env."""
    if env.step == 0:
        return False, None, "step 0 (radiation window is zero)"
    fluxes = ("ssrd", "ssr", "strd", "str")
    if not all(env.has(n) for n in fluxes):
        have = [n for n in fluxes if env.has(n)]
        return False, None, f"missing radiation fluxes (have only {have})"
    # Note: wind is NOT required here. Mean radiant temperature needs no wind;
    # the wind-dependent radiation indices (UTCI/WBGT/BGT/PMV) fail individually
    # (caught in compute_indices) if 10 m wind is absent.
    if env.has("fdir"):
        return True, "exact", "fdir present"
    if env.approximate_fdir:
        method = str(env.approximate_fdir).upper()
        return True, "approx", f"fdir estimated from ssrd via {method} (demonstration)"
    return False, None, "no fdir in source (rerun with --approximate-fdir to estimate)"


def select_specs(args, rad_ok, rad_reason):
    """Resolve which index specs to compute from the CLI flags and the source's
    radiation capability."""
    chosen = [s for s in INDEX_SPECS if getattr(args, s.dest)]
    explicit = bool(chosen)
    if not explicit:
        chosen = list(INDEX_SPECS)  # no flags -> everything possible

    runnable = []
    for spec in chosen:
        if spec.category == "radiation" and not rad_ok:
            if explicit:
                print(f"  - skipping {spec.key}: radiation unavailable ({rad_reason})")
            continue
        runnable.append(spec)
    return runnable


###########################################################################################################
# Compute + output.


def compute_indices(env, specs):
    print("\nComputed index stats:")
    results = {}
    for spec in specs:
        try:
            values = np.asarray(spec.compute(env), dtype=float)
        except Exception as exc:
            print(f"  ! {spec.key}: computation failed: {exc}")
            continue
        results[spec.key] = values
        field_stats(spec.key, values)
    return results


def _is_regular_latlon(lat, lon):
    return np.allclose(lat, lat[:, :1]) and np.allclose(lon, lon[:1, :])


def write_netcdf(results, env, source, path):
    import xarray as xr

    lat, lon = env.lat, env.lon
    if _is_regular_latlon(lat, lon):
        coords = {"latitude": lat[:, 0], "longitude": lon[0, :]}
        dims = ("latitude", "longitude")
        data_vars = {
            key: (
                dims,
                vals,
                {
                    "long_name": SPECS_BY_KEY[key].long_name,
                    "units": SPECS_BY_KEY[key].units,
                },
            )
            for key, vals in results.items()
        }
    else:
        # Curvilinear / non-regular grid: store lat/lon as 2D coordinates.
        coords = {
            "latitude": (("y", "x"), lat),
            "longitude": (("y", "x"), lon),
        }
        dims = ("y", "x")
        data_vars = {
            key: (
                dims,
                vals,
                {
                    "long_name": SPECS_BY_KEY[key].long_name,
                    "units": SPECS_BY_KEY[key].units,
                },
            )
            for key, vals in results.items()
        }

    ds = xr.Dataset(data_vars, coords=coords)
    ds = ds.assign_coords(time=np.datetime64(env.valid_time))
    ds.attrs.update(
        title="thermofeel thermal-comfort indices",
        institution="ECMWF",
        source=f"thermofeel {tf.__version__} ({source})",
        forecast_reference_time=str(env.base_time),
        valid_time=str(env.valid_time),
        forecast_step_hours=env.step,
        comment=(
            "Radiation accumulations converted to step-mean flux by dividing by "
            "step*3600 s; see script docstring."
        ),
    )
    ds.to_netcdf(path)
    print(f"  wrote NetCDF: {path}  ({len(results)} variable(s))")


def write_grib(results, env, path):
    """Encode computed indices to GRIB2 by cloning a source field as template and
    overriding paramId (WMO codes) or the local GRIB2 octets (experimental
    codes) plus the data values."""
    template = env.template
    fields = []
    experimental = []
    for key, values in results.items():
        spec = SPECS_BY_KEY[key]
        vals = np.asarray(values, dtype=float)
        # Force GRIB edition 2 FIRST (dict order is preserved and applied in
        # order): MARS/Polytope operational surface fields are edition 1, on
        # which setting a GRIB2 ECMWF-local paramId raises "Concept no match".
        # Converting to edition 2 first makes both the WMO-code and local-octet
        # paths work; for already-edition-2 templates (open data) it is a no-op.
        meta = {"metadata.edition": 2}
        if spec.paramid is not None:
            meta["metadata.paramId"] = spec.paramid
        else:
            # Unknown paramIds are rejected by ecCodes' concept database, so
            # experimental indices are encoded via the GRIB2 local-use octets.
            meta["metadata.discipline"] = LOCAL_DISCIPLINE
            meta["metadata.parameterCategory"] = LOCAL_CATEGORY
            meta["metadata.parameterNumber"] = spec.parameter_number
            experimental.append(spec)
        if np.isnan(vals).any():
            meta["metadata.bitmapPresent"] = 1
        # set(values=...) keeps the raw metadata valid; the metadata override
        # then needs sync() to rebuild the raw GRIB message.
        field = template.set(values=vals).set(meta).sync()
        fields.append(field)

    ekd.FieldList.from_fields(fields).to_target("file", path)
    print(f"  wrote GRIB: {path}  ({len(fields)} message(s))")

    if experimental:
        lines = [
            "EXPERIMENTAL GRIB CODES - not registered with WMO.",
            "The following indices have no WMO paramId and were encoded using",
            "GRIB2 local-use octets (discipline/category/number in 192-254).",
            "They read back as an ECMWF-local parameter (paramId 0 / unknown):",
        ]
        for spec in experimental:
            lines.append(
                f"  {spec.key}: local disc={LOCAL_DISCIPLINE} cat={LOCAL_CATEGORY} "
                f"num={spec.parameter_number} (nominal private paramId "
                f"{spec.nominal_paramid})"
            )
        banner(lines)


def resolve_outputs(output, fmt):
    """Map --output + --output-format to concrete {format: path} targets."""
    if output is None:
        return {}
    root, ext = os.path.splitext(output)
    if fmt == "both":
        return {"grib": root + ".grib", "netcdf": root + ".nc"}
    return {fmt: output}


###########################################################################################################
# CLI.


def command_line_options(argv=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Compute thermofeel thermal-comfort indices from ECMWF forecast "
            "fields (earthkit 1.0 stack)."
        ),
        epilog=(
            "Index selection: pass one or more --<index> flags to compute a "
            "subset, or omit them all to compute every index available for the "
            "chosen source.\n\n"
            "Examples:\n"
            "  # open data (no auth): non-radiation indices to NetCDF\n"
            "  compute-thermal-indices.py --source opendata --step 6 \\\n"
            "      --output out.nc --output-format netcdf\n\n"
            "  # open data + approximated fdir (Erbs by default): adds MRT/UTCI/WBGT\n"
            "  compute-thermal-indices.py --source opendata --approximate-fdir \\\n"
            "      --output out.nc\n\n"
            "  # open data + DISC-approximated fdir\n"
            "  compute-thermal-indices.py --source opendata --approximate-fdir disc \\\n"
            "      --output out.nc\n\n"
            "  # local GRIB file with fdir: full set to GRIB + NetCDF\n"
            "  compute-thermal-indices.py --source file --input fc.grib \\\n"
            "      --output out --output-format both"
        ),
    )
    parser.add_argument(
        "--source",
        choices=list(SOURCES),
        default="opendata",
        help="forecast data source (default: opendata)",
    )
    parser.add_argument("--input", help="input GRIB file (required for --source file)")
    parser.add_argument(
        "--step", type=int, default=6, help="forecast step in hours (default: 6)"
    )
    parser.add_argument(
        "--date",
        default=None,
        help="forecast run date (YYYYMMDD or relative int; default: latest/-1)",
    )
    parser.add_argument(
        "--time",
        default=None,
        help="forecast run time, e.g. 00/06/12/18 (default: latest for opendata, "
        "00 for mars/polytope)",
    )
    parser.add_argument(
        "--grid",
        default="0.25/0.25",
        help="output grid for mars/polytope, 'dx/dy' (default: 0.25/0.25)",
    )
    parser.add_argument(
        "--approximate-fdir",
        nargs="?",
        const="erbs",
        default=None,
        choices=["erbs", "disc"],
        metavar="{erbs,disc}",
        help="estimate fdir from ssrd to enable APPROXIMATE MRT/UTCI/WBGT/... "
        "(mainly for open data, which lacks fdir). Bare flag uses 'erbs' "
        "(Erbs et al. 1982); 'disc' selects the Maxwell (1987) DISC model. "
        "Uses thermofeel.approximations; demonstration only, not "
        "validation-grade.",
    )
    parser.add_argument("--output", default=None, help="output file path")
    parser.add_argument(
        "--output-format",
        choices=["grib", "netcdf", "both"],
        default="netcdf",
        help="output format (default: netcdf)",
    )

    selection = parser.add_argument_group(
        "index selection", "compute only the given indices (default: all available)"
    )
    for spec in INDEX_SPECS:
        selection.add_argument(
            spec.cli,
            dest=spec.dest,
            action="store_true",
            help=f"compute {spec.long_name.lower()}",
        )

    return parser.parse_args(argv)


def main(argv=None):
    args = command_line_options(argv)

    print(f"thermofeel version : {tf.__version__}")
    print(f"earthkit-data      : {ekd.__version__}")
    print(f"python             : {sys.version.split()[0]}")
    print(f"numpy              : {np.version.version}")
    print("-" * 100)

    fl = SOURCES[args.source](args)
    if len(fl) == 0:
        raise SystemExit("ERROR: the source returned no fields.")

    env = build_env(fl, args)

    rad_ok, rad_mode, rad_reason = radiation_status(env)
    print("-" * 100)
    if rad_ok and rad_mode == "approx":
        model = {"erbs": "Erbs et al. (1982)", "disc": "Maxwell (1987) DISC"}[
            env.approximate_fdir
        ]
        banner(
            [
                "APPROXIMATE RADIATION INDICES",
                f"fdir is estimated from ssrd via {model}",
                "(thermofeel.approximations).",
                "MRT / UTCI / WBGT / BGT / PMV below are a DEMONSTRATION only and",
                "are NOT validation-grade. Use a source with a real fdir field",
                "(file / polytope / mars) for quantitative work.",
            ]
        )
    elif rad_ok:
        print(f"Radiation indices: ENABLED ({rad_reason}).")
    else:
        banner(
            [
                "RADIATION INDICES SKIPPED",
                f"Reason: {rad_reason}.",
                "MRT / UTCI / WBGT / BGT / PMV need direct solar radiation (fdir).",
                "For open data, rerun with --approximate-fdir to estimate it.",
            ]
        )

    specs = select_specs(args, rad_ok, rad_reason)
    if not specs:
        raise SystemExit("ERROR: no indices to compute for this selection/source.")

    results = compute_indices(env, specs)
    if not results:
        raise SystemExit("ERROR: no indices were computed successfully.")

    print("\nComputed indices:", ", ".join(results))

    targets = resolve_outputs(args.output, args.output_format)
    if not targets:
        print("\nNo --output given: results computed but not written to disk.")
    else:
        print("\nWriting output:")
        if "netcdf" in targets:
            write_netcdf(results, env, args.source, targets["netcdf"])
        if "grib" in targets:
            write_grib(results, env, targets["grib"])

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

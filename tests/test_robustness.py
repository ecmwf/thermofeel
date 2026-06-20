# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Robustness and error-handling contracts.

These tests pin the documented behaviour at the edges of the input domain:
``NaN`` propagation (a NaN input yields a NaN output, never an exception), the
``ValueError`` contracts, and the known ``NaN`` outcomes (zero-wind globe
temperature; non-convergent Liljegren iteration). See ``plans/ROBUSTNESS.md``.
"""

import numpy as np
import pytest

import thermofeel as tmf
from thermofeel.helpers import fahrenheit_to_celsius

# Representative finite inputs (Kelvin / m s-1 / % / hPa).
T = np.array([300.0])
TD = np.array([290.0])
RH = np.array([50.0])
VA = np.array([3.0])
MRT = np.array([310.0])
NAN = np.array([np.nan])

# Each entry maps a NaN temperature through one public index.
NAN_CASES = {
    "relative_humidity_percent": lambda x: tmf.calculate_relative_humidity_percent(
        x, TD
    ),
    "dew_point_from_relative_humidity": (
        lambda x: tmf.calculate_dew_point_from_relative_humidity(RH, x)
    ),
    "humidex": lambda x: tmf.calculate_humidex(x, TD),
    "wind_chill": lambda x: tmf.calculate_wind_chill(x, VA),
    "apparent_temperature": lambda x: tmf.calculate_apparent_temperature(x, VA, RH),
    "heat_index_simplified": lambda x: tmf.calculate_heat_index_simplified(x, RH),
    "heat_index_adjusted": lambda x: tmf.calculate_heat_index_adjusted(x, TD),
    "normal_effective_temperature": (
        lambda x: tmf.calculate_normal_effective_temperature(x, VA, RH)
    ),
    "wbgt_simple": lambda x: tmf.calculate_wbgt_simple(x, RH),
    "wbgt": lambda x: tmf.calculate_wbgt(x, MRT, VA, TD),
    "utci": lambda x: tmf.calculate_utci(x, VA, MRT, td_k=TD),
}


@pytest.mark.parametrize("fn", NAN_CASES.values(), ids=list(NAN_CASES))
def test_nan_temperature_propagates(fn):
    # A NaN input must yield NaN (not raise, and not a spurious finite value).
    out = np.asarray(fn(NAN), dtype=float)
    assert out.shape == (1,)
    assert np.isnan(out).all()


def test_wbgt_liljegren_nan_element_propagates():
    # A NaN element never satisfies the convergence test, so it returns NaN,
    # while finite neighbours still converge to a real value.
    t2 = np.array([298.15, np.nan])
    rh = np.array([50.0, 50.0])
    pressure = np.array([1013.0, 1013.0])
    va = np.array([3.0, 3.0])
    ssrd = np.array([500.0, 500.0])
    fdir = np.array([0.5, 0.5])
    cossza = np.array([0.5, 0.5])
    out = tmf.calculate_wbgt_liljegren(t2, rh, pressure, va, ssrd, fdir, cossza)
    assert np.isfinite(out[0])
    assert np.isnan(out[1])


def test_bgt_zero_wind_returns_nan():
    # Documented edge: the closed-form globe temperature is NaN at exactly zero
    # wind, but finite for any positive wind speed.
    with np.errstate(invalid="ignore"):
        zero = tmf.calculate_bgt(T, MRT, np.array([0.0]))
    positive = tmf.calculate_bgt(T, MRT, np.array([0.5]))
    assert np.isnan(zero).all()
    assert np.isfinite(positive).all()


def test_utci_requires_ehpa_or_td():
    with pytest.raises(ValueError):
        tmf.calculate_utci(T, VA, MRT)


def test_wbgt_liljegren_rejects_unknown_wind_scaling():
    with pytest.raises(ValueError):
        tmf.calculate_wbgt_liljegren(
            T,
            RH,
            np.array([1013.0]),
            VA,
            np.array([500.0]),
            np.array([0.5]),
            np.array([0.5]),
            wind_scaling="bogus",
        )


def test_fahrenheit_to_celsius():
    assert tmf.fahrenheit_to_celsius(32.0) == pytest.approx(0.0)
    assert fahrenheit_to_celsius(212.0) == pytest.approx(100.0)

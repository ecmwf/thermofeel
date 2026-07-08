# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Property-based tests (Hypothesis).

These assert *invariants* that must hold across the whole valid input domain,
complementing the pinned pointwise/regression tests: physical bounds and
edge behaviour of the fdir estimators, and a couple of index invariants.
"""

import numpy as np
from hypothesis import given
from hypothesis import strategies as st

import thermofeel as tf
from thermofeel import approximations as ap

_FINITE = {"allow_nan": False, "allow_infinity": False}
ssrd_st = st.floats(min_value=0.0, max_value=1500.0, **_FINITE)
cossza_st = st.floats(min_value=0.0, max_value=1.0, **_FINITE)
doy_st = st.integers(min_value=1, max_value=366)
pressure_st = st.floats(min_value=300.0, max_value=1100.0, **_FINITE)


def _arr(x):
    return np.array([x], dtype=float)


# --- fdir estimators: physical bounds ---------------------------------------


@given(ssrd=ssrd_st, cossza=cossza_st, doy=doy_st)
def test_erbs_fdir_within_bounds(ssrd, cossza, doy):
    # direct horizontal radiation is a finite, non-negative part of the global
    # horizontal radiation: 0 <= fdir <= ssrd
    fdir = ap.approximate_fdir_erbs(_arr(ssrd), _arr(cossza), doy=_arr(doy))[0]
    assert np.isfinite(fdir)
    assert 0.0 <= fdir <= ssrd + 1e-9


@given(ssrd=ssrd_st, cossza=cossza_st, doy=doy_st, pressure=pressure_st)
def test_disc_fdir_within_bounds(ssrd, cossza, doy, pressure):
    fdir = ap.approximate_fdir_disc(
        _arr(ssrd), _arr(cossza), _arr(doy), pressure_hpa=_arr(pressure)
    )[0]
    assert np.isfinite(fdir)
    assert 0.0 <= fdir <= ssrd + 1e-9


@given(
    ssrd=ssrd_st,
    cossza=st.floats(min_value=0.0, max_value=0.065, **_FINITE),
    doy=doy_st,
)
def test_fdir_night_is_zero(ssrd, cossza, doy):
    # at or below the min cos(zenith) the Sun is below the horizon -> fdir = 0
    assert ap.approximate_fdir_erbs(_arr(ssrd), _arr(cossza), doy=_arr(doy))[0] == 0.0
    assert ap.approximate_fdir_disc(_arr(ssrd), _arr(cossza), _arr(doy))[0] == 0.0


@given(doy=doy_st)
def test_earth_sun_distance_factor_bounds(doy):
    # the eccentricity factor stays within the physical aphelion/perihelion band
    f = ap.earth_sun_distance_factor(_arr(doy))[0]
    assert 0.96 <= f <= 1.04


# --- index invariants -------------------------------------------------------


@given(pmv=st.floats(min_value=-50.0, max_value=50.0, **_FINITE))
def test_ppd_in_range_and_even(pmv):
    # PPD is bounded to [5, 100] % and is an even function of PMV
    ppd = tf.calculate_ppd(_arr(pmv))[0]
    assert 5.0 - 1e-9 <= ppd <= 100.0 + 1e-9
    assert tf.calculate_ppd(_arr(pmv))[0] == tf.calculate_ppd(_arr(-pmv))[0]


@given(
    t2=st.floats(min_value=240.0, max_value=330.0, **_FINITE),
    depression=st.floats(min_value=0.0, max_value=80.0, **_FINITE),
)
def test_relative_humidity_bounds(t2, depression):
    # with the dew point at or below the air temperature, RH is in [0, 100] %
    td = t2 - depression
    rh = tf.calculate_relative_humidity_percent(_arr(t2), _arr(td))[0]
    assert 0.0 <= rh <= 100.0 + 1e-6


if __name__ == "__main__":
    import unittest

    unittest.main()  # pragma: no cover

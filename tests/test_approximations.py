# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

import numpy as np
import pytest

from thermofeel import approximations as ap

# The Erbs and DISC reference values below are produced by pvlib
# (pvlib.irradiance.erbs / .disc), an INDEPENDENT implementation of the same
# models, matched on the solar constant, Earth-Sun distance (Spencer 1971) and
# air-mass conventions. pvlib is a validation ORACLE only -- it is not a runtime
# or test dependency; the numbers are pinned literals.


# --- Earth-Sun distance factor (Spencer 1971) -------------------------------


def test_earth_sun_distance_factor():
    f = ap.earth_sun_distance_factor(np.array([3.0, 185.0]))
    assert f[0] == pytest.approx(1.03507737, abs=1e-7)  # ~perihelion (early Jan)
    assert f[1] == pytest.approx(0.96658938, abs=1e-7)  # ~aphelion (early Jul)
    assert f[0] > 1.0 > f[1]


# --- Erbs (1982) diffuse-fraction decomposition ------------------------------

# (zdeg, ssrd, doy, expected fdir = ghi - dhi from pvlib.erbs, S0=1366.1).
# Rows span the low / mid / high clearness-index branches.
ERBS_ROWS = [
    (30.0, 800.0, 172, 603.520225),
    (45.0, 500.0, 100, 191.234892),
    (60.0, 200.0, 15, 8.177840),
    (80.0, 150.0, 300, 91.995499),
    (20.0, 950.0, 200, 786.038260),
    (85.0, 60.0, 60, 19.772949),
    (78.0, 45.0, 120, 0.651409),  # kt <= 0.22 (very diffuse) branch
]


@pytest.mark.parametrize("zdeg,ssrd,doy,expected", ERBS_ROWS)
def test_approximate_fdir_erbs_vs_pvlib(zdeg, ssrd, doy, expected):
    cossza = np.cos(np.radians([zdeg]))
    fdir = ap.approximate_fdir_erbs(
        np.array([ssrd]), cossza, doy=np.array([float(doy)]), solar_constant=1366.1
    )
    assert fdir[0] == pytest.approx(expected, rel=1e-9, abs=1e-6)


def test_erbs_high_kt_identity():
    # kt > 0.8 -> diffuse fraction floors at 0.165 -> fdir = 0.835 * ssrd.
    cossza = np.cos(np.radians([10.0]))
    fdir = ap.approximate_fdir_erbs(
        np.array([1050.0]), cossza, doy=np.array([172.0]), solar_constant=1366.1
    )
    assert fdir[0] == pytest.approx(0.835 * 1050.0, abs=1e-6)
    # a saturating input (kt clipped to 1) keeps the 0.835 fraction and never
    # exceeds ssrd.
    huge = ap.approximate_fdir_erbs(
        np.array([5000.0]),
        np.array([0.9]),
        doy=np.array([172.0]),
        solar_constant=1366.1,
    )
    assert huge[0] == pytest.approx(0.835 * 5000.0, abs=1e-6)
    assert huge[0] <= 5000.0


def test_erbs_doy_correction():
    # doy=None uses the mean-distance solar constant; supplying doy applies the
    # Earth-Sun distance correction and changes the result.
    args = (np.array([400.0]), np.array([0.7]))
    no_doy = ap.approximate_fdir_erbs(*args)
    perihelion = ap.approximate_fdir_erbs(*args, doy=np.array([3.0]))
    assert no_doy[0] == pytest.approx(76.325537, abs=1e-5)
    assert perihelion[0] == pytest.approx(67.380785, abs=1e-5)
    # nearer Sun -> higher TOA -> lower kt -> more diffuse -> less direct
    assert perihelion[0] < no_doy[0]


# --- DISC (Maxwell 1987) quasi-physical model --------------------------------

# (zdeg, ssrd, doy, pressure_hpa, expected fdir = dni*cossza from pvlib.disc,
# S0=1370, pressure-corrected Kasten-1966 air mass). Machine-precision match.
DISC_ROWS = [
    (30.0, 800.0, 172, 1013.25, 508.095649),
    (45.0, 500.0, 100, 900.00, 154.763199),
    (60.0, 200.0, 15, 1013.25, 9.977575),
    (80.0, 150.0, 300, 1013.25, 115.253655),
    (20.0, 950.0, 200, 950.00, 697.796234),
    (75.0, 90.0, 60, 1013.25, 0.0),  # Kn < 0 -> dni clipped to 0
]


@pytest.mark.parametrize("zdeg,ssrd,doy,phpa,expected", DISC_ROWS)
def test_approximate_fdir_disc_vs_pvlib(zdeg, ssrd, doy, phpa, expected):
    cossza = np.cos(np.radians([zdeg]))
    fdir = ap.approximate_fdir_disc(
        np.array([ssrd]),
        cossza,
        np.array([float(doy)]),
        pressure_hpa=np.array([phpa]),
        solar_constant=1370.0,
    )
    assert fdir[0] == pytest.approx(expected, rel=1e-9, abs=1e-6)


def test_disc_airmass_and_pressure():
    # Relative air mass ~1 at the zenith.
    am0 = ap._relative_airmass_kasten1966(np.array([1.0]))
    assert am0[0] == pytest.approx(1.0, abs=0.01)
    # both pressures give finite direct radiation at moderate sun.
    cossza = np.cos(np.radians([50.0]))
    hi = ap.approximate_fdir_disc(
        np.array([600.0]), cossza, np.array([172.0]), pressure_hpa=np.array([1013.25])
    )
    lo = ap.approximate_fdir_disc(
        np.array([600.0]), cossza, np.array([172.0]), pressure_hpa=np.array([700.0])
    )
    assert np.isfinite(hi[0]) and np.isfinite(lo[0])


# --- robustness: NaN propagation, night, dtype -------------------------------


def test_nan_propagates():
    nan = np.array([np.nan])
    # every input of each estimator propagates NaN to NaN
    assert np.isnan(
        ap.approximate_fdir_erbs(nan, np.array([0.5]), doy=np.array([100.0]))
    ).all()
    assert np.isnan(ap.approximate_fdir_erbs(np.array([500.0]), nan)).all()
    assert np.isnan(
        ap.approximate_fdir_erbs(np.array([500.0]), np.array([0.5]), doy=nan)
    ).all()
    assert np.isnan(
        ap.approximate_fdir_disc(np.array([500.0]), nan, np.array([100.0]))
    ).all()
    assert np.isnan(
        ap.approximate_fdir_disc(np.array([500.0]), np.array([0.5]), nan)
    ).all()
    assert np.isnan(
        ap.approximate_fdir_disc(
            np.array([500.0]), np.array([0.5]), np.array([100.0]), pressure_hpa=nan
        )
    ).all()


def test_night_returns_zero():
    # below min cos(zenith) the Sun is below the horizon -> fdir = 0
    night = np.array([0.02])
    assert (
        ap.approximate_fdir_erbs(np.array([100.0]), night, doy=np.array([100.0]))[0]
        == 0.0
    )
    assert (
        ap.approximate_fdir_disc(np.array([100.0]), night, np.array([100.0]))[0] == 0.0
    )


def test_integer_input_not_truncated():
    out = ap.approximate_fdir_erbs(
        np.array([800]), np.array([0.8]), doy=np.array([172])
    )
    assert out.dtype == float
    out2 = ap.approximate_fdir_disc(np.array([800]), np.array([0.8]), np.array([172]))
    assert out2.dtype == float


if __name__ == "__main__":
    import unittest

    unittest.main()  # pragma: no cover

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
from thermofeel import liljegren
from thermofeel import thermofeel as thermofeel_module
from thermofeel.helpers import fahrenheit_to_celsius

# Representative finite inputs (Kelvin / m s-1 / % / hPa).
T = np.array([300.0])
TD = np.array([290.0])
RH = np.array([50.0])
VA = np.array([3.0])
MRT = np.array([310.0])
Q = np.array([100.0])  # body-absorbed net radiation [W m-2]
NAN = np.array([np.nan])


def _solve_globe(inputs):
    return liljegren.solve_globe(*inputs)


def _solve_wetbulb(inputs):
    return liljegren.solve_wetbulb(*inputs, 1.0)


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
    "discomfort_index": lambda x: tmf.calculate_discomfort_index(x, RH),
    "summer_simmer_index": lambda x: tmf.calculate_summer_simmer_index(x, RH),
    "relative_strain_index": lambda x: tmf.calculate_relative_strain_index(x, RH),
    "apparent_temperature_radiation": (
        lambda x: tmf.calculate_apparent_temperature_radiation(x, VA, RH, Q)
    ),
    "pmv": lambda x: tmf.calculate_pmv(x, MRT, VA, rh=RH),
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


@pytest.mark.parametrize(
    ("solver", "convection", "expected"),
    [
        (_solve_globe, "h_sphere_in_air", 35.96840506),
        (_solve_wetbulb, "h_cylinder_in_air", 19.95920111),
    ],
)
def test_liljegren_solvers_skip_nan_rows_without_extra_iterations(
    monkeypatch, solver, convection, expected
):
    finite_inputs = (
        np.array([298.15]),
        np.array([0.5]),
        np.array([1013.0]),
        np.array([2.0]),
        np.array([500.0]),
        np.array([0.5]),
        np.array([0.5]),
    )
    original_convection = getattr(liljegren, convection)
    calls = 0

    def count_iterations(*args):
        nonlocal calls
        calls += 1
        return original_convection(*args)

    monkeypatch.setattr(liljegren, convection, count_iterations)

    finite = solver(finite_inputs)
    finite_calls = calls
    np.testing.assert_allclose(finite, [expected], rtol=1e-8)

    calls = 0
    mixed_inputs = tuple(
        np.array([values[0], np.nan if index == 0 else values[0]])
        for index, values in enumerate(finite_inputs)
    )
    mixed = solver(mixed_inputs)
    assert calls == finite_calls
    np.testing.assert_allclose(mixed[0], finite[0])
    assert np.isnan(mixed[1])

    calls = 0
    all_invalid_inputs = (np.array([np.nan]), *finite_inputs[1:])
    all_invalid = solver(all_invalid_inputs)
    assert calls == 0
    assert np.isnan(all_invalid).all()


def test_ppd_nan_propagates():
    # calculate_ppd is a pure map of PMV; a NaN vote (e.g. a non-converged PMV
    # element) propagates to NaN rather than raising or returning a spurious %.
    assert np.isnan(tmf.calculate_ppd(NAN)).all()


def test_relative_strain_index_singularity_and_sign_flip():
    # RSI = (Ta - 21) / (58 - e). Near saturation around 35-36 degC the vapour
    # pressure e approaches, then exceeds, 58 hPa: the denominator shrinks to zero
    # and flips sign, so RSI blows up and then turns negative. This pins the
    # documented, deliberately-unclamped out-of-domain edge (ROBUSTNESS.md R-9).
    hot_saturated = tmf.celsius_to_kelvin(np.array([35.0, 36.0]))
    rh = np.array([100.0, 100.0])
    with np.errstate(divide="ignore", invalid="ignore"):
        rsi = tmf.calculate_relative_strain_index(hot_saturated, rh)
    assert rsi[0] > 5.0  # 35 degC / 100% RH: e ~ 56 hPa -> large positive RSI
    assert rsi[1] < 0.0  # 36 degC / 100% RH: e ~ 59 hPa > 58 -> spurious negative


def test_pmv_skips_nan_rows_without_extra_iterations(monkeypatch):
    finite_inputs = (
        tmf.celsius_to_kelvin(np.array([22.0])),
        tmf.celsius_to_kelvin(np.array([22.0])),
        np.array([0.1]),
        np.array([60.0]),
    )
    original_maximum = thermofeel_module.np.maximum
    calls = 0

    def count_iterations(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original_maximum(*args, **kwargs)

    monkeypatch.setattr(thermofeel_module.np, "maximum", count_iterations)

    finite = tmf.calculate_pmv(*finite_inputs[:3], rh=finite_inputs[3])
    finite_calls = calls
    np.testing.assert_allclose(finite, [-0.75256792], rtol=1e-8)

    calls = 0
    mixed = tmf.calculate_pmv(
        tmf.celsius_to_kelvin(np.array([22.0, np.nan])),
        tmf.celsius_to_kelvin(np.array([22.0, 22.0])),
        np.array([0.1, 0.1]),
        rh=np.array([60.0, 60.0]),
    )
    assert calls == finite_calls
    np.testing.assert_allclose(mixed[0], finite[0])
    assert np.isnan(mixed[1])

    calls = 0
    all_invalid = tmf.calculate_pmv(
        tmf.celsius_to_kelvin(np.array([np.nan])),
        finite_inputs[1],
        finite_inputs[2],
        rh=finite_inputs[3],
    )
    assert calls == 0
    assert np.isnan(all_invalid).all()


def test_bgt_zero_wind_returns_mrt():
    # Calm-air limit: with no convection the globe sits at radiative equilibrium,
    # so the globe temperature equals the mean radiant temperature at va == 0
    # (the closed form is 0/0 there). Positive wind is finite as usual.
    zero = tmf.calculate_bgt(T, MRT, np.array([0.0]))
    positive = tmf.calculate_bgt(T, MRT, np.array([0.5]))
    np.testing.assert_allclose(zero, MRT)
    assert np.isfinite(positive).all()


def test_wbgt_zero_wind_is_finite():
    # calculate_wbgt depends on calculate_bgt, so the calm-air limit (bgt -> mrt)
    # makes the Stull WBGT finite at va == 0 too (previously NaN).
    out = tmf.calculate_wbgt(T, MRT, np.array([0.0]), TD)
    assert np.isfinite(out).all()


def test_bgt_negative_wind_is_nan():
    # Negative wind is invalid input; it is not masked by the calm-air limit.
    with np.errstate(invalid="ignore"):
        out = tmf.calculate_bgt(T, MRT, np.array([-1.0]))
    assert np.isnan(out).all()


def test_utci_requires_ehpa_or_td():
    with pytest.raises(ValueError):
        tmf.calculate_utci(T, VA, MRT)


def test_pmv_humidity_contract():
    # Exactly one humidity input is required (mirrors calculate_utci).
    t = tmf.celsius_to_kelvin(np.array([22.0]))
    v = np.array([0.1])
    with pytest.raises(ValueError):
        tmf.calculate_pmv(t, t, v)  # neither rh nor vapour_pressure_hpa
    with pytest.raises(ValueError):
        tmf.calculate_pmv(
            t, t, v, rh=np.array([60.0]), vapour_pressure_hpa=np.array([14.0])
        )  # both
    # Passing vapour pressure (converted from an RH row with the same ISO relation
    # the function uses internally) reproduces the rh-driven PMV.
    ta_c, rh = 22.0, 60.0
    pa_pa = rh * 10.0 * np.exp(16.6536 - 4030.183 / (ta_c + 235.0))
    vp_hpa = np.array([pa_pa / 100.0])
    pmv_rh = tmf.calculate_pmv(t, t, v, rh=np.array([rh]))
    pmv_vp = tmf.calculate_pmv(t, t, v, vapour_pressure_hpa=vp_hpa)
    np.testing.assert_allclose(pmv_vp, pmv_rh)


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


def test_multiphase_integer_input_not_truncated():
    # An integer-typed temperature array must not silently truncate the float
    # vapour pressure on assignment (regression for the np.zeros_like dtype bug).
    phase = np.array([0, 0])
    out_int = tmf.calculate_saturation_vapour_pressure_multiphase(
        np.array([290, 300]), phase
    )
    out_float = tmf.calculate_saturation_vapour_pressure_multiphase(
        np.array([290.0, 300.0]), phase
    )
    np.testing.assert_allclose(out_int, out_float)
    assert not np.allclose(out_int, np.floor(out_float))  # genuinely fractional


def test_approximate_dsrp_integer_input_not_truncated():
    # Integer-typed fdir must not truncate the division result (regression for
    # the np.copy dtype bug).
    cossza = np.array([0.3, 0.3])
    out_int = tmf.approximate_dsrp(np.array([400, 600]), cossza)
    out_float = tmf.approximate_dsrp(np.array([400.0, 600.0]), cossza)
    np.testing.assert_allclose(out_int, out_float)

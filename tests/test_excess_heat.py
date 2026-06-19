# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

import numpy as np
import pytest

import thermofeel as tmf


def test_daily_mean_temperature():
    t2_min = np.asarray([1.0, 2.0, 3.0, 4.0, 5.0])
    t2_max = np.asarray([5.0, 2.0, 5.0, 4.0, 7.0])
    dmt = tmf.excess_heat.daily_mean_temperature(t2_min, t2_max)
    np.testing.assert_allclose(dmt, [3.0, 2.0, 4.0, 4.0, 6.0])


def test_daily_mean_temperature_scalar():
    dmt = tmf.excess_heat.daily_mean_temperature(1.0, 5.0)
    assert dmt == pytest.approx(3.0)


def test_significance_index():
    dmt = np.asarray([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
    threshold = 5.0
    ehi_sig = tmf.excess_heat.significance_index(dmt, threshold)
    np.testing.assert_allclose(ehi_sig, [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0])


def test_significance_index_scalar():
    ehi_sig = tmf.excess_heat.significance_index(8.0, 5.0)
    assert ehi_sig == pytest.approx(3.0)


def test_acclimatisation_index():
    dmt = np.asarray([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
    threshold = np.asarray([2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 3.0])
    ehi_accl = tmf.excess_heat.acclimatisation_index(dmt, threshold)
    np.testing.assert_allclose(ehi_accl, [1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 6.0])


def test_acclimatisation_index_scalar():
    ehi_accl = tmf.excess_heat.acclimatisation_index(6.0, 3.0)
    assert ehi_accl == pytest.approx(3.0)


class TestExcessHeatAndColdFactors:
    @pytest.fixture
    def ehi_sig(self):
        return np.asarray([5.0, 5.0, 5.0, 0.0, 0.0, 0.0, -2.0, -2.0, -2.0])

    @pytest.fixture
    def ehi_accl(self):
        return np.asarray([3.0, 0.0, -4.0, 3.0, 0.0, -4.0, 3.0, 0.0, -4.0])

    def test_excess_heat_factor_with_clip_false(self, ehi_sig, ehi_accl):
        exhf = tmf.excess_heat.excess_heat_factor(ehi_sig, ehi_accl, clip=False)
        np.testing.assert_allclose(
            exhf, [15.0, 5.0, 5.0, 0.0, 0.0, 0.0, -6.0, -2.0, -2.0]
        )

    def test_excess_heat_factor_with_clip_true(self, ehi_sig, ehi_accl):
        exhf = tmf.excess_heat.excess_heat_factor(ehi_sig, ehi_accl, clip=True)
        np.testing.assert_allclose(exhf, [15.0, 5.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    def test_excess_heat_factor_clip_false_positive_scalar(self):
        assert tmf.excess_heat.excess_heat_factor(
            5.0, 3.0, clip=False
        ) == pytest.approx(15.0)

    def test_excess_heat_factor_clip_false_negative_scalar(self):
        assert tmf.excess_heat.excess_heat_factor(
            -2.0, 3.0, clip=False
        ) == pytest.approx(-6.0)

    def test_excess_heat_factor_clip_true_negative_scalar(self):
        assert tmf.excess_heat.excess_heat_factor(
            -2.0, 3.0, clip=True
        ) == pytest.approx(0.0)

    def test_excess_cold_factor_with_clip_false(self, ehi_sig, ehi_accl):
        excf = tmf.excess_heat.excess_cold_factor(ehi_sig, ehi_accl, clip=False)
        np.testing.assert_allclose(
            excf, [5.0, 5.0, 20.0, 0.0, 0.0, 0.0, -2.0, -2.0, -8.0]
        )

    def test_excess_cold_factor_with_clip_true(self, ehi_sig, ehi_accl):
        excf = tmf.excess_heat.excess_cold_factor(ehi_sig, ehi_accl, clip=True)
        np.testing.assert_allclose(
            excf, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, -2.0, -8.0]
        )

    def test_excess_cold_factor_clip_false_positive_scalar(self):
        assert tmf.excess_heat.excess_cold_factor(
            5.0, -4.0, clip=False
        ) == pytest.approx(20.0)

    def test_excess_cold_factor_clip_false_negative_scalar(self):
        assert tmf.excess_heat.excess_cold_factor(
            -2.0, -4.0, clip=False
        ) == pytest.approx(-8.0)

    def test_excess_cold_factor_clip_true_positive_scalar(self):
        assert tmf.excess_heat.excess_cold_factor(
            5.0, -4.0, clip=True
        ) == pytest.approx(0.0)


def test_heatwave_severity():
    exhf = np.asarray([[0.0, 1.0], [1.0, 2.0], [3.0, 1.0], [4.0, 0.0]])
    tr = np.asarray([2.5, 2.0])
    hsev = tmf.excess_heat.heatwave_severity(exhf, threshold=tr)
    ref = np.asarray([[0.0, 0.5], [0.4, 1.0], [1.2, 0.5], [1.6, 0.0]])
    np.testing.assert_allclose(hsev, ref)


def test_heatwave_severity_scalar():
    hsev = tmf.excess_heat.heatwave_severity(3.0, threshold=2.5)
    assert hsev == pytest.approx(1.2)


# The excess heat/cold factors are also exposed at the top level as
# calculate_excess_{heat,cold}_factor, wrapping the excess_heat submodule.


def test_calculate_excess_heat_factor_matches_submodule():
    ehi_sig = np.asarray([5.0, 0.0, -2.0])
    ehi_accl = np.asarray([3.0, 0.0, -4.0])
    np.testing.assert_allclose(
        tmf.calculate_excess_heat_factor(ehi_sig, ehi_accl),
        tmf.excess_heat.excess_heat_factor(ehi_sig, ehi_accl),
    )
    np.testing.assert_allclose(
        tmf.calculate_excess_heat_factor(ehi_sig, ehi_accl, clip=True),
        tmf.excess_heat.excess_heat_factor(ehi_sig, ehi_accl, clip=True),
    )


def test_calculate_excess_cold_factor_matches_submodule():
    ehi_sig = np.asarray([5.0, 0.0, -2.0])
    ehi_accl = np.asarray([-4.0, 0.0, 3.0])
    np.testing.assert_allclose(
        tmf.calculate_excess_cold_factor(ehi_sig, ehi_accl),
        tmf.excess_heat.excess_cold_factor(ehi_sig, ehi_accl),
    )
    np.testing.assert_allclose(
        tmf.calculate_excess_cold_factor(ehi_sig, ehi_accl, clip=True),
        tmf.excess_heat.excess_cold_factor(ehi_sig, ehi_accl, clip=True),
    )


def test_calculate_excess_factor_scalars():
    assert tmf.calculate_excess_heat_factor(5.0, 3.0) == pytest.approx(15.0)
    assert tmf.calculate_excess_cold_factor(5.0, -4.0) == pytest.approx(20.0)


def test_excess_factor_raw_names_not_exported_at_top_level():
    # Only the calculate_* wrappers are public at the top level; the bare
    # submodule names are deliberately not re-exported there.
    assert not hasattr(tmf, "excess_heat_factor")
    assert not hasattr(tmf, "excess_cold_factor")

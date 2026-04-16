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
    t2_min = np.asarray([1., 2., 3., 4., 5.])
    t2_max = np.asarray([5., 2., 5., 4., 7.])
    dmt = tmf.excess_heat.daily_mean_temperature(t2_min, t2_max)
    np.testing.assert_allclose(dmt, [3., 2., 4., 4., 6.])


def test_significance_index():
    dmt = np.asarray([3., 4., 5., 6., 7., 8., 9.])
    threshold = 5.
    ehi_sig = tmf.excess_heat.significance_index(dmt, threshold)
    np.testing.assert_allclose(ehi_sig, [-2., -1., 0., 1., 2., 3., 4.])


def test_acclimatisation_index():
    dmt = np.asarray([3., 4., 5., 6., 7., 8., 9.])
    threshold = np.asarray([2., 2., 3., 3., 4., 4., 3.])
    ehi_accl = tmf.excess_heat.acclimatisation_index(dmt, threshold)
    np.testing.assert_allclose(ehi_accl, [1., 2., 2., 3., 3., 4., 6.])


class TestExcessHeatAndColdFactors:

    @pytest.fixture
    def ehi_sig(self):
        return np.asarray([5.0, 5.0, 5.0, 0.0, 0.0, 0.0, -2.0, -2.0, -2.0])

    @pytest.fixture
    def ehi_accl(self):
        return np.asarray([3.0, 0.0, -4.0, 3.0, 0.0, -4.0, 3.0, 0.0, -4.0])

    def test_excess_heat_factor_with_clip_false(self, ehi_sig, ehi_accl):
        exhf = tmf.excess_heat.excess_heat_factor(ehi_sig, ehi_accl, clip=False)
        np.testing.assert_allclose(exhf, [15.0, 5.0, 5.0, 0.0, 0.0, 0.0, -6.0, -2.0, -2.0])

    def test_excess_heat_factor_with_clip_true(self, ehi_sig, ehi_accl):
        exhf = tmf.excess_heat.excess_heat_factor(ehi_sig, ehi_accl, clip=True)
        np.testing.assert_allclose(exhf, [15.0, 5.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    def test_excess_heat_factor_with_clip_false(self, ehi_sig, ehi_accl):
        excf = tmf.excess_heat.excess_cold_factor(ehi_sig, ehi_accl, clip=False)
        np.testing.assert_allclose(excf, [5.0, 5.0, 20.0, 0.0, 0.0, 0.0, -2.0, -2.0, -8.0])

    def test_excess_heat_factor_with_clip_true(self, ehi_sig, ehi_accl):
        excf = tmf.excess_heat.excess_cold_factor(ehi_sig, ehi_accl, clip=True)
        np.testing.assert_allclose(excf, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, -2.0, -8.0])


def test_heatwave_severity():
    exhf = np.asarray([[0.0, 1.0], [1.0, 2.0], [3.0, 1.0], [4.0, 0.0]])
    tr = np.asarray([2.5, 2.0])
    hsev = tmf.excess_heat.heatwave_severity(exhf, threshold=tr)
    ref = np.asarray([[0.0, 0.5], [0.4, 1.0], [1.2, 0.5], [1.6, 0.0]])
    np.testing.assert_allclose(hsev, ref)

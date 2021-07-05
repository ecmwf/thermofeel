# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import unittest

import numpy as np
from context import data_file

import thermofeel.thermofeel as tfc

# import pytest


# combining the pytest library with the functionality of np testing for use with np.array file type


class TestThermalCalculator(unittest.TestCase):
    def setUp(self):
        # variables for use in thermalindexcalculator
        t = np.genfromtxt(
            data_file("thermofeel_testcases.csv"),
            delimiter=",",
            names=True,
        )
        self.t2m = t["t2m"]
        self.ssr = t["ssr"]
        self.td = t["td"]
        self.va = t["va"]
        self.mrt = t["mrt"]
        self.ssrd = t["ssrd"]
        self.strd = t["strd"]
        self.fdir = t["fdir"]
        self.strr = t["strr"]
        self.cossza = t["cossza"]

        # indices
        tr = np.genfromtxt(
            data_file("thermofeel_test_results.csv"),
            delimiter=",",
            names=True,
        )
        self.rh = tfc.calculate_saturation_vapour_pressure(self.t2m)
        self.rhpercent = tfc.calculate_relative_humidity_percent(self.t2m, self.td)
        self.heatindex = tr["heatindex"]
        self.heatindexadjusted = tr["heatindexadjust"]
        self.utci = tr["utci"]
        self.apparenttemperature = tr["apparenttemperature"]
        self.wbgts = tr["wbgts"]
        self.wbgt = tr["wbgt"]
        self.net = tr["net"]
        self.humidex = tr["humidex"]
        self.windchill = tr["windchill"]
        self.mrtr = tr["mrt"]
        self.mrtw = tr["mrtw"]

    def assert_equal(self, result, calculation):
        self.assertequal = self.assertIsNone(
            np.testing.assert_array_almost_equal(result, calculation, decimal=6)
        )

    def assert_equal_less_precise(self, result, calculation):
        self.assertequal = self.assertIsNone(
            np.testing.assert_array_almost_equal(result, calculation, decimal=1)
        )

    def test_relative_humidity(self):
        self.assert_equal(self.rh, tfc.calculate_saturation_vapour_pressure(self.t2m))

    def test_relative_humidity_percent(self):
        self.assert_equal(
            self.rhpercent, tfc.calculate_relative_humidity_percent(self.t2m, self.td)
        )

    def test_heat_index(self):
        self.assert_equal_less_precise(
            self.heatindex, tfc.calculate_heat_index_simplified(self.t2m)
        )

    def test_mean_radiant_temperature(self):
        self.assert_equal_less_precise(
            self.mrtr,
            tfc.calculate_mean_radiant_temperature(
                self.ssrd / 3600,
                self.ssr / 3600,
                self.fdir / 3600,
                self.strd / 3600,
                self.strr / 3600,
                self.cossza / 3600,
            ),
        )

    def test_utci(self):
        rh_pc = tfc.calculate_relative_humidity_percent(self.t2m, self.td)
        ehPa = tfc.calculate_saturation_vapour_pressure(self.t2m) * rh_pc / 100.0
        self.assert_equal_less_precise(
            self.utci,
            tfc.calculate_utci(
                t2_k=self.t2m, va_ms=self.va, mrt_k=self.mrt, e_hPa=ehPa
            ),
        )

    def test_apparent_temperature(self):
        self.assert_equal_less_precise(
            self.apparenttemperature,
            tfc.calculate_apparent_temperature(self.t2m, self.va),
        )

    def test_wbgts(self):
        self.assert_equal_less_precise(self.wbgts, tfc.calculate_wbgts(self.t2m))

    def test_net(self):
        self.assert_equal_less_precise(
            self.net,
            tfc.calculate_net_effective_temperature(self.t2m, self.va, self.td),
        )

    def test_humidex(self):
        self.assert_equal(self.humidex, tfc.calculate_humidex(self.t2m, self.td))

    def test_wind_chill(self):
        self.assert_equal(self.windchill, tfc.calculate_wind_chill(self.t2m, self.va))

    def test_heat_index_adjusted(self):
        self.assert_equal(
            self.heatindexadjusted, tfc.calculate_heat_index_adjusted(self.t2m, self.td)
        )

    def test_wbgt(self):
        self.assert_equal_less_precise(
            self.wbgt, tfc.calculate_wbgt(self.t2m, self.va, self.mrt)
        )

    def test_mrt_from_wbgt(self):
        wbgt_k = tfc.celcius_to_kelvin(self.wbgt)
        mrt = tfc.calculate_mrt_from_wbgt(self.t2m, wbgt_k, self.va)
        self.assert_equal_less_precise(self.mrtw, mrt)


if __name__ == "__main__":
    unittest.main()

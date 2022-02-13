# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import os
import unittest

import numpy as np

import thermofeel as tmf


def data_file(name):
    return os.path.join(os.path.dirname(__file__), name)


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

        self.rh = np.loadtxt(data_file("rh.csv"))
        self.es = np.loadtxt(data_file("es.csv"))
        self.utci = np.loadtxt(data_file("utci.csv"))
        self.net = np.loadtxt(data_file("net.csv"))
        self.heatindexadjusted = np.loadtxt(data_file("hia.csv"))
        self.heatindex = np.loadtxt(data_file("heatindex.csv"))
        self.at = np.loadtxt(data_file("at.csv"))
        self.wbgts = np.loadtxt(data_file("wbgts.csv"))
        self.wbgt = np.loadtxt(data_file("wbgt.csv"))
        self.humidex = np.loadtxt(data_file("humidex.csv"))
        self.windchill = np.loadtxt(data_file("windchill.csv"))
        self.mrtr = np.loadtxt(data_file("mrtr.csv"))
        self.mrtw = np.loadtxt(data_file("mrtw.csv"))

    def assert_equal(self, expected, result, decimal=6):
        self.assertequal = self.assertIsNone(
            np.testing.assert_array_almost_equal(expected, result, decimal)
        )

    def test_saturation_vapour_pressure(self):
        es = tmf.calculate_saturation_vapour_pressure(self.t2m)
        # np.savetxt("es.csv", es)
        self.assert_equal(self.es, es)

    def test_relative_humidity_percent(self):
        rh = tmf.calculate_relative_humidity_percent(self.t2m, self.td)
        # np.savetxt("rh.csv", rh)
        self.assert_equal(self.rh, rh)

    def test_heat_index(self):
        heatindex = tmf.calculate_heat_index_simplified(self.t2m)
        # np.savetxt("heatindex.csv", heatindex)
        self.assert_equal(self.heatindex, heatindex)

    def test_mean_radiant_temperature(self):
        mrtr = tmf.calculate_mean_radiant_temperature(
            self.ssrd / 3600,
            self.ssr / 3600,
            self.fdir / 3600,
            self.strd / 3600,
            self.strr / 3600,
            self.cossza / 3600,
        )
        # np.savetxt("mrtr.csv", mrtr)
        self.assert_equal(self.mrtr, mrtr)

    def test_utci(self):
        rh_pc = tmf.calculate_relative_humidity_percent(self.t2m, self.td)
        ehPa = tmf.calculate_saturation_vapour_pressure(self.t2m) * rh_pc / 100.0
        utci = tmf.calculate_utci(
            t2_k=self.t2m, va_ms=self.va, mrt_k=self.mrt, ehPa=ehPa
        )
        # np.savetxt("utci.csv", utci)
        # self.assert_equal(self.utci, utci)

    def test_apparent_temperature(self):
        at = tmf.calculate_apparent_temperature(self.t2m, self.va)
        # np.savetxt("at.csv", at)
        self.assert_equal(self.at, at)

    def test_net_effective_temperature(self):
        net = tmf.calculate_net_effective_temperature(self.t2m, self.va, self.td)
        # np.savetxt("net.csv", net)
        self.assert_equal(self.net, net)

    def test_humidex(self):
        humidex = tmf.calculate_humidex(self.t2m, self.td)
        # np.savetxt("humidex.csv", humidex)
        self.assert_equal(self.humidex, humidex)

    def test_wind_chill(self):
        windchill = tmf.calculate_wind_chill(self.t2m, self.va)
        # np.savetxt("windchill.csv", windchill)
        self.assert_equal(self.windchill, windchill)

    def test_heat_index_adjusted(self):
        hia = tmf.calculate_heat_index_adjusted(self.t2m, self.td)
        #np.savetxt("hia.csv", hia)
        self.assert_equal(self.heatindexadjusted, hia)

    # def test_wbgt(self):
    #     wbgt = tmf.calculate_wbgt(self.t2m, self.va, self.mrt)
    #     # np.savetxt("wbgt.csv", wbgt)
    #     self.assert_equal(self.wbgt, wbgt)

    def test_wbgts(self):
        wbgts = tmf.calculate_wbgts(self.t2m)
        # np.savetxt("wbgts.csv", wbgts)
        self.assert_equal(self.wbgts, wbgts)

    def test_mrt_from_bgt(self):
        wbgt_k = tmf.celsius_to_kelvin(self.wbgt)
        mrtw = tmf.calculate_mrt_from_bgt(self.t2m, wbgt_k, self.va)
        # np.savetxt("mrtw.csv", mrtw)
        self.assert_equal(self.mrtw, mrtw)


if __name__ == "__main__":
    unittest.main()  # pragma: no cover

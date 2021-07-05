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

import thermofeel.thermofeel as tfc


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
        self.at = np.loadtxt(data_file("apparenttemperature.csv"))
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
        es = tfc.calculate_saturation_vapour_pressure(self.t2m)
        self.assert_equal(self.es, es)
        # np.savetxt("es.csv", es)

    def test_relative_humidity_percent(self):
        rh = tfc.calculate_relative_humidity_percent(self.t2m, self.td)
        self.assert_equal(self.rh, rh)
        # np.savetxt("rh.csv", rh)

    def test_heat_index(self):
        heatindex = tfc.calculate_heat_index_simplified(self.t2m)
        self.assert_equal(self.heatindex, heatindex)
        # np.savetxt("heatindex.csv", heatindex)

    def test_mean_radiant_temperature(self):
        mrtr = tfc.calculate_mean_radiant_temperature(
            self.ssrd / 3600,
            self.ssr / 3600,
            self.fdir / 3600,
            self.strd / 3600,
            self.strr / 3600,
            self.cossza / 3600,
        )
        self.assert_equal(self.mrtr, mrtr)
        # np.savetxt("mrtr.csv", mrtr)

    def test_utci(self):
        rh_pc = tfc.calculate_relative_humidity_percent(self.t2m, self.td)
        ehPa = tfc.calculate_saturation_vapour_pressure(self.t2m) * rh_pc / 100.0
        utci = tfc.calculate_utci(
            t2_k=self.t2m, va_ms=self.va, mrt_k=self.mrt, e_hPa=ehPa
        )
        self.assert_equal(self.utci, utci)
        # np.savetxt("utci.csv", utci)

    def test_apparent_temperature(self):
        at = tfc.calculate_apparent_temperature(self.t2m, self.va)
        self.assert_equal(self.apparenttemperature, at)
        # np.savetxt("apparenttemperature.csv", at)

    def test_net_effective_temperature(self):
        net = tfc.calculate_net_effective_temperature(self.t2m, self.va, self.td)
        self.assert_equal(self.net, net)
        # np.savetxt("net.csv", net)

    def test_humidex(self):
        humidex = tfc.calculate_humidex(self.t2m, self.td)
        self.assert_equal(self.humidex, humidex)
        # np.savetxt("humidex.csv", humidex)

    def test_wind_chill(self):
        windchill = tfc.calculate_wind_chill(self.t2m, self.va)
        self.assert_equal(self.windchill, windchill)
        # np.savetxt("windchill.csv", windchill)

    def test_heat_index_adjusted(self):
        hia = tfc.calculate_heat_index_adjusted(self.t2m, self.td)
        self.assert_equal(self.heatindexadjusted, hia)
        # np.savetxt("hia.csv", hia)

    def test_wbgt(self):
        wbgt = tfc.calculate_wbgt(self.t2m, self.va, self.mrt)
        self.assert_equal(self.wbgt, wbgt)
        # np.savetxt("wbgt.csv", wbgt)

    def test_wbgts(self):
        wbgts = tfc.calculate_wbgts(self.t2m)
        self.assert_equal(self.wbgts, wbgts)
        # np.savetxt("wbgts.csv", wbgts)

    def test_mrt_from_wbgt(self):
        wbgt_k = tfc.celcius_to_kelvin(self.wbgt)
        mrtw = tfc.calculate_mrt_from_wbgt(self.t2m, wbgt_k, self.va)
        self.assert_equal(self.mrtw, mrtw)
        # np.savetxt("mrtw.csv", mrtw)


if __name__ == "__main__":
    unittest.main()  # pragma: no cover

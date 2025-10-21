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
        self.phase = t["phase"]
        self.va_height = t["va_height"]

        self.rh = np.loadtxt(data_file("rh.csv"))
        self.es = np.loadtxt(data_file("es.csv"))
        self.ens = np.loadtxt(data_file("ens.csv"))
        self.es_multiphase = np.loadtxt(data_file("es_multiphase.csv"))
        self.va_scaled = np.loadtxt(data_file("va_scaled.csv"))
        self.utci = np.loadtxt(data_file("utci.csv"))
        self.wbgts = np.loadtxt(data_file("wbgts.csv"))
        self.wbt = np.loadtxt(data_file("wbt.csv"))
        self.bgt = np.loadtxt(data_file("bgt.csv"))
        self.wbgt = np.loadtxt(data_file("wbgt.csv"))
        self.humidex = np.loadtxt(data_file("humidex.csv"))
        self.net = np.loadtxt(data_file("net.csv"))
        self.at = np.loadtxt(data_file("at.csv"))
        self.windchill = np.loadtxt(data_file("windchill.csv"))
        self.heatindex = np.loadtxt(data_file("heatindex.csv"))
        self.heatindexadjusted = np.loadtxt(data_file("hia.csv"))
        self.mrtr = np.loadtxt(data_file("mrtr.csv"))
        self.mrt_from_bgt = np.loadtxt(data_file("mrt_from_bgt.csv"))

        self.dsrp = tmf.approximate_dsrp(self.fdir, self.cossza)

    def assert_equal(self, expected, result, decimal=6):
        self.assertequal = self.assertIsNone(
            np.testing.assert_array_almost_equal(expected, result, decimal)
        )

    def test_relative_humidity_percent(self):
        rh = tmf.calculate_relative_humidity_percent(self.t2m, self.td)
        # np.savetxt("rh.csv", rh)
        self.assert_equal(self.rh, rh)

    def test_saturation_vapour_pressure(self):
        es = tmf.calculate_saturation_vapour_pressure(self.t2m)
        # np.savetxt("es.csv", es)
        self.assert_equal(self.es, es)

    def test_saturation_vapour_pressure_multiphase(self):
        es_multiphase = tmf.calculate_saturation_vapour_pressure_multiphase(
            self.t2m, self.phase
        )
        # np.savetxt("es_multiphase.csv", es_multiphase)
        self.assert_equal(self.es_multiphase, es_multiphase)

    def test_nonsaturation_vapour_pressure(self):
        rh = tmf.calculate_relative_humidity_percent(self.t2m, self.td)
        ens = tmf.calculate_nonsaturation_vapour_pressure(self.t2m, rh)
        # np.savetxt("ens.csv", ens)
        self.assert_equal(self.ens, ens)

    def test_scale_windspeed(self):
        va_scaled = tmf.scale_windspeed(self.va, self.va_height)
        # np.savetxt("va_scaled.csv", va_scaled)
        self.assert_equal(self.va_scaled, va_scaled)

    def test_dew_point_from_relative_humidity(self):
        td = tmf.calculate_dew_point_from_relative_humidity(self.rh, self.t2m)
        self.assert_equal(self.td, td, decimal=2)

    def test_mean_radiant_temperature(self):
        mrtr = tmf.calculate_mean_radiant_temperature(
            ssrd=self.ssrd / 3600,
            ssr=self.ssr / 3600,
            dsrp=self.dsrp / 3600,
            strd=self.strd / 3600,
            fdir=self.fdir / 3600,
            strr=self.strr / 3600,
            cossza=self.cossza / 3600,
        )
        # np.savetxt("mrtr.csv", mrtr)
        self.assert_equal(self.mrtr, mrtr)  # , decimal = 3)

    def test_utci(self):
        rh_pc = tmf.calculate_relative_humidity_percent(self.t2m, self.td)
        ehPa = tmf.calculate_saturation_vapour_pressure(self.t2m) * rh_pc / 100.0
        utci = tmf.calculate_utci(
            t2_k=self.t2m, va=self.va, mrt=self.mrt, td_k=None, ehPa=ehPa
        )
        # np.savetxt("utci.csv", utci)
        self.assert_equal(self.utci, utci)

    def test_wbgt_simple(self):
        rh_pc = tmf.calculate_relative_humidity_percent(self.t2m, self.td)
        wbgts = tmf.calculate_wbgt_simple(self.t2m, rh_pc)
        # np.savetxt("wbgts.csv", wbgts)
        self.assert_equal(self.wbgts, wbgts)

    def test_wbt(self):
        rh_pc = tmf.calculate_relative_humidity_percent(self.t2m, self.td)
        wbt = tmf.calculate_wbt(self.t2m, rh_pc)
        # np.savetxt("wbt.csv", wbt)
        self.assert_equal(self.wbt, wbt)

    def test_bgt(self):
        bgt = tmf.calculate_bgt(self.t2m, self.mrt, self.va)
        # np.savetxt("bgt.csv", bgt)
        self.assert_equal(self.bgt, bgt)

    def test_wbgt(self):
        wbgt = tmf.calculate_wbgt(self.t2m, self.mrt, self.va, self.td)
        # np.savetxt("wbgt.csv", wbgt)
        self.assert_equal(self.wbgt, wbgt)

    def test_mrt_from_bgt(self):
        mrt_from_bgt = tmf.calculate_mrt_from_bgt(self.t2m, self.bgt, self.va)
        # np.savetxt("mrt_from_bgt.csv", mrt_from_bgt)
        self.assert_equal(self.mrt_from_bgt, mrt_from_bgt)

    def test_humidex(self):
        humidex = tmf.calculate_humidex(self.t2m, self.td)
        # np.savetxt("humidex.csv", humidex)
        self.assert_equal(self.humidex, humidex)

    def test_normal_effective_temperature(self):
        rh_pc = tmf.calculate_relative_humidity_percent(self.t2m, self.td)
        net = tmf.calculate_normal_effective_temperature(self.t2m, self.va, rh_pc)
        # np.savetxt("net.csv", net)
        self.assert_equal(self.net, net)

    def test_apparent_temperature(self):
        rh_pc = tmf.calculate_relative_humidity_percent(self.t2m, self.td)
        at = tmf.calculate_apparent_temperature(self.t2m, self.va, rh_pc)
        # np.savetxt("at.csv", at)
        self.assert_equal(self.at, at)

    def test_wind_chill(self):
        windchill = tmf.calculate_wind_chill(self.t2m, self.va)
        # np.savetxt("windchill.csv", windchill)
        self.assert_equal(self.windchill, windchill)

    def test_heat_index(self):
        rh_pc = tmf.calculate_relative_humidity_percent(self.t2m, self.td)
        heatindex = tmf.calculate_heat_index_simplified(self.t2m, rh_pc)
        # np.savetxt("heatindex.csv", heatindex)
        self.assert_equal(self.heatindex, heatindex)

    def test_heat_index_adjusted(self):
        hia = tmf.calculate_heat_index_adjusted(self.t2m, self.td)
        # np.savetxt("hia.csv", hia)
        self.assert_equal(self.heatindexadjusted, hia)


if __name__ == "__main__":
    unittest.main()  # pragma: no cover

# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import unittest
from math import cos, radians

import numpy as np
import pytest

import thermofeel


class TestThermalCalculator(unittest.TestCase):
    def test_relative_humidity_percent(self):
        t2mk = thermofeel.celcius_to_kelvin(30.0)
        t2dk = thermofeel.celcius_to_kelvin(28.0)

        rhpc = thermofeel.calculate_relative_humidity_percent(t2mk, t2dk)

        assert rhpc == pytest.approx(89.08526710467393, abs=1e-6)

    def test_saturation_vapour_pressure(self):
        t2mk = thermofeel.celcius_to_kelvin(25.0)
        svp = thermofeel.calculate_saturation_vapour_pressure(t2mk)
        assert svp == pytest.approx(31.699201897293, abs=1e-6)

    def test_heat_index(self):
        t2mk = thermofeel.celcius_to_kelvin(25.0)
        hi = thermofeel.calculate_heat_index_simplified(t2mk)
        # print(f"hi {hi}")
        assert hi == pytest.approx(49.321122555, abs=1e-6)

    # def test_meant_radiant_temperature(self):
    #     ssrd= 15000
    #     ssr = 14992
    #     fdir = 15002
    #     strd = 14993
    #     strr = 15001
    #     cossza = 0.4
    #     mrt = thermofeel.calculate_mean_radiant_temperature(ssrd/3600,
    #                                                         ssr/3600,
    #                                                         fdir/3600,
    #                                                         strd/3600,
    #                                                         strr/3600,
    #                                                         cossza/3600)
    #     print(mrt)

    # def test_utci(self):
    #      t2mk = 309.0
    #      va = 3
    #      mrt = 310.0
    #      e_hPa = 12
    #      utci= thermofeel.calculate_utci(t2mk,va,mrt,e_hPa)
    #      print(utci)

    def test_apparent_temperature(self):
        t2mk = thermofeel.celcius_to_kelvin(25.0)
        va = 3
        at = thermofeel.calculate_apparent_temperature(t2mk, va)
        # print(f"at {at}")
        assert at == pytest.approx(27.844072534, abs=1e-6)

    def test_wbgts(self):
        t2mk = thermofeel.celcius_to_kelvin(30.0)
        wbgts = thermofeel.calculate_wbgts(t2mk)
        # print(f"wbgts {wbgts}")
        assert wbgts == pytest.approx(22.05908916, abs=1e-6)

    def test_net(self):
        t2mk = 307
        tdk = 299
        va = 4
        net = thermofeel.calculate_net_effective_temperature(t2mk, va, tdk)
        # print(f"net {net}")
        assert net == pytest.approx(36.0638063, abs=1e-6)

    def test_humidex(self):
        t2mk = 304
        tdk = 300
        hu = thermofeel.calculate_humidex(t2mk, tdk)
        # print(f"hu {hu}")
        assert hu == pytest.approx(27.97125748479999, abs=1e-6)

    def test_windchill(self):
        t2mk = 290
        va = 10
        wc = thermofeel.calculate_wind_chill(t2mk, va)
        assert wc == pytest.approx(23.78881163916766, abs=1e-6)

    def test_heat_index_adjusted(self):
        t2mk = np.array([295])
        tdk = np.array([290])
        hia = thermofeel.calculate_heat_index_adjusted(t2mk, tdk)
        # print(f"hia {hia}")
        assert hia[0] == pytest.approx(22.00355699, abs=1e-6)

    def test_wbgt(self):
        t2mk = 301
        va = 4
        mrt = 310
        wbgt = thermofeel.calculate_wbgt(t2mk, va, mrt)
        # print(f"wbgt {wbgt}")
        assert wbgt == pytest.approx(26.769875412856436, abs=1e-6)

    def test_calculate_cos_solar_zenith_angle(self):
        # should return ~ 0.360303587797559
        cossza = thermofeel.calculate_cos_solar_zenith_angle(
            lat=48.81667, lon=2.28972, d=15, m=11, y=2006, h=10.58333
        )
        assert abs(cossza - 0.360303587797559) < 1e-6

        # London, ~ 0.8799471697555967
        cossza = thermofeel.calculate_cos_solar_zenith_angle(
            lat=51.0, lon=0.0, d=4, m=6, y=2021, h=12.0
        )
        assert abs(cossza - 0.8799471697555967) < 1e-6

        # from alternative formula
        assert abs(cossza - cos(radians(90.0 - 61.5))) < 1e-2

    def test_calculate_cos_solar_zenith_angle_integrated(self):
        cossza = thermofeel.calculate_cos_solar_zenith_angle_integrated(
            lat=48.81667, lon=2.28972, d=15, m=11, y=2006, h=10.58333, base=0, step=3
        )
        assert abs(cossza - 0.34495713937581207) < 1e-6

    def test_solar_declination_angle(self):
        sda1, tc1 = thermofeel.solar_declination_angle(jd=166, h=0)
        assert abs(sda1 - 23.32607701732299) < 1e-6
        assert abs(tc1 - -0.054061457069008334) < 1e-6

        sda2, tc2 = thermofeel.solar_declination_angle(jd=4, h=12)
        assert abs(sda2 - -22.64240042915207) < 1e-6
        assert abs(tc2 - -1.219397058249299) < 1e-6

        sda3, tc3 = thermofeel.solar_declination_angle(jd=600, h=3)
        assert abs(sda3 - 11.471993171760428) < 1e-6
        assert abs(tc3 - -0.7161824119549858) < 1e-6


if __name__ == "__main__":
    unittest.main()  # pragma: no cover

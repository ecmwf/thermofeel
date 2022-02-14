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

import thermofeel as tmf


class TestThermalCalculator(unittest.TestCase):
    def test_relative_humidity_percent(self):
        t2mk = tmf.celsius_to_kelvin(30.0)
        t2dk = tmf.celsius_to_kelvin(28.0)

        rhpc = tmf.calculate_relative_humidity_percent(t2mk, t2dk)

        assert rhpc == pytest.approx(89.08526710467393, abs=1e-6)

    def test_saturation_vapour_pressure(self):
        t2mk = tmf.celsius_to_kelvin(25.0)
        svp = tmf.calculate_saturation_vapour_pressure(t2mk)
        assert svp == pytest.approx(31.699201897293, abs=1e-6)

    def test_heat_index(self):
        t2mk = tmf.celsius_to_kelvin(25.0)
        hi = tmf.calculate_heat_index_simplified(t2mk)
        # print(f"hi {hi}")
        assert hi == pytest.approx(49.321122555, abs=1e-6)

    def test_meant_radiant_temperature(self):
        ssrd = np.array([60000])
        ssr = np.array([471818])
        fdir = np.array([374150])
        strd = np.array([1061213])
        strr = np.array([-182697])
        cossza = np.array([0.4])
        mrt = tmf.calculate_mean_radiant_temperature(
            ssrd / 3600,
            ssr / 3600,
            fdir / 3600,
            strd / 3600,
            strr / 3600,
            cossza / 3600,
        )
        # print(f"mrt {mrt}")
        assert mrt == pytest.approx(262.81089323, abs=1e-5)

    def test_utci(self):
        t2mk = np.array([309.0])
        va = np.array([3])
        mrt = np.array([310.0])
        e_hPa = np.array([12])
        utci = tmf.calculate_utci(t2mk, va, mrt, e_hPa)
        assert utci == pytest.approx(34.61530078, abs=1e-5)

    def test_mrt_from_bgt(self):
        t_k = np.array([tmf.celsius_to_kelvin(25.0)])
        bgt_k = np.array([tmf.celsius_to_kelvin(23.0)])
        # print(f"bgt_k {bgt_k}")
        va = np.array([10])
        mrt_c = tmf.calculate_mrt_from_bgt(t_k, bgt_k, va)
        # print(f"mrt_c {mrt_c}")
        bgt2_c = tmf.calculate_bgt(t_k, tmf.celsius_to_kelvin(mrt_c), va)

        assert tmf.kelvin_to_celsius(bgt_k) == pytest.approx(bgt2_c, abs=1e-6)

    def test_apparent_temperature(self):
        t2mk = tmf.celsius_to_kelvin(25.0)
        va = 3
        at = tmf.calculate_apparent_temperature(t2mk, va)
        # print(f"at {at}")
        assert at == pytest.approx(27.844072534, abs=1e-6)

    def test_wbgts(self):
        t2mk = tmf.celsius_to_kelvin(30.0)
        wbgts = tmf.calculate_wbgts(t2mk)
        # print(f"wbgts {wbgts}")
        assert wbgts == pytest.approx(22.05908916, abs=1e-6)

    def test_net(self):
        t2mk = 307
        tdk = 299
        va = 4
        net = tmf.calculate_net_effective_temperature(t2mk, va, tdk)
        # print(f"net {net}")
        assert net == pytest.approx(36.0638063, abs=1e-6)

    def test_humidex(self):
        t2mk = 304
        tdk = 300
        hu = tmf.calculate_humidex(t2mk, tdk)
        # print(f"hu {hu}")
        assert hu == pytest.approx(27.97125748479999, abs=1e-6)

    def test_windchill(self):
        t2mk = 290
        va = 10
        wc = tmf.calculate_wind_chill(t2mk, va)
        assert wc == pytest.approx(23.78881163916766, abs=1e-6)

    def test_heat_index_adjusted(self):
        t2mk = np.array([295])
        tdk = np.array([290])
        hia = tmf.calculate_heat_index_adjusted(t2mk, tdk)
        # print(f"hia {hia}")
        assert hia[0] == pytest.approx(22.00355699, abs=1e-6)

    def test_wbt(self):
        t_c = 20
        rh = 50
        wbt = tmf.calculate_wbt(t_c, rh)
        # print(f"wbt {wbt}")
        assert wbt == pytest.approx(13.6993419, abs=1e-6)  # validated with the article

    def test_bgt(self):
        t_k = np.array([278.15, 300, 300])
        va = np.array([20, 20, -10])  # negative va values are treated as 0
        mrt = np.array([278.15, 310, 310])
        bgt = tmf.calculate_bgt(t_k, va, mrt)
        # print(f"bgt {bgt}")
        assert bgt[0] == pytest.approx(4.1575627528, abs=1e-6)
        assert bgt[1] == pytest.approx(25.78402849, abs=1e-6)
        assert bgt[2] == pytest.approx(25.78400875, abs=1e-6)

    def test_wbgt(self):
        t_k = np.array([300])
        td_k = np.array([290])
        va = np.array([20])
        mrt = np.array([310])
        wbgt = tmf.calculate_wbgt(t_k, va, mrt, td_k)
        # print(f"wbgt {wbgt}")
        assert wbgt[0] == pytest.approx(22.0381904925, abs=1e-6)

        # test negative values are treated as 0
        va[0] = -10
        wbgt = tmf.calculate_wbgt(t_k, va, mrt, td_k)
        # print(f"wbgt {wbgt}")
        assert wbgt[0] == pytest.approx(22.03818628, abs=1e-6)

    def test_wbt_dj(self):
        t2m = np.array([295])
        td = np.array([290])
        p = np.array([1000])
        wbtdj = tmf.calculate_wbt_dj(t2k=t2m, p=p, tdk=td)

        assert wbtdj == pytest.approx(45.114, abs=1e-3)

    def test_calculate_cos_solar_zenith_angle(self):
        # should return ~ 0.360303587797559
        cossza = tmf.calculate_cos_solar_zenith_angle(
            lat=48.81667, lon=2.28972, d=15, m=11, y=2006, h=10.58333
        )
        # print(f"cossza {cossza}")
        assert cossza == pytest.approx(0.360303587797559, abs=1e-6)

        # London, ~ 0.8799471697555967
        cossza = tmf.calculate_cos_solar_zenith_angle(
            lat=51.0, lon=0.0, d=4, m=6, y=2021, h=12.0
        )
        # print(f"cossza {cossza}")
        assert cossza == pytest.approx(0.8799471697555967, abs=1e-6)
        # from alternative formula
        assert cossza == pytest.approx(cos(radians(90.0 - 61.5)), abs=1e-2)

    def test_calculate_cos_solar_zenith_angle_integrated(self):
        lat = 48.81667
        lon = 2.28972
        d = 15
        m = 11
        y = 2006
        h = 10.58333
        tbegin = 0
        tend = 3
        cossza = tmf.calculate_cos_solar_zenith_angle_integrated(
            lat, lon, y, m, d, h, tbegin, tend
        )
        # print(f"cossza {cossza}")
        assert cossza == pytest.approx(0.3612630470539099, abs=1e-6)

        # opposite point in the world should be dark
        lat = -lat
        lon = 180 + lon
        cossza = tmf.calculate_cos_solar_zenith_angle_integrated(
            lat, lon, y, m, d, h, tbegin, tend
        )
        # print(f"cossza {cossza}")
        assert cossza == pytest.approx(0.0, abs=1e-6)

        # lons can be > 360
        lat = -lat
        lon = 180 + lon
        cossza = tmf.calculate_cos_solar_zenith_angle_integrated(
            lat, lon, y, m, d, h, tbegin, tend
        )
        # print(f"cossza {cossza}")
        assert cossza == pytest.approx(0.3612630470539099, abs=1e-6)

        lat = 48.81667
        lon = 2.28972

        # integration with splits every 20min (3 per hour)
        cossza = tmf.calculate_cos_solar_zenith_angle_integrated(
            lat, lon, y, m, d, h, tbegin, tend, intervals_per_hour=3
        )
        # print(f"cossza {cossza}")
        assert cossza == pytest.approx(0.3612630469576353, abs=1e-7)

        # gauss integration order 2
        cossza = tmf.calculate_cos_solar_zenith_angle_integrated(
            lat, lon, y, m, d, h, tbegin, tend, integration_order=2
        )
        # print(f"cossza {cossza}")
        assert cossza == pytest.approx(0.3612623904213413, abs=1e-7)

        # gauss integration order 1
        cossza = tmf.calculate_cos_solar_zenith_angle_integrated(
            lat, lon, y, m, d, h, tbegin, tend, integration_order=1
        )
        # print(f"cossza {cossza}")
        assert cossza == pytest.approx(0.36298755581259323, abs=1e-6)

        # gauss integration order 4
        cossza = tmf.calculate_cos_solar_zenith_angle_integrated(
            lat, lon, y, m, d, h, tbegin, tend, integration_order=4
        )
        # print(f"cossza {cossza}")
        assert cossza == pytest.approx(0.36126304695749595, abs=1e-7)

    def test_solar_declination_angle(self):
        sda, tc = tmf.solar_declination_angle(jd=166, h=0)
        assert sda == pytest.approx(23.32607701732299, abs=1e-6)
        assert tc == pytest.approx(-0.054061457069008334, abs=1e-6)

        sda, tc = tmf.solar_declination_angle(jd=4, h=12)
        assert sda == pytest.approx(-22.64240042915207, abs=1e-6)
        assert tc == pytest.approx(-1.219397058249299, abs=1e-6)

        sda, tc = tmf.solar_declination_angle(jd=600, h=3)
        assert sda == pytest.approx(11.471993171760428, abs=1e-6)
        assert tc == pytest.approx(-0.7161824119549858, abs=1e-6)


if __name__ == "__main__":
    unittest.main()  # pragma: no cover

# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import unittest

import numpy as np
import pytest

import thermofeel as tmf


class TestThermalCalculator(unittest.TestCase):
    def test_relative_humidity_percent(self):
        t2_k = np.array([tmf.celsius_to_kelvin(30.0)])
        td_k = np.array([tmf.celsius_to_kelvin(28.0)])
        rhpc = np.array([tmf.calculate_relative_humidity_percent(t2_k, td_k)])
        assert rhpc == pytest.approx(89.08526710467393, abs=1e-6)

    def test_saturation_vapour_pressure(self):
        t2_k = np.array([tmf.celsius_to_kelvin(25.0)])
        svp = np.array([tmf.calculate_saturation_vapour_pressure(t2_k)])
        assert svp == pytest.approx(31.699201897293, abs=1e-6)

    def test_saturation_vapour_pressure_multiphase(self):
        t2_k = np.array([tmf.celsius_to_kelvin(-25.0)])
        phase = np.array([1])
        es = np.array(
            [tmf.calculate_saturation_vapour_pressure_multiphase(t2_k, phase)]
        )
        # assert es == pytest.approx(0.63555512, abs=1e-6) # old formula
        assert es == pytest.approx(0.63142553, abs=1e-6)

    def test_nonsaturation_vapour_pressure(self):
        t2_k = np.array([300])
        rh = np.array([87])
        svp = np.array([tmf.calculate_nonsaturation_vapour_pressure(t2_k, rh)])
        assert svp == pytest.approx(30.649976725404283, abs=1e-6)

    def test_scale_windspeed(self):
        va = np.array([7.0])
        h = np.array([2.0])
        vh = np.array([tmf.scale_windspeed(va, h)])
        assert vh == pytest.approx(5.369069989882623, abs=1e-6)

    def test_dew_point_from_relative_humidity(self):
        rh = np.array([56])
        t2_k = np.array([304.15])
        td_k = np.array([tmf.calculate_dew_point_from_relative_humidity(rh, t2_k)])
        assert td_k == pytest.approx(294.3484414118635, abs=1e-6)

    def test_mean_radiant_temperature(self):
        ssrd = np.array([60000])
        ssr = np.array([471818])
        fdir = np.array([374150])
        strd = np.array([1061213])
        strr = np.array([-182697])
        cossza = np.array([0.4])
        dsrp = np.array([tmf.approximate_dsrp(fdir, cossza)])
        mrt = np.array(
            [
                tmf.calculate_mean_radiant_temperature(
                    ssrd=ssrd / 3600,
                    ssr=ssr / 3600,
                    fdir=fdir / 3600,
                    strd=strd / 3600,
                    strr=strr / 3600,
                    cossza=cossza / 3600,
                    dsrp=dsrp / 3600,
                )
            ]
        )
        # print(f"mrt {mrt}")
        assert mrt == pytest.approx(270.85099123, abs=1e-6)

    def test_utci(self):
        # case 1
        t2_k = np.array([309.0])
        va = np.array([3])
        mrt = np.array([310.0])
        e_hPa = np.array([12])
        utci = np.array([tmf.calculate_utci(t2_k, va, mrt, td_k=None, ehPa=e_hPa)])
        assert utci == pytest.approx(307.76473586, abs=1e-5)
        # case 2
        t2_k = np.array([tmf.celsius_to_kelvin(27.0)])
        va = np.array([4])
        mrt = np.array([tmf.celsius_to_kelvin(9.2)])
        e_hPa = np.array([16.5])
        utci = np.array([tmf.calculate_utci(t2_k, va, mrt, td_k=None, ehPa=e_hPa)])
        assert tmf.kelvin_to_celsius(utci) == pytest.approx(18.93148565062157, abs=1e-5)

    def test_wbgt_simple(self):
        t2_k = np.array([tmf.celsius_to_kelvin(30.0)])
        rh = np.array([80])
        wbgts = np.array([tmf.calculate_wbgt_simple(t2_k, rh)])
        # print(f"wbgts {wbgts}")
        assert wbgts == pytest.approx(307.39508355517813, abs=1e-6)

    def test_wbt(self):
        t2_c = np.array([tmf.celsius_to_kelvin(20.0)])
        rh = np.array([50])
        wbt = np.array([tmf.calculate_wbt(t2_c, rh)])
        # print(f"wbt {wbt}")
        assert wbt == pytest.approx(286.84934189999996, abs=1e-6)

    def test_bgt(self):
        t2_k = np.array([278.15, 300, 300])
        va = np.array([20, 20, -10])  # negative va values are treated as 0
        mrt = np.array([278.15, 310, 310])
        bgt = np.array([tmf.calculate_bgt(t2_k, va, mrt)])
        # print(f"bgt {bgt}")
        assert bgt[0, 0] == pytest.approx(277.1238737724192, abs=1e-6)
        assert bgt[0, 1] == pytest.approx(298.70218703427656, abs=1e-6)
        assert bgt[0, 2] == pytest.approx(298.70216299754475, abs=1e-6)

    def test_wbgt(self):
        t2_k = np.array([300])
        td_k = np.array([290])
        va = np.array([20])
        mrt = np.array([310])
        wbgt = np.array([tmf.calculate_wbgt(t2_k, mrt, va, td_k)])
        # print(f"wbgt {wbgt}")
        assert wbgt[0] == pytest.approx(295.5769818634555, abs=1e-6)
        # # test negative values are treated as 0
        # va[0] = -10
        # wbgt = np.array([tmf.calculate_wbgt(t_k, mrt, va, td_k)])
        # # print(f"wbgt {wbgt}")
        # assert wbgt[0] == pytest.approx(295.5769818634555, abs=1e-6)

    def test_mrt_from_bgt(self):
        t2_k = np.array([tmf.celsius_to_kelvin(25.0)])
        bgt_k = np.array([tmf.celsius_to_kelvin(23.0)])
        # print(f"bgt_k {bgt_k}")
        va = np.array([10])
        mrt_c = np.array([tmf.calculate_mrt_from_bgt(t2_k, bgt_k, va)])
        assert mrt_c == pytest.approx(279.80189775556704, abs=1e-6)

    def test_humidex(self):
        t2_k = np.array([304])
        td_k = np.array([300])
        hu = np.array([tmf.calculate_humidex(t2_k, td_k)])
        # print(f"hu {hu}")
        assert hu == pytest.approx(318.4601286141123, abs=1e-6)

    def test_normal_effective_temperature(self):
        t2_k = np.array([307])
        va = np.array([4])
        rh = np.array([80])
        net = np.array([tmf.calculate_normal_effective_temperature(t2_k, va, rh)])
        # print(f"net {net}")
        assert net == pytest.approx(304.13650125, abs=1e-6)

    def test_apparent_temperature(self):
        t2_k = np.array([tmf.celsius_to_kelvin(25.0)])
        va = np.array([3])
        rh = np.array([75])
        at = np.array([tmf.calculate_apparent_temperature(t2_k, va, rh)])
        # print(f"at {at}")
        assert at == pytest.approx(299.86678322384626, abs=1e-6)

    def test_wind_chill(self):
        t2_k = np.array([270])
        va = np.array([10])
        wc_k = np.array([tmf.calculate_wind_chill(t2_k, va)])
        assert wc_k == pytest.approx(261.92338925380074, abs=1e-6)
        # reference result from from wikipedia article https://en.wikipedia.org/wiki/Wind_chill
        t2_k = np.array([tmf.celsius_to_kelvin(-20)])  # -20C to K
        va = np.array([5 / 3.6])  # 5 km/h to m/s
        wc_k = np.array([tmf.calculate_wind_chill(t2_k, va)])
        wc_c = tmf.kelvin_to_celsius(wc_k)
        assert wc_c == pytest.approx(
            -24.27850328, abs=1e-6
        )  # around ~ -24C for wind chill (not exact)
        va = np.array([30 / 3.6])  # 30 km/h to m/s
        wc_k = np.array([tmf.calculate_wind_chill(t2_k, va)])
        wc_c = tmf.kelvin_to_celsius(wc_k)
        print(f"wc_c {wc_c}")
        assert wc_c == pytest.approx(
            -32.56804448, abs=1e-6
        )  # around ~ -33C for wind chill (not exact)

    def test_heat_index_simplified(self):
        t2_k = np.array([tmf.celsius_to_kelvin(21.0)])
        rh = np.array([80])
        hi = np.array([tmf.calculate_heat_index_simplified(t2_k, rh)])
        # print(f"hi {hi}")
        assert hi == pytest.approx(294.68866082, abs=1e-6)

    def test_heat_index_adjusted(self):
        t2_k = np.array([295])
        td_k = np.array([290])
        hia = np.array([tmf.calculate_heat_index_adjusted(t2_k, td_k)])
        # print(f"hia {hia}")
        assert hia[0] == pytest.approx(295.15355699, abs=1e-6)


if __name__ == "__main__":
    unittest.main()  # pragma: no cover

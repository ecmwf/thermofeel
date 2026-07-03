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
        # signature is calculate_bgt(t2_k, mrt, va)
        t2_k = np.array([278.15, 300.0])
        mrt = np.array([278.15, 310.0])
        va = np.array([20.0, 20.0])
        bgt = np.array([tmf.calculate_bgt(t2_k, mrt, va)])
        # print(f"bgt {bgt}")
        # when mrt == t2 (no net radiation) the globe temperature equals t2
        assert bgt[0, 0] == pytest.approx(278.15, abs=1e-6)
        assert bgt[0, 1] == pytest.approx(300.877985349940, abs=1e-6)

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

    def test_wind_speed_2m_liljegren(self):
        # Hand-computed from the PyWBGT stability table and exponents:
        # va * (2/10)**urban_exp[class-1], floored at 0.13 m/s.
        # va=3, day, ssrd=800 -> col=1,row=2,class=2,exp=0.15
        # va=0.62, day, ssrd=800 -> col=1,row=0,class=1,exp=0.15
        # va=1, night, ssrd=0 -> col=5,row=0,class=5,exp=0.30
        # va=7, day, ssrd=200 -> col=2,row=4,class=4,exp=0.25
        va = np.array([3.0, 0.62, 1.0, 7.0])
        cossza = np.array([0.9, 0.9, -0.5, 0.5])
        ssrd = np.array([800.0, 800.0, 0.0, 200.0])
        v2 = tmf.calculate_wind_speed_2m_liljegren(va, cossza, ssrd)
        expected = np.array([2.356545, 0.487019, 0.617034, 4.681182])
        np.testing.assert_array_almost_equal(v2, expected, decimal=6)

    def test_wbgt_liljegren(self):
        # Reference values generated from Liljegren's reference C implementation
        # (MIT-licensed mirror github.com/mdljts/wbgt), called with the same
        # preprocessing thermofeel applies (0.62 m/s 10 m wind floor, KNMI
        # Liljegren 10->2 m stability scaling, fdir clamped to [0, 0.9] and 0
        # below 89.5 deg zenith).
        # Inputs: t2_k, rh[%], pressure[hPa], va_10m[m/s], ssrd[W/m2], fdir, cossza
        cases = [
            ((308.15, 40.0, 1013.0, 1.0, 800.0, 0.70, 0.90), 306.559802),
            ((303.15, 70.0, 1010.0, 3.0, 500.0, 0.50, 0.70), 302.954374),
            ((293.15, 60.0, 1015.0, 5.0, 200.0, 0.30, 0.50), 290.823880),
            ((298.15, 85.0, 1008.0, 2.0, 0.0, 0.00, 0.05), 296.597705),
            # va below the 0.62 m/s floor -> clamped
            ((306.15, 50.0, 1013.0, 0.2, 700.0, 0.60, 0.80), 307.358482),
            # fdir above 0.9 -> clamped
            ((305.15, 35.0, 1012.0, 2.5, 900.0, 0.97, 0.95), 301.229074),
        ]
        for args, expected_k in cases:
            arr = [np.array([v]) for v in args]
            wbgt_k = np.array([tmf.calculate_wbgt_liljegren(*arr)])
            assert wbgt_k[0] == pytest.approx(expected_k, abs=1e-4)

    def test_wbgt_liljegren_brode_option(self):
        # The "brode" wind-scaling option uses the generic scale_windspeed log
        # profile instead of the KNMI stability profile, giving a different
        # (validated) value for the same inputs.
        args = [np.array([v]) for v in (308.15, 40.0, 1013.0, 1.0, 800.0, 0.70, 0.90)]
        wbgt_brode = tmf.calculate_wbgt_liljegren(*args, wind_scaling="brode")
        assert wbgt_brode[0] == pytest.approx(306.619995, abs=1e-4)

    def test_heat_force(self):
        # Lower-closed 2 degC WBGT bands: <14 -> 0, [14,16) -> 1, ..., >=32 -> 10
        wbgt_c = np.array([10.0, 13.999, 14.0, 15.9, 16.0, 24.0, 31.999, 32.0, 35.0])
        wbgt_k = tmf.celsius_to_kelvin(wbgt_c)
        hf = tmf.calculate_heat_force(wbgt_k)
        expected = np.array([0, 0, 1, 1, 2, 6, 9, 10, 10])
        np.testing.assert_array_equal(hf, expected)

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
        assert hu == pytest.approx(318.47466703, abs=1e-6)

    def test_discomfort_index(self):
        # Reference values derived analytically from Thom's relative-humidity
        # form DI = T - 0.55 (1 - 0.01 RH)(T - 14.5) [T in degC], returned in K.
        # The last two are exact identities of the formula (independent of the
        # implementation): RH=100% -> DI equals air temperature; T=14.5 degC ->
        # DI = 14.5 degC for any RH.
        cases = [
            (30.0, 70.0, 300.5925),
            (40.0, 40.0, 304.735),
            (25.0, 50.0, 295.2625),
            (28.0, 100.0, 301.15),
            (14.5, 35.0, 287.65),
        ]
        for t_c, rh_pc, expected_k in cases:
            t2_k = np.array([tmf.celsius_to_kelvin(t_c)])
            rh = np.array([rh_pc])
            di = tmf.calculate_discomfort_index(t2_k, rh)
            assert di[0] == pytest.approx(expected_k, abs=1e-6)

    def test_summer_simmer_index(self):
        # PRIMARY gate (G2 affine-THI identity). By definition the common 1987
        # SSI is an affine transform of Thom's Fahrenheit Temperature-Humidity
        # Index: SSI_F = 1.98*THI_F - 56.83 with THI_F = Tf - (0.55 - 0.0055*RH)
        # (Tf - 58). Rebuild that whole chain independently and require the
        # function to reproduce it exactly. This proves the implemented formula
        # matches its stated definition (only the *provenance* of that formula is
        # secondary-source; the maths is pinned here).
        identity_points = [
            (33.0, 55.0),
            (18.0, 90.0),
            (37.5, 15.0),
            (10.0, 100.0),
            (28.2, 47.5),
        ]
        for t_c, rh_pc in identity_points:
            t2_k = np.array([tmf.celsius_to_kelvin(t_c)])
            rh = np.array([rh_pc])
            tf = tmf.kelvin_to_fahrenheit(t2_k)
            thi_f = tf - (0.55 - 0.0055 * rh) * (tf - 58.0)
            expected_ssi_f = 1.98 * thi_f - 56.83
            expected_k = tmf.fahrenheit_to_kelvin(expected_ssi_f)
            ssi = tmf.calculate_summer_simmer_index(t2_k, rh)
            assert ssi[0] == pytest.approx(expected_k[0], abs=1e-6)

        # Hand-derived exact evaluations. Every Tf is an integer so the whole
        # chain is checkable by hand. Notation: coeff = 0.55 - 0.0055*RH,
        # THI_F = Tf - coeff*(Tf - 58), SSI_F = 1.98*THI_F - 56.83,
        # SSI_K = (SSI_F + 459.67)*5/9.
        cases = [
            # T=30C -> Tf=86; coeff=0.55-0.275=0.275; THI_F=86-0.275*28=78.3;
            # SSI_F=1.98*78.3-56.83=98.204; SSI_K=(98.204+459.67)*5/9=309.930.
            # (Task worked example quotes SSI_F~=98.20 -> SSI_K~=309.928; the
            # unrounded value is 309.930 and both sit inside the 1e-2 gate.)
            (30.0, 50.0, 309.930),
            # T=25C -> Tf=77; coeff=0.55-0.22=0.33; THI_F=77-0.33*19=70.73;
            # SSI_F=1.98*70.73-56.83=83.2154; SSI_K=(83.2154+459.67)*5/9=301.603.
            (25.0, 40.0, 301.603),
            # T=35C -> Tf=95; coeff=0.55-0.33=0.22; THI_F=95-0.22*37=86.86;
            # SSI_F=1.98*86.86-56.83=115.1528; SSI_K=(115.1528+459.67)*5/9=319.346.
            (35.0, 60.0, 319.346),
            # T=20C -> Tf=68; coeff=0.55-0.44=0.11; THI_F=68-0.11*10=66.9;
            # SSI_F=1.98*66.9-56.83=75.632; SSI_K=(75.632+459.67)*5/9=297.390.
            (20.0, 80.0, 297.390),
            # T=15C, RH=100 -> Tf=59; coeff=0 so THI_F=Tf=59 (RH=100 identity);
            # SSI_F=1.98*59-56.83=59.99; SSI_K=(59.99+459.67)*5/9=288.700.
            (15.0, 100.0, 288.700),
            # T=40C -> Tf=104; coeff=0.55-0.11=0.44; THI_F=104-0.44*46=83.76;
            # SSI_F=1.98*83.76-56.83=109.0148; SSI_K=(109.0148+459.67)*5/9=315.936.
            (40.0, 20.0, 315.936),
        ]
        for t_c, rh_pc, expected_k in cases:
            t2_k = np.array([tmf.celsius_to_kelvin(t_c)])
            rh = np.array([rh_pc])
            ssi = tmf.calculate_summer_simmer_index(t2_k, rh)
            assert ssi[0] == pytest.approx(expected_k, abs=1e-2)

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

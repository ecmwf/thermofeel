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

    def test_pmv(self):
        # G3 - ISO 7730:2005 Annex D Table D.1: (ta, tr, v, rh, met, clo) -> PMV.
        # Every row was re-derived from the Annex D fixed-point algorithm and
        # cross-checked against the MIT oracle pythermalcomfort.pmv_ppd_iso
        # (7730-2005 model). 12 of the 13 rows reproduce the printed PMV to
        # <0.008. Row 7 (23.5, 23.5, 0.1, 40 %, 1.2 met, 1.0 clo) is a known
        # ISO Table D.1 misprint: the standard prints PMV 0.50 / PPD 10, but the
        # Annex D algorithm AND the oracle both give PMV 0.36 / PPD 7.7. The
        # formula-reconciling value is pinned (the same policy the programme
        # applied to the non-reconciling apparent-temperature dossier row), so
        # the tolerance stays uniform at 0.02 with no per-row weakening.
        ta = np.array([22, 27, 27, 23.5, 23.5, 19, 23.5, 23.5, 23, 23, 22, 27, 27])
        tr = np.array([22, 27, 27, 25.5, 25.5, 19, 23.5, 23.5, 21, 21, 22, 27, 27])
        v = np.array([0.1, 0.1, 0.3, 0.1, 0.3, 0.1, 0.1, 0.3, 0.1, 0.3, 0.1, 0.1, 0.3])
        rh = np.array([60, 60, 60, 60, 60, 40, 40, 40, 40, 40, 60, 60, 60])
        met = np.array(
            [1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.6, 1.6, 1.6]
        )
        clo = np.array(
            [0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5]
        )
        # Row 7 uses the algorithm/oracle value (0.36), not the misprint (0.50).
        pmv_expected = np.array(
            [
                -0.75,
                0.77,
                0.44,
                -0.01,
                -0.55,
                -0.60,
                0.36,
                0.12,
                0.05,
                -0.16,
                0.05,
                1.17,
                0.95,
            ]
        )
        # Row 7 PPD pinned to the reconciling 8 (not the misprinted 10).
        ppd_expected = np.array([17, 17, 9, 5, 11, 13, 8, 5, 5, 6, 5, 34, 24])
        pmv = tmf.calculate_pmv(
            tmf.celsius_to_kelvin(ta),
            tmf.celsius_to_kelvin(tr),
            v,
            rh=rh,
            met=met,
            clo=clo,
        )
        assert not np.isnan(pmv).any()  # every ISO row converges
        np.testing.assert_allclose(pmv, pmv_expected, atol=0.02)
        ppd = tmf.calculate_ppd(pmv)
        np.testing.assert_allclose(ppd, ppd_expected, atol=0.5)

    def test_ppd(self):
        # G2 identities: PPD is exactly 5% at thermal neutrality (PMV = 0) and is
        # an even function of PMV (equal warm/cold deviations give equal PPD).
        assert tmf.calculate_ppd(0.0) == 5.0
        x = np.array([-2.4, -1.3, -0.5, 0.7, 1.6, 2.9])
        np.testing.assert_array_equal(tmf.calculate_ppd(x), tmf.calculate_ppd(-x))
        # Reference points on the ISO 7730 PPD(PMV) curve (PMV +-0.5 -> 10 %,
        # +-1 -> 26 %, +-2 -> 77 %).
        pmv = np.array([-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0])
        expected = np.array([76.8, 26.1, 10.2, 5.0, 10.2, 26.1, 76.8])
        np.testing.assert_allclose(tmf.calculate_ppd(pmv), expected, atol=0.1)

    def test_pmv_humidity_contract(self):
        # Exactly one humidity input is required (mirrors calculate_utci).
        t = tmf.celsius_to_kelvin(np.array([22.0]))
        v = np.array([0.1])
        with pytest.raises(ValueError):
            tmf.calculate_pmv(t, t, v)  # neither rh nor vapour_pressure_hpa
        with pytest.raises(ValueError):
            tmf.calculate_pmv(
                t, t, v, rh=np.array([60.0]), vapour_pressure_hpa=np.array([14.0])
            )  # both
        # Passing vapour pressure (converted from an RH row with the same ISO
        # relation the function uses internally) reproduces the rh-driven PMV.
        ta_c, rh = 22.0, 60.0
        pa_pa = rh * 10.0 * np.exp(16.6536 - 4030.183 / (ta_c + 235.0))
        vp_hpa = np.array([pa_pa / 100.0])
        pmv_rh = tmf.calculate_pmv(t, t, v, rh=np.array([rh]))
        pmv_vp = tmf.calculate_pmv(t, t, v, vapour_pressure_hpa=vp_hpa)
        np.testing.assert_allclose(pmv_vp, pmv_rh)

    def test_pmv_nonconvergence_returns_nan(self):
        # A NaN element never satisfies the convergence test, so the vectorised
        # fixed point runs to its iteration cap and returns NaN for it (the
        # documented non-convergence behaviour), while a finite neighbour still
        # converges to a real value.
        t2 = tmf.celsius_to_kelvin(np.array([22.0, np.nan]))
        tr = tmf.celsius_to_kelvin(np.array([22.0, 22.0]))
        var = np.array([0.1, 0.1])
        rh = np.array([60.0, 60.0])
        pmv = tmf.calculate_pmv(t2, tr, var, rh=rh)
        assert np.isfinite(pmv[0])
        assert np.isnan(pmv[1])


if __name__ == "__main__":
    unittest.main()  # pragma: no cover

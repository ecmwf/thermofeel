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

    def test_relative_strain_index(self):
        # G2 analytic identity: at Ta = 21 degC the numerator (Ta - 21) is exactly
        # 0, so RSI == 0 for any relative humidity (independent of the vapour
        # pressure e). Encoded for several rh values.
        for rh_pc in (10.0, 30.0, 50.0, 73.74, 90.0, 100.0):
            t2_k = np.array([tmf.celsius_to_kelvin(21.0)])
            rh = np.array([rh_pc])
            rsi = tmf.calculate_relative_strain_index(t2_k, rh)
            assert rsi[0] == 0.0

        # G3 external reference: Asghari et al. (2020) Table 4, 15-yr summer
        # monthly means (DOI 10.2174/1874213002013010011). Each row is re-derived
        # here with thermofeel's non-saturation vapour pressure e(Ta, rh) via the
        # hPa closed form RSI = (Ta - 21) / (58 - e). computed vs source (delta):
        #   (34.27, 66.50) -> 0.59710 vs 0.600  (-0.00290), e = 35.7758 hPa
        #   (28.12, 22.95) -> 0.14444 vs 0.151  (-0.00656), e =  8.7074 hPa
        #   (27.54, 73.70) -> 0.21120 vs 0.210  (+0.00120), e = 27.0339 hPa
        #   (33.73, 60.30) -> 0.48003 vs 0.480  (+0.00003), e = 31.4807 hPa
        #   (27.00, 73.74) -> 0.18873 vs 0.190  (-0.00127), e = 26.2079 hPa
        # All residuals are within abs=0.02 (max |delta| = 0.0066, mixed signs so
        # no systematic bias); the spread reflects source rounding (2-3 dp) plus a
        # mean-of-RSI vs RSI-of-means (Jensen) offset, not a formula mismatch.
        cases = [
            (34.27, 66.5, 0.60),
            (28.12, 22.95, 0.151),
            (27.54, 73.7, 0.21),
            (33.73, 60.3, 0.48),
            (27.00, 73.74, 0.19),
        ]
        for t_c, rh_pc, expected in cases:
            t2_k = np.array([tmf.celsius_to_kelvin(t_c)])
            rh = np.array([rh_pc])
            rsi = tmf.calculate_relative_strain_index(t2_k, rh)
            assert rsi[0] == pytest.approx(expected, abs=0.02)

    def test_pet_reference_environment_identity(self):
        # PRIMARY gate. PET is DEFINED as the air temperature of the reference
        # environment (mrt = ta, v = 0.1 m/s, vpa = 12 hPa, clo = 0.9,
        # m_act = 80 W) reproducing the body state. Feeding the reference
        # environment to itself must therefore return the air temperature
        # exactly. This single property pins the reference environment, the
        # clothing convention AND the whole-body activity-unit convention at
        # once, so it is the highest-value test for this function.
        t_c = np.arange(-10.0, 41.0, 2.5)
        t2_k = tmf.celsius_to_kelvin(t_c)
        pet = tmf.calculate_pet(
            t2_k,
            t2_k,
            np.full_like(t_c, 0.1),
            vapour_pressure_hpa=np.full_like(t_c, 12.0),
        )
        np.testing.assert_allclose(pet, t2_k, atol=1e-4)

    def test_pet_published_reference_values(self):
        # Hoeppe (1999) Table 1. The published values, the Walther & Goestchel
        # (2018) recomputation and this implementation disagree by up to ~1.7 K
        # on the high-radiant-load rows -- the two papers already disagree with
        # each other by up to 1.3 K -- so the tolerance is per row and reflects
        # the documented spread of the model itself, not implementation error.
        # (Ta, Tmrt, v, vpa, published PET, tolerance) in degC / m s-1 / hPa
        cases = [
            (21.0, 21.0, 0.1, 12.0, 21.0, 0.05),
            (-5.0, 40.0, 0.5, 2.0, 10.0, 1.8),
            (-5.0, -5.0, 5.0, 2.0, -13.0, 0.3),
            (30.0, 60.0, 1.0, 21.0, 43.0, 1.5),
            (30.0, 30.0, 1.0, 21.0, 29.0, 0.5),
        ]
        for t_c, mrt_c, va, vpa, expected_c, tol in cases:
            pet = tmf.calculate_pet(
                np.array([tmf.celsius_to_kelvin(t_c)]),
                np.array([tmf.celsius_to_kelvin(mrt_c)]),
                np.array([va]),
                vapour_pressure_hpa=np.array([vpa]),
            )
            got_c = tmf.kelvin_to_celsius(pet[0])
            assert got_c == pytest.approx(expected_c, abs=tol)

    def test_pet_monotonicity(self):
        # Physical sanity: PET rises with air temperature, with mean radiant
        # temperature and with activity, and falls with air velocity.
        base = dict(
            t2_k=np.array([298.15]),
            mrt_k=np.array([303.15]),
            va=np.array([1.0]),
            rh=np.array([50.0]),
        )
        ref = tmf.calculate_pet(**base)
        assert tmf.calculate_pet(**{**base, "t2_k": np.array([303.15])}) > ref
        assert tmf.calculate_pet(**{**base, "mrt_k": np.array([313.15])}) > ref
        assert tmf.calculate_pet(**{**base, "m_act_w": np.array([150.0])}) > ref
        assert tmf.calculate_pet(**{**base, "va": np.array([5.0])}) < ref

    def test_pet_clothing_and_subject_variants(self):
        # Exercises every branch of the clothed-height band (clo >= 2,
        # 0.6 < clo < 2, 0.3 < clo <= 0.6, clo <= 0.3) and the female basal
        # metabolism relation. More clothing must raise PET in a cool
        # environment.
        t2_k = np.array([288.15])
        mrt_k = np.array([288.15])
        va = np.array([1.0])
        rh = np.array([50.0])
        pets = [
            float(tmf.calculate_pet(t2_k, mrt_k, va, rh=rh, clo=np.array([c]))[0])
            for c in (0.2, 0.5, 1.0, 2.5)
        ]
        assert pets == sorted(pets), "PET must increase with clothing when cool"

        female = tmf.calculate_pet(t2_k, mrt_k, va, rh=rh, sex="female")
        male = tmf.calculate_pet(t2_k, mrt_k, va, rh=rh, sex="male")
        assert np.isfinite(female).all()
        assert float(female[0]) != float(male[0])

    def test_pet_hot_and_cold_regimes(self):
        # Drives the three branches of the closed-form core-node solution:
        # set blood flow (cold), vasodilated (warm) and saturated at the
        # 90 L m-2 h-1 cap (hot). All must be finite and correctly ordered.
        t_c = np.array([-20.0, 10.0, 30.0, 45.0])
        pet = tmf.calculate_pet(
            tmf.celsius_to_kelvin(t_c),
            tmf.celsius_to_kelvin(t_c),
            np.full_like(t_c, 1.0),
            rh=np.full_like(t_c, 60.0),
        )
        assert np.isfinite(pet).all()
        assert np.all(np.diff(pet) > 0)

    def test_pet_broadcasting_and_chunking(self):
        # Scalars broadcast against arrays, 2-D shape is preserved, and an
        # array larger than the internal cache-blocking chunk gives the same
        # answer as a single small one (i.e. chunking is not observable).
        two_d = tmf.calculate_pet(
            np.full((2, 3), 300.0), np.full((2, 3), 305.0), 1.0, rh=50.0
        )
        assert two_d.shape == (2, 3)

        n = 20000  # > thermofeel.pet.CHUNK, so more than one block is processed
        big = tmf.calculate_pet(
            np.full(n, 300.0), np.full(n, 305.0), np.full(n, 1.0), rh=np.full(n, 50.0)
        )
        assert big.shape == (n,)
        np.testing.assert_allclose(big, two_d.ravel()[0], atol=1e-9)

    def test_apparent_temperature(self):
        t2_k = np.array([tmf.celsius_to_kelvin(25.0)])
        va = np.array([3])
        rh = np.array([75])
        at = np.array([tmf.calculate_apparent_temperature(t2_k, va, rh)])
        # print(f"at {at}")
        assert at == pytest.approx(299.86678322384626, abs=1e-6)

    def test_apparent_temperature_radiation(self):
        # Steadman (1994) / BoM radiation-inclusive Apparent Temperature:
        #   AT = Ta + 0.348*e - 0.70*va + 0.70*q/(va + 10) - 4.25   [Ta in degC]
        # with e the ambient vapour pressure (hPa) from the BoM approximation
        #   e = (rh/100)*6.105*exp(17.27*Ta/(237.7 + Ta)).
        # Each row's e and AT are hand-derived from those formulae (arithmetic
        # shown) and pinned as Kelvin literals; the function is never used to
        # generate its own expected value (no self-fulfilment).
        #
        # Row A: Ta=30, rh=50, va=2, q=100
        #   e       = 0.50*6.105*exp(17.27*30/267.7)   = 21.1435807 hPa
        #   0.348*e = 7.3579661 ; 0.70*va = 1.40 ; 0.70*100/12 = 5.8333333
        #   AT_c    = 30 + 7.3579661 - 1.40 + 5.8333333 - 4.25 = 37.5412994
        #   AT_K    = 310.6912994
        # Row B: Ta=35, rh=60, va=5, q=200
        #   e       = 0.60*6.105*exp(17.27*35/272.7)   = 33.6099046 hPa
        #   0.348*e = 11.6962468 ; 0.70*va = 3.50 ; 0.70*200/15 = 9.3333333
        #   AT_c    = 35 + 11.6962468 - 3.50 + 9.3333333 - 4.25 = 48.2795801
        #   AT_K    = 321.4295801
        # Row C: Ta=25, rh=40, va=0.5, q=0  (q=0 -> radiation term vanishes)
        #   e       = 0.40*6.105*exp(17.27*25/262.7)   = 12.6331850 hPa
        #   0.348*e = 4.3963484 ; 0.70*va = 0.35 ; qterm = 0
        #   AT_c    = 25 + 4.3963484 - 0.35 + 0 - 4.25 = 24.7963484
        #   AT_K    = 297.9463484
        # Row D: Ta=20, rh=80, va=3, q=300
        #   e       = 0.80*6.105*exp(17.27*20/257.7)   = 18.6581445 hPa
        #   0.348*e = 6.4930343 ; 0.70*va = 2.10 ; 0.70*300/13 = 16.1538462
        #   AT_c    = 20 + 6.4930343 - 2.10 + 16.1538462 - 4.25 = 36.2968805
        #   AT_K    = 309.4468805
        # Row E: Ta=40, rh=30, va=1, q=150
        #   e       = 0.30*6.105*exp(17.27*40/277.7)   = 22.0367568 hPa
        #   0.348*e = 7.6687914 ; 0.70*va = 0.70 ; 0.70*150/11 = 9.5454545
        #   AT_c    = 40 + 7.6687914 - 0.70 + 9.5454545 - 4.25 = 52.2642459
        #   AT_K    = 325.4142459
        cases = [
            (30.0, 50.0, 2.0, 100.0, 310.6912994),
            (35.0, 60.0, 5.0, 200.0, 321.4295801),
            (25.0, 40.0, 0.5, 0.0, 297.9463484),
            (20.0, 80.0, 3.0, 300.0, 309.4468805),
            (40.0, 30.0, 1.0, 150.0, 325.4142459),
        ]
        for t_c, rh_pc, va_ms, q_wm2, expected_k in cases:
            t2_k = np.array([tmf.celsius_to_kelvin(t_c)])
            va = np.array([va_ms])
            rh = np.array([rh_pc])
            q = np.array([q_wm2])
            at = tmf.calculate_apparent_temperature_radiation(t2_k, va, rh, q)
            assert at[0] == pytest.approx(expected_k, abs=1e-2)

        # Vapour-path check (breaks circularity): assert the shared helper
        # returns the hand-computed e for Row A independently of the AT value.
        #   e = 0.50*6.105*exp(17.27*30/267.7) = 21.1435807 hPa
        e_row_a = tmf.calculate_nonsaturation_vapour_pressure(
            np.array([tmf.celsius_to_kelvin(30.0)]), np.array([50.0])
        )
        assert e_row_a[0] == pytest.approx(21.1435807, abs=1e-6)

        # External anchor (independent MIT oracle, verified during research):
        #   Ta=23 degC (296.15 K), rh=70%, va=1, q=50 -> AT ~ 28.1 degC (301.25 K).
        # Re-derived from the BoM formula:
        #   e       = 0.70*6.105*exp(17.27*23/260.7) = 19.6104356 hPa
        #   0.348*e = 6.8244316 ; 0.70*va = 0.70 ; 0.70*50/11 = 3.1818182
        #   AT_c    = 23 + 6.8244316 - 0.70 + 3.1818182 - 4.25 = 28.0562498
        #   AT_K    = 301.2062498  (within 0.044 K of the 301.25 K anchor).
        at_anchor = tmf.calculate_apparent_temperature_radiation(
            np.array([296.15]), np.array([1.0]), np.array([70.0]), np.array([50.0])
        )
        assert at_anchor[0] == pytest.approx(301.25, abs=0.1)

        # NB: a second dossier oracle value (Ta=25 degC, rh=30%, va=0.1, q=100
        # -> 25.3 degC) does NOT reconcile with the BoM formula (hand-calc gives
        # ~30.9 degC); it traces to a different tool/convention and is
        # deliberately not pinned here.

        # Analytic identity (G2): at q=0 the radiation form differs from the
        # non-radiation calculate_apparent_temperature only by its own constants,
        #   AT_radiation - AT = (0.348 - 0.33)*e - (4.25 - 4.0) = 0.018*e - 0.25,
        # elementwise (the -0.70*va terms cancel; the q term is 0). Encoded
        # against the shipped non-radiation function.
        t2_k = np.array([tmf.celsius_to_kelvin(c) for c in (18.0, 27.0, 33.0, 41.0)])
        va = np.array([0.5, 2.0, 4.0, 6.0])
        rh = np.array([35.0, 55.0, 75.0, 95.0])
        q0 = np.zeros_like(t2_k)
        at_rad0 = tmf.calculate_apparent_temperature_radiation(t2_k, va, rh, q0)
        at_plain = tmf.calculate_apparent_temperature(t2_k, va, rh)
        e = tmf.calculate_nonsaturation_vapour_pressure(t2_k, rh)
        np.testing.assert_allclose(at_rad0 - at_plain, 0.018 * e - 0.25, atol=1e-9)

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


if __name__ == "__main__":
    unittest.main()  # pragma: no cover

"""
ThermoFeel Library
calculates different heat indexes from inputs
"""
# import statements
import numpy as np
import math
import datetime


class CalculateThermalIndex:
    """
    Heat Indices is a class that contains all the methods to calculate different heat indexes.
    it includes Relative Humidity, Universal Thermal Climate Index and Mean Radiant Temperature.
    """

    # convert kelvin to celcius
    def kelvin_to_celcius(t2m):
        t = t2m - 273.15
        return t2m

    # convert celcius to kelvin
    def celcius_to_kelvin(t2m):
        t2m = t2m + 273.15
        return t2m

    # convert from pa to hpa for e (relative humidity)
    def pa_to_hpa(rh):
        rh = rh / 10
        return rh

    # Relative Humidity
    """Relative Humidity
    :param t2m: 2m temperature
    :param percent: when set to true returns relative humidity in percentage
    """

    def calculate_rh(t2m, percent: "bool" = True) -> "returns rh":
        g = [-2.8365744e3, -6.028076559e3, 1.954263612e1, -2.737830188e-2,
             1.6261698e-5, 7.0229056e-10, -1.8680009e-13, 2.7150305]
        tk = t2m + 273.15
        ess = g[7] * np.log(tk)
        for i in range(7):
            ess = ess + g[i] * pow(tk, (i - 2))
        ess = np.exp(ess) * 0.01
        if percent == False:
            return ess
        rh = (ess / 6.105 * np.exp(17.27 * t2m / 237.7 + t2m)) * 100
        return rh

    # Heat Index
    """
    Heat Index
    :param t2m: 2m temperature
    :param rh: Relative Humidity
    """

    def calculate_heat_index(t2m, rh=None):
        if rh is None:
            rh = CalculateThermalIndex.calculate_rh(t2m, percent=False)
        t2m = CalculateThermalIndex.kelvin_to_celcius(t2m)
        rh = CalculateThermalIndex.pa_to_hpa(rh)
        hiarray = [8.784695, 1.61139411, 2.338549, 0.14611605, 1.2308094 * 10 ** -2, 2.211732 * 10 ** -3,
                   7.2546 * 10 ** -4, 3.58 * 10 ** -6]
        hi = -hiarray[0] + hiarray[1] * t2m + hiarray[2] * rh - hiarray[3] * t2m * \
             rh - hiarray[4] * rh ** 2 + hiarray[5] * t2m ** 2 * rh + hiarray[6] * \
             t2m * rh ** 2 - hiarray[7] * t2m ** 2 * rh ** 2
        return hi

    def szafunc(day, dLongitude, dLatitude):
        """
            inputs: day: datetime object
                    dLongitude: longitudes (scalar or Numpy array)
                    dLatitude: latitudes (scalar or Numpy array)
            output: solar zenith angles
        """
        dHours, dMinutes, dSeconds = day.hour, day.minute, day.second
        iYear, iMonth, iDay = day.year, day.month, day.day

        dEarthMeanRadius = 6371.01
        dAstronomicalUnit = 149597890

        ###################################################################
        # Calculate difference in days between the current Julian Day
        # and JD 2451545.0, which is noon 1 January 2000 Universal Time
        ###################################################################
        # Calculate time of the day in UT decimal hours
        dDecimalHours = dHours + (dMinutes + dSeconds / 60.) / 60.
        # Calculate current Julian Day
        liAux1 = int((iMonth - 14.) / 12.)
        liAux2 = int((1461. * (iYear + 4800. + liAux1)) / 4.) + int((367. * (iMonth - 2. - 12. * liAux1)) / 12.) - int(
            (3. * int((iYear + 4900. + liAux1) / 100.)) / 4.) + iDay - 32075.
        dJulianDate = liAux2 - 0.5 + dDecimalHours / 24.
        # Calculate difference between current Julian Day and JD 2451545.0
        dElapsedJulianDays = dJulianDate - 2451545.0

        ###################################################################
        # Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
        # ecliptic in radians but without limiting the angle to be less than 2*Pi
        # (i.e., the result may be greater than 2*Pi)
        ###################################################################
        dOmega = 2.1429 - 0.0010394594 * dElapsedJulianDays
        dMeanLongitude = 4.8950630 + 0.017202791698 * dElapsedJulianDays  # Radians
        dMeanAnomaly = 6.2400600 + 0.0172019699 * dElapsedJulianDays
        dEclipticLongitude = dMeanLongitude + 0.03341607 * np.sin(dMeanAnomaly) + 0.00034894 * np.sin(
            2. * dMeanAnomaly) - 0.0001134 - 0.0000203 * np.sin(dOmega)
        dEclipticObliquity = 0.4090928 - 6.2140e-9 * dElapsedJulianDays + 0.0000396 * np.cos(dOmega)

        ###################################################################
        # Calculate celestial coordinates ( right ascension and declination ) in radians
        # but without limiting the angle to be less than 2*Pi (i.e., the result may be
        # greater than 2*Pi)
        ###################################################################
        dSin_EclipticLongitude = np.sin(dEclipticLongitude)
        dY = np.cos(dEclipticObliquity) * dSin_EclipticLongitude
        dX = np.cos(dEclipticLongitude)
        dRightAscension = np.arctan2(dY, dX)
        if dRightAscension < 0.0:
            dRightAscension = dRightAscension + 2.0 * np.pi
        dDeclination = np.arcsin(np.sin(dEclipticObliquity) * dSin_EclipticLongitude)

        ###################################################################
        # Calculate local coordinates ( azimuth and zenith angle ) in degrees
        ###################################################################
        dGreenwichMeanSiderealTime = 6.6974243242 + 0.0657098283 * dElapsedJulianDays + dDecimalHours
        dLocalMeanSiderealTime = (dGreenwichMeanSiderealTime * 15. + dLongitude) * (np.pi / 180.)
        dHourAngle = dLocalMeanSiderealTime - dRightAscension
        dLatitudeInRadians = dLatitude * (np.pi / 180.)
        dCos_Latitude = np.cos(dLatitudeInRadians)
        dSin_Latitude = np.sin(dLatitudeInRadians)
        dCos_HourAngle = np.cos(dHourAngle)
        dZenithAngle = (
            np.arccos(dCos_Latitude * dCos_HourAngle * np.cos(dDeclination) + np.sin(dDeclination) * dSin_Latitude))
        dY = -np.sin(dHourAngle)
        dX = np.tan(dDeclination) * dCos_Latitude - dSin_Latitude * dCos_HourAngle
        dAzimuth = np.arctan2(dY, dX)
        dAzimuth[dAzimuth < 0.0] = dAzimuth[dAzimuth < 0.0] + 2.0 * np.pi
        dAzimuth = dAzimuth / (np.pi / 180.)
        # Parallax Correction
        dParallax = (dEarthMeanRadius / dAstronomicalUnit) * np.sin(dZenithAngle)
        dZenithAngle = (dZenithAngle + dParallax) / (np.pi / 180.)

        return dZenithAngle

    def calculate_solar_angles(days, months, years, lats, longs, lsmeridian, lstime):

        K = 2
        pi = 3.141592653
        if np.remainder(years.all(), 4) == 0:
            K = 1
        monthsfix = 275 * months / 9
        monthsfixp2 = months + 9 / 12
        n = np.round(monthsfix) - K * np.round(monthsfixp2) + days - 30
        d = 2 * pi * (n - 1) / 365.245

        cequtime = 0.01719 + 0.428146 * np.cos(d) - 7.352048 * np.sin(d) - \
                   3.349758 * np.cos(2 * d) - 9.362591 * np.sin(2 * d)

        longcor = 4 * (lsmeridian - longs) / 60

        locst = lstime + longcor + cequtime / 60

        solarangle = (12 - locst) * pi / 12

        solardeclin = 23.45 * np.sin(2 * pi * 284 + n / 365) * pi / 180

        solarangleeq = np.cos(solardeclin) * np.cos(solarangle) * \
                       np.cos(lats * pi / 180) + np.sin(solardeclin) * \
                       np.sin(lats * pi / 180)
        solarhangle = np.arcsin(solarangleeq)

        solarzenith = 90 - solarangleeq * 180 / pi

        return (solarzenith)

    def calculate_mean_radiant_temperature(ssrd, ssr, fdir, strd, strr, cossza):
        """
         mrt - Mean Radiant Temperature
         :param ssrd: is surface solar radiation downwards
         :param ssr: is surface net solar radiation
         :param fdir: is Total sky direct solar radiation at surface
         :param strd: is Surface thermal radiation downwards
         :param strr: is Surface net thermal radiation
         :param cossza: is cosine of solar zenith angle
         :param fp: is projected factor area
        """
        dsw = ssrd - fdir
        rsw = ssrd - ssr
        lur = strd - strr
        pi = 3.141592653
        gamma = np.arcsin(cossza) * 180 / pi
        fp = 0.308 * math.cos(pi / 180) * gamma * 0.998 - gamma * gamma / 50000
        if cossza > 0.01:
            fdir = fdir / cossza
        mrtcal = 17636684.3 * (0.5 * strd + 0.5 * lur + 0.721649485) \
                 * (0.5 * dsw + 0.5 * rsw + fp * fdir)
        mrt = CalculateThermalIndex.kelvin_to_celcius(pow(mrtcal, 0.25))
        return mrt

    def calculate_utci(t2m, va, mrt=None, rh=None):
        """
        UTCI
        :param t: is 2m temperature
        :param va: is wind speed at 10 meters
        :param mrt: is mean radiant temperature
        :param rh: is relative humidity
        Calculate UTCI with a 6th order polynomial approximation according to:
        Brode, P. et al. Deriving the operational procedure for the
        Universal Thermal Climate Index (UTCI). Int J Biometeorol (2012) 56: 48.1

        """
        t = CalculateThermalIndex.kelvin_to_celcius(t2m)
        if rh is None:
            rh = CalculateThermalIndex.calculate_rh(t2m)
        # if mrt is None:
        #    mrt = CalculateThermalIndex.calculate_mean_radiant_temperature(t, ssrd, ssr, fdir, strd, strr, cossza, fp)
        e_mrt = np.subtract(mrt, t)
        rh = rh / 10.0

        if 50 >= t2m.any() >= -50 and 17 >= va.any() >= 0 and rh.any() <= 5 and 70 >= e_mrt.any() >= -30:
            t2m2 = t2m * t2m
            t2m3 = t2m2 * t2m
            t2m4 = t2m3 * t2m
            t2m5 = t2m4 * t2m
            t2m6 = t2m5 * t2m

            va2 = va * va
            va3 = va2 * va
            va4 = va3 * va
            va5 = va4 * va
            va6 = va5 * va

            e_mrt2 = e_mrt * e_mrt
            e_mrt3 = e_mrt2 * e_mrt
            e_mrt4 = e_mrt3 * e_mrt
            e_mrt5 = e_mrt4 * e_mrt
            e_mrt6 = e_mrt5 * e_mrt

            rh2 = rh * rh
            rh3 = rh2 * rh
            rh4 = rh3 * rh
            rh5 = rh4 * rh
            rh6 = rh5 * rh

            utci_approx = t2m + 6.07562052E-01 + \
                          -2.27712343E-02 * t2m + \
                          8.06470249E-04 * t2m2 + \
                          -1.54271372E-04 * t2m3 + \
                          -3.24651735E-06 * t2m4 + \
                          7.32602852E-08 * t2m5 + \
                          1.35959073E-09 * t2m6 + \
                          -2.25836520E+00 * va + \
                          8.80326035E-02 * t2m * va + \
                          2.16844454E-03 * t2m2 * va + \
                          -1.53347087E-05 * t2m3 * va + \
                          -5.72983704E-07 * t2m4 * va + \
                          -2.55090145E-09 * t2m5 * va + \
                          -7.51269505E-01 * va2 + \
                          -4.08350271E-03 * t2m * va2 + \
                          -5.21670675E-05 * t2m2 * va2 + \
                          1.94544667E-06 * t2m3 * va2 + \
                          1.14099531E-08 * t2m4 * va2 + \
                          1.58137256E-01 * va3 + \
                          -6.57263143E-05 * t2m * va3 + \
                          2.22697524E-07 * t2m2 * va3 + \
                          -4.16117031E-08 * t2m3 * va3 + \
                          -1.27762753E-02 * va4 + \
                          9.66891875E-06 * t2m * va4 + \
                          2.52785852E-09 * t2m2 * va4 + \
                          4.56306672E-04 * va5 + \
                          -1.74202546E-07 * t2m * va5 + \
                          -5.91491269E-06 * va6 + \
                          3.98374029E-01 * e_mrt + \
                          1.83945314E-04 * t2m * e_mrt + \
                          -1.73754510E-04 * t2m2 * e_mrt + \
                          -7.60781159E-07 * t2m3 * e_mrt + \
                          3.77830287E-08 * t2m4 * e_mrt + \
                          5.43079673E-10 * t2m5 * e_mrt + \
                          -2.00518269E-02 * va * e_mrt + \
                          8.92859837E-04 * t2m * va * e_mrt + \
                          3.45433048E-06 * t2m2 * va * e_mrt + \
                          -3.77925774E-07 * t2m3 * va * e_mrt + \
                          -1.69699377E-09 * t2m4 * va * e_mrt + \
                          1.69992415E-04 * va2 * e_mrt + \
                          -4.99204314E-05 * t2m * va2 * e_mrt + \
                          2.47417178E-07 * t2m2 * va2 * e_mrt + \
                          1.07596466E-08 * t2m3 * va2 * e_mrt + \
                          8.49242932E-05 * va3 * e_mrt + \
                          1.35191328E-06 * t2m * va3 * e_mrt + \
                          -6.21531254E-09 * t2m2 * va3 * e_mrt + \
                          -4.99410301E-06 * va4 * e_mrt + \
                          -1.89489258E-08 * t * va4 * e_mrt + \
                          8.15300114E-08 * va5 * e_mrt + \
                          7.55043090E-04 * e_mrt2 + \
                          -5.65095215E-05 * t2m * e_mrt2 + \
                          -4.52166564E-07 * t2m * e_mrt2 + \
                          2.46688878E-08 * t2m3 * e_mrt2 + \
                          2.42674348E-10 * t2m4 * e_mrt2 + \
                          1.54547250E-04 * va * e_mrt2 + \
                          5.24110970E-06 * t2m * va * e_mrt2 + \
                          -8.75874982E-08 * t2m2 * va * e_mrt2 + \
                          -1.50743064E-09 * t2m3 * va * e_mrt2 + \
                          -1.56236307E-05 * va2 * e_mrt2 + \
                          -1.33895614E-07 * t2m * va2 * e_mrt2 + \
                          2.49709824E-09 * t2m2 * va2 * e_mrt2 + \
                          6.51711721E-07 * va3 * e_mrt2 + \
                          1.94960053E-09 * t2m * va3 * e_mrt2 + \
                          -1.00361113E-08 * va4 * e_mrt2 + \
                          -1.21206673E-05 * e_mrt3 + \
                          -2.18203660E-07 * t2m * e_mrt3 + \
                          7.51269482E-09 * t2m2 * e_mrt3 + \
                          9.79063848E-11 * t2m3 * e_mrt3 + \
                          1.25006734E-06 * va * e_mrt3 + \
                          -1.81584736E-09 * t2m * va * e_mrt3 + \
                          -3.52197671E-10 * t2m2 * va * e_mrt3 + \
                          -3.36514630E-08 * va2 * e_mrt3 + \
                          1.35908359E-10 * t2m * va2 * e_mrt3 + \
                          4.17032620E-10 * va3 * e_mrt3 + \
                          -1.30369025E-09 * e_mrt4 + \
                          4.13908461E-10 * t2m * e_mrt4 + \
                          9.22652254E-12 * t2m2 * e_mrt4 + \
                          -5.08220384E-09 * va * e_mrt4 + \
                          -2.24730961E-11 * t2m * va * e_mrt4 + \
                          1.17139133E-10 * va2 * e_mrt4 + \
                          6.62154879E-10 * e_mrt5 + \
                          4.03863260E-13 * t2m * e_mrt5 + \
                          1.95087203E-12 * va * e_mrt5 + \
                          -4.73602469E-12 * e_mrt6 + \
                          5.12733497E+00 * rh + \
                          -3.12788561E-01 * t2m * rh + \
                          -1.96701861E-02 * t2m2 * rh + \
                          9.99690870E-04 * t2m3 * rh + \
                          9.51738512E-06 * t2m4 * rh + \
                          -4.66426341E-07 * t2m5 * rh + \
                          5.48050612E-01 * va * rh + \
                          -3.30552823E-03 * t2m * va * rh + \
                          -1.64119440E-03 * t2m2 * va * rh + \
                          -5.16670694E-06 * t2m3 * va * rh + \
                          9.52692432E-07 * t2m4 * va * rh + \
                          -4.29223622E-02 * va2 * rh + \
                          5.00845667E-03 * t2m * va2 * rh + \
                          1.00601257E-06 * t2m2 * va2 * rh + \
                          -1.81748644E-06 * t2m3 * va2 * rh + \
                          -1.25813502E-03 * va3 * rh + \
                          -1.79330391E-04 * t2m * va3 * rh + \
                          2.34994441E-06 * t2m2 * va3 * rh + \
                          1.29735808E-04 * va4 * rh + \
                          1.29064870E-06 * t2m * va4 * rh + \
                          -2.28558686E-06 * va5 * rh + \
                          -3.69476348E-02 * e_mrt * rh + \
                          1.62325322E-03 * t2m * e_mrt * rh + \
                          -3.14279680E-05 * t2m2 * e_mrt * rh + \
                          2.59835559E-06 * t2m3 * e_mrt * rh + \
                          -4.77136523E-08 * t2m4 * e_mrt * rh + \
                          8.64203390E-03 * va * e_mrt * rh + \
                          -6.87405181E-04 * t2m * va * e_mrt * rh + \
                          -9.13863872E-06 * t2m2 * va * e_mrt * rh + \
                          5.15916806E-07 * t2m3 * va * e_mrt * rh + \
                          -3.59217476E-05 * va2 * e_mrt * rh + \
                          3.28696511E-05 * t2m * va2 * e_mrt * rh + \
                          -7.10542454E-07 * t2m2 * va2 * e_mrt * rh + \
                          -1.24382300E-05 * va3 * e_mrt * rh + \
                          -7.38584400E-09 * t2m * va3 * e_mrt * rh + \
                          2.20609296E-07 * va4 * e_mrt * rh + \
                          -7.32469180E-04 * e_mrt2 * rh + \
                          -1.87381964E-05 * t2m * e_mrt2 * rh + \
                          4.80925239E-06 * t2m2 * e_mrt2 * rh + \
                          -8.75492040E-08 * t2m3 * e_mrt2 * rh + \
                          2.77862930E-05 * va * e_mrt2 * rh + \
                          -5.06004592E-06 * t2m * va * e_mrt2 * rh + \
                          1.14325367E-07 * t2m2 * va * e_mrt2 * rh + \
                          2.53016723E-06 * va2 * e_mrt2 * rh + \
                          -1.72857035E-08 * t2m * va2 * e_mrt2 * rh + \
                          -3.95079398E-08 * va3 * e_mrt2 * rh + \
                          -3.59413173E-07 * e_mrt3 * rh + \
                          7.04388046E-07 * t2m * e_mrt3 * rh + \
                          -1.89309167E-08 * t2m2 * e_mrt3 * rh + \
                          -4.79768731E-07 * va * e_mrt3 * rh + \
                          7.96079978E-09 * t2m * va * e_mrt3 * rh + \
                          1.62897058E-09 * va2 * e_mrt3 * rh + \
                          3.94367674E-08 * e_mrt4 * rh + \
                          -1.18566247E-09 * t2m * e_mrt4 * rh + \
                          3.34678041E-10 * va * e_mrt4 * rh + \
                          -1.15606447E-10 * e_mrt5 * rh + \
                          -2.80626406E+00 * rh2 + \
                          5.48712484E-01 * t2m * rh2 + \
                          -3.99428410E-03 * t2m2 * rh2 + \
                          -9.54009191E-04 * t2m3 * rh2 + \
                          1.93090978E-05 * t2m4 * rh2 + \
                          -3.08806365E-01 * va * rh2 + \
                          1.16952364E-02 * t2m * va * rh2 + \
                          4.95271903E-04 * t2m2 * va * rh2 + \
                          -1.90710882E-05 * t2m3 * va * rh2 + \
                          2.10787756E-03 * va2 * rh2 + \
                          -6.98445738E-04 * t2m * va2 * rh2 + \
                          2.30109073E-05 * t2m2 * va2 * rh2 + \
                          4.17856590E-04 * va3 * rh2 + \
                          -1.27043871E-05 * t2m * va3 * rh2 + \
                          -3.04620472E-06 * va4 * rh2 + \
                          5.14507424E-02 * e_mrt * rh2 + \
                          -4.32510997E-03 * t2m * e_mrt * rh2 + \
                          8.99281156E-05 * t2m2 * e_mrt * rh2 + \
                          -7.14663943E-07 * t2m3 * e_mrt * rh2 + \
                          -2.66016305E-04 * va * e_mrt * rh2 + \
                          2.63789586E-04 * t2m * va * e_mrt * rh2 + \
                          -7.01199003E-06 * t2m2 * va * e_mrt * rh2 + \
                          -1.06823306E-04 * va2 * e_mrt * rh2 + \
                          3.61341136E-06 * t2m * va2 * e_mrt * rh2 + \
                          2.29748967E-07 * va3 * e_mrt * rh2 + \
                          3.04788893E-04 * e_mrt2 * rh2 + \
                          -6.42070836E-05 * t2m * e_mrt2 * rh2 + \
                          1.16257971E-06 * t2m2 * e_mrt2 * rh2 + \
                          7.68023384E-06 * va * e_mrt2 * rh2 + \
                          -5.47446896E-07 * t2m * va * e_mrt2 * rh2 + \
                          -3.59937910E-08 * va2 * e_mrt2 * rh2 + \
                          -4.36497725E-06 * e_mrt3 * rh2 + \
                          1.68737969E-07 * t2m * e_mrt3 * rh2 + \
                          2.67489271E-08 * va * e_mrt3 * rh2 + \
                          3.23926897E-09 * e_mrt4 * rh2 + \
                          -3.53874123E-02 * rh3 + \
                          -2.21201190E-01 * t2m * rh3 + \
                          1.55126038E-02 * t2m2 * rh3 + \
                          -2.63917279E-04 * t2m3 * rh3 + \
                          4.53433455E-02 * va * rh3 + \
                          -4.32943862E-03 * t2m * va * rh3 + \
                          1.45389826E-04 * t2m2 * va * rh3 + \
                          2.17508610E-04 * va2 * rh3 + \
                          -6.66724702E-05 * t2m * va2 * rh3 + \
                          3.33217140E-05 * va3 * rh3 + \
                          -2.26921615E-03 * e_mrt * rh3 + \
                          3.80261982E-04 * t2m * e_mrt * rh3 + \
                          -5.45314314E-09 * t2m2 * e_mrt * rh3 + \
                          -7.96355448E-04 * va * e_mrt * rh3 + \
                          2.53458034E-05 * t2m * va * e_mrt * rh3 + \
                          -6.31223658E-06 * va2 * e_mrt * rh3 + \
                          3.02122035E-04 * e_mrt2 * rh3 + \
                          -4.77403547E-06 * t2m * e_mrt2 * rh3 + \
                          1.73825715E-06 * va * e_mrt2 * rh3 + \
                          -4.09087898E-07 * e_mrt3 * rh3 + \
                          6.14155345E-01 * rh4 + \
                          -6.16755931E-02 * t2m * rh4 + \
                          1.33374846E-03 * t2m2 * rh4 + \
                          3.55375387E-03 * va * rh4 + \
                          -5.13027851E-04 * t2m * va * rh4 + \
                          1.02449757E-04 * va2 * rh4 + \
                          -1.48526421E-03 * e_mrt * rh4 + \
                          -4.11469183E-05 * t2m * e_mrt * rh4 + \
                          -6.80434415E-06 * va * e_mrt * rh4 + \
                          -9.77675906E-06 * e_mrt2 * rh4 + \
                          8.82773108E-02 * rh5 + \
                          -3.01859306E-03 * t2m * rh5 + \
                          1.04452989E-03 * va * rh5 + \
                          2.47090539E-04 * e_mrt * rh5 + \
                          1.48348065E-03 * rh6
            return utci_approx

    def calculate_wbgts(t2m):
        """wgbts - Wet Bulb Globe Temperature Simple
        :param t: 2m temperature
        :param rh: relative humidity
        """
        rh = CalculateThermalIndex.calculate_rh(t2m, percent=False)
        rh = CalculateThermalIndex.pa_to_hpa(rh)
        t2m = CalculateThermalIndex.kelvin_to_celcius(t2m)
        td = 243.15 * np.log(rh / 0.6112) / (17.67 - np.log(rh / 0.6112))
        wbgt = 0.066 * t2m + 4098 * rh * td / (td + 273.3) ** 2 / 0.066 + 4098 * rh / (td + 273.3) ** 2
        return (wbgt)

        # wbgts = 0.567 * t + 0.393 * rh + 3.94

    def calculate_wbgt(t2m, mrt, va):
        f = (1.1e8 * va ** 0.6) / (0.98 * 0.15 ** 0.4)
        a = f / 2
        b = -f * t2m - mrt ** 4
        rt1 = 3 ** (1 / 3)
        rt2 = (np.sqrt(3) * np.sqrt(27 * a ** 4 - 16 * b ** 3) + 9 * a ** 2)
        va = va
        wbgt_quartic = - 1 / 2 * np.sqrt((2 * 2 ** (2 / 3) * b) / (rt1 * rt2 ** (1 / 3)) +
                                         (2 ** (1 / 3) * rt2 ** (1 / 3)) / 3 ** (2 / 3)) + 1 / 2 * np.sqrt(
            (4 * a) / np.sqrt((2 * 2 ** (2 / 3) * b) /
                              (rt1 * rt2 ** (1 / 3)) + (2 ** (1 / 3) * rt2 ** (1 / 3)) / 3 ** (2 / 3)) - (
                        2 ** (1 / 3) * rt2 ** (1 / 3)) /
            3 ** (2 / 3) - (2 * 2 ** (2 / 3) * b) / (rt1 * rt2 ** (1 / 3)))
        return (wbgt_quartic)

    def calculate_mrt_from_wbgt(t2m, wbgt, va):
        f = (1.1e8 * va ** 0.6) / (0.98 * 0.15 ** 0.4)
        wbgt4 = wbgt ** 4
        dit = wbgt - t2m
        mrtcalculation = wbgt4 + f * dit
        mrtcalc2 = pow(mrtcalculation, 0.25)
        return (mrtcalc2)

    def calculate_humidex(t2m, td):
        """
        PROVISIONAL
        humidex - heat index used by the Canadian Meteorological Serivce
        :param t: 2m temperature in degrees celcius
        :param td: dew point temperature in kelvin
        """
        e = 6.11 * np.exp(5417.7530 * 1 / 273.16 - 1 / td)
        h = 0.5555 * e - 10.0
        humidex = t2m + h
        return humidex

    def calculate_net_effective_temperature(t2m, rh, va):
        ditermeq = 1 / 1.76 + 1.4 * va ** 0.75
        diterm = 0.68 - 0.0014 * rh + ditermeq
        net = 37 - (37 - t2m / diterm) - 0.29 * t2m * (1 - 0.01 * rh)
        return net


# The main method to run the code
if __name__ == '__main__':
    main()

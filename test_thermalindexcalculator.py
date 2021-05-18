import unittest
import ThermoFeel
import numpy as np
import pandas as pd
import ThermoFeel as tfc

# combining the pytest library with the functionality of np testing for use with np.array file type


class TestThermalCalculator(unittest.TestCase):
    def setUp(self):
        #variables for use in thermalindexcalculator
        t = pd.read_csv('thermofeeltestcases.csv', delimiter=',')
        self.t2m = t[['t2m']].to_numpy()
        self.td = t[['td']].to_numpy()
        self.va = t[['va']].to_numpy()
        self.mrt = t[['mrt']].to_numpy()
        self.lat = t[['lat']].to_numpy()
        self.lon = t[['lon']].to_numpy()
        self.y = t[['y']].to_numpy()
        self.m = t[['m']].to_numpy()
        self.d = t[['d']].to_numpy()
        self.h = t[['h']].to_numpy()
        self.base = 6
        self.step = 3

        #variables to raise errors

        self.varnone = np.array([None,None])

        #indices
        tr = pd.read_csv('thermofeeltestresults.csv', delimiter=',')
        self.rh = tfc.calculate_relative_humidity(self.t2m)
        self.rhpercent = tfc.calculate_relative_humidity_percent(self.t2m, self.td)
        self.heatindex = tr[['heatindex']].to_numpy()
        self.utci = tr[['utci']].to_numpy()
        self.apparenttemperature = tr[['apparenttemperature']].to_numpy()
        self.wbgts = tr[['wbgts']].to_numpy()
        self.wbgt = tr[['wbgt']].to_numpy()
        self.net = tr[['net']].to_numpy()
        self.humidex = tr[['humidex']].to_numpy()
        self.windchill = tr[['windchill']].to_numpy()

    def assert_equal(self,result,calculation):
        self.assertequal = self.assertIsNone(np.testing.assert_array_almost_equal(result, calculation))

    def assert_error(self,var):
        self.asserterror = self.assertIsNone(np.testing.assert_raises(ValueError,var))

    def test_relative_humidity(self):
        self.assert_equal(self.rh,tfc.calculate_relative_humidity(self.t2m))

    def test_relative_humidity_percent(self):
        self.assert_equal(self.rhpercent, tfc.calculate_relative_humidity_percent(self.t2m, self.td))

    def test_heat_index(self):
        self.assert_equal(self.heatindex,tfc.calculate_heat_index(self.t2m))

    def test_solar_zenith_angle(self):
        print('testing solar zenith angle')
        # self.assertIsNone(np.testing.assert_array_almost_equal(self.solarzenithangle,
        #                                                 ThermoFeel.ThermalIndexCalculator.
        #                                                 calculate_solar_zenith_angle(self.t2m)))
        pass

    def test_mean_radiant_temperature(self):
        # print('testing mean radiant temperature')
        # self.assertIsNone(np.testing.assert_array_almost_equal(self.mrt,
        #                                                ThermoFeel.ThermalIndexCalculator.
        #                                                calculate_mean_radiant_temperature(self.t2m)))
        pass
    def test_validate(self):
        self.assert_error(self.varnone)

    def test_utci(self):
        self.assert_equal(self.utci, tfc.calculate_utci(self.t2m, self.va, self.mrt))

    def test_apparent_temperature(self):
        self.assert_equal(self.apparenttemperature, tfc.calculate_apparent_temperature(self.t2m, self.rh, self.va))

    def test_wbgts(self):
        self.assert_equal(self.wbgts, tfc.calculate_wbgts(self.t2m))

    def test_wbgt(self):
        self.assert_equal(self.wbgt, tfc.calculate_wbgt(self.t2m,self.va,self.mrt))

    def test_net(self):
        self.assert_equal(self.net, tfc.calculate_net_effective_temperature(self.t2m,self.rh,self.va))

    def test_humidex(self):
        self.assert_equal(self.humidex, tfc.calculate_humidex(self.t2m, self.td))

    def test_wind_chill(self):
        self.assert_equal(self.windchill, tfc.calculate_wind_chill(self.t2m, self.va))


if __name__ == '__main__':
    unittest.main()

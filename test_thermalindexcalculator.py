import unittest
import ThermoFeel
import numpy as np

class TestThermalCalculator(unittest.TestCase):
    def setUp(self):
        self.t2m =np.array([310,300])
        self.td = np.array([280,290])
        self.va = np.array([2,0.02])
        self.rh = ThermoFeel.ThermalIndexCalculator.calculate_relative_humidity(self.t2m)
        self.rhpercent = ThermoFeel.ThermalIndexCalculator.calculate_relative_humidity_percent(self.t2m,self.td)
        self.heatindex = np.array([50.59517795,34.48123685])
        self.utci = np.array([31.68332968, 27.18415706])
        self.apparenttemperature = np.array([31.45, 22.836])
        self.wbgts = np.array([24.27395, 18.60395])
        self.wbgt = np.array([34.4907777, 35.225786])
        self.net = np.array([39.39678865, 31.01867697])
        self.humidex = np.array([31.81680821, 23.11586174])
        self.windchill = np.array([61.68250738, 51.50823866])
        self.lat = np.array([15,20])
        self.lon = np.array([30,40])
        self.y = np.array([2000,2001])
        self.m = np.array([12,1])
        self.d = np.array([12,20])
        self.h =np.array([22,13])
        self.base = 6
        self.step = 3
    def tearDown(self):
        pass
    def test_relative_humidity(self):
        print('test relative humidity')
        self.assertIsNone(np.testing.assert_array_almost_equal(self.rh, ThermoFeel.ThermalIndexCalculator.
                                                        calculate_relative_humidity(self.t2m)))
    def test_relative_humidity_percent(self):
        print('test relative humidity')
        self.assertIsNone(np.testing.assert_array_almost_equal(self.rhpercent, ThermoFeel.ThermalIndexCalculator.
                                                        calculate_relative_humidity_percent(self.t2m,self.td)))
    def test_heat_index(self):
        print('testing heat index')
        self.assertIsNone(np.testing.assert_array_almost_equal(self.heatindex,
                                                        ThermoFeel.ThermalIndexCalculator.calculate_heat_index
                                                        (self.t2m)))
    def test_solar_zenith_angle(self):
       # print('testing solar zenith angle')
        #self.assertIsNone(np.testing.assert_array_almost_equal(self.solarzenithangle,
       #                                                 ThermoFeel.ThermalIndexCalculator.
       #                                                 calculate_solar_zenith_angle(self.t2m)))
        pass
    def test_mean_radiant_temperature(self):
            #print('testing mean radiant temperature')
           # self.assertIsNone(np.testing.assert_array_almost_equal(self.mrt,
            #                                                ThermoFeel.ThermalIndexCalculator.
            #                                                calculate_mean_radiant_temperature(self.t2m)))
            pass
    def test_utci(self):
        print('testing utci')
        self.assertIsNone(np.testing.assert_almost_equal(self.utci, ThermoFeel.ThermalIndexCalculator.
                                                        calculate_utci(self.t2m, self.va, self.mrt)))

    def test_apparent_temperature(self):
        print('testing apparent temperature')
        self.assertIsNone(np.testing.assert_almost_equal(self.apparenttemperature, ThermoFeel.ThermalIndexCalculator.
                                                        calculate_apparent_temperature(self.t2m, self.rh, self.va)))

    def test_wbgts(self):
        print('testing wet bulb globe temperature simple')
        self.assertIsNone(np.testing.assert_array_almost_equal(self.wbgts, ThermoFeel.ThermalIndexCalculator.
                                                        calculate_wbgts(self.t2m)))

    def test_wbgt(self):
        print('testing wet bulb globe temperature')
        self.assertIsNone(np.testing.assert_array_almost_equal(self.wbgt, ThermoFeel.ThermalIndexCalculator.
                                                        calculate_wbgt(self.t2m)))
    def test_net(self):
        print('testing net effective temperature')
        self.assertIsNone(np.testing.assert_almost_equal(self.net, ThermoFeel.ThermalIndexCalculator.
                                          calculate_net_effective_temperature(self.t2m)))
    def test_humidex(self):
        print('testing humidex')
        self.assertIsNone(np.testing.assert_array_almost_equal(self.humidex, ThermoFeel.ThermalIndexCalculator.
                                                        calculate_humidex(self.t2m,self.td)))
    def test_wind_chill(self):
        print('testing wind chill')
        self.assertIsNone(np.testing.assert_array_almost_equal(self.windchill, ThermoFeel.ThermalIndexCalculator.
                                                        calculate_wind_chill(self.t2m,self.va)))


if __name__ == '__main__':
    unittest.main()
import unittest
import ThermoFeel
import numpy as np

class TestThermalCalculator(unittest.TestCase):
    def setUp(self):
        self.t2m =np.array[300,320,330]
        self.td = np.array[280,300,310]
        self.rh = ThermoFeel.ThermalIndexCalculator.calculate_relative_humidity(self.t2m)
        self.rhpercent = ThermoFeel.ThermalIndexCalculator.calculate_relative_humidity_percent(self.t2m,self.td)
        self.lat = np.array[15,20,40]
        self.lon = np.array[30,40,20]
        self.y = np.array[2000,2001,2010]
        self.m = np.array[12,1,3]
        self.d = np.array[12,20,23]
        self.h =np.array[22,13,15]
        self.base = 6
        self.step = 3
    def tearDown(self):
        pass
    def test_relative_humidity(self):
        result = ThermoFeel.ThermalIndexCalculator.calculate_relative_humidity(np.array[320,330])

if __name__ == '__main__':
    unittest.main()
import unittest
import ThermoFeel
import numpy as np

class TestThermalCalculator(unittest.TestCase):

    def test_relative_humidity(self):
        result = ThermoFeel.ThermalIndexCalculator.calculate_relative_humidity(np.array[320,330])
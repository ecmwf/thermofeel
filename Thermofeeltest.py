
"""
ThermoFeel Unit/Intergration Tests
"""
import ThermoFeel
def test():
    from ThermoFeel import ThermalIndexCalculator
    import numpy as np

    t2m = np.array([330,350])
    td= np.array([260,270])
    va = np.array([2,2])
    mrt = np.array([330,350])
    rh= ThermalIndexCalculator.calculate_relative_humidity(t2m)
    print(ThermalIndexCalculator.calculate_apparent_temperature(t2m,rh,va))
    #print(ThermalIndexCalculator.calculate_relative_humidity_percent(t2m,td))
    print(ThermalIndexCalculator.calculate_wbgts(t2m))
    print(ThermalIndexCalculator.calculate_utci(t2m=t2m,va=va,mrt=mrt,rh=None))
    #print(ThermalIndexCalculator.calculate_heat_index(t2m=t2m,rh=rh))
    #print(ThermalIndexCalculator.calculate_humidex(t2m=t2m,td=td))
    print(ThermalIndexCalculator.calculate_wbgt(t2m,mrt=mrt,va=va))
    #print(ThermalIndexCalculator.calculate_rh(t2m))
    #print(ThermalIndexCalculator.calculate_net_effective_temperature(t2m, rh, va))
    #tc = ThermalIndexCalculator.calculate_solar_zenith_angle(30,40,2010,12,23,10,6,6)
    #print(tc)
test()
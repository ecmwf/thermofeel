
"""
ThermoFeel Unit/Intergration Tests
"""
import ThermoFeel
def test():
    from ThermoFeel import ThermalIndexCalculator
    import numpy as np

    t2m = np.array([300,300])
    td= np.array([260,260])
    va = np.array([1,1])
    mrt = np.array([300,290])
    rh= ThermalIndexCalculator.calculate_relative_humidity(t2m)

    #print(ThermalIndexCalculator.calculate_wbgts(t2m))
    #print(ThermalIndexCalculator.calculate_utci(t2m=t2m,va=va,mrt=mrt,rh=None))
    #print(ThermalIndexCalculator.calculate_heat_index(t2m=t2m,rh=None))
    #print(ThermalIndexCalculator.calculate_humidex(t2m=t2m,td=td))
    #print(ThermalIndexCalculator.calculate_wbgt(t,mrt=mrt,va=va))
    #print(ThermalIndexCalculator.calculate_rh(t2m))
    tc = ThermalIndexCalculator.calculate_solar_zenith_angle(30,40,2010,12,23,10,6,2)
    print(tc)
test()
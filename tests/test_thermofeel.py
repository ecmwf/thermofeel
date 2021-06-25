
"""
thermofeel Unit/Intergration Tests
"""
import thermofeel
# import voidfunction
from thermofeel import *
def test():
    import numpy as np
    import pandas as pd
    t = pd.read_csv('thermofeeltestcases.csv', delimiter=',')
    t2m = t[['t2m']].to_numpy()
    #t2m = 280
    td = t[['td']].to_numpy()
    va = t[['va']].to_numpy()
    mrt = t[['mrt']].to_numpy()
    ssrd = t[['ssrd']].to_numpy()
    strd = t[['strd']].to_numpy()
    fdir = t[['fdir']].to_numpy()
    strr = t[['strr']].to_numpy()
    cossza = t[['cossza']].to_numpy()
    ssr = t[['ssr']].to_numpy()
    #print(calculate_solar_zenith_angle(np.array([40,40]),np.array([60,60]),np.array([2000,2000]),np.array([7,7]),np.array([12,12]),np.array([21,21])))
    #print(calculate_solar_zenith_angle_f(50,50,2000,3,21,14,12,3))
    #np.savetxt("mrt2.csv",calculate_mean_radiant_temperature(ssrd/3600,ssr/3600,fdir/3600,strd/3600,strr/3600,cossza/3600))
    #print(thermofeel.ThermalIndexCalculator.calculate_relative_humidity(t2m))
    #print(thermofeel.ThermalIndexCalculator.calculate_relative_humidity(t2m))
    #print(thermofeel.ThermalIndexCalculator.calculate_wind_chill(t2m,va))
    #t2m = np.array([330,350])
    #td= np.array([260,270])
    #va = np.array([2,2])
    #mrt = np.array([330,350])
    #print(t2m-273.15,"t2m")
    #np.savetxt("windchill3.csv",ThermalIndexCalculator.calculate_wind_chill(t2m,va))
    #np.savetxt("rh.csv",ThermalIndexCalculator.calculate_relative_humidity(t2m))
    np.savetxt("apparenttemperature.csv",calculate_apparent_temperature(t2m,va))
    # np.savetxt("RelativeHumid.csv",np.ravel(ThermalIndexCalculator.calculate_relative_humidity_percent(t2m,td)))
    #np.savetxt("wbgts2.csv",np.ravel(ThermalIndexCalculator.calculate_wbgts(t2m)))
    #np.savetxt("utci3.csv",calculate_utci(t2m=t2m,va=va,mrt=mrt,rh=None))
    # np.savetxt("hi.csv",ThermalIndexCalculator.calculate_heat_index(t2m=t2m,rh=rh))
    # np.savetxt("humi.csv",ThermalIndexCalculator.calculate_humidex(t2m=t2m,td=td))
    #np.savetxt("wbgt2.csv",calculate_wbgt(t2m,mrt=mrt,va=va))
    # np.savetxt("rh.csv",ThermalIndexCalculator.calculate_rh(t2m))
    #np.savetxt("NET.csv",calculate_net_effective_temperature(t2m,va))
    #tc = ThermalIndexCalculator.calculate_solar_zenith_angle(30,40,2010,12,23,10,6,6)
    #print(tc)
test()
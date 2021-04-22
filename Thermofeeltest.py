
"""
ThermoFeel Unit/Intergration Tests
"""
import ThermoFeel
def test():
    from ThermoFeel import CalculateThermalIndex
    from netCDF4 import Dataset
    import numpy as np
    temperature = 'Nov2mTemp.nc'
    windspeed = 'NovWindspeed.nc'
    mrt = 'mrtNov3.nc'
    utci = 'UTCINov.nc'

    temperature = Dataset(temperature,mode='r')
    windspeed = Dataset(windspeed,mode='r')
    mrt = Dataset(mrt,mode='r')
    utci = Dataset(utci,mode='r')
    
    lons = temperature.variables['longitude'][:]
    lats = temperature.variables['latitude'][:]
    t = temperature.variables['t2m'][:]
    mrt = mrt.variables['mrt'][:]
    va = windspeed.variables['v10'][:]
    utci = utci.variables['utci'][:]
    #print(CalculateThermalIndex.calculate_utci(t,va,va,rh=None))
    print(CalculateThermalIndex.calculate_wbgt(t,mrt=mrt,va=va))
    #print(CalculateThermalIndex.calculate_rh(t))
test()
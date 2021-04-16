
"""
ThermoFeel Unit/Intergration Tests
"""
import ThermoFeel
def test():
    from ThermoFeel import HeatIndicies
    from netCDF4 import Dataset
    import numpy as np
    temperature = 'filename'
    windspeed = 'filename'
    mrt = 'filename'
    temperature = Dataset(temperature,mode='r')
    windspeed = Dataset(windspeed,mode='r')
    mrt = Dataset(mrt,mode='r')
    
    lons = temperature.variables['longitude'][:]
    lats = temperature.variables['latitude'][:]
    T = temperature.variables['t2m'][:]
    MRT = mrt.variables['mrt'][:]
    va = windspeed.variables['v10'][:]
    print(HeatIndicies.utci(T,va,MRT,RH=None))

test()
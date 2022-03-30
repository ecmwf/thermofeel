# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


# To run this code you need to use pip install thermofeel
# import statements
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset  # , date2num, num2date

from thermofeel import (
    calculate_apparent_temperature,
    calculate_bgt,
    calculate_cos_solar_zenith_angle,
    calculate_cos_solar_zenith_angle_integrated,
    calculate_heat_index_adjusted,
    calculate_heat_index_simplified,
    calculate_humidex,
    calculate_mean_radiant_temperature,
    calculate_mrt_from_bgt,
    calculate_net_effective_temperature,
    calculate_relative_humidity_percent,
    calculate_saturation_vapour_pressure,
    calculate_utci,
    calculate_wbgt,
    calculate_wbgts,
    calculate_wbt,
    calculate_wind_chill,
    kelvin_to_celsius,
)
from thermofeel.thermofeel import approximate_dsrp

# read in two netcdf files containing all the variables
# to calculate the thermal indexes
my_example_nc_file2 = "radiation.nc"
fh2 = Dataset(my_example_nc_file2, mode="r")
lons = fh2.variables["longitude"][:]
lats = fh2.variables["latitude"][:]
ssrd = fh2.variables["ssrd"][0]
ssr = fh2.variables["ssr"][0]
fdir = fh2.variables["fdir"][0]
strd = fh2.variables["strd"][0]
strr = fh2.variables["str"][0]

lon_mg, lat_mg = np.meshgrid(lons, lats)

my_example_nc_file3 = "utcicomponents1.nc"
fh3 = Dataset(my_example_nc_file3, mode="r")
windspeed = np.sqrt(fh3.variables["u10"][0] ** 2 + fh3.variables["v10"][0] ** 2)
t2m = fh3.variables["t2m"][0]
td = fh3.variables["d2m"][0]

# calculate all indexes from thermofeel
rhp = calculate_relative_humidity_percent(t2k=t2m, tdk=td)
svp = calculate_saturation_vapour_pressure(tk=t2m)
cosszainstant = calculate_cos_solar_zenith_angle(
    h=6, lat=lat_mg, lon=lon_mg, y=2020, m=6, d=8
)
cosszaintegrated = calculate_cos_solar_zenith_angle_integrated(
    lat=lat_mg, lon=lon_mg, y=2020, m=6, d=8, h=6, tbegin=0, tend=6
)

approx_dsrp = approximate_dsrp(fdir=fdir, cossza=cosszainstant)

mrtinstant = calculate_mean_radiant_temperature(
    ssrd=ssrd,
    ssr=ssr,
    dsrp=approx_dsrp,
    fdir=fdir,
    strd=strd,
    strr=strr,
    cossza=cosszainstant,
)
mrtintegrate = calculate_mean_radiant_temperature(
    ssrd=ssrd, ssr=ssr, fdir=fdir, strd=strd, strr=strr, cossza=cosszaintegrated
)
utci = calculate_utci(t2_k=t2m, va_ms=windspeed, mrt_k=mrtintegrate, td_k=td)
wbgts = calculate_wbgts(t2m=t2m)
wbt = calculate_wbt(tc=kelvin_to_celsius(t2m), rh=rhp)
bgt = calculate_bgt(t_k=t2m, mrt=mrtintegrate, va=windspeed)
wbgt = calculate_wbgt(t_k=t2m, mrt=mrtintegrate, va=windspeed, td=td)
mrtbg = calculate_mrt_from_bgt(t2m=t2m, bgt=bgt, va=windspeed)
humidex = calculate_humidex(t2m=t2m, td=td)
net = calculate_net_effective_temperature(t2m=t2m, va=windspeed, td=td)
aptmp = calculate_apparent_temperature(t2m=t2m, va=windspeed)
windchill = calculate_wind_chill(t2m=t2m, va=windspeed)
hisimple = calculate_heat_index_simplified(t2m=t2m)
hia = calculate_heat_index_adjusted(t2m=t2m, td=td)

# to plot a single figure the example is for humidex index
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
filled_c = plt.pcolormesh(
    lon_mg, lat_mg, (humidex), transform=ccrs.PlateCarree(), cmap="RdBu_r"
)
fig.colorbar(filled_c, orientation="horizontal")
plt.show()

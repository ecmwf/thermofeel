import sys
import datetime
import numpy as np
import pandas as pd

import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade

def plot(lats, lons, vals):
    
    assert lats.shape == lons.shape
    assert lats.shape == vals.shape

    print(lats.shape)

    ax = plt.axes(projection=ccrs.PlateCarree())
    plt.contourf(lons, lats, vals, 60,
                 transform=ccrs.PlateCarree(), cmap='nipy_spectral')
    ax.coastlines()
    ax.set_global()
    plt.savefig(sys.argv[2])
    # plt.show()

def main():
    npz = np.load(sys.argv[1])
    # print(sorted(npz.files))    
    plot(npz['lats'], npz['lons'], npz['values'])

if __name__ == "__main__":
    sys.exit(main())

# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np


def plot(lats, lons, vals):

    assert lats.shape == lons.shape
    assert lats.shape == vals.shape

    print(lats.shape)

    ax = plt.axes(projection=ccrs.PlateCarree())
    plt.contourf(
        lons, lats, vals, 60, transform=ccrs.PlateCarree(), cmap="nipy_spectral"
    )
    ax.coastlines()
    ax.set_global()
    plt.savefig(sys.argv[2])
    # plt.show()


def main():
    npz = np.load(sys.argv[1])
    # print(sorted(npz.files))
    plot(npz["lats"], npz["lons"], npz["values"])


if __name__ == "__main__":
    sys.exit(main())

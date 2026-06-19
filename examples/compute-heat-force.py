# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Self-contained example: compute the Liljegren Wet Bulb Globe Temperature and
# the KNMI 0-10 heat-force scale from standard meteorological variables.

import numpy as np

import thermofeel as thermofeel


def main():
    # Inputs (SI / library conventions), as NumPy arrays so the same call works
    # on a single point or on a full gridded field.
    t2_k = np.array([298.15, 303.15, 308.15])  # 2m air temperature [K]
    rh = np.array([60.0, 70.0, 40.0])  # relative humidity [%]
    pressure = np.array([1013.0, 1010.0, 1013.0])  # surface pressure [hPa]
    va = np.array([3.0, 2.0, 1.0])  # 10m wind speed [m/s]
    ssrd = np.array([400.0, 600.0, 800.0])  # instantaneous SSRD [W/m2]
    fdir = np.array([0.5, 0.6, 0.7])  # direct-beam fraction [0-1]
    cossza = np.array([0.6, 0.8, 0.9])  # cosine of solar zenith angle

    # The cosine of the solar zenith angle can be obtained from earthkit-meteo:
    #   from earthkit.meteo import solar
    #   cossza = solar.cos_solar_zenith_angle(...)

    wbgt_k = thermofeel.calculate_wbgt_liljegren(
        t2_k, rh, pressure, va, ssrd, fdir, cossza
    )
    wbgt_c = thermofeel.kelvin_to_celsius(wbgt_k)
    heat_force = thermofeel.calculate_heat_force(wbgt_k)

    print(f"thermofeel version: {thermofeel.__version__}")
    print(f"{'WBGT [K]':>12} {'WBGT [C]':>12} {'heat force':>12}")
    for wk, wc, hf in zip(wbgt_k, wbgt_c, heat_force):
        print(f"{wk:12.3f} {wc:12.3f} {int(hf):12d}")


if __name__ == "__main__":
    main()

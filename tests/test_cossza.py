# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from math import cos, radians

import thermofeel


def test_calculate_cos_solar_zenith_angle():

    # should return ~ 0.360303587797559
    cossza = thermofeel.calculate_cos_solar_zenith_angle(
        lat=48.81667, lon=2.28972, d=15, m=11, y=2006, h=10.58333
    )

    assert abs(cossza - 0.360303587797559) < 1e-6

    # London, ~ 0.8799471697555967
    cossza = thermofeel.calculate_cos_solar_zenith_angle(
        lat=51.0, lon=0.0, d=4, m=6, y=2021, h=12.0
    )

    assert abs(cossza - 0.8799471697555967) < 1e-6

    # from alternative formula
    assert abs(cossza - cos(radians(90.0 - 61.5))) < 1e-2


def test_calculate_cos_solar_zenith_angle_integrated():

    cossza = thermofeel.calculate_cos_solar_zenith_angle_integrated(
        lat=48.81667, lon=2.28972, d=15, m=11, y=2006, h=10.58333, base=0, step=3
    )

    assert abs(cossza - 0.34495713937581207) < 1e-6


if __name__ == "__main__":
    test_calculate_cos_solar_zenith_angle()
    test_calculate_cos_solar_zenith_angle_integrated()

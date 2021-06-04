import sys
import eccodes
import datetime

from math import cos

from thermofeel.thermofeel import *

def main():
    
    # should return ~ 0.360303587797559
    cossza = calculate_cos_solar_zenith_angle(lat=48.81667, lon=2.28972, d=15, m=11, y=2006, h=10.58333)
    print(cossza)

    # London, ~ 0.8788171126619654
    cossza = calculate_cos_solar_zenith_angle(lat=51.0, lon=0.0, d=4, m=6, y=2021, h=12.0)
    
    print(cossza, cos((90-61.5) * to_radians))

if __name__ == "__main__":
    sys.exit(main())

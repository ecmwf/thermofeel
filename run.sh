#!/bin/bash

# reference
python tonumpy.py cossza11.grib cossza11.npz
# calculations with thermofeel
python cossza.py cossza11.grib cossza.npz cossza-int.npz 

# make plots
python plot.py cossza11.npz cossza11.png
python plot.py cossza.npz cossza.png
python plot.py cossza-int.npz cossza-int.png


# Examples

Runnable examples live in the
[`examples/`](https://github.com/ecmwf/thermofeel/tree/main/examples) directory
of the repository:

- `compute-heat-force.py` — Liljegren WBGT and the KNMI 0–10 heat-force scale
  (self-contained).
- `compute-obs.py` — computing indices from station observations.
- `thermofeeljupyterexamples.ipynb` — a Jupyter notebook walkthrough.

## Heat force from standard variables

The functions are vectorised over NumPy, so the same call works on a single
point or on a full gridded field.

```python
import numpy as np

import thermofeel as thermofeel

t2_k = np.array([298.15, 303.15, 308.15])  # 2 m air temperature [K]
rh = np.array([60.0, 70.0, 40.0])  # relative humidity [%]
pressure = np.array([1013.0, 1010.0, 1013.0])  # surface pressure [hPa]
va = np.array([3.0, 2.0, 1.0])  # 10 m wind speed [m/s]
ssrd = np.array([400.0, 600.0, 800.0])  # instantaneous SSRD [W/m2]
fdir = np.array([0.5, 0.6, 0.7])  # direct-beam fraction [0-1]
cossza = np.array([0.6, 0.8, 0.9])  # cosine of solar zenith angle

wbgt_k = thermofeel.calculate_wbgt_liljegren(
    t2_k, rh, pressure, va, ssrd, fdir, cossza
)
heat_force = thermofeel.calculate_heat_force(wbgt_k)
```

The cosine of the solar zenith angle can be obtained from
[earthkit-meteo](https://github.com/ecmwf/earthkit-meteo):

```python
from earthkit.meteo import solar

cossza = solar.cos_solar_zenith_angle(...)
```

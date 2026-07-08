# Approximations

Some indices need an input that a dataset does not provide. The
`thermofeel.approximations` submodule offers clearly-labelled **estimators** for
those inputs. They are deliberately kept out of the top level — you call them as
`thermofeel.approximations.<name>` so it is always visible that an approximation
is in play — and they are **estimation grade, not observations**: use a dataset
that carries the real field for quantitative work.

## Estimating direct solar radiation (`fdir`)

The [mean radiant temperature](mrt.md) — and hence UTCI, WBGT and PMV — needs
the direct (beam) solar radiation on a horizontal surface, ECMWF `fdir`. Some
datasets (notably [ECMWF open data](https://www.ecmwf.int/en/forecasts/datasets/open-data))
provide the global horizontal radiation `ssrd` but **not** `fdir`. The
direct/diffuse split cannot be recovered exactly from `ssrd` (two skies with the
same `ssrd` can have very different beam fractions), so these estimators use an
empirical decomposition driven by the clearness index. Expect errors of order
0.1–0.25 of the global radiation, larger in broken cloud and at low sun; the
resulting MRT/UTCI are a demonstration, not a validation-grade product.

Two models are provided:

| Function | Model | Notes |
| --- | --- | --- |
| `approximate_fdir_erbs(ssrd, cossza, *, doy=None, ...)` | Erbs et al. (1982) | diffuse-fraction split from the clearness index; `doy` applies the Earth–Sun distance correction |
| `approximate_fdir_disc(ssrd, cossza, doy, *, pressure_hpa=1013.25, ...)` | DISC — Maxwell (1987) | quasi-physical direct-normal model using the pressure-corrected air mass; generally more accurate |

You need the global horizontal radiation as a flux in W/m² (de-accumulate ECMWF
`ssrd` and divide by the period), the cosine of the solar zenith angle (e.g. from
[earthkit-meteo](https://github.com/ecmwf/earthkit-meteo)), and the day of the
year. Both return the estimated direct horizontal radiation in the same unit as
`ssrd`.

```python
from thermofeel import approximations

fdir = approximations.approximate_fdir_erbs(ssrd, cossza, doy=day_of_year)
# or, generally more accurate:
fdir = approximations.approximate_fdir_disc(ssrd, cossza, day_of_year, pressure_hpa=surface_pressure)
```

The runnable example `examples/compute-thermal-indices.py` wires this in behind
its `--approximate-fdir[=erbs|disc]` flag, so the open-data source can produce
(approximate) MRT / UTCI / WBGT.

More information:

Erbs, D.G., Klein, S.A. & Duffie, J.A. (1982). Estimation of the diffuse
radiation fraction for hourly, daily and monthly-average global radiation. Solar
Energy 28(4), 293–302. <https://doi.org/10.1016/0038-092X(82)90302-4>

Maxwell, E.L. (1987). A Quasi-Physical Model for Converting Hourly Global
Horizontal to Direct Normal Insolation. SERI/TR-215-3087, Solar Energy Research
Institute, Golden, CO.

Spencer, J.W. (1971). Fourier series representation of the position of the Sun.
Search 2(5), 172. (Earth–Sun distance factor.)

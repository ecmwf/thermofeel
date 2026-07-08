# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""
Approximations used to *estimate* an input that a thermal index needs but that a
dataset does not provide.

These are clearly-labelled empirical **estimates**, deliberately kept separate
from the validated index formulas in ``thermofeel.thermofeel`` and NOT
re-exported at the top level: calling code carries the ``.approximations.``
marker so it is always visible that an approximation is in play. Their outputs
are demonstration / estimation grade, not observations; use a dataset that
provides the real field for quantitative work.

Currently provides estimators of the direct (beam) solar radiation on a
horizontal surface -- ECMWF ``fdir`` -- from the global horizontal radiation
``ssrd`` and the solar geometry, for datasets (e.g. ECMWF open data) that carry
``ssrd`` but not ``fdir``. ``fdir`` is required by
``calculate_mean_radiant_temperature`` and hence by UTCI / WBGT / PMV. The
direct/diffuse split cannot be reconstructed exactly from ``ssrd`` (two skies
with the same ``ssrd`` can have very different beam fractions), so these use an
empirical decomposition driven by the clearness index; expect errors of order
0.1-0.25 of the global radiation, larger in broken cloud and at low sun.

References
----------
- Erbs, D.G., Klein, S.A. & Duffie, J.A. (1982). Estimation of the diffuse
  radiation fraction for hourly, daily and monthly-average global radiation.
  Solar Energy 28(4), 293-302. https://doi.org/10.1016/0038-092X(82)90302-4
- Maxwell, E.L. (1987). A Quasi-Physical Model for Converting Hourly Global
  Horizontal to Direct Normal Insolation. SERI/TR-215-3087, Solar Energy
  Research Institute, Golden, CO.
- Spencer, J.W. (1971). Fourier series representation of the position of the
  Sun. Search 2(5), 172. (Earth-Sun distance factor.)
- Kasten, F. (1966). A new table and approximation formula for the relative
  optical air mass. Arch. Meteorol. Geophys. Bioklimatol. B 14, 206-223.
"""

import numpy as np
from numpy.typing import ArrayLike

# Mean total solar irradiance at 1 AU [W m-2] (Kopp & Lean, 2011).
SOLAR_CONSTANT = 1361.0


def earth_sun_distance_factor(doy: ArrayLike) -> np.ndarray:
    """
    Earth-Sun distance (eccentricity) factor (mean-distance normal irradiance
    scaling), via the Spencer (1971) Fourier series.

        :param doy: (float array) day of the year (1..365/366)
        returns the dimensionless factor (a0/a)^2, ranging ~0.967 (aphelion,
        early July) to ~1.035 (perihelion, early January)

    Multiply the solar constant by this factor to obtain the actual
    extraterrestrial normal irradiance for the day.

    Reference: Spencer (1971); see also Iqbal, An Introduction to Solar
    Radiation (1983).
    """
    doy = np.asarray(doy, dtype=float)
    b = 2.0 * np.pi * (doy - 1.0) / 365.0
    return (
        1.00011
        + 0.034221 * np.cos(b)
        + 0.00128 * np.sin(b)
        + 0.000719 * np.cos(2.0 * b)
        + 7.7e-05 * np.sin(2.0 * b)
    )


def _relative_airmass_kasten1966(cossza: np.ndarray) -> np.ndarray:
    """Kasten (1966) relative optical air mass from cos(solar zenith angle).

    ``cossza`` is clipped to a small positive minimum so the formula stays
    finite below the horizon; those elements are irrelevant to the callers (the
    air mass is capped and the result masked to night)."""
    cz = np.clip(cossza, 1.0e-3, 1.0)
    zdeg = np.degrees(np.arccos(cz))
    return 1.0 / (cz + 0.15 * (93.885 - zdeg) ** (-1.253))


def approximate_fdir_erbs(
    ssrd: ArrayLike,
    cossza: ArrayLike,
    *,
    doy: ArrayLike | None = None,
    solar_constant: float = SOLAR_CONSTANT,
    min_cossza: float = 0.065,
) -> np.ndarray:
    """
    Estimate direct horizontal solar radiation (ECMWF ``fdir``) from global
    horizontal ``ssrd`` and cos(solar zenith angle) via the Erbs et al. (1982)
    diffuse-fraction decomposition.

        :param ssrd: (float array) surface global horizontal solar radiation
            (ECMWF ssrd), a mean flux [W m-2] (or any energy unit; the output
            follows the input unit)
        :param cossza: (float array) cosine of the solar zenith angle
        :param doy: (optional, float array) day of the year; if given, the
            extraterrestrial irradiance is corrected for the Earth-Sun distance
            (``earth_sun_distance_factor``); if omitted the mean-distance solar
            constant is used
        :param solar_constant: (float) solar constant [W m-2] (default 1361)
        :param min_cossza: (float) below this cos(zenith) the Sun is treated as
            below the horizon; ``fdir`` returns 0 there (default 0.065)
        returns the estimated direct horizontal radiation, in the same unit as
        ``ssrd``, clipped to [0, ssrd] (and to 0 where ``ssrd`` is negative)

    The Erbs correlation splits global horizontal radiation into its direct and
    diffuse parts from the clearness index ``kt = ssrd / TOA``. It is an
    empirical hourly fit; applying it to longer step-mean fluxes is a further
    approximation. NaN inputs propagate to NaN. Estimation grade only -- use a
    real ``fdir`` field for quantitative work.

    Reference: Erbs, Klein & Duffie (1982), Solar Energy 28(4):293-302,
    https://doi.org/10.1016/0038-092X(82)90302-4
    """
    ssrd = np.asarray(ssrd, dtype=float)
    cossza = np.asarray(cossza, dtype=float)
    e0 = earth_sun_distance_factor(doy) if doy is not None else 1.0
    ssrd, cossza, e0 = np.broadcast_arrays(ssrd, cossza, np.asarray(e0, dtype=float))
    finite = np.isfinite(ssrd) & np.isfinite(cossza) & np.isfinite(e0)

    toa = solar_constant * e0 * np.maximum(cossza, min_cossza)
    day = (cossza > min_cossza) & (ssrd > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        kt = np.clip(np.where(day, ssrd / toa, 0.0), 0.0, 1.0)
    kd = np.where(
        kt <= 0.22,
        1.0 - 0.09 * kt,
        np.where(
            kt <= 0.80,
            0.9511 - 0.1604 * kt + 4.388 * kt**2 - 16.638 * kt**3 + 12.336 * kt**4,
            0.165,
        ),
    )
    fdir = np.clip(np.where(day, ssrd * (1.0 - kd), 0.0), 0.0, np.maximum(ssrd, 0.0))
    return np.where(finite, fdir, np.nan)


def approximate_fdir_disc(
    ssrd: ArrayLike,
    cossza: ArrayLike,
    doy: ArrayLike,
    *,
    pressure_hpa: ArrayLike = 1013.25,
    solar_constant: float = SOLAR_CONSTANT,
    min_cossza: float = 0.065,
    max_airmass: float = 12.0,
) -> np.ndarray:
    """
    Estimate direct horizontal solar radiation (ECMWF ``fdir``) from global
    horizontal ``ssrd`` and the solar geometry via the DISC model (Maxwell,
    1987) -- a quasi-physical direct-normal estimator that generally improves on
    the simpler Erbs split.

        :param ssrd: (float array) surface global horizontal solar radiation
            (ECMWF ssrd), a mean flux [W m-2]
        :param cossza: (float array) cosine of the solar zenith angle
        :param doy: (float array) day of the year (required; sets the
            Earth-Sun distance correction of the extraterrestrial irradiance)
        :param pressure_hpa: (float array) surface pressure [hPa] for the
            pressure-corrected air mass (default 1013.25)
        :param solar_constant: (float) solar constant [W m-2] (default 1361)
        :param min_cossza: (float) below this cos(zenith) the Sun is treated as
            below the horizon; ``fdir`` returns 0 there (default 0.065)
        :param max_airmass: (float) air mass is capped at this value in the Kn
            fit, per the original model (default 12)
        returns the estimated direct horizontal radiation [W m-2], clipped to
        [0, ssrd] (and to 0 where ``ssrd`` is negative)

    DISC converts the global clearness index to a direct-beam clearness index
    ``Kn`` using empirical polynomials in ``kt`` and the (pressure-corrected
    Kasten 1966) air mass, then ``dni = Kn * I0`` and ``fdir = dni * cossza``.
    NaN inputs propagate to NaN. Estimation grade only.

    Reference: Maxwell, E.L. (1987), A Quasi-Physical Model for Converting Hourly
    Global Horizontal to Direct Normal Insolation, SERI/TR-215-3087. The
    coefficients are those of the DISC model as implemented in pvlib
    (``pvlib.irradiance.disc``).
    """
    ssrd = np.asarray(ssrd, dtype=float)
    cossza = np.asarray(cossza, dtype=float)
    pressure_hpa = np.asarray(pressure_hpa, dtype=float)
    doy = np.asarray(doy, dtype=float)
    ssrd, cossza, pressure_hpa, doy = np.broadcast_arrays(
        ssrd, cossza, pressure_hpa, doy
    )
    finite = (
        np.isfinite(ssrd)
        & np.isfinite(cossza)
        & np.isfinite(pressure_hpa)
        & np.isfinite(doy)
    )

    i0 = solar_constant * earth_sun_distance_factor(doy)  # extraterrestrial normal
    day = (cossza > min_cossza) & (ssrd > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        kt = np.clip(
            np.where(day, ssrd / (i0 * np.maximum(cossza, min_cossza)), 0.0), 0.0, 1.0
        )
        am = _relative_airmass_kasten1966(cossza) * pressure_hpa / 1013.25
    am = np.minimum(am, max_airmass)

    cloudy = kt <= 0.6
    a = np.where(
        cloudy,
        0.512 + kt * (-1.56 + kt * (2.286 - 2.222 * kt)),
        -5.743 + kt * (21.77 + kt * (-27.49 + 11.56 * kt)),
    )
    b = np.where(
        cloudy,
        0.37 + 0.962 * kt,
        41.4 + kt * (-118.5 + kt * (66.05 + 31.9 * kt)),
    )
    c = np.where(
        cloudy,
        -0.28 + kt * (0.932 - 2.048 * kt),
        -47.01 + kt * (184.2 + kt * (-222.0 + 73.81 * kt)),
    )
    with np.errstate(over="ignore"):
        delta_kn = a + b * np.exp(c * am)
    knc = 0.866 + am * (-0.122 + am * (0.0121 + am * (-0.000653 + 1.4e-05 * am)))
    kn = knc - delta_kn
    dni = np.where(day & (kn > 0), kn * i0, 0.0)
    fdir = np.clip(dni * cossza, 0.0, np.maximum(ssrd, 0.0))
    return np.where(finite, fdir, np.nan)

# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

import numpy as np


def daily_mean_temperature(t2_min, t2_max):
    """
    Daily mean temperature computed from min and max

        :param t2_min: 2-metre daily minimum temperature.
        :param t2_max: 2-metre daily maximum temperature.

    Reference: Nairn and Fawcett (2014)
    https://doi.org/10.3390/ijerph120100227
    """
    return 0.5 * (t2_min + t2_max)


def significance_index(dmt, threshold):
    """
    Significance index

        :param dmt: 3-day running mean daily mean temperature (day i, i+1, i+2).
        :param threshold: climatological threshold, e.g., 95th (heat extremes)
            or 5th (cold extremes) percentile of daily mean temperature over a
            reference period.

    Reference: Nairn and Fawcett (2014), equation (1)
    https://doi.org/10.3390/ijerph120100227
    """
    return dmt - threshold


def acclimatisation_index(dmt, threshold):
    """
    Acclimatisation index

        :param dmt: 3-day running mean daily mean temperature (days i, i+1, i+2).
        :param: threshold: 30-day running mean daily mean temperature (days i-30, ..., i-1).

    Reference: Nairn and Fawcett (2014), equation (2)
    https://doi.org/10.3390/ijerph120100227
    """
    return dmt - threshold


def excess_heat_factor(ehi_sig, ehi_accl, clip=False):
    """
    Excess heat factor

        :param ehi_sig: Significance index, e.g., computed with respect to 95th
            percentile of daily mean temperature over a reference period.
        :param ehi_accl: Acclimatisation index.
        :param clip: Whether to clip the lower value range at zero. Disabled by default.

    Reference: Nairn and Fawcett (2014), equation (3)
    https://doi.org/10.3390/ijerph120100227
    """
    if clip:
        ehi_sig = np.maximum(0, ehi_sig)
    return ehi_sig * np.maximum(1.0, ehi_accl)


def heatwave_severity(exhf, threshold):
    """
    Heatwave severity index

        :param exhf: Excess heat factor.
        :param threshold: Excess heat factor threshold, the 85th percentile
            of all positive excess heat factor values in a reference period.

    Reference: Nairn (2018), equation (4)
    https://doi.org/10.3390/ijerph15112494
    """
    return exhf / threshold


def excess_cold_factor(ehi_sig, ehi_accl, clip=False):
    """
    Excess cold factor

        :param ehi_sig: Significance index, e.g., computed with respect to 5th
            percentile of daily mean temperature over a reference period..
        :param ehi_accl: Acclimatisation index.
        :param clip: Whether to clip the upper value range at zero. Disabled by default.

    Reference: Nairn (2013), equation (7)
    https://www.cawcr.gov.au/technical-reports/CTR_060.pdf
    """
    if clip:
        ehi_sig = np.minimum(0, ehi_sig)
    return -ehi_sig * np.minimum(-1.0, ehi_accl)

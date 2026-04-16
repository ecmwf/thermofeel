Excess Heat and Cold Factors
============================

The **Excess Heat Factor** (EXHF) and **Excess Cold Factor** (EXCF) are thermal
indices that quantify the hazard of extreme heat and cold events.

Each index combines two components: a *significance index* that measures how
much the current temperature exceeds a climatological threshold, and an
*acclimatisation index* that measures how much the current temperature departs
from the recent past. EXHF or EXCF indicate events that are both unusually hot
or cold relative to climatology *and* to which the population has had little
time to acclimatise, and are therefore associated with elevated health risk.

The indices and their reference definitions typically use a day defined from
9 am to 9 am local time, so that the daily maximum typically precedes the daily
minimum. This accounts for the greater physiological significance of a hot night
following a hot day compared to the other way around.

References
----------

* Nairn, J., & Fawcett, R. (2014). The Excess Heat Factor: A Metric for Heatwave Intensity and Its Use in Classifying Heatwave Severity. *Int J Environ Res Public Health*, 12(1), 227–253. https://doi.org/10.3390/ijerph120100227
* Nairn, J. (2013). Defining heatwaves: heatwave defined as a heat-impact event servicing all community and business sectors in Australia. CAWCR Technical Report No. 060. https://www.cawcr.gov.au/technical-reports/CTR_060.pdf
* Nairn, J. (2018). Heatwave Severity. *Int J Environ Res Public Health*, 15(11), 2494. https://doi.org/10.3390/ijerph15112494


How To Use
----------

The indices are computed in a pipeline of steps, as described in the subsections
below. Each step builds on the output of the previous one.

.. note::

    The functions defined here only implement the basic formulas of the indices,
    but not the required aggregations of the input data. The intermediate
    aggregated quantities required by the pipeline, including running means
    of daily mean temperature and climatological percentile thresholds, can, for
    example, be computed with `earthkit.transforms <https://earthkit-transforms.readthedocs.io>`_.

**Daily Mean Temperature**

The first step is to compute the daily mean temperature from the 2-metre daily
minimum and maximum temperatures. Both inputs must be in the same unit (K or °C);
the output is in the same unit.

.. code-block:: python

    from thermofeel.excess_heat import daily_mean_temperature

    dmt = daily_mean_temperature(t2_min, t2_max)

**Significance Index**

The significance index measures how far the 3-day running mean of daily mean
temperature exceeds a climatological threshold. For heat events the threshold is
typically the 95th percentile of daily mean temperature over a reference period;
for cold events the 5th percentile is used instead. The output unit matches the
input unit (K or °C).

.. code-block:: python

    from thermofeel.excess_heat import significance_index

    # heat events: threshold = 95th-percentile climatology of dmt
    ehi_sig_heat = significance_index(dmt_3day_mean, threshold_95th)

    # cold events: threshold = 5th-percentile climatology of dmt
    ehi_sig_cold = significance_index(dmt_3day_mean, threshold_5th)

**Acclimatisation Index**

The acclimatisation index measures how far the 3-day running mean of daily mean
temperature exceeds the 30-day running mean of daily mean temperature in the
preceding period. The output unit matches the input unit (K or °C).

.. code-block:: python

    from thermofeel.excess_heat import acclimatisation_index

    ehi_accl = acclimatisation_index(dmt_3day_mean, dmt_30day_mean)

**Excess Heat Factor**

The excess heat factor is the product of the significance index and the
acclimatisation index (clamped from below at 1). It is returned in K² (or °C²
if the inputs are in °C). Optionally, set ``clip=True`` to restrict the result
to non-negative values, retaining only heat events.

.. code-block:: python

    from thermofeel.excess_heat import excess_heat_factor

    exhf = excess_heat_factor(ehi_sig_heat, ehi_accl)
    exhf_clipped = excess_heat_factor(ehi_sig_heat, ehi_accl, clip=True)

**Heatwave Severity**

The heatwave severity index normalises the excess heat factor by a threshold —
typically the 85th percentile of all *positive* EXHF values over a reference
period — to produce a dimensionless measure of event severity relative to a
historical baseline.

.. code-block:: python

    from thermofeel.excess_heat import heatwave_severity

    hsev = heatwave_severity(exhf, threshold_exhf_85th)

**Excess Cold Factor**

The excess cold factor uses the same formulation as the excess heat factor but
is symmetric about zero to capture cold extremes. Use the significance index
computed against the 5th-percentile threshold. Optionally set ``clip=True`` to
restrict the result to non-negative values, retaining only cold events.

.. code-block:: python

    from thermofeel.excess_heat import excess_cold_factor

    excf = excess_cold_factor(ehi_sig_cold, ehi_accl)
    excf_clipped = excess_cold_factor(ehi_sig_cold, ehi_accl, clip=True)

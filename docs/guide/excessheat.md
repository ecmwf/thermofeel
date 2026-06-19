# Excess Heat and Cold Factors

The **Excess Heat Factor** (EXHF) and **Excess Cold Factor** (EXCF) are thermal
indices that quantify the hazard of extreme heat and cold events.

Each index combines two components: a *significance index* that measures how much
the current temperature exceeds a climatological threshold, and an
*acclimatisation index* that measures how much the current temperature departs
from the recent past. EXHF or EXCF indicate events that are both unusually hot or
cold relative to climatology *and* to which the population has had little time to
acclimatise, and are therefore associated with elevated health risk.

The indices and their reference definitions typically use a day defined from 9 am
to 9 am local time, so that the daily maximum typically precedes the daily
minimum. This accounts for the greater physiological significance of a hot night
following a hot day compared to the other way around.

More information:

- Nairn, J. & Fawcett, R. The Excess Heat Factor: A Metric for Heatwave Intensity and Its Use in Classifying Heatwave Severity. Int J Environ Res Public Health 12(1), 227–253 (2014). <https://doi.org/10.3390/ijerph120100227>
- Nairn, J. Defining heatwaves: heatwave defined as a heat-impact event servicing all community and business sectors in Australia. CAWCR Technical Report No. 060 (2013). <https://www.cawcr.gov.au/technical-reports/CTR_060.pdf>
- Nairn, J. Heatwave Severity. Int J Environ Res Public Health 15(11), 2494 (2018). <https://doi.org/10.3390/ijerph15112494>

## How to use

!!! note

    The functions defined here only implement the basic formulas of the indices,
    not the required aggregations of the input data. The intermediate aggregated
    quantities required by the pipeline — including running means of daily mean
    temperature and climatological percentile thresholds — can for example be
    computed with
    [earthkit-transforms](https://earthkit-transforms.readthedocs.io).

### Daily mean temperature

The excess heat and cold factor computations are based on 2 m daily mean
temperature. Nairn and Fawcett (2014) approximate the daily mean temperature with
the average of the daily minimum and maximum temperatures, which are more readily
available as inputs from long-running historical records.

```python
from thermofeel.excess_heat import daily_mean_temperature

dmt = daily_mean_temperature(t2_min, t2_max)
```

Both inputs must be in the same unit (K or °C); the output is in the same unit.

### Significance index

The significance index measures how far the 3-day running mean of daily mean
temperature exceeds a climatological threshold. For heat events the threshold is
typically the 95th percentile of daily mean temperature over a reference period;
for cold events the 5th percentile is used instead.

```python
from thermofeel.excess_heat import significance_index

# heat events: threshold = 95th-percentile climatology of dmt
ehi_sig_heat = significance_index(dmt_3day_mean, threshold_95th)

# cold events: threshold = 5th-percentile climatology of dmt
ehi_sig_cold = significance_index(dmt_3day_mean, threshold_5th)
```

The output unit matches the input unit.

!!! note

    The use of a static climatological threshold in the significance index means
    that indicated hot and cold extremes will almost exclusively occur in the
    summer or winter seasons, respectively. A time-varying climatological
    threshold allows for the detection of warm and cold *spells* (relative
    anomalies) instead.

### Acclimatisation index

The acclimatisation index measures how far the 3-day running mean of daily mean
temperature (days i, i+1, i+2) exceeds the 30-day running mean of daily mean
temperature in the preceding period (days i-30, ..., i-1).

```python
from thermofeel.excess_heat import acclimatisation_index

ehi_accl = acclimatisation_index(dmt_3day_mean, dmt_30day_mean)
```

The output unit matches the input unit.

### Excess heat factor

The excess heat factor combines the significance index and the acclimatisation
index. Optionally, set `clip=True` to restrict the EXHF to non-negative values.

```python
from thermofeel.excess_heat import excess_heat_factor

exhf = excess_heat_factor(ehi_sig_heat, ehi_accl)
exhf_clipped = excess_heat_factor(ehi_sig_heat, ehi_accl, clip=True)
```

The output unit is the input unit squared, e.g. K².

### Heatwave severity

The heatwave severity index normalises the excess heat factor by a threshold to
produce a dimensionless measure of event severity relative to a historical
baseline. The 85th percentile of all *positive* EXHF values over a reference
period is typically used as the threshold value.

```python
from thermofeel.excess_heat import heatwave_severity

hsev = heatwave_severity(exhf, threshold_exhf_85th)
```

The severity index is nondimensional.

### Excess cold factor

The excess cold factor uses the same formulation as the excess heat factor but is
symmetric about zero to capture cold extremes. Use the significance index computed
against the 5th-percentile threshold. Optionally set `clip=True` to restrict the
result to cold events only.

```python
from thermofeel.excess_heat import excess_cold_factor

excf = excess_cold_factor(ehi_sig_cold, ehi_accl)
excf_clipped = excess_cold_factor(ehi_sig_cold, ehi_accl, clip=True)
```

The output unit is the input unit squared, e.g. K².

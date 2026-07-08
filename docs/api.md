# API reference

The thermal-index functions are available at the top level of the package, for
example `thermofeel.calculate_utci(...)` or
`thermofeel.calculate_excess_heat_factor(...)`. The supporting functions for the
excess heat/cold factors (daily mean temperature, significance and
acclimatisation indices, heatwave severity) live in the `thermofeel.excess_heat`
submodule. This page is generated from the source docstrings.

All functions are vectorised over NumPy and take SI inputs; **pass arrays and
wrap scalars in an array** (e.g. `np.array([300.0])`). See
[Calling convention](guide/overview.md#calling-convention) for details.

## Thermal indices and supporting quantities

::: thermofeel.thermofeel
    options:
      show_root_heading: false
      show_root_toc_entry: false
      heading_level: 3

## Excess heat and cold factors

The Excess Heat Factor and Excess Cold Factor are available at the top level as
`thermofeel.calculate_excess_heat_factor(...)` and
`thermofeel.calculate_excess_cold_factor(...)` (listed above). The
`thermofeel.excess_heat` submodule below additionally provides the supporting
daily mean temperature, significance and acclimatisation indices, and heatwave
severity.

::: thermofeel.excess_heat
    options:
      show_root_heading: false
      show_root_toc_entry: false
      heading_level: 3

## Approximations

The `thermofeel.approximations` submodule provides clearly-labelled **estimators**
(not validated index formulas), kept out of the top-level namespace so calling
code always shows the `.approximations.` marker. It currently estimates the
direct solar radiation `fdir` (needed for the mean radiant temperature, and hence
UTCI / WBGT / PMV) from the global horizontal radiation `ssrd` and solar geometry,
for datasets that lack `fdir`. These are demonstration/estimation grade — see the
[Approximations guide](guide/approximations.md).

::: thermofeel.approximations
    options:
      show_root_heading: false
      show_root_toc_entry: false
      heading_level: 3

## Unit converters

::: thermofeel.helpers
    options:
      show_root_heading: false
      show_root_toc_entry: false
      heading_level: 3

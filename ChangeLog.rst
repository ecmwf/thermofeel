ChangeLog
=========

Version 2.1.0
-------------

 * Fix computaiton of normal effective temperature
 * Fix documentation links to new repository URL
 
Version 2.0.0
-------------

**Standardisation**
 * I/O variables converted into International System of Units (SI) (e.g., K)
 * Docstrings standardised in their descriptions (e.g., inputs as float arrays) and references
 * Variable names standardised across function declarations (i.e., t2m vs t2k vs t_k → now t2_k)
 * Use of temperature converter functions
 * Input variables are limited to those explicitly needed in the function formula (see e.g., ``calculate_normal_effective_temperature`` function)

**New Functions**
 * ``scale_windspeed`` function to scale 10m wind speed to height h with h < 10m
 * ``calculate_nonsaturation_vapour_pressure`` to calculate vapour pressure at a given relative humidity and air temperature
 * ``calculate_wbgt_simple`` function renamed
 * ``calculate_normal_effective_temperature`` function renamed

**Improvements**
 * thermofeel library docstring lists computed variables in alphabetical order
 * the cosine of the solar zenith angle now computed via the `earthkit-meteo <https://github.com/ecmwf/earthkit-meteo>`_ library
 * ``calculate_bgt`` function is calculated via a 4x faster formula
 * ``calculate_saturation_vapour_pressure_multiphase`` formulas replaced with those used in the IFS
 * changeable threshold in ``approximate_dsrp`` function (set to 0.1 by default)
 * invalidity outside input variables range specified in ``calculate_wind_chill`` function docstring

**Bug Fixes**
 * fixed ``approximate_dsrp`` to avoid fdir being overwritten with dsrp when calculating MRT
 * fixed ``calculate_wbgt_simple`` constant value and vapour pressure calculated from non-saturated formula
 * fixed ``calculate_bgt`` wind speed at 1.1m in input
 * fixed ``calculate_mrt_from_bgt`` wind speed at 1.1m in input
 * fixed ``calculate_normal_effective_temperature`` wind speed at 1.2m in input
 * fixed ``calculate_apparent_temperature`` wind speed at 10m in input; vapour pressure calculated from non-saturated formula
 * fixed ``calculate_wind_chill`` wind speed in km/h and operation symbols in main formula
 * fixed ``calculate_heat_index_simplified`` a wrong sign and missing constant value in hiarray; hi set to 2m air temperature when the latter is below 20°C
 * fixed ``fahrenheit_to_kelvin`` converter function 

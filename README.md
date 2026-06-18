![thermofeel logo](https://raw.githubusercontent.com/ecmwf/thermofeel/main/thermofeel.png)

[![license](https://img.shields.io/github/license/ecmwf/thermofeel)](https://www.apache.org/licenses/LICENSE-2.0.html)
[![tag release](https://img.shields.io/github/v/release/ecmwf/thermofeel?sort=semver)](https://github.com/ecmwf/thermofeel)
[![docs](https://readthedocs.org/projects/thermofeel/badge/?version=latest)](https://thermofeel.readthedocs.io/en/latest/?badge=latest)
[![ci](https://img.shields.io/github/actions/workflow/status/ecmwf/thermofeel/ci.yml)](https://github.com/ecmwf/thermofeel/actions)

# thermofeel

**thermofeel** (pronounced *thermo-feel*)

A library to calculate human thermal comfort indices.

Currently calculates the thermal indices:

- Universal Thermal Climate Index
- Apparent Temperature
- Heat Index Adjusted
- Heat Index Simplified
- Humidex
- Normal Effective Temperature
- Wet Bulb Globe Temperature
- Wet Bulb Globe Temperature Simple
- Wet Bulb Globe Temperature (Liljegren method)
- Heat Force (KNMI 0–10 heat-stress scale)
- Wind Chill

In support of the above indices, it also calculates:

- Globe Temperature
- Mean Radiant Temperature
- Mean Radiant Temperature from Globe Temperature
- Relative Humidity Percentage
- Saturation vapour pressure
- Wet Bulb Temperature

## PyPI

[![pypi status](https://img.shields.io/pypi/status/thermofeel)](https://pypi.org/project/thermofeel)
[![pypi release](https://img.shields.io/pypi/v/thermofeel?color=green)](https://pypi.org/project/thermofeel)
[![pypi downloads](https://img.shields.io/pypi/dm/thermofeel)](https://pypi.org/project/thermofeel)
[![code size](https://img.shields.io/github/languages/code-size/ecmwf/thermofeel?color=green)](https://github.com/ecmwf/thermofeel)

Install with:

```bash
pip install thermofeel
```

## System dependencies

thermofeel core functions depend on:

- numpy

Optionally, thermofeel depends on:

- pytest — for unit testing

## Release notes

Thermofeel 2.0 brings a number of changes to the underlying code but most
importantly to the API.

Consequently, downstream packages using thermofeel 1.\* will require code changes
to migrate to version 2.0 and beyond.

The main changes are:

- standardisation of input and output variables
- standardisation of variable names
- removal of dependency on numba for code acceleration
- removal of solar zenith angle calculation (now provided by earthkit-meteo)
- several bug fixes and improvements

Please consult the [ChangeLog](https://github.com/ecmwf/thermofeel/blob/main/CHANGELOG.md)
for more details.

## Contributing

The main repository is hosted on
[GitHub](https://github.com/ecmwf/thermofeel). Testing, bug reports and
contributions are highly welcomed and appreciated.

Please see the
[Contributing](https://github.com/ecmwf/thermofeel/blob/main/CONTRIBUTING.md)
document for the best way to help.

Current developers:

- Claudia Di Napoli — [ECMWF](https://ecmwf.int)
- Tiago Quintino — [ECMWF](https://ecmwf.int)

See also the
[contributors](https://github.com/ecmwf/thermofeel/contributors) for a more
complete list.

## License

Copyright 2021 European Centre for Medium-Range Weather Forecasts (ECMWF)

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at

> http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation nor
does it submit to any jurisdiction.

## Citing

In publications, please use our paper in SoftwareX as the main citation for
**thermofeel**:

> Brimicombe, C., Di Napoli, C., Quintino, T., Pappenberger, F., Cornforth, R., & Cloke, H. L. (2022).
> Thermofeel: A python thermal comfort indices library. SoftwareX, 18, 101005.
> <https://doi.org/10.1016/j.softx.2022.101005>

To cite **thermofeel** the code currently please use:

> Brimicombe, C., Di Napoli, C., Quintino, T., Pappenberger, F., Cornforth, R., & Cloke, H. L. (2021).
> *thermofeel: a python thermal comfort indices library* <https://doi.org/10.21957/mp6v-fd16>

## Acknowledgements

Past and current funding and support for **thermofeel** is listed in the
adjoining
[Acknowledgements](https://github.com/ecmwf/thermofeel/blob/main/ACKNOWLEDGEMENTS.md).

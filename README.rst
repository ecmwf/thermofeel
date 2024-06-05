.. image:: https://raw.githubusercontent.com/ecmwf/thermofeel/master/thermofeel.png
  :width: 600
  :alt: thermofeel logo

|license| |tag_release| |docs| |ci|

**thermofeel** (pronounced *thermo-feel*)

A library to calculate human thermal comfort indexes.

Currently calculates the thermal indexes:
  * Universal Thermal Climate Index
  * Apparent Temperature
  * Heat Index Adjusted
  * Heat Index Simplified
  * Humidex
  * Normal Effective Temperature
  * Wet Bulb Globe Temperature
  * Wet Bulb Globe Temperature Simple
  * Wind Chill

In support of the above indexes, it also calculates:
  * Globe Temperature
  * Mean Radiant Temperature
  * Mean Radiant Temperature from Globe Temperature
  * Relative Humidity Percentage
  * Saturation vapour pressure
  * Wet Bulb Temperature

PyPi
====

|pypi_status|  |pypi_release| |pypi_downloads| |code_size|

Install with::

    $ pip install thermofeel

System dependencies
===================

thermofeel core functions depend on:
 * numpy
 * earthkit-meteo > 0.0.1 - for solar zenith angle calculation

Optionally, thermofeel depends on:
 * pytest - for unit testing


Release notes
=============

Thermofeel 2.0 brings a number of changes to the underlying code but most importantly to the API.

Consequently, downstream packages using thermofeel 1.* will require code changes to migrate to version 2.0 and beyond.

The main changes are:
 * standardisation of input and output variables
 * standardisation of variable names
 * removal of dependency on numba for code acceleration
 * removal of solar zenith angle calculation (now provided by earthkit-meteo)
 * several bug fixes and improvements

Please consult ChangeLog_ for more details.

.. _ChangeLog: https://github.com/ecmwf/thermofeel/blob/master/ChangeLog.rst


Contributing
============

The main repository is hosted on `GitHub <https://github.com/ecmwf/thermofeel>`_. Testing, bug reports and contributions are highly welcomed and appreciated.

Please see the Contributing_ document for the best way to help.

.. _Contributing: https://github.com/ecmwf/thermofeel/blob/master/CONTRIBUTING.rst

Current developers:

- Claudia Di Napoli - `ECMWF <https://ecmwf.int>`_
- Tiago Quintino - `ECMWF <https://ecmwf.int>`_

See also the `contributors <https://github.com/ecmwf/thermofeel/contributors>`_ for a more complete list.

License
=======

Copyright 2021 European Centre for Medium-Range Weather Forecasts (ECMWF)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation nor
does it submit to any jurisdiction.

Citing
======


In publications, please use our paper in SoftwareX as the main citation for **thermofeel**:

Brimicombe, C., Di Napoli, C., Quintino, T., Pappenberger, F., Cornforth, R., & Cloke, H. L. (2022). 
Thermofeel: A python thermal comfort indices library. SoftwareX, 18, 101005. 
https://doi.org/10.1016/j.softx.2022.101005


To cite **thermofeel** the code currently please use:

Brimicombe, C., Di Napoli, C., Quintino, T., Pappenberger, F., Cornforth, R., & Cloke, H. L. (2021).
*thermofeel: a python thermal comfort indices library* https://doi.org/10.21957/mp6v-fd16


Acknowledgements
================
Past and current funding and support for **thermofeel** is listed in the adjoning Acknowledgements_


.. _Acknowledgements: https://github.com/ecmwf/thermofeel/blob/master/ACKNOWLEDGEMENTS.rst


.. |last_commit| image:: https://img.shields.io/github/last-commit/ecmwf/thermofeel
    :target: https://github.com/ecmwf/thermofeel

.. |commits_since_release| image:: https://img.shields.io/github/commits-since/ecmwf/thermofeel/latest?sort=semver
    :target: https://github.com/ecmwf/thermofeel

.. |license| image:: https://img.shields.io/github/license/ecmwf/thermofeel
    :target: https://www.apache.org/licenses/LICENSE-2.0.html

.. |pypi_release| image:: https://img.shields.io/pypi/v/thermofeel?color=green
    :target: https://pypi.org/project/thermofeel

.. |pypi_status| image:: https://img.shields.io/pypi/status/thermofeel
    :target: https://pypi.org/project/thermofeel

.. |tag_release| image:: https://img.shields.io/github/v/release/ecmwf/thermofeel?sort=semver
    :target: https://github.com/ecmwf/thermofeel

.. |codecov| image:: https://codecov.io/gh/ecmwf/thermofeel/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/ecmwf/thermofeel

.. |ci| image:: https://img.shields.io/github/actions/workflow/status/ecmwf/thermofeel/ci.yml
  :target: https://github.com/ecmwf/thermofeel/actions

.. |pypi_downloads| image:: https://img.shields.io/pypi/dm/thermofeel
  :target: https://pypi.org/project/thermofeel

.. |code_size| image:: https://img.shields.io/github/languages/code-size/ecmwf/thermofeel?color=green
  :target: https://github.com/ecmwf/thermofeel
  
.. |docs| image:: https://readthedocs.org/projects/thermofeel/badge/?version=latest
  :target: https://thermofeel.readthedocs.io/en/latest/?badge=latest


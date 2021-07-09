.. image:: https://raw.githubusercontent.com/ecmwf-projects/thermofeel/master/thermofeel.png
  :width: 600
  :alt: thermofeel logo

|license| |tag_release| |commits_since_release| |last_commit|

**thermofeel** (pronounced *thermo-feel*)

A library to calculate human thermal comfort indexes.

Currently calculates the thermal indexes:
 * Universal Thermal Climate Index
 * Mean Radiant Temperature
 * Mean Radiant Temperature from Wet Bulb Globe Temperature
 * Heat Index Simplified
 * Heat Index Adjusted
 * Humidex
 * Apparent Temperature
 * Wind Chill
 * Net Effective Temperature
 * Wet Bulb Globe Temperature Simple
 * Wet Bulb Globe Temperature
 
In support of the above indexes, it also calculates:
 * Solar Declination Angle
 * Solar Zenith Angle
 * Relative Humidity Percentage
 * Saturation vapour pressure


PyPi
====

|pypi_status|  |pypi_release| |pypi_downloads| |code_size| 

Install with::

    $ pip install thermofeel

Testing
=======

|python-publish| |test-coverage| |codecov| 

System dependencies
===================

thermofeel core functions depend on:
 * numpy

Contributing
============

The main repository is hosted on GitHub, testing, bug reports and contributions are highly welcomed and appreciated:

https://github.com/ecmwf-projects/thermofeel

Please see the Contributing_ document for the best way to help.

.. _Contributing: https://github.com/ecmwf-projects/thermofeel/blob/master/CONTRIBUTING.rst

Main contributors:

- Chloe Brimicombe - `ECMWF <https://ecmwf.int>`_
- Tiago Quintino - `ECMWF <https://ecmwf.int>`_

See also the `contributors <https://github.com/ecmwf-projects/thermofeel/contributors>`_ for a more complete list.


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

..
  In publications, please use our paper in SoftwareX as the main citation for **thermofeel**. 
  
To cite **thermofeel** the code currently please use: 
Brimicombe C,Di Napoli C, Quintino T,Pappenberger F, Cornforth R and Cloke H,2021 
*thermofeel: a python thermal comfort indices library* https://doi.org/10.21957/mp6v-fd16

For referring to the latest release of **thermofeel** please use this DOI: https://doi.org/10.21957/mp6v-fd16



Acknowledgements
================
Past and current funding and support for **thermofeel** is listed in the adjoning Acknowledgements_


.. _Acknowledgements: https://github.com/ecmwf-projects/thermofeel/blob/master/ACKNOWLEDGEMENTS.rst


.. |last_commit| image:: https://img.shields.io/github/last-commit/ecmwf-projects/thermofeel
    :target: https://github.com/ecmwf-projects/thermofeel

.. |commits_since_release| image:: https://img.shields.io/github/commits-since/ecmwf-projects/thermofeel/latest?sort=semver
    :target: https://github.com/ecmwf-projects/thermofeel

.. |license| image:: https://img.shields.io/github/license/ecmwf-projects/thermofeel
    :target: https://www.apache.org/licenses/LICENSE-2.0.html

.. |pypi_release| image:: https://img.shields.io/pypi/v/thermofeel?color=green
    :target: https://pypi.org/project/thermofeel

.. |pypi_status| image:: https://img.shields.io/pypi/status/thermofeel
    :target: https://pypi.org/project/thermofeel

.. |tag_release| image:: https://img.shields.io/github/v/release/ecmwf-projects/thermofeel?sort=semver
    :target: https://github.com/ecmwf-projects/thermofeel

.. |codecov| image:: https://codecov.io/gh/ecmwf-projects/thermofeel/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/ecmwf-projects/thermofeel

.. |python-publish| image:: https://github.com/ecmwf-projects/thermofeel/actions/workflows/python-publish.yml/badge.svg
  :target: https://github.com/ecmwf-projects/thermofeel/actions

.. |test-coverage| image:: https://github.com/ecmwf-projects/thermofeel/actions/workflows/test-coverage.yml/badge.svg
  :target: https://github.com/ecmwf-projects/thermofeel/actions

.. |pypi_downloads| image:: https://img.shields.io/pypi/dm/thermofeel
  :target: https://pypi.org/project/thermofeel

.. |code_size| image:: https://img.shields.io/github/languages/code-size/ecmwf-projects/thermofeel?color=green
  :target: https://github.com/ecmwf-projects/thermofeel




  

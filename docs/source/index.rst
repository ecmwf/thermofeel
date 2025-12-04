.. thermofeel documentation master file, created by
   sphinx-quickstart on Mon Jun  7 11:08:20 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to thermofeel's documentation!
======================================
thermofeel is a python library developed by the European Centre for Medium-Range Weather Forecasts (ECMWF) to integrate thermal (heat and cold) indexes \
into our internal systems and then ongoing maintenance.

It allows users to:

* integrate our methods into their Numerical Weather Prediction (NWP) services, and
* calculate thermal indices from their own data.

This allows for transparency and reproducibility, and allows users to expand \
the scope of these methods.


Indices and tables
==================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   guide/overview
   guide/apparenttemperature
   guide/heatindex
   guide/humidex
   guide/net
   guide/utci
   guide/wbgt
   guide/windchill
   guide/relativehumidity
   guide/mrt

License
-------

*thermofeel* is available under the open source `Apache License`__. In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.

__ http://www.apache.org/licenses/LICENSE-2.0.html

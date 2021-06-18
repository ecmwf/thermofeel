# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import re

version_file = "thermofeel/_version.py"
version_line = open(version_file, "rt").read()
version_regex = r"^__version__ = ['\"]([^'\"]*)['\"]"
res = re.search(version_regex, version_line, re.M)
if res:
    version_string = res.group(1)
else:
    raise RuntimeError("Cannot parse version string from {file}.".format(file=version_file))

# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Lightweight, offline smoke checks for the scripts in ``examples/``.

The examples are not part of the shipped library, so these are best-effort: the
compile check runs everywhere; the ``--help`` check for the earthkit-based
example is skipped unless the optional ``examples`` dependencies are installed;
none of them touches the network.
"""

import pathlib
import subprocess
import sys

import pytest

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"

pytestmark = pytest.mark.skipif(
    not EXAMPLES.is_dir(), reason="examples/ not available (installed package)"
)

EXAMPLE_SCRIPTS = sorted(EXAMPLES.glob("*.py")) if EXAMPLES.is_dir() else []


@pytest.mark.parametrize("script", EXAMPLE_SCRIPTS, ids=lambda p: p.name)
def test_example_compiles(script):
    # every example must be syntactically valid Python (parses without importing)
    import py_compile

    py_compile.compile(str(script), doraise=True)


def test_compute_heat_force_runs():
    # the self-contained example needs only numpy + thermofeel, so it runs
    # end to end offline in the test environment.
    result = subprocess.run(
        [sys.executable, str(EXAMPLES / "compute-heat-force.py")],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stderr
    assert "heat force" in result.stdout.lower()


def test_compute_thermal_indices_help():
    # offline --help smoke test for the earthkit-based example; needs the
    # optional [examples] dependencies, so it is skipped in the numpy-only test
    # environment. --help never touches the network.
    pytest.importorskip("earthkit.data")
    pytest.importorskip("earthkit.meteo")
    pytest.importorskip("xarray")
    result = subprocess.run(
        [sys.executable, str(EXAMPLES / "compute-thermal-indices.py"), "--help"],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
    out = result.stdout.lower()
    assert "usage" in out
    assert "--source" in out
    assert "--approximate-fdir" in out

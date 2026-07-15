# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Matplotlib defaults and figure-saving helper for the campaign plots."""

from pathlib import Path


def use_style():
    """Apply consistent, compact defaults for all campaign figures."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 110,
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "legend.fontsize": 8,
            "figure.constrained_layout.use": True,
        }
    )
    return plt


def save_png(fig, path):
    """Save a figure as a compact PNG and report its size."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    kb = path.stat().st_size / 1024.0
    print(f"  plot: {path} ({kb:.0f} KB)")
    if kb > 400:
        print(f"  !! warning: {path.name} is large ({kb:.0f} KB)")

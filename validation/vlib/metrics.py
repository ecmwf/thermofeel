# (C) Copyright 2026- ECMWF and individual contributors.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Standard verification metrics used across the validation campaign.

Conventions (as in the solar / verification literature):

- error = predicted - observed
- MBE   = mean(error)                      (bias; same unit as the variable)
- MAE   = mean(|error|)
- RMSE  = sqrt(mean(error^2))
- rMBE, rRMSE = MBE, RMSE normalised by |mean(observed)|, in percent
  (the absolute value keeps the relative metrics' signs meaningful for
  quantities whose mean can be negative, e.g. PMV)
- R2    = squared Pearson correlation coefficient
"""

import numpy as np


def error_stats(pred, obs):
    """Return the standard error statistics of ``pred`` against ``obs``.

    NaNs are dropped pairwise. Returns a dict with N, mean_obs, mean_pred,
    MBE, rMBE_pct, MAE, RMSE, rRMSE_pct and R2.
    """
    pred = np.asarray(pred, dtype=float).ravel()
    obs = np.asarray(obs, dtype=float).ravel()
    if pred.shape != obs.shape:
        raise ValueError(f"shape mismatch: {pred.shape} vs {obs.shape}")
    ok = np.isfinite(pred) & np.isfinite(obs)
    pred, obs = pred[ok], obs[ok]
    n = pred.size
    if n == 0:
        raise ValueError("no valid (finite) pairs to score")
    err = pred - obs
    mean_obs = float(np.mean(obs))
    mbe = float(np.mean(err))
    mae = float(np.mean(np.abs(err)))
    rmse = float(np.sqrt(np.mean(err**2)))
    if np.std(pred) > 0 and np.std(obs) > 0:
        r2 = float(np.corrcoef(pred, obs)[0, 1] ** 2)
    else:
        r2 = float("nan")
    rel = 100.0 / abs(mean_obs) if mean_obs != 0 else float("nan")
    return {
        "N": int(n),
        "mean_obs": mean_obs,
        "mean_pred": float(np.mean(pred)),
        "MBE": mbe,
        "rMBE_pct": mbe * rel,
        "MAE": mae,
        "RMSE": rmse,
        "rRMSE_pct": rmse * rel,
        "R2": r2,
    }


def stats_frame(rows):
    """Build a pandas DataFrame from ``(label_dict, stats_dict)`` pairs."""
    import pandas as pd

    records = []
    for labels, stats in rows:
        rec = dict(labels)
        rec.update(stats)
        records.append(rec)
    return pd.DataFrame.from_records(records)


def write_results(df, path, float_format="%.4g"):
    """Write a results table as CSV (versioned) and echo it to stdout."""
    from pathlib import Path

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, float_format=float_format)
    print(f"\n== {path} ==")
    print(df.to_string(index=False, float_format=lambda v: float_format % v))

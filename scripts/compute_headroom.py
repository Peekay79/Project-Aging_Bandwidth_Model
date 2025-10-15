"""
Compute headroom as Capacity − Load at sample level and aggregate.
"""
from __future__ import annotations
import numpy as np
import pandas as pd


def compute_headroom(
    capacity_per_sample: pd.DataFrame,
    load_per_sample: pd.DataFrame,
) -> pd.DataFrame:
    """Compute headroom at sample level and aggregate.

    Parameters
    - capacity_per_sample: ['tissue','age_months','sample_id','capacity_sum','capacity_z']
    - load_per_sample: ['tissue','age_months','sample_id','load_sum','load_z']

    Returns
    - DataFrame with ['tissue','age_months','sample_id','headroom_z'] and group means
    """
    c = capacity_per_sample.rename(columns={"capacity_z": "capacity"})
    l = load_per_sample.rename(columns={"load_z": "load"})
    merged = c.merge(l[["tissue", "age_months", "sample_id", "load"]],
                     on=["tissue", "age_months", "sample_id"], how="outer")
    merged["capacity"] = merged["capacity"].fillna(0.0)
    merged["load"] = merged["load"].fillna(0.0)
    merged["headroom_z"] = merged["capacity"] - merged["load"]
    return merged


def aggregate_headroom(headroom_per_sample: pd.DataFrame) -> pd.DataFrame:
    """Aggregate sample-level headroom to per-tissue × per-age with mean and 95% CI.

    Returns DataFrame with ['tissue','age_months','headroom_mean','headroom_ci_lo','headroom_ci_hi','n']
    """
    def _agg(g: pd.DataFrame) -> pd.Series:
        n = g.shape[0]
        m = g["headroom_z"].mean()
        sd = g["headroom_z"].std(ddof=0)
        se = 0.0 if n == 0 else sd / max(n, 1) ** 0.5
        return pd.Series({
            "headroom_mean": m,
            "headroom_ci_lo": m - 1.96 * se,
            "headroom_ci_hi": m + 1.96 * se,
            "n": n,
        })

    return headroom_per_sample.groupby(["tissue", "age_months"], as_index=False).apply(_agg)


if __name__ == "__main__":
    # Minimal smoke test
    c = pd.DataFrame({"tissue": ["Brain"], "age": [3], "capacity": [10.0]})
    l = pd.DataFrame({"tissue": ["Brain"], "age": [3], "load": [4.0]})
    print(compute_headroom(c, l))

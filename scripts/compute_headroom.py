"""
Compute headroom as capacity minus load.
"""
from __future__ import annotations
import pandas as pd


def compute_headroom(capacity_df: pd.DataFrame, load_df: pd.DataFrame) -> pd.DataFrame:
    """Merge capacity and load, compute headroom per tissue×age.

    Parameters
    - capacity_df: DataFrame with ['tissue', 'age', 'capacity']
    - load_df: DataFrame with ['tissue', 'age', 'load']

    Returns
    - DataFrame with ['tissue', 'age', 'capacity', 'load', 'headroom']
    """
    merged = capacity_df.merge(load_df, on=["tissue", "age"], how="outer")
    merged["capacity"] = merged["capacity"].fillna(0.0)
    merged["load"] = merged["load"].fillna(0.0)
    merged["headroom"] = merged["capacity"] - merged["load"]
    return merged


if __name__ == "__main__":
    # Minimal smoke test
    c = pd.DataFrame({"tissue": ["Brain"], "age": [3], "capacity": [10.0]})
    l = pd.DataFrame({"tissue": ["Brain"], "age": [3], "load": [4.0]})
    print(compute_headroom(c, l))

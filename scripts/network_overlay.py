"""
Stub for network overlay on headroom results.

Adds a simple within-tissue z-score of headroom as an example
"network"-style normalization.
"""
from __future__ import annotations
import pandas as pd


def overlay_network_metrics(headroom_df: pd.DataFrame) -> pd.DataFrame:
    """Add a simple per-tissue z-score of headroom as 'headroom_z'.

    This is a placeholder demonstrating how a network overlay step
    could enrich the headroom results with additional metrics.
    """
    df = headroom_df.copy()
    def zscore(s: pd.Series) -> pd.Series:
        mu = s.mean()
        sd = s.std(ddof=0)
        if sd == 0:
            return (s - mu)
        return (s - mu) / sd

    df["headroom_z"] = df.groupby("tissue")["headroom"].transform(zscore)
    return df


if __name__ == "__main__":
    demo = pd.DataFrame({
        "tissue": ["Brain"] * 3 + ["Liver"] * 3,
        "age": [3, 12, 24] * 2,
        "headroom": [10.0, 9.0, 8.0, 7.0, 7.5, 6.5],
    })
    print(overlay_network_metrics(demo))

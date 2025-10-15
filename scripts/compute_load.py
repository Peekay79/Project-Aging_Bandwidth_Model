"""
Compute proteostasis load at the sample level and aggregate.

Definition (documented in docs/ASSUMPTIONS.md):
- Load per sample = z-scored sum of normalized abundances for proteins
  annotated as aggregation-prone and generic stress markers when available.
  Fallback: include the top 10% most disordered proteins via proxy column
  `high_disorder_proxy` in annotations, or by within-tissue missingness rank.
"""
from __future__ import annotations
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd


def compute_load(
    proteome_long_df: pd.DataFrame,
    annotations_df: pd.DataFrame,
    value_col: str = "protein_abundance_z",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute load at sample level and aggregate to tissue×age.

    Parameters
    - proteome_long_df: tidy DataFrame with ['tissue','age_months','sample_id','gene_symbol', value_col]
    - annotations_df: DataFrame with ['gene_symbol','is_aggregation_prone','high_disorder_proxy', ...]
    - value_col: abundance column to use (default z-scored abundances)

    Returns
    - per_sample: ['tissue','age_months','sample_id','load_sum','load_z']
    - per_group: ['tissue','age_months','load_mean','load_ci_lo','load_ci_hi','n']
    """
    needed = {"tissue", "age_months", "sample_id", "gene_symbol", value_col}
    missing = needed - set(proteome_long_df.columns)
    if missing:
        raise ValueError(f"proteome_long_df missing columns: {sorted(missing)}")

    ann = annotations_df.copy()
    if "is_aggregation_prone" not in ann.columns:
        ann["is_aggregation_prone"] = False
    if "high_disorder_proxy" not in ann.columns:
        ann["high_disorder_proxy"] = False
    ann = ann[["gene_symbol", "is_aggregation_prone", "high_disorder_proxy"]]

    merged = proteome_long_df.merge(ann, on="gene_symbol", how="left")
    merged["is_aggregation_prone"] = merged["is_aggregation_prone"].fillna(False)
    merged["high_disorder_proxy"] = merged["high_disorder_proxy"].fillna(False)

    # Missingness-based fallback marker: compute per tissue gene missingness
    # Rank genes by fraction missing; top 10% included as fallback
    def _mark_top_missing(g: pd.DataFrame) -> pd.DataFrame:
        miss_frac = (
            g.groupby("gene_symbol")[value_col]
            .apply(lambda s: s.isna().mean())
            .rename("missing_frac")
        )
        if len(miss_frac) == 0:
            g["fallback_missing"] = False
            return g
        cutoff = np.quantile(miss_frac.values, 0.90)
        high = miss_frac[miss_frac >= cutoff].index
        g = g.copy()
        g["fallback_missing"] = g["gene_symbol"].isin(set(high))
        return g

    merged = merged.groupby("tissue", group_keys=False).apply(_mark_top_missing)

    merged["is_load_marker"] = (
        merged["is_aggregation_prone"]
        | merged["high_disorder_proxy"]
        | merged["fallback_missing"]
    )

    load_only = merged.loc[merged["is_load_marker"]]
    per_sample_sum = (
        load_only.groupby(["tissue", "age_months", "sample_id"], as_index=False)[value_col]
        .sum()
        .rename(columns={value_col: "load_sum"})
    )

    # z-score within tissue
    def _z(s: pd.Series) -> pd.Series:
        mu = s.mean()
        sd = s.std(ddof=0)
        if sd == 0 or np.isnan(sd):
            return s * 0.0
        return (s - mu) / sd
    per_sample_sum["load_z"] = per_sample_sum.groupby("tissue")["load_sum"].transform(_z)

    def _agg(g: pd.DataFrame) -> pd.Series:
        n = g.shape[0]
        m = g["load_z"].mean()
        sd = g["load_z"].std(ddof=0)
        se = 0.0 if n == 0 else sd / max(n, 1) ** 0.5
        return pd.Series({
            "load_mean": m,
            "load_ci_lo": m - 1.96 * se,
            "load_ci_hi": m + 1.96 * se,
            "n": n,
        })

    per_group = per_sample_sum.groupby(["tissue", "age_months"], as_index=False).apply(_agg)
    return per_sample_sum, per_group


if __name__ == "__main__":
    base_dir = Path(__file__).resolve().parents[1]
    tidy_path = base_dir / "data" / "processed" / "mouse_proteome_long.csv"
    ann_path = base_dir / "data" / "protein_annotations_seed.csv"
    if tidy_path.exists() and ann_path.exists():
        p = pd.read_csv(tidy_path)
        a = pd.read_csv(ann_path)
        per_sample, per_group = compute_load(p, a)
        print(per_group.head())
    else:
        print("Data not found. This module is intended to be imported by the pipeline.")

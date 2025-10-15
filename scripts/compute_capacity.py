"""
Compute proteostasis capacity at the sample level and aggregate.

Definition (documented in docs/ASSUMPTIONS.md):
- Capacity per sample = z-scored sum of normalized abundances for
  proteins annotated as chaperone OR UPS OR autophagy. We expect the
  input abundance to already be normalized per prepare_data.py, and we
  apply a z-score across samples within each tissue to the summed value.

Returns both per-sample capacity and a per-tissue×age summary.
"""
from __future__ import annotations
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd


def compute_capacity(
    proteome_long_df: pd.DataFrame,
    annotations_df: pd.DataFrame,
    qc_columns: Optional[List[str]] = None,
    value_col: str = "protein_abundance_z",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute capacity at sample level and aggregate to tissue×age.

    Parameters
    - proteome_long_df: tidy DataFrame with ['tissue','age_months','sample_id','gene_symbol', value_col]
    - annotations_df: DataFrame with ['gene_symbol', 'is_chaperone','is_UPS','is_autophagy', ...]
    - qc_columns: Optional list of QC annotation columns to consider
    - value_col: abundance column to use (default: z-scored abundances)

    Returns
    - per_sample: ['tissue','age_months','sample_id','capacity_sum','capacity_z']
    - per_group: ['tissue','age_months','capacity_mean','capacity_ci_lo','capacity_ci_hi','n']
    """
    if qc_columns is None:
        qc_columns = ["is_chaperone", "is_UPS", "is_autophagy"]

    # Minimal columns sanity
    needed = {"tissue", "age_months", "sample_id", "gene_symbol", value_col}
    missing = needed - set(proteome_long_df.columns)
    if missing:
        raise ValueError(f"proteome_long_df missing columns: {sorted(missing)}")

    ann = annotations_df.copy()
    ann_cols = ["gene_symbol"] + [c for c in qc_columns if c in ann.columns]
    ann = ann[ann_cols]
    for c in qc_columns:
        if c not in ann.columns:
            ann[c] = False

    merged = proteome_long_df.merge(ann, on="gene_symbol", how="left")
    for c in qc_columns:
        merged[c] = merged[c].fillna(False)

    merged["is_qc"] = merged[qc_columns].any(axis=1)
    qc_only = merged.loc[merged["is_qc"]]

    per_sample_sum = (
        qc_only.groupby(["tissue", "age_months", "sample_id"], as_index=False)[value_col]
        .sum()
        .rename(columns={value_col: "capacity_sum"})
    )

    # z-score capacity within tissue
    def _z(s: pd.Series) -> pd.Series:
        mu = s.mean()
        sd = s.std(ddof=0)
        if sd == 0 or np.isnan(sd):
            return s * 0.0
        return (s - mu) / sd
    per_sample_sum["capacity_z"] = per_sample_sum.groupby("tissue")["capacity_sum"].transform(_z)

    # Aggregate to per tissue × age with CI (normal approx)
    def _agg(g: pd.DataFrame) -> pd.Series:
        n = g.shape[0]
        m = g["capacity_z"].mean()
        sd = g["capacity_z"].std(ddof=0)
        se = 0.0 if n == 0 else sd / max(n, 1) ** 0.5
        return pd.Series({
            "capacity_mean": m,
            "capacity_ci_lo": m - 1.96 * se,
            "capacity_ci_hi": m + 1.96 * se,
            "n": n,
        })

    per_group = per_sample_sum.groupby(["tissue", "age_months"], as_index=False).apply(_agg)
    return per_sample_sum, per_group


if __name__ == "__main__":
    # Minimal demo if run directly (expects processed tidy CSV and seed annotations)
    base_dir = Path(__file__).resolve().parents[1]
    data_dir = base_dir / "data"
    tidy_path = data_dir / "processed" / "mouse_proteome_long.csv"
    ann_path = data_dir / "protein_annotations_seed.csv"
    if tidy_path.exists() and ann_path.exists():
        p = pd.read_csv(tidy_path)
        a = pd.read_csv(ann_path)
        per_sample, per_group = compute_capacity(p, a)
        print(per_group.head())
    else:
        print("Data not found. This module is intended to be imported by the pipeline.")

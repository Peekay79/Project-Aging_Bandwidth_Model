"""
Compute proteostasis capacity per tissue and age.

Capacity is defined as the summed abundance of proteins involved in
protein quality control (chaperones, UPS, autophagy) for each
(tissue, age) group.
"""
from __future__ import annotations
from pathlib import Path
from typing import List, Optional

import pandas as pd


def compute_capacity(
    proteome_df: pd.DataFrame,
    annotations_df: pd.DataFrame,
    qc_columns: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Compute capacity by summing QC protein abundances per tissue×age.

    Parameters
    - proteome_df: DataFrame with columns ['tissue', 'age', 'protein_id', 'abundance']
    - annotations_df: DataFrame with columns ['protein_id', 'is_chaperone', 'is_UPS', 'is_autophagy', ...]
    - qc_columns: Optional list of QC annotation columns to consider (defaults to chaperone/UPS/autophagy)

    Returns
    - DataFrame with columns ['tissue', 'age', 'capacity']
    """
    if qc_columns is None:
        qc_columns = ["is_chaperone", "is_UPS", "is_autophagy"]

    cols_to_keep = ["protein_id"] + [c for c in qc_columns if c in annotations_df.columns]
    merged = proteome_df.merge(annotations_df[cols_to_keep], on="protein_id", how="left")

    # Treat missing annotations as False
    for c in qc_columns:
        if c not in merged:
            merged[c] = False
        merged[c] = merged[c].fillna(False)

    merged["is_qc"] = merged[qc_columns].any(axis=1)
    capacity = (
        merged.loc[merged["is_qc"]]
        .groupby(["tissue", "age"], as_index=False)["abundance"].sum()
        .rename(columns={"abundance": "capacity"})
    )
    return capacity


if __name__ == "__main__":
    # Minimal demo if run directly (expects CSVs in ../data)
    base_dir = Path(__file__).resolve().parents[1]
    data_dir = base_dir / "data"
    proteome_path = data_dir / "mouse_proteome_mock.csv"
    ann_path = data_dir / "protein_annotations_mock.csv"
    if proteome_path.exists() and ann_path.exists():
        p = pd.read_csv(proteome_path)
        a = pd.read_csv(ann_path)
        cap = compute_capacity(p, a)
        print(cap.head())
    else:
        print("Data not found. This module is intended to be imported by the pipeline.")

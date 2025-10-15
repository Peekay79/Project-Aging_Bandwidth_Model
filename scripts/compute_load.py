"""
Compute proteostasis load per tissue and age.

Load is defined as the summed abundance of aggregation-prone proteins
for each (tissue, age) group.
"""
from __future__ import annotations
from pathlib import Path

import pandas as pd


def compute_load(proteome_df: pd.DataFrame, annotations_df: pd.DataFrame) -> pd.DataFrame:
    """Compute load by summing aggregation-prone protein abundances per tissue×age.

    Parameters
    - proteome_df: DataFrame with columns ['tissue', 'age', 'protein_id', 'abundance']
    - annotations_df: DataFrame with columns including ['protein_id', 'is_aggregation_prone']

    Returns
    - DataFrame with columns ['tissue', 'age', 'load']
    """
    cols_to_keep = ["protein_id", "is_aggregation_prone"]
    cols_to_keep = [c for c in cols_to_keep if c in annotations_df.columns]
    merged = proteome_df.merge(annotations_df[cols_to_keep], on="protein_id", how="left")
    if "is_aggregation_prone" not in merged:
        merged["is_aggregation_prone"] = False
    merged["is_aggregation_prone"] = merged["is_aggregation_prone"].fillna(False)

    load = (
        merged.loc[merged["is_aggregation_prone"]]
        .groupby(["tissue", "age"], as_index=False)["abundance"].sum()
        .rename(columns={"abundance": "load"})
    )
    return load


if __name__ == "__main__":
    base_dir = Path(__file__).resolve().parents[1]
    data_dir = base_dir / "data"
    proteome_path = data_dir / "mouse_proteome_mock.csv"
    ann_path = data_dir / "protein_annotations_mock.csv"
    if proteome_path.exists() and ann_path.exists():
        p = pd.read_csv(proteome_path)
        a = pd.read_csv(ann_path)
        ld = compute_load(p, a)
        print(ld.head())
    else:
        print("Data not found. This module is intended to be imported by the pipeline.")

"""
Network overlay utilities.

Phase 3 requirements:
- Load data/ppi_edges.tsv, build a graph
- Compute degree and betweenness centrality for proteins present
- For each tissue, list top 50 central proteins whose abundance increases with age
  and contribute positively to Load or negatively to Capacity.
"""
from __future__ import annotations
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
import networkx as nx


def overlay_network_metrics(headroom_per_sample: pd.DataFrame) -> pd.DataFrame:
    """Pass-through for now; headroom already z-scored.

    Placeholder: In future, we could weight by network centrality.
    """
    return headroom_per_sample


def compute_network_hotspots(
    proteome_long_df: pd.DataFrame,
    annotations_df: pd.DataFrame,
    outputs_dir: Path,
) -> Path:
    """Compute network hotspots per tissue and write outputs/network_hotspots.csv.

    Criteria:
    - Proteins present in the data
    - Positive age slope in abundance within tissue
    - Contribute to Load (aggregation/stress) or reduce Capacity (negative member)
      — approximated by annotation flags
    - Ranked by graph centrality (degree, betweenness)
    """
    ppi_path = Path(__file__).resolve().parents[1] / "data" / "ppi_edges.tsv"
    if not ppi_path.exists():
        raise FileNotFoundError("data/ppi_edges.tsv not found. Run scripts/download_data.py first.")

    edges = pd.read_csv(ppi_path, sep="\t")
    G = nx.Graph()
    for _, r in edges.iterrows():
        G.add_edge(str(r["proteinA"]), str(r["proteinB"]), weight=float(r.get("score", 0)))

    # Centrality (limited to nodes we observe)
    observed = set(proteome_long_df["gene_symbol"].astype(str).unique())
    sub_nodes = [n for n in G.nodes() if n in observed]
    SG = G.subgraph(sub_nodes).copy()
    degree = nx.degree_centrality(SG)
    # Sample limited betweenness for scalability
    k_sample = min(1000, max(50, len(SG))) if len(SG) > 50 else None
    if k_sample:
        between = nx.betweenness_centrality(SG, k=k_sample, seed=42)
    else:
        between = nx.betweenness_centrality(SG)

    ann = annotations_df.copy()
    load_like_cols = [c for c in ["is_aggregation_prone", "high_disorder_proxy"] if c in ann.columns]
    cap_like_cols = [c for c in ["is_chaperone", "is_UPS", "is_autophagy"] if c in ann.columns]
    for c in load_like_cols + cap_like_cols:
        if c not in ann.columns:
            ann[c] = False
    ann = ann[["gene_symbol"] + load_like_cols + cap_like_cols]

    def _slope(x, y) -> float:
        x = pd.Series(x).astype(float)
        y = pd.Series(y).astype(float)
        if x.nunique() < 2:
            return 0.0
        try:
            m, b = np.polyfit(x, y, deg=1)
            return float(m)
        except Exception:
            return 0.0

    trends = (
        proteome_long_df.groupby(["tissue", "gene_symbol"], as_index=False)
        .apply(lambda g: pd.Series({"age_slope": _slope(g["age_months"].values, g["protein_abundance_z"].values)}))
    )
    df = trends.merge(ann, on="gene_symbol", how="left")
    for c in load_like_cols + cap_like_cols:
        df[c] = df[c].fillna(False)

    def _centrality_row(gene: str) -> Tuple[float, float]:
        return degree.get(gene, 0.0), between.get(gene, 0.0)

    df["degree_cent"], df["betweenness_cent"] = zip(*df["gene_symbol"].map(_centrality_row))

    # Filter: increasing with age and load-like or anti-capacity
    df = df[df["age_slope"] > 0]
    if load_like_cols:
        df = df[df[load_like_cols].any(axis=1)]

    df["centrality_score"] = df["degree_cent"] * 0.5 + df["betweenness_cent"] * 0.5
    out = (
        df.sort_values(["tissue", "centrality_score"], ascending=[True, False])
        .groupby("tissue")
        .head(50)
    )

    outputs_dir.mkdir(parents=True, exist_ok=True)
    out_path = outputs_dir / "network_hotspots.csv"
    cols = [
        "tissue", "gene_symbol", "age_slope", "centrality_score",
        "degree_cent", "betweenness_cent"
    ] + load_like_cols + cap_like_cols
    out[cols].to_csv(out_path, index=False)
    return out_path


if __name__ == "__main__":
    base = Path(__file__).resolve().parents[1]
    tidy = base / "data" / "processed" / "mouse_proteome_long.csv"
    ann = base / "data" / "protein_annotations_seed.csv"
    if not tidy.exists():
        raise SystemExit("Run scripts/prepare_data.py first")
    df = pd.read_csv(tidy)
    if ann.exists():
        adf = pd.read_csv(ann)
    else:
        # Fallback to mock annotations if seed missing
        ann_mock = base / "data" / "protein_annotations_mock.csv"
        if not ann_mock.exists():
            raise SystemExit("Missing annotations: data/protein_annotations_seed.csv")
        adf = pd.read_csv(ann_mock)
        if "gene_symbol" not in adf.columns and "protein_id" in adf.columns:
            adf = adf.rename(columns={"protein_id": "gene_symbol"})
    out = compute_network_hotspots(df, adf, base / "outputs")
    print(out)

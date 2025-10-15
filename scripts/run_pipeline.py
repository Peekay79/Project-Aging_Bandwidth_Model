"""
End-to-end pipeline: generate mock data, compute capacity/load/headroom,
plot results, and write a short Markdown summary.

Run with:
    python scripts/run_pipeline.py
"""
from __future__ import annotations
from pathlib import Path
from typing import Tuple, List

import numpy as np
import pandas as pd

from compute_capacity import compute_capacity
from compute_load import compute_load
from compute_headroom import compute_headroom
from plot_headroom import plot_headroom_lines
from network_overlay import overlay_network_metrics


# ------------------------------
# Configuration
# ------------------------------
TISSUES: List[str] = [
    "Brain",
    "Heart",
    "Liver",
    "Kidney",
    "Lung",
    "Muscle",
    "Spleen",
    "Intestine",
]
AGES: List[int] = [3, 12, 18, 24]  # months
NUM_PROTEINS: int = 1000
RNG_SEED: int = 42


# ------------------------------
# Data generation
# ------------------------------

def _generate_mock_proteome(base_dir: Path) -> Path:
    data_dir = base_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    proteome_path = data_dir / "mouse_proteome_mock.csv"

    rng = np.random.default_rng(RNG_SEED)

    protein_ids = [f"P{idx:04d}" for idx in range(1, NUM_PROTEINS + 1)]

    records = []
    for tissue in TISSUES:
        for age in AGES:
            # base abundance distribution per protein
            # Log-normal for positivity and dynamic range
            abundance = rng.lognormal(mean=2.0, sigma=0.6, size=NUM_PROTEINS)
            # Introduce subtle age trend to make plots interesting
            age_factor = 1.0 - (age - min(AGES)) / (max(AGES) - min(AGES)) * 0.15
            abundance = abundance * age_factor
            for pid, ab in zip(protein_ids, abundance):
                records.append({
                    "tissue": tissue,
                    "age": age,
                    "protein_id": pid,
                    "abundance": float(ab),
                })

    df = pd.DataFrame.from_records(records)
    df.to_csv(proteome_path, index=False)
    return proteome_path


def _generate_mock_annotations(base_dir: Path) -> Path:
    data_dir = base_dir / "data"
    ann_path = data_dir / "protein_annotations_mock.csv"
    rng = np.random.default_rng(RNG_SEED)

    protein_ids = [f"P{idx:04d}" for idx in range(1, NUM_PROTEINS + 1)]

    # Assign flags with reasonable sparsity
    is_chaperone = rng.random(NUM_PROTEINS) < 0.04
    is_UPS = rng.random(NUM_PROTEINS) < 0.06
    is_autophagy = rng.random(NUM_PROTEINS) < 0.03
    is_aggregation_prone = rng.random(NUM_PROTEINS) < 0.10

    ann = pd.DataFrame({
        "protein_id": protein_ids,
        "is_chaperone": is_chaperone,
        "is_UPS": is_UPS,
        "is_autophagy": is_autophagy,
        "is_aggregation_prone": is_aggregation_prone,
    })
    ann.to_csv(ann_path, index=False)
    return ann_path


def ensure_mock_data(base_dir: Path) -> Tuple[Path, Path]:
    proteome_path = base_dir / "data" / "mouse_proteome_mock.csv"
    ann_path = base_dir / "data" / "protein_annotations_mock.csv"

    if not proteome_path.exists():
        _generate_mock_proteome(base_dir)
    if not ann_path.exists():
        _generate_mock_annotations(base_dir)

    return proteome_path, ann_path


# ------------------------------
# Pipeline
# ------------------------------

def main() -> None:
    base_dir = Path(__file__).resolve().parents[1]
    outputs_dir = base_dir / "outputs"
    docs_dir = base_dir / "docs"
    outputs_dir.mkdir(parents=True, exist_ok=True)
    docs_dir.mkdir(parents=True, exist_ok=True)

    # 1) Ensure mock data
    proteome_csv, ann_csv = ensure_mock_data(base_dir)

    # 2) Load data
    proteome_df = pd.read_csv(proteome_csv)
    ann_df = pd.read_csv(ann_csv)

    # 3) Compute metrics
    capacity_df = compute_capacity(proteome_df, ann_df)
    load_df = compute_load(proteome_df, ann_df)
    headroom_df = compute_headroom(capacity_df, load_df)

    # 4) Optional overlay
    headroom_df = overlay_network_metrics(headroom_df)

    # 5) Plot
    plot_paths = plot_headroom_lines(headroom_df, outputs_dir)

    # 6) Write short summary Markdown
    summary_md = _write_summary(headroom_df, docs_dir)

    print("Pipeline finished successfully.")
    print(f"Plots saved to: {outputs_dir}")
    print(f"Summary written to: {summary_md}")


def _write_summary(headroom_df: pd.DataFrame, docs_dir: Path) -> Path:
    # Compute simple trend: slope of headroom vs age per tissue
    lines: List[str] = []
    lines.append("# Headroom trends (mock data)\n")
    lines.append("This report summarizes headroom (capacity − load) trends across ages.\n")

    ages_sorted = sorted(AGES)
    for tissue, sub in headroom_df.groupby("tissue"):
        sub = sub.sort_values("age")
        x = sub["age"].values.astype(float)
        y = sub["headroom"].values.astype(float)
        if len(np.unique(x)) >= 2:
            slope, intercept = np.polyfit(x, y, deg=1)
            direction = "decreases" if slope < 0 else ("increases" if slope > 0 else "is flat")
            lines.append(f"- {tissue}: headroom {direction} with age (slope={slope:.3f}).")
        else:
            lines.append(f"- {tissue}: insufficient age points to estimate trend.")

    out_path = docs_dir / "summary.md"
    with out_path.open("w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return out_path


if __name__ == "__main__":
    main()

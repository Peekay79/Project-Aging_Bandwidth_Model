"""
Plot headroom vs. age per tissue and combined.
"""
from __future__ import annotations
from pathlib import Path
from typing import List

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import pandas as pd


def plot_headroom_lines(headroom_df: pd.DataFrame, output_dir: Path) -> List[Path]:
    """Create and save line plots of headroom vs age.

    - Per-tissue plots: one file per tissue
    - Combined plot: all tissues together

    Returns list of saved file paths.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Ensure correct dtypes
    df = headroom_df.copy()
    df["age"] = pd.to_numeric(df["age"])  # expected months

    saved: List[Path] = []

    # Per-tissue plots
    for tissue, sub in df.groupby("tissue"):
        sub = sub.sort_values("age")
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(sub["age"], sub["headroom"], marker="o", linewidth=2)
        ax.set_title(f"Headroom vs Age — {tissue}")
        ax.set_xlabel("Age (months)")
        ax.set_ylabel("Headroom (capacity − load)")
        ax.grid(True, alpha=0.3)
        fp = output_dir / f"headroom_{tissue.replace(' ', '_')}.png"
        fig.tight_layout()
        fig.savefig(fp, dpi=150)
        plt.close(fig)
        saved.append(fp)

    # Combined plot
    fig, ax = plt.subplots(figsize=(8, 5))
    for tissue, sub in df.groupby("tissue"):
        sub = sub.sort_values("age")
        ax.plot(sub["age"], sub["headroom"], marker="o", linewidth=2, label=tissue)
    ax.set_title("Headroom vs Age — All Tissues")
    ax.set_xlabel("Age (months)")
    ax.set_ylabel("Headroom (capacity − load)")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2, fontsize=8)
    combined_fp = output_dir / "headroom_all_tissues.png"
    fig.tight_layout()
    fig.savefig(combined_fp, dpi=150)
    plt.close(fig)
    saved.append(combined_fp)

    return saved


if __name__ == "__main__":
    # Small demo with fake data
    demo = pd.DataFrame(
        {
            "tissue": ["Brain"] * 4 + ["Liver"] * 4,
            "age": [3, 12, 18, 24] * 2,
            "headroom": [10, 9, 7, 6, 8, 7, 7, 5],
        }
    )
    out = Path(__file__).resolve().parents[1] / "outputs"
    paths = plot_headroom_lines(demo, out)
    print("Saved:", ", ".join(str(p) for p in paths))

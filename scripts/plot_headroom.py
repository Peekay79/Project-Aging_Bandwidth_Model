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
        fig.savefig(output_dir / f"headroom_{tissue.replace(' ', '_')}.svg")
        plt.close(fig)
        saved.append(fp)

    # Combined plot (PNG + SVG)
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
    fig.savefig(output_dir / "headroom_all_tissues.svg")
    plt.close(fig)
    saved.append(combined_fp)

    # Small multiples grid
    tissues = sorted(df["tissue"].unique().tolist())
    n = len(tissues)
    ncols = 4 if n >= 8 else (3 if n >= 6 else 2)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 3*nrows), sharex=True, sharey=True)
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]
    for i, tissue in enumerate(tissues):
        ax = axes[i]
        sub = df[df["tissue"] == tissue].sort_values("age")
        ax.plot(sub["age"], sub["headroom"], marker="o", linewidth=2)
        ax.set_title(tissue, fontsize=10)
        ax.grid(True, alpha=0.3)
    # Hide any unused axes
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    fig.suptitle("Headroom vs Age — Small Multiples", fontsize=14)
    for ax in axes[-ncols:]:
        ax.set_xlabel("Age (months)")
    for ax in axes[::ncols]:
        ax.set_ylabel("Headroom")
    grid_fp = output_dir / "headroom_all_tissues.png"
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(grid_fp, dpi=150)
    plt.close(fig)
    saved.append(grid_fp)

    return saved


def plot_inflection_bars(inflection_csv: Path, output_dir: Path) -> List[Path]:
    """Plot a bar chart of inflection ages across tissues if CSV exists.

    Expects columns: ['tissue','model','inflection_age','aic']
    """
    if not inflection_csv.exists():
        return []
    df = pd.read_csv(inflection_csv)
    df = df.dropna(subset=["inflection_age"]).copy()
    if df.empty:
        return []
    df = df.sort_values("inflection_age")
    fig, ax = plt.subplots(figsize=(8, max(4, 0.3 * len(df))))
    ax.barh(df["tissue"], df["inflection_age"], color="#4C72B0")
    ax.set_xlabel("Inflection age (months)")
    ax.set_title("Estimated inflection ages by tissue (AIC-selected model)")
    fig.tight_layout()
    out_png = output_dir / "inflection_ages_bar.png"
    out_svg = output_dir / "inflection_ages_bar.svg"
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_svg)
    plt.close(fig)
    return [out_png, out_svg]


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

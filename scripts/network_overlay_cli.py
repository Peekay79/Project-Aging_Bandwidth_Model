#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import pandas as pd
from network_overlay import compute_network_hotspots

if __name__ == "__main__":
    base = Path(__file__).resolve().parents[1]
    tidy = base / "data" / "processed" / "mouse_proteome_long.csv"
    ann = base / "data" / "protein_annotations_seed.csv"
    if not tidy.exists():
        raise SystemExit("Run scripts/prepare_data.py first")
    df = pd.read_csv(tidy)
    adf = pd.read_csv(ann)
    out = compute_network_hotspots(df, adf, base / "outputs")
    print(out)

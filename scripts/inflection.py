#!/usr/bin/env python3
"""
Fit inflection models for headroom vs age per tissue and export results.

Models per tissue:
 a) Piecewise linear with one breakpoint (segmented regression)
 b) Logistic/sigmoid

Model selection by AIC. Export estimated inflection age or breakpoint per tissue to
outputs/inflection_points.csv
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

BASE_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = BASE_DIR / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _aic(n: int, rss: float, k: int) -> float:
    return float(n * np.log(rss / n if n > 0 else 1.0) + 2 * k)


def _piecewise(x, x0, k1, k2, b):
    return np.where(x < x0, k1 * x + b, k2 * x + (k1 - k2) * x0 + b)


def _logistic(x, L, x0, k, b):
    # L / (1 + exp(-k (x - x0))) + b
    return L / (1 + np.exp(-k * (x - x0))) + b


def fit_piecewise(x, y) -> Tuple[float, Dict[str, float]]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(np.unique(x)) < 3:
        return np.nan, {"rss": np.inf, "params": {}}
    x0_init = np.median(x)
    k1_init = 0.0
    k2_init = 0.0
    b_init = float(np.mean(y))
    try:
        params, _ = curve_fit(_piecewise, x, y, p0=[x0_init, k1_init, k2_init, b_init], maxfev=10000)
        preds = _piecewise(x, *params)
        rss = float(np.sum((y - preds) ** 2))
        return float(params[0]), {"rss": rss, "params": {"x0": params[0], "k1": params[1], "k2": params[2], "b": params[3]}}
    except Exception:
        return np.nan, {"rss": np.inf, "params": {}}


def fit_logistic(x, y) -> Tuple[float, Dict[str, float]]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(np.unique(x)) < 3:
        return np.nan, {"rss": np.inf, "params": {}}
    L_init = float(np.max(y) - np.min(y) + 1e-6)
    x0_init = float(np.median(x))
    k_init = 0.1
    b_init = float(np.min(y))
    try:
        params, _ = curve_fit(_logistic, x, y, p0=[L_init, x0_init, k_init, b_init], maxfev=10000)
        preds = _logistic(x, *params)
        rss = float(np.sum((y - preds) ** 2))
        # inflection at x0
        return float(params[1]), {"rss": rss, "params": {"L": params[0], "x0": params[1], "k": params[2], "b": params[3]}}
    except Exception:
        return np.nan, {"rss": np.inf, "params": {}}


def main() -> Path:
    tidy_path = BASE_DIR / "data" / "processed" / "mouse_proteome_long.csv"
    if not tidy_path.exists():
        raise FileNotFoundError("Run scripts/prepare_data.py first")

    # Recompute headroom aggregates to ensure up to date
    from compute_capacity import compute_capacity
    from compute_load import compute_load
    from compute_headroom import compute_headroom, aggregate_headroom

    ann_path = BASE_DIR / "data" / "protein_annotations_seed.csv"
    if not ann_path.exists():
        raise FileNotFoundError("data/protein_annotations_seed.csv not found")

    proteome_df = pd.read_csv(tidy_path)
    ann_df = pd.read_csv(ann_path)

    cap_sample, _ = compute_capacity(proteome_df, ann_df)
    load_sample, _ = compute_load(proteome_df, ann_df)
    headroom_sample = compute_headroom(cap_sample, load_sample)
    headroom_group = aggregate_headroom(headroom_sample).rename(columns={"headroom_mean": "headroom"})

    records: List[Dict[str, float]] = []
    for tissue, sub in headroom_group.groupby("tissue"):
        sub = sub.sort_values("age_months")
        x = sub["age_months"].values
        y = sub["headroom"].values
        n = len(x)
        if n < 3:
            records.append({"tissue": tissue, "model": "NA", "inflection_age": np.nan, "aic": np.nan})
            continue
        # Fit models
        bp, m1 = fit_piecewise(x, y)
        x0, m2 = fit_logistic(x, y)
        aic1 = _aic(n, m1["rss"], k=4)
        aic2 = _aic(n, m2["rss"], k=4)
        if aic1 <= aic2:
            records.append({"tissue": tissue, "model": "piecewise", "inflection_age": bp, "aic": aic1})
        else:
            records.append({"tissue": tissue, "model": "logistic", "inflection_age": x0, "aic": aic2})

    out = pd.DataFrame.from_records(records)
    out_path = OUT_DIR / "inflection_points.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")
    return out_path


if __name__ == "__main__":
    main()

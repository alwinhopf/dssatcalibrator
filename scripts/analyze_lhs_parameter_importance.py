"""Summarize parameter importance and contraction in an evaluated LHS design."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.feature_selection import mutual_info_regression

from dssatcalibrator.config import load_config
from dssatcalibrator.spaces import ParameterSpace


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("design")
    parser.add_argument("output")
    parser.add_argument("--top-fraction", type=float, default=0.1)
    parser.add_argument("--trees", type=int, default=300)
    args = parser.parse_args()

    cfg = load_config(args.config)
    space = ParameterSpace.from_config(cfg)
    design = pd.read_csv(args.design)
    missing = [name for name in space.names if name not in design]
    if missing:
        raise ValueError(f"Design is missing active parameters: {missing}")
    valid = design[np.isfinite(pd.to_numeric(design["score"], errors="coerce"))].copy()
    if len(valid) < 20:
        raise ValueError("At least 20 finite evaluated candidates are required.")

    x = valid[space.names].astype(float)
    y = valid["score"].astype(float).to_numpy()
    x_scaled = (x.to_numpy() - space.low) / np.maximum(space.high - space.low, 1e-12)
    mutual = mutual_info_regression(x_scaled, y, random_state=1348)
    forest = ExtraTreesRegressor(
        n_estimators=int(args.trees),
        min_samples_leaf=max(2, len(valid) // 1000),
        max_features=1.0,
        random_state=1348,
        n_jobs=-1,
    ).fit(x_scaled, y)

    fraction = float(args.top_fraction)
    if not 0 < fraction < 1:
        raise ValueError("--top-fraction must be between zero and one.")
    threshold = float(np.quantile(y, fraction))
    top = valid[valid["score"] <= threshold]
    rows = []
    for index, (name, spec) in enumerate(zip(space.names, space.specs)):
        rho, p_value = spearmanr(x[name], y)
        prior_width = float(space.high[index] - space.low[index])
        top_width = float(top[name].quantile(0.95) - top[name].quantile(0.05))
        rows.append({
            "parameter": name,
            "group": spec["group"],
            "scope": spec.get("scope", "global"),
            "cultivar": spec.get("cultivar", ""),
            "n_candidates": len(valid),
            "top_fraction": fraction,
            "spearman_rho": float(rho),
            "spearman_abs": abs(float(rho)),
            "spearman_p": float(p_value),
            "mutual_information": float(mutual[index]),
            "extra_trees_importance": float(forest.feature_importances_[index]),
            "prior_width": prior_width,
            "top_p05_p95_width": top_width,
            "top_width_fraction": top_width / prior_width if prior_width > 0 else np.nan,
            "top_median": float(top[name].median()),
            "all_median": float(valid[name].median()),
        })

    out = pd.DataFrame(rows)
    for column in ("spearman_abs", "mutual_information", "extra_trees_importance"):
        out[f"{column}_rank"] = out[column].rank(ascending=False, method="average")
    out["importance_mean_rank"] = out[
        ["spearman_abs_rank", "mutual_information_rank", "extra_trees_importance_rank"]
    ].mean(axis=1)
    out = out.sort_values(
        ["importance_mean_rank", "top_width_fraction"], ascending=[True, True]
    ).reset_index(drop=True)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    print(out.head(30).to_string(index=False))


if __name__ == "__main__":
    main()

"""Select high-performing, separated starting points from an evaluated LHS."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dssatcalibrator.config import load_config
from dssatcalibrator.spaces import ParameterSpace


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("design")
    parser.add_argument("output_dir")
    parser.add_argument("--count", type=int, default=3)
    parser.add_argument("--pool", type=int, default=50)
    args = parser.parse_args()

    space = ParameterSpace.from_config(load_config(args.config))
    design = pd.read_csv(args.design)
    design["score"] = pd.to_numeric(design["score"], errors="coerce")
    pool = design[np.isfinite(design["score"])].nsmallest(args.pool, "score").copy()
    if pool.empty:
        raise ValueError("No finite LHS candidates are available.")
    x = pool[space.names].to_numpy(float)
    scale = np.maximum(space.high - space.low, 1e-12)
    z = (x - space.low) / scale
    selected = [0]
    while len(selected) < min(args.count, len(pool)):
        distance = np.min(
            np.stack([np.linalg.norm(z - z[index], axis=1) for index in selected]),
            axis=0,
        )
        distance[selected] = -1.0
        selected.append(int(np.argmax(distance)))

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for rank, index in enumerate(selected, start=1):
        row = pool.iloc[index]
        theta = {name: float(row[name]) for name in space.names}
        path = outdir / f"start_{rank}.json"
        path.write_text(json.dumps(theta, indent=2) + "\n", encoding="utf-8")
        rows.append({
            "start": rank,
            "sample_id": int(row["sample_id"]),
            "lhs_score": float(row["score"]),
            "path": str(path),
        })
    pd.DataFrame(rows).to_csv(outdir / "selected_starts.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()

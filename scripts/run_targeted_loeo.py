"""Run and persist a lightweight leave-one-environment-out calibration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from dssatcalibrator.config import load_config
from dssatcalibrator.cv import run_cross_validation


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    cfg = load_config(args.config)
    cfg["cross_validation"] = {
        "enabled": True,
        "strategy": "leave_one_out",
        "final_theta": "best_fold",
        "overfit_threshold": 1.5,
        "report": True,
    }
    result = run_cross_validation(cfg, progress=True)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    result.report_df.to_csv(outdir / "cv_report.csv", index=False)
    (outdir / "cv_summary.json").write_text(
        json.dumps(result.summary, indent=2) + "\n", encoding="utf-8"
    )
    (outdir / "recommended_theta.json").write_text(
        json.dumps(result.final_theta, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result.summary, indent=2))


if __name__ == "__main__":
    main()

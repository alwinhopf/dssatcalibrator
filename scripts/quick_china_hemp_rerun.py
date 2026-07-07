from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from dssatcalibrator import orchestrator, viz  # noqa: E402
from dssatcalibrator.config import load_config  # noqa: E402


DATE_TAG = "20260704"
RESULTS_DIR = Path("results/china_hemp_calibration")
FIGURES_DIR = Path("figures/china_hemp_calibration")


def quick_common(cfg: dict, name: str) -> dict:
    cfg = deepcopy(cfg)
    cfg["calibrator"]["name"] = name
    cfg["calibrator"]["results_dir"] = str(RESULTS_DIR)
    cfg["calibrator"]["figures_dir"] = str(FIGURES_DIR)
    cfg["calibrator"]["workdir"] = "results/_workdir_quick"
    cfg["calibrator"]["cache_spawns"] = False
    cfg["calibrator"]["keep_run_dirs"] = False
    cfg["calibrator"]["num_cores"] = 0
    cfg["calibrator"]["batch_size"] = 18
    cfg["calibrator"]["spawn_timeout"] = max(60, int(cfg["calibrator"].get("spawn_timeout", 60)))
    cfg["method"].setdefault("sensitivity", {})["active"] = False
    cfg["method"].setdefault("select", {})["active"] = False
    cfg["method"].setdefault("surrogate", {})["active"] = False
    cfg["method"].setdefault("multiobjective", {})["engine"] = "none"
    return cfg


def configure_stage1() -> dict:
    cfg = load_config("calibration_china_hemp/stage1_emergence_anthesis_cultivar.yaml")
    cfg = quick_common(cfg, f"china_hemp_quick_stage1_phenology_{DATE_TAG}")
    cfg["method"]["preset"] = "B"
    cfg["method"].setdefault("bayesian", {})["engine"] = "none"
    cfg["method"].setdefault("optimizer", {}).update({
        "engine": "cmaes",
        "maxiter": 2,
        "popsize": 6,
        "restarts": 1,
    })
    return cfg


def configure_stage2(stage1_theta: dict) -> dict:
    cfg = load_config("calibration_china_hemp/stage2_biomass_after_phenology_cultivar.yaml")
    cfg = quick_common(cfg, f"china_hemp_quick_stage2_joint_{DATE_TAG}")
    cfg["method"]["preset"] = "B"
    cfg["method"].setdefault("bayesian", {})["engine"] = "none"
    cfg["method"].setdefault("optimizer", {}).update({
        "engine": "cmaes",
        "maxiter": 2,
        "popsize": 8,
        "restarts": 1,
    })

    # Keep the quick growth run honest: allow the key phenology / node-timing
    # coefficients to move, while still fitting biomass, LAI, height, width, and L#SD.
    cultivar_pheno = {"CSDL", "PPSEN", "EM-FL", "FL-SH", "FL-SD", "SD-PM"}
    ecotype_phases = {"PL-EM", "EM-V1", "V1-JU", "JU-R0"}
    for pname in cultivar_pheno:
        spec = cfg["parameters"]["genetic_cultivar"][pname]
        spec["active"] = True
        spec["fixed"] = False
        spec["scope"] = "cultivar"
        spec.pop("cultivars", None)
        # Stage 1 only has IB0008 phenology event data; use it as a warm start.
        val = stage1_theta.get(f"{pname}__IB0008", stage1_theta.get(pname))
        if val is not None:
            spec["start"] = float(val)
    for pname in ecotype_phases:
        spec = cfg["parameters"]["genetic_ecotype"][pname]
        spec["active"] = True
        spec["fixed"] = False
        spec["scope"] = "cultivar"
        spec.pop("cultivars", None)
        val = stage1_theta.get(f"{pname}__IB0008", stage1_theta.get(pname))
        if val is not None:
            spec["start"] = float(val)
    return cfg


def run_and_report(cfg: dict):
    name = cfg["calibrator"]["name"]
    print(f"\n=== Running {name} ===", flush=True)
    result = orchestrator.calibrate(cfg, progress=True)
    spawns = orchestrator.spawn_results_for(cfg, result.best_theta, result.experiments)
    outdir = RESULTS_DIR / name
    figdir = FIGURES_DIR / name
    paths = viz.make_report(result, outdir, best_spawns=spawns, figdir=figdir)

    (outdir / "quick_best_theta.json").write_text(json.dumps(result.best_theta, indent=2, sort_keys=True))
    if result.best.residuals is not None and not result.best.residuals.empty:
        resid = result.best.residuals.copy()
        resid.to_csv(outdir / "quick_residuals.csv", index=False)
        summary = resid.groupby(["exp_id", "user_var"], as_index=False).agg(
            n=("resid", "size"),
            rmse=("resid", lambda x: float((x.pow(2).mean()) ** 0.5)),
            mbe=("resid", "mean"),
            last_date=("date", "max"),
        )
        summary.to_csv(outdir / "quick_residual_summary_by_exp_var.csv", index=False)
        print(summary.round(3).to_string(index=False), flush=True)

    print("\nBest theta:", flush=True)
    for k, v in sorted(result.best_theta.items()):
        print(f"  {k}: {v:.5g}", flush=True)
    print("\nFit summary:", flush=True)
    print(viz.summary_fit_table(result).round(3).to_string(index=False), flush=True)
    print(f"Tables:  {(RESULTS_DIR / name).resolve()}", flush=True)
    print(f"Figures: {(FIGURES_DIR / name).resolve()}", flush=True)
    print(f"Experiment panels: {len(paths.get('experiment_diagnostics', []))}", flush=True)
    return result


def main() -> None:
    stage1 = run_and_report(configure_stage1())
    stage2 = run_and_report(configure_stage2(stage1.best_theta))
    print("\n=== Quick rerun complete ===", flush=True)
    print(json.dumps({
        "stage1": stage1.cfg["calibrator"]["name"],
        "stage2": stage2.cfg["calibrator"]["name"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()

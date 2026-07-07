from __future__ import annotations

import json
import os
import sys
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from dssatcalibrator import orchestrator, viz  # noqa: E402
from dssatcalibrator.config import load_config  # noqa: E402


DATE_TAG = os.environ.get("DSSATCAL_RUN_TAG", "20260704_allcores")
RESULTS_DIR = Path("results/china_hemp_calibration")
FIGURES_DIR = Path("figures/china_hemp_calibration")
CULTIVARS = ["IB0002", "IB0008"]

STAGE1_PHENO_PARAMS = {
    "genetic_cultivar": ["CSDL", "PPSEN", "EM-FL", "FL-SH", "FL-SD", "SD-PM"],
    "genetic_ecotype": ["THVAR", "PL-EM", "EM-V1", "V1-JU", "JU-R0", "R1PPO", "OPTBI", "SLOBI"],
}


def common(cfg: dict, name: str) -> dict:
    cfg = deepcopy(cfg)
    cfg["calibrator"]["name"] = name
    cfg["calibrator"]["results_dir"] = str(RESULTS_DIR)
    cfg["calibrator"]["figures_dir"] = str(FIGURES_DIR)
    cfg["calibrator"]["workdir"] = "results/_workdir_long_protected_allcores"
    cfg["calibrator"]["cache_spawns"] = False
    cfg["calibrator"]["keep_run_dirs"] = False
    cfg["calibrator"]["num_cores"] = int(os.environ.get("DSSATCAL_LONG_CORES", os.cpu_count() or 1))
    cfg["calibrator"]["batch_size"] = 24
    cfg["calibrator"]["spawn_timeout"] = max(90, int(cfg["calibrator"].get("spawn_timeout", 60)))
    cfg["method"].setdefault("sensitivity", {})["active"] = False
    cfg["method"].setdefault("select", {})["active"] = False
    cfg["method"].setdefault("surrogate", {})["active"] = False
    cfg["method"].setdefault("multiobjective", {})["engine"] = "none"
    cfg["method"]["preset"] = "B"
    cfg["method"].setdefault("bayesian", {})["engine"] = "none"
    cfg["crops"][0]["calibration_cultivars"] = list(CULTIVARS)
    return cfg


def widen_vegetative_bounds(cfg: dict) -> None:
    cultivar = cfg["parameters"]["genetic_cultivar"]
    ecotype = cfg["parameters"]["genetic_ecotype"]
    cultivar["EM-FL"]["max"] = max(float(cultivar["EM-FL"]["max"]), 110.0)
    ecotype["V1-JU"]["max"] = max(float(ecotype["V1-JU"]["max"]), 30.0)
    ecotype["JU-R0"]["max"] = max(float(ecotype["JU-R0"]["max"]), 35.0)


def configure_stage1() -> dict:
    cfg = load_config("calibration_china_hemp/stage1_emergence_anthesis_cultivar.yaml")
    cfg = common(cfg, f"china_hemp_long_stage1_phenology_node_{DATE_TAG}")
    widen_vegetative_bounds(cfg)
    cfg["engine"].setdefault("timeseries_outputs", {})["node_stage"] = "L#SD"
    cfg["objective"]["obs_autocorr"] = True
    cfg["objective"].setdefault("weights", {}).update({
        "emergence": 1.5,
        "anthesis": 2.0,
        "node_stage": 1.0,
    })
    cfg["objective"].setdefault("error_model", {})["node_stage"] = {
        "type": "absolute",
        "value": 1.0,
    }
    cfg["method"].setdefault("optimizer", {}).update({
        "engine": "cmaes",
        "maxiter": int(os.environ.get("DSSATCAL_STAGE1_MAXITER", 10)),
        "popsize": int(os.environ.get("DSSATCAL_STAGE1_POPSIZE", 12)),
        "restarts": 1,
    })
    return cfg


def _stage1_value(theta: dict, name: str, cultivar: str, fallback: float) -> float:
    return float(theta.get(f"{name}__{cultivar}", theta.get(name, fallback)))


def protect_stage1_phenology(cfg: dict, stage1_theta: dict) -> dict:
    """Freeze stage-1 phenology values per cultivar in the stage-2 genotype writes."""
    for group, names in STAGE1_PHENO_PARAMS.items():
        params = cfg["parameters"].get(group, {})
        for name in names:
            if name not in params:
                continue
            spec = params[name]
            spec["active"] = False
            spec["fixed"] = True
            spec["scope"] = "cultivar"
            spec.pop("cultivars", None)
            fallback = float(spec.get("start", 0.5 * (float(spec["min"]) + float(spec["max"]))))
            starts = {cultivar: _stage1_value(stage1_theta, name, cultivar, fallback) for cultivar in CULTIVARS}
            spec["start_by_cultivar"] = starts
            # Keep the generic start valid for config consumers that do not expand
            # cultivar-specific fixed specs.
            spec["start"] = starts.get("IB0008", next(iter(starts.values())))
    return cfg


def configure_stage2(stage1_theta: dict) -> dict:
    cfg = load_config("calibration_china_hemp/stage2_biomass_after_phenology_cultivar.yaml")
    cfg = common(cfg, f"china_hemp_long_stage2_biomass_protected_{DATE_TAG}")
    widen_vegetative_bounds(cfg)
    cfg["objective"]["ignore_zero_observations"] = ["width"]
    cfg["method"].setdefault("optimizer", {}).update({
        "engine": "cmaes",
        "maxiter": int(os.environ.get("DSSATCAL_STAGE2_MAXITER", 12)),
        "popsize": int(os.environ.get("DSSATCAL_STAGE2_POPSIZE", 14)),
        "restarts": 1,
    })
    return protect_stage1_phenology(cfg, stage1_theta)


def write_run_tables(result, outdir: Path, *, label: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{label}_best_theta.json").write_text(json.dumps(result.best_theta, indent=2, sort_keys=True))
    if result.best.residuals is None or result.best.residuals.empty:
        return
    resid = result.best.residuals.copy()
    resid.to_csv(outdir / f"{label}_residuals.csv", index=False)
    summary = resid.groupby(["exp_id", "user_var"], as_index=False).agg(
        n=("resid", "size"),
        rmse=("resid", lambda x: float((x.pow(2).mean()) ** 0.5)),
        mbe=("resid", "mean"),
        last_date=("date", "max"),
    )
    summary.to_csv(outdir / f"{label}_residual_summary_by_exp_var.csv", index=False)
    print(summary.round(3).to_string(index=False), flush=True)


def run_and_report(cfg: dict, *, label: str, extra_json: dict | None = None):
    name = cfg["calibrator"]["name"]
    print(f"\n=== Running {name} ===", flush=True)
    result = orchestrator.calibrate(cfg, progress=True)
    spawns = orchestrator.spawn_results_for(cfg, result.best_theta, result.experiments)
    outdir = RESULTS_DIR / name
    figdir = FIGURES_DIR / name
    paths = viz.make_report(result, outdir, best_spawns=spawns, figdir=figdir)
    write_run_tables(result, outdir, label=label)
    if extra_json:
        for filename, payload in extra_json.items():
            (outdir / filename).write_text(json.dumps(payload, indent=2, sort_keys=True))

    print("\nBest theta:", flush=True)
    for k, v in sorted(result.best_theta.items()):
        print(f"  {k}: {v:.5g}", flush=True)
    print("\nFit summary:", flush=True)
    print(viz.summary_fit_table(result).round(3).to_string(index=False), flush=True)
    print(f"Tables:  {outdir.resolve()}", flush=True)
    print(f"Figures: {figdir.resolve()}", flush=True)
    print(f"Experiment panels: {len(paths.get('experiment_diagnostics', []))}", flush=True)
    return result


def protected_values(cfg: dict) -> dict:
    out = {}
    for group, names in STAGE1_PHENO_PARAMS.items():
        for name in names:
            spec = cfg["parameters"].get(group, {}).get(name)
            if isinstance(spec, dict) and spec.get("start_by_cultivar"):
                out[f"{group}.{name}"] = spec["start_by_cultivar"]
    return out


def main() -> None:
    part = os.environ.get("DSSATCAL_RUN_PART", "both").lower()
    if part == "stage1":
        stage1 = run_and_report(configure_stage1(), label="stage1")
        print("\n=== Stage 1 complete ===", flush=True)
        print(json.dumps({"stage1": stage1.cfg["calibrator"]["name"]}, indent=2), flush=True)
        return
    if part == "stage2":
        theta_path = Path(os.environ["DSSATCAL_STAGE1_THETA"])
        stage1_theta = json.loads(theta_path.read_text())
        stage2_cfg = configure_stage2(stage1_theta)
        stage2 = run_and_report(
            stage2_cfg,
            label="stage2",
            extra_json={
                "protected_stage1_phenology.json": protected_values(stage2_cfg),
                "stage1_best_theta_used_for_protection.json": stage1_theta,
            },
        )
        print("\n=== Stage 2 complete ===", flush=True)
        print(json.dumps({"stage2": stage2.cfg["calibrator"]["name"]}, indent=2), flush=True)
        return

    stage1_cfg = configure_stage1()
    stage1 = run_and_report(stage1_cfg, label="stage1")
    stage2_cfg = configure_stage2(stage1.best_theta)
    stage2 = run_and_report(
        stage2_cfg,
        label="stage2",
        extra_json={
            "protected_stage1_phenology.json": protected_values(stage2_cfg),
            "stage1_best_theta_used_for_protection.json": stage1.best_theta,
        },
    )
    print("\n=== Protected long rerun complete ===", flush=True)
    print(json.dumps({
        "stage1": stage1.cfg["calibrator"]["name"],
        "stage2": stage2.cfg["calibrator"]["name"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()

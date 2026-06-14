"""Command-line entry point for a DSSAT calibration run.

    python run_calibration.py config_hemp.yaml --n 300
    python run_calibration.py config_hemp.yaml --n 50 --experiments YUKU2101 YUFE2201
    python run_calibration.py config_hemp.yaml --validate        # leave-one-env-out

Writes figures + CSV summaries to results/<name>/ (or --outdir).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from dssatcalibrator import orchestrator, viz
from dssatcalibrator.config import load_config


def main():
    ap = argparse.ArgumentParser(description="Run a DSSAT calibration.")
    ap.add_argument("config", help="path to the calibration YAML")
    ap.add_argument("--n", type=int, default=None, help="number of samples (overrides config)")
    ap.add_argument("--experiments", nargs="*", default=None, help="subset of experiments")
    ap.add_argument("--engine", default=None, help="sample engine: lhs|sobol|montecarlo|grid")
    ap.add_argument("--bayesian-engine", default=None, help="bayesian engine: glue|smc_pf")
    ap.add_argument("--n-particles", type=int, default=None,
                    help="SMC particle count (smc_pf only; overrides config)")
    ap.add_argument("--outdir", default=None, help="output directory for the report")
    ap.add_argument("--validate", action="store_true", help="run leave-one-environment-out validation")
    ap.add_argument("--no-progress", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    if args.n is not None:
        cfg["method"].setdefault("sample", {})["n"] = args.n
    if args.engine:
        cfg["method"].setdefault("sample", {})["engine"] = args.engine
    if args.bayesian_engine:
        cfg["method"].setdefault("bayesian", {})["engine"] = args.bayesian_engine
    if args.n_particles is not None:
        cfg["method"].setdefault("bayesian", {})["n_particles"] = args.n_particles
    if args.experiments:
        cfg["experiments"] = args.experiments

    name = cfg["calibrator"]["name"]
    outdir = Path(args.outdir or f"{cfg['calibrator'].get('results_dir', 'results')}/{name}")
    figdir = Path(cfg["calibrator"].get("figures_dir", "figures")) / name

    if args.validate:
        df = orchestrator.validate_loeo(cfg, progress=not args.no_progress)
        outdir.mkdir(parents=True, exist_ok=True)
        df.round(3).to_csv(outdir / "validation_loeo.csv", index=False)
        print(f"\nLeave-one-environment-out validation -> {outdir/'validation_loeo.csv'}")
        print(df.round(2).to_string(index=False))
        return

    bayes = cfg["method"].get("bayesian", {}).get("engine", "glue")
    if bayes == "smc_pf":
        size = f"smc_pf, n_particles={cfg['method'].get('bayesian', {}).get('n_particles', 200)}"
    else:
        size = f"glue, n={cfg['method']['sample']['n']}"
    print(f"Calibrating '{cfg['calibrator']['name']}' (preset {cfg['method'].get('preset')}; {size})")
    result = orchestrator.calibrate(cfg, progress=not args.no_progress)
    best_spawns = orchestrator.spawn_results_for(cfg, result.best_theta, result.experiments)
    paths = viz.make_report(result, outdir, best_spawns=best_spawns, figdir=figdir)

    print("\n=== Best-fit parameters ===")
    for k, v in result.best_theta.items():
        print(f"  {k:8s} {v:.4f}")
    print("\n=== Fit summary (best) ===")
    print(viz.summary_fit_table(result).to_string(index=False))
    print(f"\nData tables -> {outdir.resolve()}")
    print(f"Figures     -> {figdir.resolve()}")
    for k, p in paths.items():
        print(f"  {k}: {Path(p).relative_to(Path.cwd()) if Path(p).is_absolute() else p}")

    # clean up temporary work directories if configured
    if not cfg["calibrator"].get("keep_run_dirs", False):
        workdir = Path(cfg["calibrator"].get("workdir", "results/_workdir"))
        if workdir.exists():
            import shutil
            print(f"\nCleaning up temporary work files under {workdir.resolve()}...")
            shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()

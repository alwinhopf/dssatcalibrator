"""Command-line entry point for a DSSAT calibration run.

    python run_calibration.py config_hemp.yaml --n 300
    python run_calibration.py config_hemp.yaml --n 50 --experiments YUKU2101 YUFE2201
    python run_calibration.py config_hemp.yaml --validate        # leave-one-env-out

Writes figures + CSV summaries to results/<name>/ (or --outdir).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from dssatcalibrator import orchestrator, viz
from dssatcalibrator.config import load_config


def _write_assimilation(res: dict, outdir: Path) -> None:
    """Write assimilation results: a full JSON dump plus, for recalibration, a
    tidy per-checkpoint parameter trace CSV."""
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "assimilation.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2, default=str)
    if res.get("mode") == "recalibration" and res.get("trace"):
        rows = [{"date": t["date"], **(t.get("theta") or {})} for t in res["trace"]]
        pd.DataFrame(rows).to_csv(outdir / "assimilation_trace.csv", index=False)
        print(f"  recalibration trace -> {(outdir / 'assimilation_trace.csv')}")


def _write_forecast(forecasts: dict, outdir: Path) -> None:
    """Write per-variable, per-experiment forecast tables (LAI percentiles)."""
    if not forecasts:
        return
    outdir.mkdir(parents=True, exist_ok=True)
    for var, per_exp in forecasts.items():
        frames = []
        for exp, df in (per_exp or {}).items():
            if df is None or df.empty:
                continue
            d = df.copy()
            d.insert(0, "exp_id", exp)
            frames.append(d)
        if frames:
            out = pd.concat(frames, ignore_index=True)
            path = outdir / f"forecast_{var}.csv"
            out.round(4).to_csv(path, index=False)
            print(f"  forecast ({var}) -> {path}")


def main():
    ap = argparse.ArgumentParser(description="Run a DSSAT calibration.")
    ap.add_argument("config", help="path to the calibration YAML")
    ap.add_argument("--n", type=int, default=None, help="number of samples (overrides config)")
    ap.add_argument("--experiments", nargs="*", default=None, help="subset of experiments")
    ap.add_argument("--engine", default=None, help="sample engine: lhs|sobol|montecarlo|grid")
    ap.add_argument("--preset", default=None, help="pipeline preset: A|B|C|D (overrides config)")
    ap.add_argument("--bayesian-engine", default=None, help="estimator: glue|smc_pf|mcmc|none")
    ap.add_argument("--optimizer", default=None, help="optimizer engine: nelder_mead|diffevo")
    ap.add_argument("--sensitivity", default=None,
                    help="turn on screening: morris|sobol (auto-keeps influential params)")
    ap.add_argument("--select", default=None, help="stepwise selection: bic|aicc")
    ap.add_argument("--surrogate", default=None, help="emulator acceleration: gp|rf")
    ap.add_argument("--n-particles", type=int, default=None,
                    help="SMC particle count (smc_pf only; overrides config)")
    ap.add_argument("--outdir", default=None, help="output directory for the report")
    ap.add_argument("--validate", action="store_true", help="run leave-one-environment-out validation")
    ap.add_argument("--assimilate", action="store_true",
                    help="run in-season assimilation (mode from config; recalibration is the coupled path)")
    ap.add_argument("--assim-mode", default=None,
                    help="override assimilation.mode: recalibration | enkf | forcing")
    ap.add_argument("--combined", action="store_true",
                    help="calibrate, then run in-season assimilation seeded with the result")
    ap.add_argument("--nowcast", default=None, metavar="YYYY-MM-DD",
                    help="operational nowcast: (re)calibrate on data up to this date, persist, forecast")
    ap.add_argument("--forecast", action="store_true",
                    help="activate the in-season LAI forecast (forecast.active=true)")
    ap.add_argument("--cv-scheme", default=None,
                    help="cross-validation scheme for --validate: loeo|year|site|random")
    ap.add_argument("--diagnostics", action="store_true",
                    help="write identifiability + structural-adequacy tables after a run")
    ap.add_argument("--no-progress", action="store_true")
    ap.add_argument("--combine", nargs="+", default=None,
                    help="list of result directories to combine (loads design.csv from each)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    if args.n is not None:
        cfg["method"].setdefault("sample", {})["n"] = args.n
    if args.engine:
        cfg["method"].setdefault("sample", {})["engine"] = args.engine
    if args.preset:
        cfg["method"]["preset"] = args.preset
    if args.bayesian_engine:
        cfg["method"].setdefault("bayesian", {})["engine"] = args.bayesian_engine
    if args.optimizer:
        cfg["method"].setdefault("optimizer", {})["engine"] = args.optimizer
        cfg["method"].setdefault("bayesian", {})["engine"] = "none"
    if args.sensitivity:
        cfg["method"].setdefault("sensitivity", {}).update(
            {"engine": args.sensitivity, "active": True, "auto_activate": True})
    if args.select:
        cfg["method"].setdefault("select", {}).update(
            {"engine": f"stepwise_{args.select}", "active": True})
    if args.surrogate:
        cfg["method"].setdefault("surrogate", {}).update(
            {"engine": args.surrogate, "active": True})
    if args.n_particles is not None:
        cfg["method"].setdefault("bayesian", {})["n_particles"] = args.n_particles
    if args.experiments:
        cfg["experiments"] = args.experiments
    if args.assim_mode:
        cfg.setdefault("assimilation", {})["mode"] = args.assim_mode
    if args.forecast:
        cfg.setdefault("forecast", {})["active"] = True

    name = cfg["calibrator"]["name"]
    outdir = Path(args.outdir or f"{cfg['calibrator'].get('results_dir', 'results')}/{name}")
    figdir = Path(cfg["calibrator"].get("figures_dir", "figures")) / name

    if args.nowcast:
        res = orchestrator.nowcast(cfg, args.nowcast, progress=not args.no_progress)
        _write_forecast(res.get("forecast", {}), outdir)
        print(f"\nNowcast as of {res['as_of']} -> {outdir.resolve()}")
        print("=== Best-fit parameters ===")
        for k, v in res["best_theta"].items():
            print(f"  {k:8s} {v:.4f}")
        return

    if args.combined:
        res = orchestrator.combined_mode(cfg, progress=not args.no_progress)
        _write_assimilation(res["assimilation"], outdir)
        print("\n=== Calibrated (base) parameters ===")
        for k, v in res["calibration"].best_theta.items():
            print(f"  {k:8s} {v:.4f}")
        print(f"\nAssimilation outputs -> {outdir.resolve()}")
        return

    if args.assimilate or cfg.get("assimilation", {}).get("active", False):
        res = orchestrator.assimilate(cfg, progress=not args.no_progress)
        _write_assimilation(res, outdir)
        print(f"\nAssimilation ({res.get('mode')}) outputs -> {outdir.resolve()}")
        return

    if args.validate:
        scheme = args.cv_scheme or cfg.get("method", {}).get("validation", {}).get("scheme", "loeo")
        df = orchestrator.validate_cv(cfg, scheme=scheme, progress=not args.no_progress)
        outdir.mkdir(parents=True, exist_ok=True)
        df.round(3).to_csv(outdir / f"validation_{scheme}.csv", index=False)
        print(f"\nCross-validation ({scheme}) -> {outdir/('validation_'+scheme+'.csv')}")
        print(df.round(2).to_string(index=False))
        return

    if args.combine:
        print(f"Combining {len(args.combine)} runs: {args.combine}")
        result = orchestrator.combine_runs(cfg, args.combine)
    else:
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
    if args.diagnostics or cfg.get("diagnostics", {}).get("active", False):
        from dssatcalibrator import diagnostics
        outdir.mkdir(parents=True, exist_ok=True)
        ident = diagnostics.identifiability(result)
        struct = diagnostics.structural_adequacy(result)
        ident.round(4).to_csv(outdir / "identifiability.csv", index=False)
        struct.round(4).to_csv(outdir / "structural_adequacy.csv", index=False)
        print("\n=== Identifiability (posterior vs prior width) ===")
        print(ident.round(3).to_string(index=False))

    if cfg.get("forecast", {}).get("active", False):
        from dssatcalibrator import forecast as fc
        fcs = {v: fc.forecast_lai(cfg, result, variable=v)
               for v in cfg["forecast"].get("variables", ["LAID"])}
        _write_forecast(fcs, outdir)

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

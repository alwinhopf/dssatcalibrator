"""Orchestration: evaluate a sampled design and drive a calibration engine.

Flow (preset C / GLUE):
    config -> ParameterSpace -> sample design
           -> spawn (sample x experiment) in parallel -> parse -> score
           -> GLUE post-process (weights, behavioural set, best theta)

Also provides ``evaluate_theta`` (one theta across all experiments) used by the
NSGA-II engine and by leave-one-environment-out validation.
"""

from __future__ import annotations

import gc
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import objective as obj
from .config import active_parameters, crop_for, resolve_exe
from .observations import Observations
from .runner import resolve_cores, run_many
from .samplers import sample
from .spaces import ParameterSpace
from .spawn import parse_treatments, spawn_and_run


@dataclass
class CalibrationResult:
    cfg: dict
    space: ParameterSpace
    obs: Observations
    experiments: list[str]
    design: pd.DataFrame
    obj_results: dict                       # sample_id -> ObjectiveResult
    best_theta: dict
    best: obj.ObjectiveResult
    glue: object = None
    nsga2: object = None
    extras: dict = field(default_factory=dict)


def _setup(cfg: dict):
    space = ParameterSpace.from_config(cfg)
    crop = crop_for(cfg, (cfg.get("crops") or [{}])[0].get("code", "HM"))
    exe = resolve_exe(cfg)
    specs = active_parameters(cfg)
    hemp_dir = Path(cfg["source"]["hemp_dir"])
    run_root = Path(cfg["calibrator"]["workdir"]) / cfg["calibrator"]["name"]
    run_root.mkdir(parents=True, exist_ok=True)
    src = cfg.get("source", {}).get("observations", "dssat")
    if src == "dssat":
        obs = Observations.from_dssat(hemp_dir, cfg.get("experiments", []), crop_ext=crop["code"])
    else:
        obs = Observations.from_csv(src)
    experiments = [e for e in cfg.get("experiments", []) if e in set(obs.experiments())]
    treatments = {e: parse_treatments(hemp_dir / f"{e}.{crop['filex_ext']}") for e in experiments}
    return space, crop, exe, specs, run_root, obs, experiments, treatments


def _score_theta(theta: dict, experiments, *, cfg, crop, specs, run_root, treatments,
                 exe, obs, n_workers) -> obj.ObjectiveResult:
    jobs, keys = [], []
    for exp in experiments:
        jobs.append(dict(theta=dict(theta), exp_id=exp, cfg=cfg, crop=crop, param_specs=specs,
                         run_root=run_root, treatments=treatments[exp], exe=exe))
        keys.append(exp)
    results = run_many(jobs, n_workers=n_workers)
    rmap = {k: r for k, r in zip(keys, results)}
    return obj.score(rmap, obs.table, cfg)


def evaluate_thetas(cfg: dict, thetas: list[dict], *, setup=None, n_workers=None,
                    progress=False) -> tuple[list, tuple]:
    """Score many parameter sets in ONE parallel batch — the shared engine primitive.

    Every engine that needs to evaluate a list of candidate parameter sets
    (optimizers, sensitivity screening, MCMC, the surrogate, NSGA-II) calls this.
    All ``(theta x experiment)`` DSSAT runs are flattened into a *single* parallel
    batch, so every core is used no matter how few experiments there are — this is
    what fixed the old "NSGA-II evaluates its population serially" bottleneck.

    Parameters
    ----------
    thetas
        list of parameter dicts ``{name: value}`` to evaluate.
    setup
        the tuple returned by :func:`_setup`; pass it in to avoid re-reading the
        experiment files on every call (cheap but not free).

    Returns ``(results, setup)`` where ``results[i]`` is the
    :class:`objective.ObjectiveResult` for ``thetas[i]`` and ``setup`` is reusable.
    """
    if setup is None:
        setup = _setup(cfg)
    space, crop, exe, specs, run_root, obs, experiments, treatments = setup
    if n_workers is None:
        n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))

    jobs, idx = [], []
    for ti, theta in enumerate(thetas):
        for exp in experiments:
            jobs.append(dict(theta=dict(theta), exp_id=exp, cfg=cfg, crop=crop,
                             param_specs=specs, run_root=run_root,
                             treatments=treatments[exp], exe=exe))
            idx.append((ti, exp))

    total = len(jobs)
    done = {"n": 0}

    def _cb(_r):
        done["n"] += 1
        if progress and (done["n"] % max(1, total // 20) == 0 or done["n"] == total):
            print(f"  spawns {done['n']}/{total}", flush=True)

    results = run_many(jobs, n_workers=n_workers, on_done=_cb if progress else None)

    per_theta: list[dict] = [{} for _ in thetas]
    for (ti, exp), res in zip(idx, results):
        per_theta[ti][exp] = res
    return [obj.score(rmap, obs.table, cfg) for rmap in per_theta], setup


def evaluate_design(cfg: dict, samples: pd.DataFrame, *, progress=True):
    """Run every (sample x experiment) spawn in parallel, score per sample, in memory-efficient batches."""
    space, crop, exe, specs, run_root, obs, experiments, treatments = _setup(cfg)
    n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))

    batch_size = cfg.get("calibrator", {}).get("batch_size", 50)
    sample_ids = list(samples.index)
    
    rows = []
    obj_results = {}
    best_score = float("inf")
    best_obj_res = None
    best_sid = None

    total_jobs = len(samples) * len(experiments)
    done = {"n": 0}

    def _cb(_res):
        done["n"] += 1
        if progress and (done["n"] % max(1, total_jobs // 20) == 0 or done["n"] == total_jobs):
            print(f"  spawns {done['n']}/{total_jobs}", flush=True)

    for b_start in range(0, len(sample_ids), batch_size):
        batch_ids = sample_ids[b_start : b_start + batch_size]
        jobs, idx = [], []
        
        for sid in batch_ids:
            row = samples.loc[sid]
            theta = space.to_theta(row.to_numpy())
            for exp in experiments:
                jobs.append(dict(theta=dict(theta), exp_id=exp, cfg=cfg, crop=crop, param_specs=specs,
                                 run_root=run_root, treatments=treatments[exp], exe=exe))
                idx.append((sid, exp))

        # Run this batch
        results = run_many(jobs, n_workers=n_workers, on_done=_cb if progress else None)

        # Score this batch
        per_sample: dict = {}
        for (sid, exp), res in zip(idx, results):
            per_sample.setdefault(sid, {})[exp] = res

        for sid, rmap in per_sample.items():
            o = obj.score(rmap, obs.table, cfg)
            
            # Check if this is the best so far
            if o.score < best_score:
                best_score = o.score
                best_obj_res = o
                best_sid = sid

            rec = {"sample_id": sid, **space.to_theta(samples.loc[sid].to_numpy()),
                   "score": o.score, "loglik": o.loglik, "n_obs": len(o.residuals)}
            rows.append(rec)
            
        # Free memory of this batch
        del jobs, idx, results, per_sample
        gc.collect()

    # To satisfy lookups in calibrate() and validate_loeo(), make sure the best
    # objective result is in obj_results
    if best_sid is not None:
        obj_results[best_sid] = best_obj_res

    design = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    return design, obj_results, (space, obs, experiments)


# ---------------------------------------------------------------------------
# Preset resolution — make method.preset actually do something (CONCEPT.md §14a)
# ---------------------------------------------------------------------------
# A preset pre-fills the recommended engine for each stage; any explicit stage
# block in the user's config overrides it. The OPTIONAL stages (sensitivity,
# select, surrogate) only run when switched on (`active: true` / `engine` set), so
# a bare preset just picks the main estimator and you opt into screening/selection.
PRESETS: dict[str, dict] = {
    "A": {"sample": {"engine": "lhs"}, "bayesian": {"engine": "smc_pf"},
          "sensitivity": {"engine": "morris"}},          # screen -> map -> particle filter
    "B": {"optimizer": {"engine": "diffevo"}, "bayesian": {"engine": "none"},
          "sensitivity": {"engine": "morris"}},          # screen -> optimise (best-fit point)
    "C": {"sample": {"engine": "lhs"}, "bayesian": {"engine": "glue"}},   # map -> GLUE (default)
    "D": {"sample": {"engine": "sobol"}, "bayesian": {"engine": "mcmc"},
          "sensitivity": {"engine": "morris"}},          # screen -> map -> full MCMC posterior
}


def _resolve_method(cfg: dict) -> dict:
    """Apply ``method.preset`` defaults without clobbering explicit stage blocks."""
    from copy import deepcopy
    method = deepcopy(cfg.get("method", {}) or {})
    preset = str(method.get("preset", "C")).upper()
    if preset in PRESETS:
        for stage, defaults in PRESETS[preset].items():
            method[stage] = {**defaults, **(method.get(stage) or {})}   # user wins
    method.setdefault("sample", {"engine": "lhs", "n": 200})
    method.setdefault("bayesian", {"engine": "glue"})
    return method


def _stage_on(block: dict | None) -> bool:
    """Is an OPTIONAL stage (sensitivity / select / surrogate) switched on?

    Optional stages are off unless you explicitly set ``active: true``. A preset may
    pre-fill the recommended *engine* for the stage, but it never turns the stage on
    by itself — so a bare preset just picks the main estimator and you opt into
    screening/selection/surrogate deliberately (and predictably)."""
    return bool(block and block.get("active", False))


def _apply_active_subset(cfg: dict, keep: list[str]) -> dict:
    """Return a copy of ``cfg`` with only ``keep`` parameters left ``active: true``.

    Used by the sensitivity and selection stages to hand the chosen parameter
    subset to the downstream estimator — pure config, no code paths change.
    """
    from copy import deepcopy
    keep = set(keep)
    out = deepcopy(cfg)
    for _group, params in (out.get("parameters") or {}).items():
        if not isinstance(params, dict):
            continue
        for name, spec in params.items():
            if isinstance(spec, dict) and spec.get("active", False):
                spec["active"] = name in keep
    return out


def _results_scorer(cfg: dict, setup, n_workers):
    """Build a ``score_results(list_of_theta) -> list[ObjectiveResult]`` closure
    bound to a fixed setup, for the screening / selection / optimiser engines."""
    def score_results(thetas):
        res, _ = evaluate_thetas(cfg, list(thetas), setup=setup, n_workers=n_workers)
        return res
    return score_results


def calibrate(cfg: dict, *, progress=True) -> CalibrationResult:
    """Run the configured calibration pipeline.

    Steps (each optional one only runs when switched on):
      0. resolve ``method.preset`` -> per-stage engines
      1. [sensitivity] screen parameters; optionally keep only influential ones
      2. [select]      AgMIP stepwise BIC/AICc -> keep the parameters that earn their place
      3. estimate with the chosen engine: glue | smc_pf | mcmc | optimizer | surrogate
      4. [multiobjective] NSGA-II Pareto front add-on
    """
    from .engines import (run_glue, run_mcmc, run_nsga2, run_optimizer,
                          run_sensitivity, run_smc_pf, run_surrogate, stepwise_select)

    method = _resolve_method(cfg)
    seed = int(cfg["calibrator"].get("seed", 42))
    n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))
    extras: dict = {}

    # Work on a (possibly param-pruned) config; estimators read active params from it.
    work_cfg = cfg
    setup = _setup(work_cfg)
    space = setup[0]

    # --- 1. sensitivity screening (optional) --------------------------------------
    sens_block = method.get("sensitivity")
    if _stage_on(sens_block):
        if progress:
            print(f"[1/4] Sensitivity screening ({sens_block.get('engine','morris')})...", flush=True)
        scorer = _results_scorer(work_cfg, setup, n_workers)
        sens = run_sensitivity(space, scorer, method=sens_block.get("engine", "morris"),
                               trajectories=sens_block.get("trajectories", 10),
                               n_base=sens_block.get("n_base", 256), seed=seed, progress=progress)
        extras["sensitivity"] = sens.ranking
        if sens_block.get("auto_activate", False):
            keep = sens.influential(keep=sens_block.get("keep"),
                                    rel_threshold=float(sens_block.get("rel_threshold", 0.1)))
            if progress:
                print(f"      auto-activate keeps {keep}", flush=True)
            work_cfg = _apply_active_subset(work_cfg, keep)
            setup = _setup(work_cfg)
            space = setup[0]

    # --- 2. AgMIP stepwise selection (optional) -----------------------------------
    sel_block = method.get("select")
    if _stage_on(sel_block):
        crit = "aicc" if "aicc" in str(sel_block.get("engine", "")).lower() else "bic"
        if progress:
            print(f"[2/4] Stepwise parameter selection ({crit.upper()})...", flush=True)
        scorer = _results_scorer(work_cfg, setup, n_workers)
        sel = stepwise_select(space, scorer, criterion=crit,
                              optimizer=sel_block.get("optimizer", "nelder_mead"),
                              optimizer_restarts=int(sel_block.get("restarts", 2)),
                              maxiter=sel_block.get("maxiter"), seed=seed, progress=progress)
        extras["selection"] = sel
        work_cfg = _apply_active_subset(work_cfg, sel.selected)
        setup = _setup(work_cfg)
        space = setup[0]

    # --- optional AgMIP iterative-reweighted (WLS) weighting ----------------------
    # Run one quick pass, set each variable's weight to 1/residual-variance, refit.
    if str(cfg.get("objective", {}).get("weighting", "")).lower() == "agmip_wls":
        work_cfg = _agmip_reweight(work_cfg, setup, n_workers, seed, progress)

    # --- 3. main estimator --------------------------------------------------------
    bayes = str(method.get("bayesian", {}).get("engine", "glue")).lower()
    opt = str(method.get("optimizer", {}).get("engine", "none")).lower()
    surrogate_on = _stage_on(method.get("surrogate"))
    _, _crop, _exe, _specs, _run_root, obs, experiments, _trts = setup

    if progress:
        print(f"[3/4] Estimating with "
              f"{'surrogate+' if surrogate_on else ''}"
              f"{bayes if bayes not in ('none','') else opt}...", flush=True)

    if surrogate_on:
        scorer = _results_scorer(work_cfg, setup, n_workers)
        sur = run_surrogate(work_cfg, space, scorer, progress=progress)
        glue = run_glue(sur.design, space.names, work_cfg, space=space)
        best = sur.obj_results[glue.best_sample_id]
        extras["surrogate_info"] = sur.info
        result = CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                                   design=glue.design, obj_results=sur.obj_results,
                                   best_theta=glue.best_theta, best=best, glue=glue, extras=extras)

    elif bayes == "smc_pf":
        smc = run_smc_pf(work_cfg, progress=progress)
        extras.update(initial_design=smc.initial_design, engine="smc_pf")
        result = CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                                   design=smc.design, obj_results=smc.obj_results,
                                   best_theta=smc.best_theta, best=smc.best, glue=smc, extras=extras)

    elif bayes == "mcmc":
        scorer = _results_scorer(work_cfg, setup, n_workers)
        mc = run_mcmc(work_cfg, scorer, space, progress=progress)
        extras.update(initial_design=mc.initial_design, engine="mcmc",
                      mcmc_chain=mc.chain, acceptance=mc.acceptance)
        result = CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                                   design=mc.design, obj_results=mc.obj_results,
                                   best_theta=mc.best_theta, best=mc.best, glue=mc, extras=extras)

    elif bayes in ("none", "") and opt in ("nelder_mead", "diffevo", "neldermead", "nm", "de"):
        ocfg = method.get("optimizer", {})
        scorer = _results_scorer(work_cfg, setup, n_workers)
        score_batch = lambda ths: [r.score for r in scorer(ths)]   # noqa: E731
        ores = run_optimizer(space, score_batch, method=opt, seed=seed,
                             maxiter=ocfg.get("maxiter"), popsize=int(ocfg.get("popsize", 15)),
                             restarts=int(ocfg.get("restarts", 1)), progress=progress)
        best = scorer([ores.best_theta])[0]
        # Represent the optimiser's trace as a one-row "design" so reporting works.
        design = pd.DataFrame([{"sample_id": 0, **ores.best_theta, "score": best.score,
                                "loglik": best.loglik, "n_obs": len(best.residuals), "weight": 1.0}])
        extras.update(engine="optimizer", optimizer_history=ores.history, n_eval=ores.n_eval)
        result = CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                                   design=design, obj_results={0: best},
                                   best_theta=ores.best_theta, best=best, glue=None, extras=extras)

    else:  # glue (default) — sample a design, evaluate, GLUE post-process
        n = int(method.get("sample", {}).get("n", 200))
        engine = method.get("sample", {}).get("engine", "lhs")
        samples = sample(space, n=n, engine=engine, seed=seed, include_start=True)
        design, obj_results, (space, obs, experiments) = evaluate_design(work_cfg, samples, progress=progress)
        glue = run_glue(design, space.names, work_cfg, space=space)
        best = obj_results[glue.best_sample_id]
        extras.update(engine="glue")
        result = CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                                   design=glue.design, obj_results=obj_results,
                                   best_theta=glue.best_theta, best=best, glue=glue, extras=extras)

    # --- 4. multi-objective Pareto front add-on (optional) ------------------------
    if method.get("multiobjective", {}).get("engine") == "nsga2":
        mo = method["multiobjective"]
        obj_vars = mo.get("variables") or sorted({uv for uv in result.best.per_var})
        if progress:
            print(f"[4/4] NSGA-II Pareto front over {obj_vars}...", flush=True)

        def eval_batch(thetas):
            res, _ = evaluate_thetas(work_cfg, thetas, setup=setup, n_workers=n_workers)
            return [{v: r.per_var.get(v, {}).get("nRMSE_pct", 1e6) for v in obj_vars} for r in res]

        result.nsga2 = run_nsga2(eval_batch, space, obj_vars,
                                 pop_size=int(mo.get("pop_size", 16)),
                                 n_gen=int(mo.get("n_gen", 5)), seed=seed)
    return result


def _agmip_reweight(cfg: dict, setup, n_workers, seed, progress) -> dict:
    """One-step AgMIP iterative-reweighted least squares (``weighting: agmip_wls``).

    Quick OLS-ish pass (a small LHS design), then set each variable group's weight
    to ``1 / variance(residuals)`` so noisy variables count less — the AgMIP default.
    Returns a config copy with the derived ``objective.weights`` baked in.
    """
    from copy import deepcopy
    space = setup[0]
    samples = sample(space, n=int(cfg.get("calibrator", {}).get("wls_probe_n", 40)),
                     engine="lhs", seed=seed, include_start=True)
    thetas = [space.to_theta(samples.iloc[i].to_numpy()) for i in range(len(samples))]
    results, _ = evaluate_thetas(cfg, thetas, setup=setup, n_workers=n_workers)
    best = min(results, key=lambda r: r.score if np.isfinite(r.score) else float("inf"))
    out = deepcopy(cfg)
    weights = dict(out.get("objective", {}).get("weights", {}) or {})
    if not best.residuals.empty:
        for uv, g in best.residuals.groupby("user_var"):
            var = float(np.var(g["resid"])) if len(g) > 1 else float(g["resid"].iloc[0] ** 2)
            weights[uv] = 1.0 / var if var > 0 else 1.0
    out.setdefault("objective", {})["weights"] = weights
    out["objective"]["weighting"] = "unified"   # apply the derived weights on the unified surface
    if progress:
        print(f"      agmip_wls weights -> { {k: round(v,3) for k,v in weights.items()} }", flush=True)
    return out


def spawn_results_for(cfg: dict, theta: dict, experiments=None) -> dict:
    """Return {exp_id: SpawnResult} for a theta (cached spawns are instant).

    Used to recover full PlantGro curves for the best fit when plotting.
    """
    space, crop, exe, specs, run_root, obs, exps, treatments = _setup(cfg)
    experiments = experiments or exps
    out = {}
    for e in experiments:
        out[e] = spawn_and_run(dict(theta), exp_id=e, cfg=cfg, crop=crop, param_specs=specs,
                               run_root=run_root, treatments=treatments[e], exe=exe)
    return out


def combine_runs(cfg: dict, run_dirs: list[str | Path]) -> CalibrationResult:
    """Combine results from multiple completed calibration directories."""
    from .engines import run_glue

    space, crop, exe, specs, run_root, obs, experiments, treatments = _setup(cfg)
    n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))

    dfs = []
    for rdir in run_dirs:
        p = Path(rdir) / "design.csv"
        if not p.exists():
            raise FileNotFoundError(f"No design.csv found in {rdir}")
        dfs.append(pd.read_csv(p))

    combined = pd.concat(dfs, ignore_index=True)

    # Drop duplicates by parameter values
    param_cols = [n for n in space.names if n in combined.columns]
    if not param_cols:
         raise ValueError("No active parameters found in design.csv files matching config")
    
    # We round parameters slightly to avoid float representation differences when checking duplicates
    combined["_dup_key"] = combined[param_cols].apply(lambda r: tuple(np.round(r.to_numpy(), 6)), axis=1)
    combined = combined.drop_duplicates(subset=["_dup_key"])
    combined = combined.drop(columns=["_dup_key"])

    combined = combined.reset_index(drop=True)
    combined["sample_id"] = combined.index

    # Run GLUE on the aggregated dataset
    glue = run_glue(combined, space.names, cfg, space=space)

    # Evaluate the overall best theta across all experiments
    print(f"Evaluating combined best theta across all {len(experiments)} experiments...")
    best = _score_theta(glue.best_theta, experiments, cfg=cfg, crop=crop, specs=specs,
                       run_root=run_root, treatments=treatments, exe=exe, obs=obs,
                       n_workers=n_workers)

    obj_results = {glue.best_sample_id: best}

    return CalibrationResult(
        cfg=cfg,
        space=space,
        obs=obs,
        experiments=experiments,
        design=glue.design,
        obj_results=obj_results,
        best_theta=glue.best_theta,
        best=best,
        glue=glue
    )


def validate_loeo(cfg: dict, *, progress=False) -> pd.DataFrame:
    """Leave-one-environment-out: calibrate on n-1 experiments, evaluate on the held-out one.

    Returns a tidy table of calibration vs evaluation fit per held-out experiment.
    """
    space, crop, exe, specs, run_root, obs, experiments, treatments = _setup(cfg)
    n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))
    method = cfg.get("method", {})
    n = int(method.get("sample", {}).get("n", 100))
    seed = int(cfg["calibrator"].get("seed", 42))

    rows = []
    for held in experiments:
        train = [e for e in experiments if e != held]
        if not train:
            continue
        cfg_train = {**cfg, "experiments": train}
        samples = sample(space, n=n, engine=method.get("sample", {}).get("engine", "lhs"),
                         seed=seed, include_start=True)
        design, obj_results, _ = evaluate_design(cfg_train, samples, progress=progress)
        from .engines import run_glue
        glue = run_glue(design, space.names, cfg_train, space=space)
        best_theta = glue.best_theta
        cal = obj_results[glue.best_sample_id]
        # evaluate best on the held-out experiment
        ev = _score_theta(best_theta, [held], cfg={**cfg, "experiments": [held]}, crop=crop,
                          specs=specs, run_root=run_root, treatments=treatments, exe=exe,
                          obs=obs, n_workers=n_workers)
        for uv, m in cal.per_var.items():
            rows.append({"held_out": held, "split": "calibration", "variable": uv, **m})
        for uv, m in ev.per_var.items():
            rows.append({"held_out": held, "split": "evaluation", "variable": uv, **m})
    return pd.DataFrame(rows)

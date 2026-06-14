"""Orchestration: evaluate a sampled design and drive a calibration engine.

Flow (preset C / GLUE):
    config -> ParameterSpace -> sample design
           -> spawn (sample x experiment) in parallel -> parse -> score
           -> GLUE post-process (weights, behavioural set, best theta)

Also provides ``evaluate_theta`` (one theta across all experiments) used by the
NSGA-II engine and by leave-one-environment-out validation.
"""

from __future__ import annotations

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


def evaluate_design(cfg: dict, samples: pd.DataFrame, *, progress=True):
    """Run every (sample x experiment) spawn in parallel, score per sample."""
    space, crop, exe, specs, run_root, obs, experiments, treatments = _setup(cfg)
    n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))

    jobs, idx = [], []
    for sid, row in samples.iterrows():
        theta = space.to_theta(row.to_numpy())
        for exp in experiments:
            jobs.append(dict(theta=dict(theta), exp_id=exp, cfg=cfg, crop=crop, param_specs=specs,
                             run_root=run_root, treatments=treatments[exp], exe=exe))
            idx.append((sid, exp))

    done = {"n": 0}
    total = len(jobs)

    def _cb(_res):
        done["n"] += 1
        if progress and (done["n"] % max(1, total // 20) == 0 or done["n"] == total):
            print(f"  spawns {done['n']}/{total}", flush=True)

    results = run_many(jobs, n_workers=n_workers, on_done=_cb if progress else None)

    per_sample: dict = {}
    for (sid, exp), res in zip(idx, results):
        per_sample.setdefault(sid, {})[exp] = res

    rows, obj_results = [], {}
    for sid, rmap in per_sample.items():
        o = obj.score(rmap, obs.table, cfg)
        obj_results[sid] = o
        rec = {"sample_id": sid, **space.to_theta(samples.loc[sid].to_numpy()),
               "score": o.score, "loglik": o.loglik, "n_obs": len(o.residuals)}
        rows.append(rec)
    design = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    return design, obj_results, (space, obs, experiments)


def calibrate(cfg: dict, *, progress=True) -> CalibrationResult:
    """Run the configured calibration (GLUE by default; NSGA-II if requested)."""
    from .engines import run_glue, run_nsga2, run_smc_pf

    space = ParameterSpace.from_config(cfg)
    method = cfg.get("method", {})
    seed = int(cfg["calibrator"].get("seed", 42))
    bayesian_engine = method.get("bayesian", {}).get("engine", "glue")

    if bayesian_engine == "smc_pf":
        smc = run_smc_pf(cfg, progress=progress)
        best = smc.best
        _, _, _, _, _, obs, experiments, _ = _setup(cfg)
        result = CalibrationResult(cfg=cfg, space=space, obs=obs, experiments=experiments,
                                   design=smc.design, obj_results=smc.obj_results,
                                   best_theta=smc.best_theta, best=best, glue=smc,
                                   extras={"initial_design": smc.initial_design,
                                           "engine": "smc_pf"})
    else:
        n = int(method.get("sample", {}).get("n", 200))
        engine = method.get("sample", {}).get("engine", "lhs")
        samples = sample(space, n=n, engine=engine, seed=seed, include_start=True)

        design, obj_results, (space, obs, experiments) = evaluate_design(cfg, samples, progress=progress)

        glue = run_glue(design, space.names, cfg)
        best = obj_results[glue.best_sample_id]

        result = CalibrationResult(cfg=cfg, space=space, obs=obs, experiments=experiments,
                                   design=glue.design, obj_results=obj_results,
                                   best_theta=glue.best_theta, best=best, glue=glue)

    # optional multi-objective Pareto front (per-variable trade-offs)
    if method.get("multiobjective", {}).get("engine") == "nsga2":
        mo = method["multiobjective"]
        obj_vars = mo.get("variables") or sorted({uv for uv in best.per_var})
        space2, crop, exe, specs, run_root, obs2, exps, trts = _setup(cfg)
        n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))

        def eval_theta(theta):
            r = _score_theta(theta, exps, cfg=cfg, crop=crop, specs=specs, run_root=run_root,
                             treatments=trts, exe=exe, obs=obs2, n_workers=n_workers)
            return {v: r.per_var.get(v, {}).get("nRMSE_pct", 1e6) for v in obj_vars}

        result.nsga2 = run_nsga2(eval_theta, space2, obj_vars,
                                 pop_size=int(mo.get("pop_size", 16)),
                                 n_gen=int(mo.get("n_gen", 5)), seed=seed)
    return result


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
        glue = run_glue(design, space.names, cfg_train)
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

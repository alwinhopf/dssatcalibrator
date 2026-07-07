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
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _warn_unmatched(obs, cfg: dict) -> None:
    """Warn (once) if observations carry variables the objective can't score."""
    try:
        unmatched = obj.unmatched_variables(getattr(obs, "table", obs), cfg)
    except Exception:
        return
    if unmatched:
        logger.warning(
            "Observations include variable(s) %s that are not mapped in "
            "engine.timeseries_outputs / engine.scalar_outputs (and may not be a "
            "DSSAT output). These rows will be IGNORED when scoring. Map them to a "
            "DSSAT output column to use them.", unmatched,
        )

from . import objective as obj
from .config import crop_for, fixed_parameters, resolve_exe
from .observations import Observations
from .runner import resolve_cores, run_many
from .samplers import sample
from .spaces import ParameterSpace, expand_parameter_specs
from .spawn import SpawnResult, parse_treatments, spawn_and_run, theta_hash


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
    fixed_specs = expand_parameter_specs(cfg, fixed_parameters(cfg))
    active_names = {s["name"] for s in space.specs}
    specs = space.specs + [s for s in fixed_specs if s["name"] not in active_names]
    hemp_dir = Path(cfg["source"]["hemp_dir"])
    run_root = Path(cfg["calibrator"]["workdir"]) / cfg["calibrator"]["name"]
    run_root.mkdir(parents=True, exist_ok=True)
    
    src_block = cfg.get("source", {})
    if "table" in src_block:
        obs = Observations(src_block["table"])
    elif cfg.get("observation_sources"):
        obs = Observations.from_sources(cfg, cfg.get("experiments", []))
    else:
        src = src_block.get("observations", "dssat")
        if src == "dssat":
            obs = Observations.from_dssat(hemp_dir, cfg.get("experiments", []), crop_ext=crop["code"])
        else:
            obs = Observations.from_csv(src)
            
    experiments = [e for e in cfg.get("experiments", []) if e in set(obs.experiments())]
    treatments = {e: parse_treatments(hemp_dir / f"{e}.{crop['filex_ext']}") for e in experiments}
    selected_by_exp = cfg.get("calibration_treatments_by_experiment")
    selected_treatments = cfg.get("calibration_treatments")
    if selected_by_exp:
        filtered = {}
        for e, available in treatments.items():
            requested = selected_by_exp.get(e, selected_by_exp.get(str(e)))
            if requested is None:
                filtered[e] = list(available)
                continue
            if isinstance(requested, (str, int)):
                requested = [requested]
            selected_set = {int(t) for t in requested}
            filtered[e] = [t for t in available if int(t) in selected_set]
        treatments = filtered
    elif selected_treatments is not None:
        selected_set = {int(t) for t in selected_treatments}
        treatments = {
            e: [t for t in available if int(t) in selected_set]
            for e, available in treatments.items()
        }

    # Optional: take planting dates from ingested farm-management rows and set them
    # as a FileX PDATE override (an input, not a calibrated parameter). Opt-in via
    # management.use_source_planting_date; resolved once and stashed on the cfg.
    if (cfg.get("management_options", {}).get("use_source_planting_date")
            and "_planting_dates" not in cfg):
        cfg["_planting_dates"] = obs.planting_dates()

    return space, crop, exe, specs, run_root, obs, experiments, treatments


def _observed_treatments(obs: Observations, exp_id: str, available: list[int]) -> list[int]:
    table = getattr(obs, "table", pd.DataFrame())
    if table is None or table.empty or "treatment" not in table.columns:
        return list(available)
    rows = table[table["exp_id"].astype(str) == str(exp_id)]
    vals = sorted({int(v) for v in pd.to_numeric(rows["treatment"], errors="coerce").dropna()})
    if not vals:
        return list(available)
    available_set = set(map(int, available))
    selected = [v for v in vals if v in available_set]
    return selected or list(available)


def _treatment_units(experiments: list[str], treatments: dict, obs: Observations) -> list[tuple[str, int]]:
    units = []
    for exp in experiments:
        for trt in _observed_treatments(obs, exp, treatments[exp]):
            units.append((exp, int(trt)))
    return units


def _obs_table_for_units(obs: Observations, units: list[tuple[str, int]]) -> pd.DataFrame:
    table = getattr(obs, "table", pd.DataFrame())
    if table is None or table.empty or "exp_id" not in table.columns or "treatment" not in table.columns:
        return table
    allowed = {(str(exp), int(trt)) for exp, trt in units}
    trt = pd.to_numeric(table["treatment"], errors="coerce")
    mask = [
        (str(exp), int(t)) in allowed if pd.notna(t) else False
        for exp, t in zip(table["exp_id"], trt)
    ]
    return table.loc[mask].copy()


def _concat_frames(frames: list[pd.DataFrame]) -> pd.DataFrame:
    frames = [f for f in frames if f is not None and not f.empty]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def _merge_treatment_results(exp_id: str, results: list[SpawnResult]) -> SpawnResult:
    ok = [r for r in results if getattr(r, "status", "") in {"success", "cached"}]
    if not ok:
        first = results[0] if results else None
        return SpawnResult(
            status="error",
            run_dir=getattr(first, "run_dir", Path("")),
            theta=getattr(first, "theta", {}),
            message="; ".join(str(getattr(r, "message", "")) for r in results if getattr(r, "message", "")),
        )
    output_keys = sorted({k for r in ok for k in (getattr(r, "outputs", {}) or {}).keys()})
    outputs = {
        key: _concat_frames([(getattr(r, "outputs", {}) or {}).get(key) for r in ok])
        for key in output_keys
    }
    bad = [r for r in results if getattr(r, "status", "") not in {"success", "cached"}]
    message = ""
    if bad:
        message = "Partial treatment failures: " + "; ".join(
            str(getattr(r, "message", "")) for r in bad if getattr(r, "message", "")
        )
    status = "success" if any(getattr(r, "status", "") == "success" for r in ok) else "cached"
    return SpawnResult(
        status=status,
        run_dir=getattr(ok[0], "run_dir", Path("")).parent,
        theta=getattr(ok[0], "theta", {}),
        plantgro=_concat_frames([getattr(r, "plantgro", pd.DataFrame()) for r in ok]),
        evaluate=_concat_frames([getattr(r, "evaluate", pd.DataFrame()) for r in ok]),
        outputs=outputs,
        message=message,
    )


def _merge_result_map(raw: dict[str, list[SpawnResult]]) -> dict[str, SpawnResult]:
    return {exp: _merge_treatment_results(exp, results) for exp, results in raw.items()}


def _score_theta(theta: dict, experiments, *, cfg, crop, specs, run_root, treatments,
                 exe, obs, n_workers) -> obj.ObjectiveResult:
    jobs, keys = [], []
    spawn_timeout = cfg.get("calibrator", {}).get("spawn_timeout")
    for exp, trt in _treatment_units(list(experiments), treatments, obs):
        jobs.append(dict(theta=dict(theta), exp_id=exp, cfg=cfg, crop=crop, param_specs=specs,
                         run_root=run_root, treatments=[trt], exe=exe,
                         **({"timeout": int(spawn_timeout)} if spawn_timeout else {})))
        keys.append((exp, trt))
    results = run_many(jobs, n_workers=n_workers)
    raw: dict[str, list[SpawnResult]] = {}
    for (exp, _trt), res in zip(keys, results):
        raw.setdefault(exp, []).append(res)
    rmap = _merge_result_map(raw)
    return obj.score(rmap, _obs_table_for_units(obs, keys), cfg)


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
    spawn_timeout = cfg.get("calibrator", {}).get("spawn_timeout")

    jobs, idx = [], []
    units = _treatment_units(experiments, treatments, obs)
    for ti, theta in enumerate(thetas):
        for exp, trt in units:
            jobs.append(dict(theta=dict(theta), exp_id=exp, cfg=cfg, crop=crop,
                             param_specs=specs, run_root=run_root,
                             treatments=[trt], exe=exe,
                             **({"timeout": int(spawn_timeout)} if spawn_timeout else {})))
            idx.append((ti, exp, trt))

    total = len(jobs)
    done = {"n": 0}

    def _cb(_r):
        done["n"] += 1
        if progress and (done["n"] % max(1, total // 20) == 0 or done["n"] == total):
            print(f"  spawns {done['n']}/{total}", flush=True)

    results = run_many(jobs, n_workers=n_workers, on_done=_cb if progress else None)

    per_theta: list[dict] = [{} for _ in thetas]
    for (ti, exp, _trt), res in zip(idx, results):
        per_theta[ti].setdefault(exp, []).append(res)
    obs_table = _obs_table_for_units(obs, units)
    return [obj.score(_merge_result_map(rmap), obs_table, cfg) for rmap in per_theta], setup


def evaluate_design(cfg: dict, samples: pd.DataFrame, *, progress=True):
    """Run every (sample x experiment) spawn in parallel, score per sample, in memory-efficient batches."""
    space, crop, exe, specs, run_root, obs, experiments, treatments = _setup(cfg)
    n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))

    batch_size = cfg.get("calibrator", {}).get("batch_size", 50)
    spawn_timeout = cfg.get("calibrator", {}).get("spawn_timeout")
    sample_ids = list(samples.index)
    
    rows = []
    spawn_manifest_rows = []
    obj_results = {}
    best_score = float("inf")
    best_obj_res = None
    best_sid = None

    units = _treatment_units(experiments, treatments, obs)
    total_jobs = len(samples) * len(units)
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
            for exp, trt in units:
                jobs.append(dict(theta=dict(theta), exp_id=exp, cfg=cfg, crop=crop, param_specs=specs,
                                 run_root=run_root, treatments=[trt], exe=exe,
                                 **({"timeout": int(spawn_timeout)} if spawn_timeout else {})))
                idx.append((sid, exp, trt))

        # Run this batch
        results = run_many(jobs, n_workers=n_workers, on_done=_cb if progress else None)

        # Score this batch
        per_sample: dict = {}
        for (sid, exp, trt), res, job in zip(idx, results, jobs):
            per_sample.setdefault(sid, {}).setdefault(exp, []).append(res)
            theta = dict(getattr(res, "theta", {}) or job["theta"])
            theta_jsonable = {k: (v.item() if hasattr(v, "item") else v) for k, v in theta.items()}
            spawn_manifest_rows.append({
                "sample_id": sid,
                "exp_id": exp,
                "treatment": trt,
                "theta_hash": theta_hash(theta) if theta else "",
                "status": getattr(res, "status", ""),
                "message": getattr(res, "message", ""),
                "run_dir": str(getattr(res, "run_dir", "")),
                "theta_json": json.dumps(theta_jsonable, sort_keys=True),
                **{f"theta_{k}": v for k, v in theta.items()},
            })

        obs_table = _obs_table_for_units(obs, units)
        for sid, rawmap in per_sample.items():
            o = obj.score(_merge_result_map(rawmap), obs_table, cfg)
            
            # Check if this is the best so far
            if o.score < best_score:
                best_score = o.score
                best_obj_res = o
                best_sid = sid

            rec = {"sample_id": sid, **space.to_theta(samples.loc[sid].to_numpy()),
                   "score": o.score, "loglik": o.loglik, "n_obs": len(o.residuals)}
            if getattr(o, "per_exp_var", None) is not None and not o.per_exp_var.empty:
                for _, metric_row in o.per_exp_var.iterrows():
                    user_var = str(metric_row["user_var"])
                    exp_id = str(metric_row["exp_id"])
                    prefix = f"{user_var}__{exp_id}"
                    for metric_name in ("RMSE", "nRMSE_pct", "MBE", "d", "EF", "R2"):
                        if metric_name in metric_row:
                            rec[f"{prefix}__{metric_name}"] = metric_row[metric_name]
                    if "MBE" in metric_row:
                        rec[f"{prefix}__abs_MBE"] = abs(float(metric_row["MBE"]))
                for user_var, metric_group in o.per_exp_var.groupby("user_var"):
                    abs_mbe_cols = [
                        f"{user_var}__{exp_id}__abs_MBE"
                        for exp_id in metric_group["exp_id"].astype(str)
                        if f"{user_var}__{exp_id}__abs_MBE" in rec
                    ]
                    if abs_mbe_cols:
                        rec[f"{user_var}__mean_abs_MBE"] = float(np.mean([rec[c] for c in abs_mbe_cols]))
                        rec[f"{user_var}__max_abs_MBE"] = float(np.max([rec[c] for c in abs_mbe_cols]))
            rows.append(rec)
            
        # Free memory of this batch
        del jobs, idx, results, per_sample
        gc.collect()

    # To satisfy lookups in calibrate() and validate_loeo(), make sure the best
    # objective result is in obj_results
    if best_sid is not None:
        obj_results[best_sid] = best_obj_res

    design = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    design.attrs["spawn_manifest"] = pd.DataFrame(spawn_manifest_rows)
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


def _apply_staging(cfg: dict) -> dict:
    """Freeze whole parameter groups / named parameters (``method.staging``).

    Supports the AgMIP-style staged workflow and sparse-data discipline: e.g.
    by mid-season freeze the seed/yield group so only phenology+canopy params are
    estimated. ``freeze_groups`` deactivates entire groups; ``freeze_params`` named
    ones. Pure config — returns a copy with those parameters set ``active: false``.
    """
    st = (cfg.get("method", {}) or {}).get("staging", {}) or {}
    fg = set(st.get("freeze_groups", []) or [])
    fp = set(st.get("freeze_params", []) or [])
    if not fg and not fp:
        return cfg
    from copy import deepcopy
    out = deepcopy(cfg)
    frozen = []
    for group, params in (out.get("parameters") or {}).items():
        if not isinstance(params, dict):
            continue
        for name, spec in params.items():
            if isinstance(spec, dict) and spec.get("active", False) and (group in fg or name in fp):
                spec["active"] = False
                frozen.append(name)
    if frozen:
        logger.info("Staging: froze %d parameter(s): %s", len(frozen), frozen)
    return out


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


# ---------------------------------------------------------------------------
# Estimator registry — the main "stage 3" engines (CONCEPT.md §14a)
# ---------------------------------------------------------------------------
# Each estimator has one signature and returns a CalibrationResult, so adding a
# new main estimator is a function + one ``@_register_estimator`` line (no edit
# to ``calibrate``). ``_resolve_estimator`` maps the configured method block to
# a registry key; the optimiser aliases all resolve to the one "optimizer" entry.

ESTIMATOR_REGISTRY: dict[str, callable] = {}
_OPTIMIZER_ALIASES = ("nelder_mead", "diffevo", "neldermead", "nm", "de",
                      "cmaes", "cma_es", "cma")


def _register_estimator(name):
    def deco(fn):
        ESTIMATOR_REGISTRY[name] = fn
        return fn
    return deco


def _resolve_estimator(method: dict) -> str:
    """Pick the estimator registry key from the resolved ``method`` block.

    Precedence: surrogate accelerator > explicit bayesian engine > optimiser
    (when ``bayesian.engine`` is none) > GLUE default.
    """
    if _stage_on(method.get("surrogate")):
        return "surrogate"
    bayes = str(method.get("bayesian", {}).get("engine", "glue")).lower()
    if bayes in ESTIMATOR_REGISTRY and bayes not in ("none", ""):
        return bayes
    opt = str(method.get("optimizer", {}).get("engine", "none")).lower()
    if bayes in ("none", "") and opt in _OPTIMIZER_ALIASES:
        return "optimizer"
    return "glue"


@_register_estimator("surrogate")
def _estimate_surrogate(work_cfg, space, setup, method, *, seed, n_workers, extras, progress):
    from .engines import run_glue, run_surrogate
    _, _crop, _exe, _specs, _run_root, obs, experiments, _trts = setup
    scorer = _results_scorer(work_cfg, setup, n_workers)
    sur = run_surrogate(work_cfg, space, scorer, progress=progress)
    glue = run_glue(sur.design, space.names, work_cfg, space=space)
    best = sur.obj_results[glue.best_sample_id]
    extras["surrogate_info"] = sur.info
    return CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                             design=glue.design, obj_results=sur.obj_results,
                             best_theta=glue.best_theta, best=best, glue=glue, extras=extras)


@_register_estimator("smc_pf")
def _estimate_smc_pf(work_cfg, space, setup, method, *, seed, n_workers, extras, progress):
    from .engines import run_smc_pf
    _, _crop, _exe, _specs, _run_root, obs, experiments, _trts = setup
    smc = run_smc_pf(work_cfg, progress=progress)
    extras.update(initial_design=smc.initial_design, engine="smc_pf")
    return CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                             design=smc.design, obj_results=smc.obj_results,
                             best_theta=smc.best_theta, best=smc.best, glue=smc, extras=extras)


@_register_estimator("mcmc")
def _estimate_mcmc(work_cfg, space, setup, method, *, seed, n_workers, extras, progress):
    from .engines import run_mcmc
    _, _crop, _exe, _specs, _run_root, obs, experiments, _trts = setup
    scorer = _results_scorer(work_cfg, setup, n_workers)
    mc = run_mcmc(work_cfg, scorer, space, progress=progress)
    extras.update(initial_design=mc.initial_design, engine="mcmc",
                  mcmc_chain=mc.chain, acceptance=mc.acceptance)
    return CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                             design=mc.design, obj_results=mc.obj_results,
                             best_theta=mc.best_theta, best=mc.best, glue=mc, extras=extras)


@_register_estimator("optimizer")
def _estimate_optimizer(work_cfg, space, setup, method, *, seed, n_workers, extras, progress):
    from .engines import run_optimizer
    _, _crop, _exe, _specs, _run_root, obs, experiments, _trts = setup
    ocfg = method.get("optimizer", {})
    opt = str(ocfg.get("engine", "none")).lower()
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
    return CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                             design=design, obj_results={0: best},
                             best_theta=ores.best_theta, best=best, glue=None, extras=extras)


@_register_estimator("dream")
def _estimate_dream(work_cfg, space, setup, method, *, seed, n_workers, extras, progress):
    from .engines import run_dream
    _, _crop, _exe, _specs, _run_root, obs, experiments, _trts = setup
    scorer = _results_scorer(work_cfg, setup, n_workers)
    mc = run_dream(work_cfg, scorer, space, progress=progress)
    extras.update(initial_design=mc.initial_design, engine="dream",
                  mcmc_chain=mc.chain, acceptance=mc.acceptance)
    return CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                             design=mc.design, obj_results=mc.obj_results,
                             best_theta=mc.best_theta, best=mc.best, glue=mc, extras=extras)


@_register_estimator("es_mda")
def _estimate_es_mda(work_cfg, space, setup, method, *, seed, n_workers, extras, progress):
    from .engines import run_es_mda
    _, _crop, _exe, _specs, _run_root, obs, experiments, _trts = setup
    scorer = _results_scorer(work_cfg, setup, n_workers)
    es = run_es_mda(work_cfg, scorer, space, progress=progress)
    extras.update(initial_design=es.initial_design, engine="es_mda")
    return CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                             design=es.design, obj_results=es.obj_results,
                             best_theta=es.best_theta, best=es.best, glue=es, extras=extras)


@_register_estimator("bayesopt")
def _estimate_bayesopt(work_cfg, space, setup, method, *, seed, n_workers, extras, progress):
    from .engines import run_bayesopt
    _, _crop, _exe, _specs, _run_root, obs, experiments, _trts = setup
    scorer = _results_scorer(work_cfg, setup, n_workers)
    bo = run_bayesopt(work_cfg, scorer, space, progress=progress)
    extras.update(engine="bayesopt", bayesopt_info=bo.info)
    return CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                             design=bo.design, obj_results=bo.obj_results,
                             best_theta=bo.best_theta, best=bo.best, glue=None, extras=extras)


@_register_estimator("history")
def _estimate_history(work_cfg, space, setup, method, *, seed, n_workers, extras, progress):
    from .engines import run_history_matching
    _, _crop, _exe, _specs, _run_root, obs, experiments, _trts = setup
    scorer = _results_scorer(work_cfg, setup, n_workers)
    hm = run_history_matching(work_cfg, scorer, space, progress=progress)
    extras.update(engine="history", history_waves=hm.waves)
    return CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                             design=hm.design, obj_results=hm.obj_results,
                             best_theta=hm.best_theta, best=hm.best, glue=hm, extras=extras)


@_register_estimator("abc_smc")
def _estimate_abc_smc(work_cfg, space, setup, method, *, seed, n_workers, extras, progress):
    from .engines import run_abc_smc
    _, _crop, _exe, _specs, _run_root, obs, experiments, _trts = setup
    scorer = _results_scorer(work_cfg, setup, n_workers)
    abc = run_abc_smc(work_cfg, scorer, space, progress=progress)
    extras.update(engine="abc_smc", thresholds=abc.thresholds,
                  initial_design=abc.initial_design)
    return CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                             design=abc.design, obj_results=abc.obj_results,
                             best_theta=abc.best_theta, best=abc.best, glue=abc, extras=extras)


@_register_estimator("glue")
def _estimate_glue(work_cfg, space, setup, method, *, seed, n_workers, extras, progress):
    from .engines import run_glue
    sample_cfg = method.get("sample", {})
    n = int(sample_cfg.get("n", 200))
    engine = sample_cfg.get("engine", "lhs")
    include_start = bool(sample_cfg.get("include_start", True))
    samples = sample(space, n=n, engine=engine, seed=seed, include_start=include_start)
    design, obj_results, (space, obs, experiments) = evaluate_design(work_cfg, samples, progress=progress)
    spawn_manifest = design.attrs.get("spawn_manifest")
    glue = run_glue(design, space.names, work_cfg, space=space)
    best = obj_results[glue.best_sample_id]
    extras.update(engine="glue")
    if spawn_manifest is not None:
        extras["spawn_manifest"] = spawn_manifest
    return CalibrationResult(cfg=work_cfg, space=space, obs=obs, experiments=experiments,
                             design=glue.design, obj_results=obj_results,
                             best_theta=glue.best_theta, best=best, glue=glue, extras=extras)


def calibrate(cfg: dict, *, progress=True) -> CalibrationResult:
    """Run the configured calibration pipeline.

    Steps (each optional one only runs when switched on):
      0. resolve ``method.preset`` -> per-stage engines
      1. [sensitivity] screen parameters; optionally keep only influential ones
      2. [select]      AgMIP stepwise BIC/AICc -> keep the parameters that earn their place
      3. estimate with the chosen engine: glue | smc_pf | mcmc | optimizer | surrogate
      4. [multiobjective] NSGA-II Pareto front add-on
    """
    from .engines import run_nsga2, run_sensitivity, stepwise_select
    from .sparse import apply_sparse_config, calibrate_staged

    if not cfg.get("_sparse_applied", False):
        cfg = apply_sparse_config(cfg)
    if (cfg.get("method", {}) or {}).get("staged", {}).get("active", False):
        return calibrate_staged(cfg, progress=progress)

    cfg = _apply_staging(cfg)
    method = _resolve_method(cfg)
    seed = int(cfg["calibrator"].get("seed", 42))
    n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))
    extras: dict = {}

    # Work on a (possibly param-pruned) config; estimators read active params from it.
    work_cfg = cfg
    setup = _setup(work_cfg)
    space = setup[0]
    _warn_unmatched(setup[5], work_cfg)

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

    # --- 3. main estimator (registry dispatch) ------------------------------------
    estimator = _resolve_estimator(method)
    if progress:
        bayes = str(method.get("bayesian", {}).get("engine", "glue")).lower()
        opt = str(method.get("optimizer", {}).get("engine", "none")).lower()
        main = bayes if bayes not in ("none", "") else opt
        label = f"surrogate+{main}" if estimator == "surrogate" else \
            (opt if estimator == "optimizer" else main)
        print(f"[3/4] Estimating with {label}...", flush=True)

    result = ESTIMATOR_REGISTRY[estimator](
        work_cfg, space, setup, method,
        seed=seed, n_workers=n_workers, extras=extras, progress=progress)
    # The GLUE estimator rebuilds the space from its own setup; keep the local
    # ``space`` aligned with the result so the NSGA-II add-on uses the same one.
    space = result.space

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

    Used to recover full PlantGro curves for the best fit when plotting. Runs
    multi-treatment experiments as separate treatment spawns, then merges the
    parsed outputs back by experiment for scoring/report compatibility.
    """
    space, crop, exe, specs, run_root, obs, exps, treatments = _setup(cfg)
    experiments = experiments or exps
    raw: dict[str, list[SpawnResult]] = {}
    spawn_timeout = cfg.get("calibrator", {}).get("spawn_timeout")
    for e, trt in _treatment_units(list(experiments), treatments, obs):
        raw.setdefault(e, []).append(
            spawn_and_run(dict(theta), exp_id=e, cfg=cfg, crop=crop, param_specs=specs,
                          run_root=run_root, treatments=[trt], exe=exe,
                          **({"timeout": int(spawn_timeout)} if spawn_timeout else {}))
        )
    return _merge_result_map(raw)


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


def _year_key(exp: str) -> str:
    # DSSAT experiment codes are SSII<YY><NN> (site/inst, 2-digit year, 2-digit
    # number), so the year is the FIRST two of the trailing digits, e.g.
    # "YUKU2101" -> "21".
    digits = "".join(ch for ch in exp if ch.isdigit())
    return digits[:2] if len(digits) >= 2 else exp


def _site_key(exp: str) -> str:
    return "".join(ch for ch in exp if not ch.isdigit()) or exp


def _make_folds(experiments: list[str], scheme: str, seed: int) -> list[tuple[str, list[str]]]:
    """Return ``[(fold_label, held_experiments)]`` for a CV scheme.

    ``loeo`` holds out one experiment at a time; ``year``/``site`` hold out a whole
    year/site group (more honest for transfer); ``random`` holds out random k-folds
    (``method.validation.k``, default 5).
    """
    exps = list(experiments)
    if scheme in ("loeo", "none", "", None):
        return [(e, [e]) for e in exps]
    if scheme == "year":
        groups: dict[str, list[str]] = {}
        for e in exps:
            groups.setdefault(_year_key(e), []).append(e)
        return [(f"year_{k}", v) for k, v in groups.items()]
    if scheme == "site":
        groups = {}
        for e in exps:
            groups.setdefault(_site_key(e), []).append(e)
        return [(f"site_{k}", v) for k, v in groups.items()]
    if scheme == "random":
        rng = np.random.default_rng(seed)
        shuffled = list(exps)
        rng.shuffle(shuffled)
        k = min(len(shuffled), 5)
        return [(f"fold_{i}", shuffled[i::k]) for i in range(k) if shuffled[i::k]]
    raise ValueError(f"Unknown validation scheme '{scheme}'. "
                     "Use loeo | year | site | random.")


def validate_cv(cfg: dict, *, scheme: str | None = None, progress=False) -> pd.DataFrame:
    """Generalised cross-validation: calibrate on train experiments, evaluate on held-out.

    ``scheme`` (or ``method.validation.scheme``): loeo | year | site | random.
    Returns a tidy table of calibration-vs-evaluation fit per fold, with a ``fold``
    label and the held-out experiment. With only a few site-years, prefer ``year``
    or ``site`` folds and read the gap between splits as your overfit signal.
    """
    space, crop, exe, specs, run_root, obs, experiments, treatments = _setup(cfg)
    n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))
    method = cfg.get("method", {})
    n = int(method.get("sample", {}).get("n", 100))
    seed = int(cfg["calibrator"].get("seed", 42))
    scheme = scheme or method.get("validation", {}).get("scheme", "loeo")

    rows = []
    for label, held in _make_folds(experiments, scheme, seed):
        train = [e for e in experiments if e not in set(held)]
        if not train or not held:
            continue
        cfg_train = {**cfg, "experiments": train}
        samples = sample(space, n=n, engine=method.get("sample", {}).get("engine", "lhs"),
                         seed=seed, include_start=True)
        design, obj_results, _ = evaluate_design(cfg_train, samples, progress=progress)
        from .engines import run_glue
        glue = run_glue(design, space.names, cfg_train, space=space)
        best_theta = glue.best_theta
        cal = obj_results[glue.best_sample_id]
        ev = _score_theta(best_theta, held, cfg={**cfg, "experiments": held}, crop=crop,
                          specs=specs, run_root=run_root, treatments=treatments, exe=exe,
                          obs=obs, n_workers=n_workers)
        for uv, m in cal.per_var.items():
            rows.append({"fold": label, "held_out": ",".join(held), "split": "calibration",
                         "variable": uv, **m})
        for uv, m in ev.per_var.items():
            rows.append({"fold": label, "held_out": ",".join(held), "split": "evaluation",
                         "variable": uv, **m})
    return pd.DataFrame(rows)


def validate_loeo(cfg: dict, *, progress=False) -> pd.DataFrame:
    """Leave-one-environment-out CV (back-compat wrapper over :func:`validate_cv`)."""
    return validate_cv(cfg, scheme="loeo", progress=progress)


def nowcast(cfg: dict, as_of_date, *, progress=True) -> dict:
    """Operational in-season run: (re)calibrate on data up to ``as_of_date``, persist
    the calibration, and forecast the target variable(s) forward.

    Persists ``results/<name>/nowcast_state.json`` (the latest best theta) and
    warm-starts the next call from it — so between cloud-free satellite passes you
    reuse the stored calibration, and a new observation triggers a refresh. When
    ``forecast.active`` is set, returns per-variable, per-experiment forecast tables
    (LAI percentiles, optionally anchored to the last observation).
    """
    import json
    from copy import deepcopy
    from . import forecast as fc

    name = cfg["calibrator"]["name"]
    out_dir = Path(cfg["calibrator"].get("results_dir", "results")) / name
    out_dir.mkdir(parents=True, exist_ok=True)
    state_path = out_dir / "nowcast_state.json"

    prev_theta = None
    if state_path.exists():
        try:
            prev_theta = json.loads(state_path.read_text()).get("theta")
        except Exception:
            prev_theta = None

    obs_all = _setup(cfg)[5].table
    ts = pd.Timestamp(as_of_date)
    filtered = obs_all[obs_all["date"].isna() | (obs_all["date"] <= ts)].copy()
    if progress:
        print(f"Nowcast as of {ts.date()}: {len(filtered)} observations in scope", flush=True)

    work = deepcopy(cfg)
    work["source"] = {**work.get("source", {}), "table": filtered}
    if prev_theta:                                   # warm start from the last calibration
        for _g, params in (work.get("parameters") or {}).items():
            if not isinstance(params, dict):
                continue
            for nm, spec in params.items():
                if isinstance(spec, dict) and nm in prev_theta:
                    spec["start"] = float(prev_theta[nm])
    recal_n = cfg.get("assimilation", {}).get("recalibration", {}).get("recal_sample_size")
    if recal_n:
        work.setdefault("method", {}).setdefault("sample", {})["n"] = int(recal_n)

    result = calibrate(work, progress=progress)
    state_path.write_text(json.dumps({"as_of": str(ts.date()), "theta": result.best_theta},
                                     indent=2, default=str))

    # last LAI observation per experiment, for the anchor-continuity correction
    last_obs = {}
    lai = filtered[(filtered["variable"] == "LAID") & filtered["date"].notna()]
    for exp, g in lai.groupby("exp_id"):
        r = g.sort_values("date").iloc[-1]
        last_obs[exp] = (pd.Timestamp(r["date"]), float(r["value"]))

    forecasts = {}
    if cfg.get("forecast", {}).get("active", False):
        for var in cfg.get("forecast", {}).get("variables", ["LAID"]):
            forecasts[var] = fc.forecast_lai(work, result, last_obs=last_obs, variable=var)

    return {"as_of": str(ts.date()), "best_theta": result.best_theta,
            "result": result, "forecast": forecasts, "last_obs": last_obs}


def assimilate(cfg: dict, *, progress=True) -> dict:
    """Run in-season data assimilation according to configuration.
    
    Supported modes:
    - "recalibration": Mid-season parameter re-estimation.
    - "forcing": Direct state replacement.
    - "enkf": Ensemble Kalman Filter.
    """
    from .engines import InSeasonRecalibrator, ForcingAssimilator, EnsembleKalmanFilter

    assim_cfg = cfg.get("assimilation", {})
    mode = assim_cfg.get("mode", "recalibration")

    # Validate the mode FIRST (before the expensive _setup): enkf / forcing are
    # uncoupled prototypes (state is never re-injected into a running DSSAT
    # simulation), so refuse to run unless the user explicitly opts in.
    if mode not in ("recalibration", "enkf", "forcing"):
        raise ValueError(f"Unknown assimilation mode '{mode}'. "
                         "Expected: recalibration | enkf | forcing.")
    if mode in ("enkf", "forcing") and not assim_cfg.get("allow_uncoupled", False):
        raise NotImplementedError(
            f"Assimilation mode '{mode}' is an UNCOUPLED prototype: the updated state "
            "is not fed back into a running DSSAT simulation, so its output is "
            "illustrative only. Set assimilation.allow_uncoupled: true to run it "
            "anyway, or use mode: recalibration (the coupled in-season path)."
        )
    if mode in ("enkf", "forcing"):
        logger.warning("Running UNCOUPLED assimilation mode '%s' — output is "
                       "illustrative only (not coupled to DSSAT state).", mode)

    setup = _setup(cfg)
    space, crop, exe, specs, run_root, obs, experiments, treatments = setup
    _warn_unmatched(obs, cfg)

    if progress:
        print(f"Starting in-season data assimilation (mode: {mode})...", flush=True)

    results = {}

    if mode == "recalibration":
        obs_df = obs.table
        valid_dates = obs_df[obs_df["date"].notna()]["date"].dt.date.unique()
        valid_dates = sorted(valid_dates)
        
        if not valid_dates:
            logger.warning("No time-series observations found for recalibration.")
            return {}
            
        freq = assim_cfg.get("recalibration", {}).get("update_frequency", "on_observation")
        if freq == "weekly" and len(valid_dates) > 1:
            resampled = pd.date_range(start=min(valid_dates), end=max(valid_dates), freq="W").date
            valid_dates = [d for d in resampled if d in set(valid_dates)]
        elif freq == "biweekly" and len(valid_dates) > 1:
            resampled = pd.date_range(start=min(valid_dates), end=max(valid_dates), freq="2W").date
            valid_dates = [d for d in resampled if d in set(valid_dates)]
            
        engine = InSeasonRecalibrator(cfg)
        warm_start = assim_cfg.get("recalibration", {}).get("warm_start", True)
        trace = []
        prev_theta = None
        for d in valid_dates:
            if progress:
                print(f"  Recalibrating at checkpoint date: {d}", flush=True)
            best_theta = engine.recalibrate(obs_df, d,
                                            warm_start_theta=prev_theta if warm_start else None)
            trace.append({"date": d, "theta": best_theta})
            prev_theta = best_theta

        results = {
            "mode": mode,
            "trace": trace,
            "final_theta": trace[-1]["theta"] if trace else None
        }
        
    elif mode == "forcing":
        engine = ForcingAssimilator(cfg)
        obs_df = obs.table
        dated_obs = obs_df[obs_df["date"].notna()].sort_values("date")
        
        state_history = []
        current_state = {}
        
        for _, r in dated_obs.iterrows():
            current_state = engine.apply(current_state, {
                "variable": r["variable"],
                "value": r["value"],
                "confidence": r["weight"] if not pd.isna(r["weight"]) else 1.0
            })
            state_history.append({
                "date": r["date"],
                "variable": r["variable"],
                "state": current_state.copy()
            })
            
        results = {
            "mode": mode,
            "state_history": state_history,
            "final_state": current_state
        }
        
    elif mode == "enkf":
        engine = EnsembleKalmanFilter(cfg)
        obs_df = obs.table
        n_ens = engine.n_ensemble
        n_vars = len(engine.state_vars)
        ensemble = np.random.default_rng(42).normal(loc=1.0, scale=0.2, size=(n_ens, n_vars))
        
        dated_obs = obs_df[obs_df["date"].notna()].sort_values("date")
        filter_history = []
        
        for _, r in dated_obs.iterrows():
            var = r["variable"]
            if var in engine.state_vars:
                obs_val = r["value"]
                obs_sig = r["sigma"] if not pd.isna(r["sigma"]) else 0.1
                
                ensemble = engine.assimilate(ensemble, var, obs_val, obs_sig)
                
                filter_history.append({
                    "date": r["date"],
                    "variable": var,
                    "mean_state": ensemble.mean(axis=0).tolist(),
                    "std_state": ensemble.std(axis=0).tolist()
                })
                
        results = {
            "mode": mode,
            "filter_history": filter_history,
            "final_ensemble_mean": ensemble.mean(axis=0).tolist()
        }
        
    return results


def combined_mode(cfg: dict, *, progress=True) -> dict:
    """Run combined mode: parameter calibration followed by in-season state assimilation."""
    if progress:
        print("=== Step 1: Base Parameter Calibration ===", flush=True)
    cal_result = calibrate(cfg, progress=progress)
    
    if progress:
        print("=== Step 2: In-Season State Assimilation ===", flush=True)
    from copy import deepcopy
    cfg_assim = deepcopy(cfg)

    for group, params in cfg_assim.get("parameters", {}).items():
        if not isinstance(params, dict):
            continue
        for name, spec in params.items():
            if isinstance(spec, dict) and name in cal_result.best_theta:
                spec["start"] = float(cal_result.best_theta[name])
                
    assim_results = assimilate(cfg_assim, progress=progress)
    
    return {
        "calibration": cal_result,
        "assimilation": assim_results
    }


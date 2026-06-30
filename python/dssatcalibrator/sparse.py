"""Sparse-data cultivar/species calibration helpers.

These tools are intentionally modular.  They do not replace the simple GLUE or
optimizer paths; they wrap them with safer priors, staged estimation, cautious
history matching, and measurement recommendations for cases where observations
are scarce.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass

import numpy as np
import pandas as pd

_PHENOLOGY_HINT = {"CSDL", "PPSEN", "EM-FL", "FL-SD", "SD-PM", "FL-SH", "FL-LF",
                   "PL-EM", "PLEM", "P1", "P2", "P3", "P4", "P5", "P1V", "P1D",
                   "PHINT", "EM-V1", "PHTHRS", "ADAP", "MDAP"}


def _iter_specs(cfg: dict):
    for group, params in (cfg.get("parameters") or {}).items():
        if not isinstance(params, dict):
            continue
        for name, spec in params.items():
            if isinstance(spec, dict):
                yield group, name, spec


def _center_bounds(spec: dict, center: float, half_width: float) -> tuple[float, float]:
    lo0, hi0 = float(spec["min"]), float(spec["max"])
    lo, hi = max(lo0, center - half_width), min(hi0, center + half_width)
    if lo >= hi:
        lo, hi = lo0, hi0
    return lo, hi


def _set_normal_prior(spec: dict, center: float, sd: float, *, overwrite: bool) -> None:
    if overwrite or not spec.get("prior"):
        spec["prior"] = {"dist": "normal", "mean": float(center), "sd": float(max(sd, 1e-9))}


def apply_delta_from_analog(cfg: dict) -> dict:
    """Constrain active parameters around an analog cultivar/species.

    This implements the practical form of delta calibration:
    ``theta_new = theta_analog + delta``.  DSSAT still receives native
    coefficients, but the bounds, starts, and priors are centered on the analog.
    """
    out = deepcopy(cfg)
    dc = (out.get("sparse", {}) or {}).get("delta_from_analog", {}) or {}
    if not dc.get("active", False):
        return out
    analog = dc.get("theta") or dc.get("analog_theta") or {}
    scale = float(dc.get("relative_width", dc.get("scale", 0.20)))
    min_width = float(dc.get("min_width", 1e-6))
    widths = dc.get("widths", {}) or {}
    overwrite = bool(dc.get("overwrite_prior", True))
    for _group, name, spec in _iter_specs(out):
        if name not in analog:
            continue
        center = float(analog[name])
        half_width = float(widths.get(name, max(abs(center) * scale, min_width)))
        lo, hi = _center_bounds(spec, center, half_width)
        spec["min"], spec["max"] = float(lo), float(hi)
        spec["start"] = float(np.clip(center, lo, hi))
        spec["analog"] = center
        _set_normal_prior(spec, center, half_width / 2.0, overwrite=overwrite)
    return out


def apply_hierarchical_priors(cfg: dict) -> dict:
    """Apply empirical-Bayes partial-pooling priors from analog populations.

    The current DSSAT runner has one parameter vector per run, so this is not a
    full multi-cultivar hierarchical sampler.  It is the sparse-data workhorse:
    active coefficients are shrunk toward analog-population means with declared
    between-cultivar standard deviations.
    """
    out = deepcopy(cfg)
    hc = (out.get("sparse", {}) or {}).get("hierarchical_priors", {}) or {}
    if not hc.get("active", False):
        return out
    params = deepcopy(hc.get("parameters", {}) or {})
    population = hc.get("population", []) or []
    if population:
        names = sorted({k for row in population if isinstance(row, dict) for k in row})
        for name in names:
            vals = np.array([float(row[name]) for row in population
                             if isinstance(row, dict) and name in row], dtype=float)
            if len(vals):
                params.setdefault(name, {})
                params[name].setdefault("mean", float(np.mean(vals)))
                params[name].setdefault("sd", float(np.std(vals, ddof=1) if len(vals) > 1 else 0.1 * abs(np.mean(vals) or 1.0)))
    bounds_sd = float(hc.get("bounds_sd", 3.0))
    overwrite = bool(hc.get("overwrite_prior", False))
    for _group, name, spec in _iter_specs(out):
        if name not in params:
            continue
        pc = params[name] or {}
        center = float(pc.get("mean", pc.get("mu", spec.get("start", 0.5 * (float(spec["min"]) + float(spec["max"]))))))
        sd = float(pc.get("sd", pc.get("sigma", max(0.1 * abs(center), 1.0))))
        lo, hi = _center_bounds(spec, center, bounds_sd * sd)
        if bool(hc.get("shrink_bounds", False)):
            spec["min"], spec["max"] = float(lo), float(hi)
        if bool(hc.get("set_start", True)):
            spec["start"] = float(np.clip(center, float(spec["min"]), float(spec["max"])))
        spec["hierarchical_mean"] = center
        spec["hierarchical_sd"] = sd
        _set_normal_prior(spec, center, sd, overwrite=overwrite)
    return out


def _trait_center(name: str, spec: dict, traits: dict):
    start = float(spec.get("start", 0.5 * (float(spec["min"]) + float(spec["max"]))))
    if "maturity_days" in traits and name.upper() in _PHENOLOGY_HINT:
        # Conservative default: later-maturity material nudges phenology thermal
        # duration upward but keeps the declared bounds authoritative.
        delta = float(traits["maturity_days"]) - 120.0
        strength = 0.35 if name.upper() in {"SD-PM", "P5", "PHTHRS"} else 0.15
        return start + strength * delta
    if "photoperiod_sensitivity" in traits and name.upper() in {"PPSEN", "CSDL", "P1D"}:
        return start * (1.0 + 0.25 * float(traits["photoperiod_sensitivity"]))
    return None


def apply_trait_priors(cfg: dict) -> dict:
    """Convert known traits into parameter priors.

    Prefer explicit rules in ``sparse.trait_priors.rules``.  Built-in rules are
    deliberately weak and only affect common phenology/photoperiod coefficients.
    """
    out = deepcopy(cfg)
    tc = (out.get("sparse", {}) or {}).get("trait_priors", {}) or {}
    if not tc.get("active", False):
        return out
    traits = tc.get("traits") or (out.get("sparse", {}) or {}).get("traits", {}) or {}
    rules = tc.get("rules", {}) or {}
    overwrite = bool(tc.get("overwrite_prior", False))
    for _group, name, spec in _iter_specs(out):
        center = None
        sd = None
        rule = rules.get(name)
        if rule:
            if "value" in rule:
                center = float(rule["value"])
            elif "trait" in rule and rule["trait"] in traits:
                center = float(rule.get("intercept", 0.0)) + float(rule.get("slope", 1.0)) * float(traits[rule["trait"]])
            sd = rule.get("sd")
        if center is None:
            center = _trait_center(name, spec, traits)
        if center is None:
            continue
        center = float(np.clip(center, float(spec["min"]), float(spec["max"])))
        span = float(spec["max"]) - float(spec["min"])
        sd = float(sd if sd is not None else max(0.20 * span, 1e-6))
        if bool(tc.get("set_start", True)):
            spec["start"] = center
        spec["trait_prior_center"] = center
        _set_normal_prior(spec, center, sd, overwrite=overwrite)
    return out


def apply_sparse_config(cfg: dict) -> dict:
    """Apply all opt-in sparse-data config transforms exactly once."""
    if cfg.get("_sparse_applied", False):
        return cfg
    out = deepcopy(cfg)
    out = apply_delta_from_analog(out)
    out = apply_hierarchical_priors(out)
    out = apply_trait_priors(out)
    out["_sparse_applied"] = True
    return out


def make_quick_dirty_config(cfg: dict, *, n: int | None = None,
                            engine: str = "lhs", estimator: str = "glue") -> dict:
    """Return a small, robust calibration config for first-pass fitting."""
    out = apply_sparse_config(cfg)
    out = deepcopy(out)
    out.setdefault("method", {})
    out["method"]["preset"] = "C"
    out["method"]["sample"] = {"engine": engine, "n": int(n or 40)}
    out["method"]["bayesian"] = {"engine": estimator, "behavioural_quantile": 0.25}
    out["method"].pop("staged", None)
    out.setdefault("objective", {}).setdefault("likelihood", {"type": "huber", "delta": 2.0})
    return out


def calibrate_quick(cfg: dict, *, n: int | None = None, progress: bool = True):
    """Run the quick first-pass calibration."""
    from .orchestrator import calibrate
    return calibrate(make_quick_dirty_config(cfg, n=n), progress=progress)


def _active_names(cfg: dict) -> list[str]:
    return [name for _group, name, spec in _iter_specs(cfg) if spec.get("active", False)]


def _stage_keep(cfg: dict, stage: dict) -> list[str]:
    active = set(_active_names(cfg))
    keep = set(stage.get("params", []) or [])
    groups = set(stage.get("groups", []) or [])
    roles = set(stage.get("roles", []) or [])
    for group, name, spec in _iter_specs(cfg):
        if name not in active:
            continue
        if group in groups or str(spec.get("role", "")) in roles:
            keep.add(name)
    lname = str(stage.get("name", "")).lower()
    if not keep and "phen" in lname:
        keep = {name for _g, name, spec in _iter_specs(cfg)
                if name in active and (name.upper() in _PHENOLOGY_HINT or spec.get("role") == "obligatory")}
    if not keep:
        keep = active
    return [n for n in _active_names(cfg) if n in keep]


def _apply_active_subset_local(cfg: dict, keep: list[str]) -> dict:
    out = deepcopy(cfg)
    keep = set(keep)
    for _group, name, spec in _iter_specs(out):
        if spec.get("active", False):
            spec["active"] = name in keep
    return out


def _focus_objective(cfg: dict, variables: list[str] | None) -> dict:
    if not variables:
        return cfg
    out = deepcopy(cfg)
    known = set((out.get("engine", {}) or {}).get("timeseries_outputs", {}).keys())
    known |= set((out.get("engine", {}) or {}).get("scalar_outputs", {}).keys())
    weights = dict((out.get("objective", {}) or {}).get("weights", {}) or {})
    for v in known:
        weights[v] = 1.0 if v in set(variables) else 0.0
    out.setdefault("objective", {})["weights"] = weights
    return out


def default_stages(cfg: dict) -> list[dict]:
    return [
        {"name": "phenology", "roles": ["obligatory"],
         "variables": ["phenology", "anthesis", "maturity"]},
        {"name": "canopy_growth", "roles": ["candidate"],
         "variables": ["lai", "biomass"]},
        {"name": "yield", "groups": ["genetic_cultivar", "management"],
         "variables": ["grain_yield", "yield"]},
    ]


def build_staged_configs(cfg: dict) -> list[tuple[str, dict]]:
    """Build per-stage configs for staged sparse-data calibration."""
    base = apply_sparse_config(cfg)
    scfg = (base.get("method", {}) or {}).get("staged", {}) or {}
    stages = scfg.get("stages") or default_stages(base)
    out = []
    for stage in stages:
        keep = _stage_keep(base, stage)
        st = _apply_active_subset_local(base, keep)
        st = _focus_objective(st, stage.get("variables"))
        st.setdefault("method", {})
        st["method"]["staged"] = {"active": False}
        if stage.get("n") is not None:
            st.setdefault("method", {}).setdefault("sample", {})["n"] = int(stage["n"])
        if stage.get("engine"):
            st.setdefault("method", {}).setdefault("bayesian", {})["engine"] = stage["engine"]
        out.append((str(stage.get("name", f"stage_{len(out)+1}")), st))
    return out


def calibrate_staged(cfg: dict, *, progress: bool = True):
    """Run sequential staged calibration, warm-starting each stage from the last."""
    from .orchestrator import calibrate
    stages = build_staged_configs(cfg)
    results = []
    prev_theta = None
    for name, stage_cfg in stages:
        stage_cfg = deepcopy(stage_cfg)
        if prev_theta:
            for _group, pname, spec in _iter_specs(stage_cfg):
                if pname in prev_theta:
                    spec["start"] = float(prev_theta[pname])
        if progress:
            print(f"[staged] {name}: {len(_active_names(stage_cfg))} active parameter(s)", flush=True)
        res = calibrate(stage_cfg, progress=progress)
        results.append({"stage": name, "result": res, "active": _active_names(stage_cfg)})
        prev_theta = {**(prev_theta or {}), **res.best_theta}
    final = results[-1]["result"]
    final.extras.setdefault("stages", results)
    final.best_theta = prev_theta or final.best_theta
    return final


def apply_identifiability_gate(cfg: dict, result, *, max_sd_ratio: float | None = None,
                               protect_roles: tuple[str, ...] = ("obligatory",)):
    """Freeze weakly identified active parameters for a follow-up run."""
    from .diagnostics import identifiability
    out = deepcopy(cfg)
    gate = (out.get("sparse", {}) or {}).get("identifiability_gate", {}) or {}
    max_sd_ratio = float(max_sd_ratio if max_sd_ratio is not None
                         else gate.get("max_sd_ratio", 0.8))
    diag = identifiability(result)
    weak = set(diag.loc[diag["sd_ratio"] >= max_sd_ratio, "parameter"].tolist()) if not diag.empty else set()
    frozen = []
    for _group, name, spec in _iter_specs(out):
        if name in weak and spec.get("role") not in protect_roles:
            spec["active"] = False
            frozen.append(name)
    return {"cfg": out, "diagnostics": diag, "frozen": frozen}


@dataclass
class ObservationRecommendation:
    table: pd.DataFrame


def recommend_observations(result, candidates: list[dict] | None = None) -> ObservationRecommendation:
    """Rank candidate next observations by a posterior-spread proxy.

    If per-sample residuals are available, variables whose simulations still vary
    strongly across plausible parameter sets receive high utility.  This is a
    lightweight expected-information proxy, not a full decision-theoretic design.
    """
    cfg = getattr(result, "cfg", {}) or {}
    if candidates is None:
        vars_ = sorted(set((cfg.get("engine", {}) or {}).get("timeseries_outputs", {}).keys()) |
                       set((cfg.get("engine", {}) or {}).get("scalar_outputs", {}).keys()))
        if not vars_ and getattr(result, "obj_results", None):
            vars_ = sorted({
                uv
                for ores in result.obj_results.values()
                for uv in getattr(ores, "residuals", pd.DataFrame()).get("user_var", pd.Series(dtype=str)).dropna().unique()
            })
        candidates = [{"variable": v} for v in vars_]
    weights = result.design.get("weight", pd.Series(np.ones(len(result.design)) / max(len(result.design), 1))).to_numpy(dtype=float)
    rows = []
    for cand in candidates:
        uv = cand.get("variable")
        sims, sigmas = [], []
        for sid, ores in (result.obj_results or {}).items():
            resid = getattr(ores, "residuals", None)
            if resid is None or resid.empty or uv not in set(resid["user_var"]):
                continue
            g = resid[resid["user_var"] == uv]
            sims.append(float(g["sim"].mean()))
            sigmas.append(float(g["sigma"].mean()))
        if len(sims) >= 2:
            w = weights[:len(sims)]
            w = w / w.sum() if w.sum() > 0 else np.ones(len(sims)) / len(sims)
            mean = float(np.sum(w * np.asarray(sims)))
            var = float(np.sum(w * (np.asarray(sims) - mean) ** 2))
            sigma = float(np.nanmean(sigmas)) if sigmas else 1.0
            utility = var / max(sigma ** 2, 1e-12)
        else:
            utility = 1.0
        rows.append({**cand, "utility": float(utility), "available_samples": len(sims)})
    table = pd.DataFrame(rows).sort_values("utility", ascending=False).reset_index(drop=True)
    return ObservationRecommendation(table=table)

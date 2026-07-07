"""Parameter search space built from the config's active parameters."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import active_parameters


_GLOBAL_SCOPES = {"", "global", "shared", "pooled", "pool"}
_EXPERIMENT_SCOPES = {
    "experiment", "experiments", "per_experiment", "per-experiment",
    "experiment_specific", "experiment-specific", "local",
}
_CULTIVAR_SCOPES = {
    "cultivar", "cultivars", "per_cultivar", "per-cultivar",
    "cultivar_specific", "cultivar-specific",
}


def _default_scope_for_group(cfg: dict, group: str | None) -> str:
    defaults = cfg.get("parameter_defaults", {}) or {}
    by_group = defaults.get("scope_by_group", {}) or {}
    if group in by_group:
        return str(by_group[group])
    return str(defaults.get("scope", "global"))


def _scope_of(cfg: dict, spec: dict) -> str:
    default = _default_scope_for_group(cfg, spec.get("group"))
    raw = str(spec.get("scope", spec.get("pooling", default))).lower()
    if raw in _EXPERIMENT_SCOPES:
        return "experiment"
    if raw in _CULTIVAR_SCOPES:
        return "cultivar"
    if raw in _GLOBAL_SCOPES:
        return "global"
    raise ValueError(
        f"Unknown parameter scope {raw!r} for {spec.get('group')}.{spec.get('name')}; "
        "use 'global', 'experiment', or 'cultivar'."
    )


def _scoped_name(base: str, scope_id: str | None = None) -> str:
    return f"{base}__{scope_id}" if scope_id else base


def _as_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(v) for v in value]


def _cultivars_for_scope(cfg: dict, spec: dict) -> list[str]:
    explicit = _as_list(
        spec.get("cultivars", spec.get("cultivar_anchors", spec.get("cultivar_codes")))
    )
    if explicit:
        return explicit
    cultivars: list[str] = []
    for crop in cfg.get("crops") or []:
        if crop.get("calibration_cultivars"):
            cultivars.extend(_as_list(crop.get("calibration_cultivars")))
            continue
        cultivars.extend(_as_list(crop.get("cultivar_anchors")))
        if crop.get("cultivar_anchor"):
            cultivars.append(str(crop["cultivar_anchor"]))
    out = []
    seen = set()
    for cultivar in cultivars:
        if cultivar and cultivar not in seen:
            seen.add(cultivar)
            out.append(cultivar)
    return out


def _apply_cultivar_overrides(rec: dict, cultivar: str) -> dict:
    """Apply per-cultivar start/bound overrides to an expanded spec."""
    aliases = {
        "start": ("start_by_cultivar", "starts_by_cultivar"),
        "min": ("min_by_cultivar", "mins_by_cultivar"),
        "max": ("max_by_cultivar", "maxs_by_cultivar"),
    }
    for field, keys in aliases.items():
        for key in keys:
            values = rec.get(key)
            if isinstance(values, dict) and cultivar in values:
                rec[field] = values[cultivar]
                break
    return rec


def expand_parameter_specs(cfg: dict, specs: list[dict]) -> list[dict]:
    """Expand active specs into optimizer dimensions.

    By default a parameter is global: one theta value is shared by all
    experiments. With ``scope: experiment`` the parameter becomes one optimizer
    dimension per configured experiment, while retaining ``base_name`` so the
    spawn layer writes the original DSSAT coefficient name into that experiment's
    local CUL/ECO/SPE copy.
    """
    experiments = list(cfg.get("experiments") or [])
    out: list[dict] = []
    for spec in specs:
        base = spec["name"]
        scope = _scope_of(cfg, spec)
        if scope == "experiment":
            if not experiments:
                raise ValueError(f"Parameter {spec.get('group')}.{base} has scope=experiment but no experiments are configured.")
            for exp_id in experiments:
                rec = dict(spec)
                rec["base_name"] = base
                rec["name"] = _scoped_name(base, str(exp_id))
                rec["scope"] = "experiment"
                rec["exp_id"] = str(exp_id)
                out.append(rec)
        elif scope == "cultivar":
            cultivars = _cultivars_for_scope(cfg, spec)
            if not cultivars:
                raise ValueError(
                    f"Parameter {spec.get('group')}.{base} has scope=cultivar "
                    "but no cultivars are configured."
                )
            for cultivar in cultivars:
                rec = dict(spec)
                rec["base_name"] = base
                rec["name"] = _scoped_name(base, str(cultivar))
                rec["scope"] = "cultivar"
                rec["cultivar"] = str(cultivar)
                rec = _apply_cultivar_overrides(rec, str(cultivar))
                out.append(rec)
        else:
            rec = dict(spec)
            rec["base_name"] = base
            rec["scope"] = "global"
            out.append(rec)
    return out


@dataclass
class ParameterSpace:
    names: list[str]
    low: np.ndarray
    high: np.ndarray
    start: np.ndarray
    specs: list[dict]

    @property
    def ndim(self) -> int:
        return len(self.names)

    def to_theta(self, vector) -> dict[str, float]:
        """Map a 1-D parameter vector (in native units) to a named theta dict."""
        return {n: float(v) for n, v in zip(self.names, vector)}

    def clip(self, vector) -> np.ndarray:
        return np.clip(np.asarray(vector, float), self.low, self.high)

    @classmethod
    def from_config(cls, cfg: dict) -> "ParameterSpace":
        specs = expand_parameter_specs(cfg, active_parameters(cfg))
        if not specs:
            raise ValueError("No active parameters in config (set active: true on some).")
        names = [s["name"] for s in specs]
        low = np.array([float(s["min"]) for s in specs])
        high = np.array([float(s["max"]) for s in specs])
        start = np.array([float(s.get("start", 0.5 * (lo + hi)))
                          for s, lo, hi in zip(specs, low, high)])
        return cls(names=names, low=low, high=high, start=np.clip(start, low, high), specs=specs)

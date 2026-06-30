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


def _scope_of(spec: dict) -> str:
    raw = str(spec.get("scope", spec.get("pooling", "global"))).lower()
    if raw in _EXPERIMENT_SCOPES:
        return "experiment"
    if raw in _GLOBAL_SCOPES:
        return "global"
    raise ValueError(
        f"Unknown parameter scope {raw!r} for {spec.get('group')}.{spec.get('name')}; "
        "use 'global' or 'experiment'."
    )


def _scoped_name(base: str, exp_id: str | None = None) -> str:
    return f"{base}__{exp_id}" if exp_id else base


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
        scope = _scope_of(spec)
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

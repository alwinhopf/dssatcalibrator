"""Configuration loading: ``environment > YAML > built-in default``.

A run is described by one YAML file (see ``config_hemp.yaml``). This module
loads it, merges defaults, applies a handful of environment overrides, and
offers small helpers to enumerate the *active* parameters and resolve paths.

Design choice: the config is kept as plain nested dicts (not bespoke classes) so
it round-trips to YAML/JSON unchanged and the manifest of a run is just the dict.
"""

from __future__ import annotations

import os
from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

# Minimal built-in defaults; the YAML supplies the rest.
DEFAULTS: dict[str, Any] = {
    "calibrator": {
        "name": "run",
        "seed": 42,
        "workdir": "results/_workdir",   # DSSAT spawn scratch lives under results/
        "results_dir": "results",
        "figures_dir": "figures",        # all PNGs collected here, by run name
        "dssat_exe": "",
        "dssat_dir": "C:/DSSAT48",
        "num_cores": 0,
        "cache_spawns": True,
        "keep_run_dirs": False,
    },
    "method": {
        "preset": "C",
        "sample": {"engine": "lhs", "n": 200},
        "validation": {"scheme": "none"},
    },
    "objective": {"weighting": "unified", "weights": {}, "error_model": {}},
    "parameters": {},
    "crops": [],
    "experiments": [],
}

# Environment overrides: ENV name -> (section, key)
_ENV_OVERRIDES = {
    "DSSATCAL_DSSAT_EXE": ("calibrator", "dssat_exe"),
    "DSSATCAL_DSSAT_DIR": ("calibrator", "dssat_dir"),
    "DSSATCAL_NUM_CORES": ("calibrator", "num_cores"),
    "DSSATCAL_WORKDIR": ("calibrator", "workdir"),
    "DSSATCAL_RESULTS_DIR": ("calibrator", "results_dir"),
}


def _deep_merge(base: dict, over: dict) -> dict:
    out = deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = deepcopy(v)
    return out


def load_config(path: str | Path) -> dict:
    """Load a config YAML, merge defaults, and apply env overrides."""
    path = Path(path)
    with open(path, "r", encoding="utf-8") as fh:
        user = yaml.safe_load(fh) or {}
    cfg = _deep_merge(DEFAULTS, user)

    for env, (sec, key) in _ENV_OVERRIDES.items():
        val = os.environ.get(env)
        if val is not None and val != "":
            cfg.setdefault(sec, {})
            cfg[sec][key] = int(val) if key == "num_cores" else val

    cfg["_config_path"] = str(path.resolve())
    return cfg


def active_parameters(cfg: dict) -> list[dict]:
    """Flatten the ``parameters`` block to the list of ACTIVE parameter specs.

    Each returned dict carries: ``group``, ``name``, ``min``, ``max``, ``start``,
    plus any extra spec keys (``prior``, ``role``, ``type``, ``dssat`` ...).
    """
    out = []
    for group, params in (cfg.get("parameters") or {}).items():
        if not isinstance(params, dict):
            continue
        for name, spec in params.items():
            if not isinstance(spec, dict) or not spec.get("active", False):
                continue
            rec = {"group": group, "name": name}
            rec.update(spec)
            out.append(rec)
    return out


def all_parameters(cfg: dict) -> list[dict]:
    """Every declared parameter (active or not) — for documentation / start values."""
    out = []
    for group, params in (cfg.get("parameters") or {}).items():
        if not isinstance(params, dict):
            continue
        for name, spec in params.items():
            if isinstance(spec, dict):
                rec = {"group": group, "name": name}
                rec.update(spec)
                out.append(rec)
    return out


def crop_for(cfg: dict, code: str) -> dict:
    """Return the crop block matching a 2-letter DSSAT code (first crop if one)."""
    crops = cfg.get("crops") or []
    for c in crops:
        if c.get("code") == code:
            return c
    return crops[0] if crops else {}


def resolve_exe(cfg: dict) -> Path:
    """Resolve the DSSAT executable path (explicit, else dir default)."""
    exe = cfg["calibrator"].get("dssat_exe", "")
    if exe:
        return Path(exe)
    return Path(cfg["calibrator"]["dssat_dir"]) / "DSCSM048.EXE"

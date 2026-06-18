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
    # Execution backend: native keeps the historic local subprocess wrapper;
    # dssatengine delegates DSSAT spawning and DSSBatch writing to the shared engine.
    "execution": {"backend": "native"},
    # Shared FileX/genotype template directory for synthesized experiments and
    # new-crop scaffolding. Empty means DSSAT_TEMPLATE_DIR, then the sibling
    # DSSAT_Gridded_Run_Tutorial/dssat_templates directory if present.
    "templates": {"template_dir": ""},
    # Parameter-file gating safety (dssatcal model): CUL free, ECO gated, SPE blocked.
    # Set species: free only for new-species adaptation from an analog template.
    "gating": {"cultivar": "free", "ecotype": "gated", "species": "blocked"},
    # Optional: take planting date (and, future, other events) from ingested
    # farm-management rows and set them as FileX inputs instead of calibrating.
    "management_options": {"use_source_planting_date": False},
    # Optional weather-driver layer (default 'file' = use DSSAT's own .WTH).
    "weather": {"provider": "file", "gap_fill": "none", "horizon": 0,
                "cache_dir": "weather_cache"},
    # Optional soil acquisition for synthesized/new-site experiments. The real
    # experiment path keeps using DSSAT's installed Soil/ by default.
    "soil": {"provider": "file", "source": "ssurgo", "cache_dir": "soil_cache"},
    # Optional in-season LAI forecast/nowcast from the calibrated posterior.
    "forecast": {"active": False, "variables": ["LAID"], "n_ensemble": 0,
                 "anchor_continuity": True, "decay_days": 21},
    # Optional identifiability + structural-adequacy diagnostics after a run.
    "diagnostics": {"active": False},
    "observation_sources": {},
    "fusion": {
        "conflict_resolution": "keep_all",
        "source_priority": []
    },
    "assimilation": {
        # `active: true` makes ``run_calibration.py`` run the assimilation step
        # (equivalently, pass ``--assimilate``).
        "active": False,
        # Only "recalibration" is coupled to DSSAT: it re-estimates parameters from
        # the observations received so far by re-running the calibration pipeline.
        # "enkf" and "forcing" are UNCOUPLED prototypes — the updated state is never
        # re-injected into a running DSSAT simulation, so their output is
        # illustrative only. They refuse to run unless `allow_uncoupled: true`.
        "mode": "recalibration",          # recalibration | enkf | forcing
        "allow_uncoupled": False,
        "recalibration": {
            "engine": "glue",             # estimator used at each checkpoint
            "recal_sample_size": 100,     # design size per checkpoint (kept small for speed)
            "warm_start": True,           # seed each checkpoint with the previous best theta
            "update_frequency": "on_observation",   # on_observation | weekly | biweekly
        },
        "enkf": {                         # PROTOTYPE / uncoupled (see allow_uncoupled)
            "n_ensemble": 50,
            "inflation": 1.05,
            "state_variables": ["LAID", "CWAD"],
        },
        "forcing": {                      # PROTOTYPE / uncoupled (see allow_uncoupled)
            "min_confidence": 0.8,
            "smoothing": True,
        },
    },
}

# Environment overrides: ENV name -> (section, key)
_ENV_OVERRIDES = {
    "DSSATCAL_DSSAT_EXE": ("calibrator", "dssat_exe"),
    "DSSATCAL_DSSAT_DIR": ("calibrator", "dssat_dir"),
    "DSSATCAL_NUM_CORES": ("calibrator", "num_cores"),
    "DSSATCAL_WORKDIR": ("calibrator", "workdir"),
    "DSSATCAL_RESULTS_DIR": ("calibrator", "results_dir"),
    "DSSAT_TEMPLATE_DIR": ("templates", "template_dir"),
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


def resolve_dssat_paths(cfg: dict) -> dict[str, Path]:
    """Resolve the DSSAT48 install layout used by every spawn.

    The calibrator treats ``dssat_dir`` as the single root for the binary,
    ``Genotype/``, ``Weather/`` and ``Soil/``. Callers can still override the
    executable with ``calibrator.dssat_exe`` when using a custom compiled model.
    """
    root = Path(cfg["calibrator"]["dssat_dir"])
    return {
        "root": root,
        "exe": resolve_exe(cfg),
        "genotype": root / "Genotype",
        "weather": root / "Weather",
        "soil": root / "Soil",
    }


def _workspace_template_dir() -> Path:
    return (Path(__file__).resolve().parents[2] /
            "DSSAT_Gridded_Run_Tutorial" / "dssat_templates")


def resolve_template_dir(cfg: dict | None = None, *, required: bool = False) -> Path | None:
    """Resolve the shared ``dssat_templates`` directory.

    Precedence is ``DSSAT_TEMPLATE_DIR`` / ``templates.template_dir`` /
    top-level ``template_dir`` / the sibling gridded tutorial template folder.
    Returns ``None`` when no template directory is configured or discoverable,
    unless ``required`` is true.
    """
    cfg = cfg or {}
    configured = (
        os.environ.get("DSSAT_TEMPLATE_DIR")
        or (cfg.get("templates") or {}).get("template_dir")
        or cfg.get("template_dir")
    )
    if configured:
        return Path(configured)

    default = _workspace_template_dir()
    if default.exists():
        return default
    if required:
        raise FileNotFoundError(
            "Shared DSSAT template directory not found. Set DSSAT_TEMPLATE_DIR "
            "or templates.template_dir."
        )
    return None

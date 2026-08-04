"""Build focused, staged global-hemp calibration configurations."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import yaml

from dssatcalibrator.config import load_config, validate_config
from dssatcalibrator.spaces import ParameterSpace


ALL_CULTIVARS = ("IB0001", "IB0002", "IB0007", "IB0008")


def _load_theta(path: str | Path) -> dict[str, float]:
    return {str(k): float(v) for k, v in json.loads(Path(path).read_text()).items()}


def _set_start(spec: dict, name: str, value: float) -> None:
    if "__" in name:
        _, cultivar = name.rsplit("__", 1)
        spec.setdefault("start_by_cultivar", {})[cultivar] = float(value)
    else:
        spec["start"] = float(value)


def _set_bounds(spec: dict, name: str, lo: float, hi: float) -> None:
    if "__" in name:
        _, cultivar = name.rsplit("__", 1)
        spec.setdefault("min_by_cultivar", {})[cultivar] = float(lo)
        spec.setdefault("max_by_cultivar", {})[cultivar] = float(hi)
    else:
        spec["min"] = float(lo)
        spec["max"] = float(hi)


def _parameter_lookup(cfg: dict) -> tuple[ParameterSpace, dict[str, dict]]:
    space = ParameterSpace.from_config(cfg)
    return space, {spec["name"]: spec for spec in space.specs}


def _activate(cfg: dict, selected: set[str]) -> None:
    selected_by_base: dict[tuple[str, str], list[str | None]] = {}
    for group, params in cfg["parameters"].items():
        for base, spec in params.items():
            spec["active"] = False
            spec["fixed"] = True
            spec.pop("cultivars", None)
            spec.pop("fixed_cultivars", None)
            scoped = str(spec.get("scope", "")).lower() in {
                "cultivar", "cultivars", "per_cultivar", "per-cultivar",
                "cultivar_specific", "cultivar-specific",
            }
            names = [f"{base}__{c}" for c in ALL_CULTIVARS] if scoped else [base]
            hits = [name for name in names if name in selected]
            if hits:
                selected_by_base[(group, base)] = [
                    name.rsplit("__", 1)[1] if "__" in name else None for name in hits
                ]

    for (group, base), cultivars in selected_by_base.items():
        spec = cfg["parameters"][group][base]
        spec["active"] = True
        if cultivars != [None]:
            active_cultivars = [str(v) for v in cultivars]
            spec["cultivars"] = active_cultivars
            spec["fixed_cultivars"] = [
                cultivar for cultivar in ALL_CULTIVARS
                if cultivar not in active_cultivars
            ]


def _localized_bounds(
    name: str,
    center: float,
    lo: float,
    hi: float,
    step: float | None,
    stage: str,
) -> tuple[float, float]:
    futura_structural = {
        "CSDL__IB0007": (15.3, 16.9),
        "PPSEN__IB0007": (0.001, 0.15),
        "EM-FL__IB0007": (10.68, 30.68),
    }
    if stage == "phenology" and name in futura_structural:
        target_lo, target_hi = futura_structural[name]
        return min(target_lo, center), max(target_hi, center)
    if step is not None and math.isfinite(step) and step > 0:
        half = 2.0 * step
    else:
        half = (0.18 if stage == "phenology" else 0.20) * (hi - lo)
    new_lo = max(lo, center - half)
    new_hi = min(hi, center + half)
    new_lo = min(new_lo, center)
    new_hi = max(new_hi, center)
    if new_hi <= new_lo:
        epsilon = step if step is not None and step > 0 else max(abs(center) * 1e-6, 1e-6)
        return center - epsilon, center + epsilon
    return new_lo, new_hi


def build(args: argparse.Namespace) -> dict:
    cfg = load_config(args.base)
    reference = load_config(args.reference)
    effective = _load_theta(args.best)
    for override in args.override:
        effective.update(_load_theta(override))

    ref_space, ref_specs = _parameter_lookup(reference)
    _, base_specs = _parameter_lookup(cfg)
    missing = sorted(set(base_specs) - set(effective))
    if missing:
        raise ValueError(f"Effective theta is missing {len(missing)} parameters: {missing}")

    recommendations = list(__import__("csv").DictReader(
        Path(args.recommendations).open(encoding="utf-8")
    ))
    selected = {
        row["parameter"] for row in recommendations
        if row["stage"] == args.stage
    }
    if args.stage == "final":
        selected = set(base_specs)

    for name, expanded in base_specs.items():
        base = expanded.get("base_name", name)
        spec = cfg["parameters"][expanded["group"]][base]
        _set_start(spec, name, effective[name])
        if "__" in name:
            _, cultivar = name.rsplit("__", 1)
            scoped_lo = float((spec.get("min_by_cultivar") or {}).get(cultivar, spec["min"]))
            scoped_hi = float((spec.get("max_by_cultivar") or {}).get(cultivar, spec["max"]))
            _set_bounds(
                spec, name,
                min(scoped_lo, effective[name]),
                max(scoped_hi, effective[name]),
            )
        else:
            spec["min"] = min(float(spec["min"]), effective[name])
            spec["max"] = max(float(spec["max"]), effective[name])

    _activate(cfg, selected)
    for name in selected:
        expanded = base_specs[name]
        ref = ref_specs[name]
        ref_index = ref_space.names.index(name)
        step = float(ref_space.step[ref_index])
        lo, hi = _localized_bounds(
            name, effective[name], float(ref["min"]), float(ref["max"]),
            step if math.isfinite(step) else None, args.stage,
        )
        spec = cfg["parameters"][expanded["group"]][expanded.get("base_name", name)]
        _set_bounds(spec, name, lo, hi)

    cfg["calibrator"]["name"] = args.name
    cfg["calibrator"]["workdir"] = f"results/_workdir_{args.name}"
    cfg["calibrator"]["cache_spawns"] = False
    cfg["calibrator"]["cache_evaluations"] = False
    cfg["calibrator"]["keep_run_dirs"] = False
    cfg["calibrator"]["num_cores"] = 0
    cfg["method"] = {
        "preset": "custom",
        "sample": {
            "engine": "lhs",
            "n": max(0, int(args.samples) - 1),
            "include_start": True,
        },
        "sensitivity": {"engine": "none", "active": False},
        "select": {"engine": "none", "active": False},
        "surrogate": {"engine": "none", "active": False},
        "bayesian": {"engine": "glue"},
        "optimizer": {"engine": "none"},
        "validation": {"scheme": "none"},
        "multiobjective": {"engine": "none"},
    }
    if args.stage == "phenology":
        cfg["engine"]["timeseries_outputs"] = {}
        cfg["engine"]["scalar_outputs"] = {"anthesis": "ADAP"}
        cfg["objective"]["weights"] = {"anthesis": 1.0}
        cfg["objective"]["error_model"] = {
            "anthesis": {"type": "absolute", "value": 1.0}
        }
        cfg["objective"]["max_bias_penalty"] = {
            "variable": "anthesis", "lambda": 0.5, "sigma": 1.0,
            "power": 1.0, "tolerance": 0.0,
        }

    cfg.pop("_config_path", None)
    validate_config(cfg)
    active = ParameterSpace.from_config(cfg)
    if set(active.names) != selected:
        raise RuntimeError(
            f"Active-space mismatch: expected {sorted(selected)}, got {active.names}"
        )
    return cfg


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base")
    parser.add_argument("reference")
    parser.add_argument("best")
    parser.add_argument("recommendations")
    parser.add_argument("output")
    parser.add_argument("--stage", choices=("phenology", "growth", "final"), required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--samples", type=int, required=True)
    parser.add_argument("--override", action="append", default=[])
    args = parser.parse_args()
    cfg = build(args)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    space = ParameterSpace.from_config(cfg)
    print(f"Wrote {args.stage} config with {space.ndim} active parameters -> {output}")


if __name__ == "__main__":
    main()

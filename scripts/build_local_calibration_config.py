"""Build a bounded local LHS or DE config around a completed best theta."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from dssatcalibrator.config import load_config
from dssatcalibrator.spaces import ParameterSpace


def _set_scoped(spec: dict, field: str, value: float, cultivar: str | None) -> None:
    if cultivar is None:
        spec[field] = float(value)
        return
    key = f"{field}_by_cultivar"
    spec.setdefault(key, {})[cultivar] = float(value)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("best_json")
    ap.add_argument("output")
    ap.add_argument("--name", required=True)
    ap.add_argument("--fraction", type=float, default=0.15,
                    help="local half-width as a fraction of each original full range")
    ap.add_argument("--engine", choices=("lhs", "de"), default="de")
    ap.add_argument("--samples", type=int, default=5000)
    ap.add_argument("--de-popsize", type=int, default=3)
    ap.add_argument("--de-maxiter", type=int, default=7)
    args = ap.parse_args()

    cfg = load_config(args.base)
    with open(args.best_json, encoding="utf-8") as fh:
        best = {str(k): float(v) for k, v in json.load(fh).items()}
    space = ParameterSpace.from_config(cfg)
    unknown = sorted(set(space.names) - set(best))
    if unknown:
        raise ValueError(f"Best-theta file is missing active parameters: {unknown}")

    params = cfg["parameters"]
    frac = max(float(args.fraction), 1e-6)
    for expanded, old_lo, old_hi in zip(space.specs, space.low, space.high):
        name = expanded["name"]
        center = best[name]
        half = frac * float(old_hi - old_lo)
        lo = max(float(old_lo), center - half)
        hi = min(float(old_hi), center + half)
        if hi <= lo:
            lo, hi = float(old_lo), float(old_hi)
        base_name = expanded.get("base_name", name)
        base_spec = params[expanded["group"]][base_name]
        cultivar = expanded.get("cultivar") if expanded.get("scope") == "cultivar" else None
        _set_scoped(base_spec, "start", center, cultivar)
        _set_scoped(base_spec, "min", lo, cultivar)
        _set_scoped(base_spec, "max", hi, cultivar)

    cfg["calibrator"]["name"] = args.name
    cfg["calibrator"]["workdir"] = f"results/_workdir_{args.name}"
    cfg["calibrator"]["cache_spawns"] = False
    cfg["calibrator"]["cache_evaluations"] = False
    method = cfg.setdefault("method", {})
    method["preset"] = "custom"
    method["sensitivity"] = {"engine": "none", "active": False}
    method["select"] = {"engine": "none", "active": False}
    method["surrogate"] = {"engine": "none", "active": False}
    method["multiobjective"] = {"engine": "none"}
    if args.engine == "de":
        method["sample"] = {"engine": "none", "active": False}
        method["bayesian"] = {"engine": "none"}
        method["optimizer"] = {
            "engine": "diffevo",
            "popsize": int(args.de_popsize),
            "maxiter": int(args.de_maxiter),
            "tol": 0.001,
            "eval_batch_size": 16,
        }
    else:
        method["sample"] = {
            "engine": "lhs", "n": int(args.samples) - 1, "include_start": True,
        }
        method["bayesian"] = {"engine": "glue"}
        method["optimizer"] = {"engine": "none"}

    cfg.pop("_config_path", None)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    print(f"Wrote {args.engine.upper()} config with {space.ndim} parameters -> {output}")


if __name__ == "__main__":
    main()

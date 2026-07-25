"""Build the joint anthesis plus growth/yield LHS config from a phenology best."""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path

import yaml

from dssatcalibrator.config import load_config
from dssatcalibrator.spaces import ParameterSpace
from dssatcalibrator.writers import read_cultivar_values, read_ecotype_values


CUL_GROWTH = (
    "FL-SH", "FL-SD", "SD-PM", "FL-LF", "LFMAX", "SLAVR", "SIZLF",
    "XFRT", "WTPSD", "SFDUR", "SDPDV", "PODUR", "THRSH",
)
ECO_GROWTH = ("PM09", "LNGSH", "FL-VS", "TRIFL", "RWDTH", "RHGHT", "R1PPO")
LATE_SPE = ("TB3", "TO1_3", "TO2_3", "TM3")


def _set_scoped(spec: dict, field: str, value: float, cultivar: str | None) -> None:
    if cultivar is None:
        spec[field] = float(value)
    else:
        spec.setdefault(f"{field}_by_cultivar", {})[cultivar] = float(value)


def _localize_existing(cfg: dict, best: dict[str, float], fraction: float) -> None:
    space = ParameterSpace.from_config(cfg)
    for expanded, old_lo, old_hi in zip(space.specs, space.low, space.high):
        name = expanded["name"]
        if name not in best:
            continue
        center = float(best[name])
        half = fraction * float(old_hi - old_lo)
        lo = max(float(old_lo), center - half)
        hi = min(float(old_hi), center + half)
        base_name = expanded.get("base_name", name)
        spec = cfg["parameters"][expanded["group"]][base_name]
        cultivar = expanded.get("cultivar") if expanded.get("scope") == "cultivar" else None
        _set_scoped(spec, "start", center, cultivar)
        _set_scoped(spec, "min", lo, cultivar)
        _set_scoped(spec, "max", hi, cultivar)


def _activate_growth_specs(cfg: dict, donor: dict, fraction: float) -> None:
    crop = cfg["crops"][0]
    genotype = Path(cfg["calibrator"]["dssat_dir"]) / "Genotype"
    stem = crop["genotype_stem"]
    cul_path = genotype / f"{stem}.CUL"
    eco_path = genotype / f"{stem}.ECO"
    cultivars = [str(v) for v in crop["calibration_cultivars"]]
    eco_map = {str(k): str(v) for k, v in crop["cultivar_ecotypes"].items()}

    for group, names in (("genetic_cultivar", CUL_GROWTH),
                         ("genetic_ecotype", ECO_GROWTH)):
        cfg["parameters"].setdefault(group, {})
        for name in names:
            if name not in cfg["parameters"][group]:
                cfg["parameters"][group][name] = deepcopy(donor["parameters"][group][name])
            spec = cfg["parameters"][group][name]
            spec["active"] = True
            spec["scope"] = "cultivar"
            donor_spec = donor["parameters"][group][name]
            global_lo = float(donor_spec["min"])
            global_hi = float(donor_spec["max"])
            for cultivar in cultivars:
                values = (read_cultivar_values(cul_path, cultivar) if group == "genetic_cultivar"
                          else read_ecotype_values(eco_path, eco_map[cultivar]))
                center = float(values[name])
                half = fraction * (global_hi - global_lo)
                lo = max(global_lo, center - half)
                hi = min(global_hi, center + half)
                _set_scoped(spec, "start", center, cultivar)
                _set_scoped(spec, "min", lo, cultivar)
                _set_scoped(spec, "max", hi, cultivar)

    for name in LATE_SPE:
        cfg["parameters"]["genetic_species"][name]["active"] = True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("phenology_config")
    ap.add_argument("phenology_best_json")
    ap.add_argument("output")
    ap.add_argument("--name", required=True)
    ap.add_argument("--samples", type=int, default=5000)
    ap.add_argument("--phenology-fraction", type=float, default=0.10)
    ap.add_argument("--growth-fraction", type=float, default=0.25)
    ap.add_argument("--donor", default="config_hemp.yaml")
    args = ap.parse_args()

    cfg = load_config(args.phenology_config)
    donor = load_config(args.donor)
    with open(args.phenology_best_json, encoding="utf-8") as fh:
        best = {str(k): float(v) for k, v in json.load(fh).items()}

    _localize_existing(cfg, best, max(float(args.phenology_fraction), 1e-6))
    _activate_growth_specs(cfg, donor, max(float(args.growth_fraction), 1e-6))

    cfg["calibrator"]["name"] = args.name
    cfg["calibrator"]["workdir"] = f"results/_workdir_{args.name}"
    cfg["calibrator"]["cache_spawns"] = False
    cfg["calibrator"]["cache_evaluations"] = False
    cfg.pop("filex_overrides", None)  # Growth/yield runs need the experiment's water and N settings.

    cfg["engine"] = {
        "run_mode": "experiment",
        "timeseries_outputs": {
            # FLWAD is deliberately excluded: CROPGRO-Hemp does not emit it,
            # and UKAB2101 values are duplicated from Italian LWAD records.
            "biomass": "CWAD", "leaf": "LWAD", "stem": "SWAD", "root": "RWAD",
            "grain": "GWAD", "LAI": "LAID", "height": "CHTD",
            "width": "CWID", "leaf_number": "L#SD",
        },
        "scalar_outputs": {"anthesis": "ADAP"},
    }
    cfg["objective"] = {
        "weighting": "unified",
        "score_metric": "rmse",
        "obs_autocorr": True,
        "ignore_zero_observations": ["grain"],
        "weights": {
            "anthesis": 3.0, "biomass": 1.0, "grain": 1.0, "LAI": 0.75,
            "leaf": 0.5, "stem": 0.5, "root": 0.5, "height": 0.5,
            "width": 0.25, "leaf_number": 0.25,
        },
        "error_model": {
            "anthesis": {"type": "absolute", "value": 3.0},
            "biomass": {"type": "relative", "value": 0.20},
            "grain": {"type": "relative", "value": 0.20},
            "LAI": {"type": "relative", "value": 0.20},
            "leaf": {"type": "relative", "value": 0.20},
            "stem": {"type": "relative", "value": 0.20},
            "root": {"type": "relative", "value": 0.25},
            "height": {"type": "relative", "value": 0.15},
            "width": {"type": "relative", "value": 0.15},
            "leaf_number": {"type": "relative", "value": 0.15},
        },
        "max_bias_penalty": {
            "variable": "anthesis", "lambda": 0.5, "sigma": 3.0,
            "power": 1.0, "tolerance": 0.0,
        },
    }
    cfg["method"] = {
        "preset": "custom",
        "sample": {"engine": "lhs", "n": int(args.samples) - 1, "include_start": True},
        "sensitivity": {"engine": "none", "active": False},
        "select": {"engine": "none", "active": False},
        "surrogate": {"engine": "none", "active": False},
        "bayesian": {"engine": "glue"},
        "optimizer": {"engine": "none"},
        "validation": {"scheme": "none"},
        "multiobjective": {"engine": "none"},
    }
    cfg.setdefault("diagnostics", {})["keep_all_objectives"] = False
    cfg.pop("_config_path", None)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    dims = ParameterSpace.from_config(cfg).ndim
    print(f"Wrote joint {dims}-parameter, {args.samples}-sample config -> {output}")


if __name__ == "__main__":
    main()

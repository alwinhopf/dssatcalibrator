"""Build the second-generation global hemp growth calibration config."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from build_focused_global_config import build as build_focused
from dssatcalibrator.config import validate_config


WEIGHTS = {
    "anthesis": 3.0,
    "biomass": 1.5,
    "grain": 1.5,
    "LAI": 1.0,
    "leaf": 1.0,
    "stem": 1.0,
    "root": 1.0,
    "height": 0.5,
    "width": 0.5,
    "leaf_number": 0.5,
}


def _potential_overrides() -> dict:
    return {
        "all": [
            {
                "section": "SIMULATION CONTROLS",
                "header_prefix": "@N OPTIONS",
                "field": "WATER",
                "value": "N",
                "type": "code",
                "required": True,
            },
            {
                "section": "SIMULATION CONTROLS",
                "header_prefix": "@N OPTIONS",
                "field": "NITRO",
                "value": "N",
                "type": "code",
                "required": True,
            },
        ]
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base")
    parser.add_argument("reference")
    parser.add_argument("best")
    parser.add_argument("output")
    parser.add_argument("--mode", choices=("diagnostic", "search"), required=True)
    parser.add_argument("--stress", choices=("source", "potential"), default="source")
    parser.add_argument("--name", required=True)
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--override", action="append", default=[])
    args = parser.parse_args()

    recommendations = (
        Path(__file__).resolve().parents[1]
        / "calibration_global_hemp"
        / "targeted_growth_v2_parameters.csv"
    )
    focused_args = argparse.Namespace(
        base=args.base,
        reference=args.reference,
        best=args.best,
        recommendations=str(recommendations),
        output=args.output,
        stage="final" if args.mode == "diagnostic" else "growth",
        name=args.name,
        samples=args.samples,
        override=args.override,
    )
    cfg = build_focused(focused_args)
    cfg["objective"]["weighting"] = "site_variable"
    cfg["objective"]["weights"] = dict(WEIGHTS)
    cfg["objective"]["obs_autocorr"] = True
    cfg["objective"]["missing_simulation_policy"] = "reject"
    cfg["objective"]["timeseries_after_simulation_policy"] = "carry_forward"
    cfg["objective"]["max_bias_penalty"] = {
        "variable": "anthesis",
        "lambda": 4.0,
        "sigma": 1.0,
        "power": 2.0,
        "tolerance": 7.0,
    }
    if args.stress == "potential":
        cfg["filex_overrides"] = _potential_overrides()
    else:
        cfg.pop("filex_overrides", None)
    cfg.setdefault("metadata", {})["growth_v2"] = {
        "stress_mode": args.stress,
        "objective": "equal experiment-variable weighting",
        "anthesis_hard_tolerance_days": 7.0,
    }
    validate_config(cfg)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    print(f"Wrote {args.mode} config ({args.stress}) -> {output}")


if __name__ == "__main__":
    main()

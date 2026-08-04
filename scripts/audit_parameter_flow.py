"""Verify that every active genotype parameter resolves to intended DSSAT inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from dssatcalibrator.config import load_config, resolve_dssat_paths
from dssatcalibrator.spawn import _spec_applies, parse_cultivars
from dssatcalibrator.spaces import ParameterSpace
from dssatcalibrator.writers import (
    cultivar_field_map,
    ecotype_field_map,
    read_cultivar_values,
    read_ecotype_values,
)


def _species_support(path: Path, spec: dict) -> tuple[bool, str]:
    key = str(spec.get("spe_key", spec.get("base_name", spec["name"])))
    index = int(spec.get("spe_index", spec.get("token_index", 0)))
    matches = [
        line for line in path.read_text(errors="replace").splitlines()
        if key in line and line.strip() and not line.lstrip().startswith(("*", "@", "!"))
    ]
    return len(matches) == 1, f"{key}[{index}]"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("output")
    args = parser.parse_args()

    cfg = load_config(args.config)
    space = ParameterSpace.from_config(cfg)
    crop = cfg["crops"][0]
    genotype = resolve_dssat_paths(cfg)["genotype"]
    hemp = Path(cfg["source"]["hemp_dir"])
    stem = crop["genotype_stem"]
    cultivar_ecotypes = {
        str(key): str(value)
        for key, value in (crop.get("cultivar_ecotypes") or {}).items()
    }
    experiment_cultivars = {
        exp: parse_cultivars(hemp / f"{exp}.{crop['filex_ext']}")
        for exp in cfg["experiments"]
    }
    cul_path = genotype / f"{stem}.CUL"
    eco_path = genotype / f"{stem}.ECO"
    spe_path = genotype / f"{stem}.SPE"
    cul_fields = cultivar_field_map(cul_path)
    eco_fields = ecotype_field_map(eco_path)

    rows = []
    for spec in space.specs:
        name = spec["name"]
        base = str(spec.get("base_name", name))
        group = str(spec["group"])
        cultivar = str(spec.get("cultivar", ""))
        applies = [
            exp for exp, cultivars in experiment_cultivars.items()
            if _spec_applies(spec, exp, cultivars)
        ]
        supported = False
        target_row = ""
        current_value = None
        if group == "genetic_cultivar":
            target_row = cultivar or str(crop.get("cultivar_anchor", ""))
            supported = base in cul_fields
            if supported and target_row:
                current_value = read_cultivar_values(cul_path, target_row).get(base)
        elif group == "genetic_ecotype":
            target_row = cultivar_ecotypes.get(
                cultivar, str(crop.get("ecotype", ""))
            )
            supported = bool(target_row) and base in eco_fields
            if supported:
                current_value = read_ecotype_values(eco_path, target_row).get(base)
        elif group == "genetic_species":
            supported, target_row = _species_support(spe_path, spec)
        rows.append({
            "parameter": name,
            "group": group,
            "scope": spec.get("scope", "global"),
            "cultivar": cultivar,
            "target_file": f"{stem}." + {
                "genetic_cultivar": "CUL",
                "genetic_ecotype": "ECO",
                "genetic_species": "SPE",
            }.get(group, ""),
            "target_row_or_key": target_row,
            "coefficient": base,
            "source_support": supported,
            "source_value": current_value,
            "n_applicable_experiments": len(applies),
            "applicable_experiments": ";".join(applies),
        })

    audit = pd.DataFrame(rows)
    bad = audit[
        (~audit["source_support"])
        | (audit["n_applicable_experiments"] == 0)
        | audit["parameter"].duplicated(keep=False)
    ]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output, index=False)
    if not bad.empty:
        raise RuntimeError(
            f"{len(bad)} parameter-flow checks failed; inspect {output}."
        )
    print(
        f"Verified {len(audit)} unique active parameters; all resolve to a "
        "supported source field and at least one configured experiment."
    )


if __name__ == "__main__":
    main()

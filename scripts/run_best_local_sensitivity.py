"""Audit local parameter sensitivity around a saved DSSAT calibration result.

The design contains the saved best point plus one lower and one upper
perturbation per active parameter. Perturbations are scaled to a reference
configuration's bounds so a narrow optimizer-polish configuration can still be
tested over a scientifically meaningful neighborhood.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path
import re
import shutil
import tempfile

import numpy as np
import pandas as pd

from dssatcalibrator.config import load_config, resolve_dssat_paths
from dssatcalibrator.objective import _group_loss, _weighted_loss
from dssatcalibrator.orchestrator import evaluate_design
from dssatcalibrator.spaces import ParameterSpace
from dssatcalibrator.writers import (
    edit_cultivar,
    edit_ecotype,
    edit_species,
    read_cultivar_values,
    read_ecotype_values,
)


def build_oat_design(
    space: ParameterSpace,
    best: dict[str, float],
    reference_space: ParameterSpace,
    fraction: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the best point and symmetric one-at-a-time perturbations."""
    if space.names != reference_space.names:
        missing = sorted(set(space.names) - set(reference_space.names))
        extra = sorted(set(reference_space.names) - set(space.names))
        raise ValueError(
            "Calibration and reference parameter spaces differ; "
            f"missing from reference={missing}, extra in reference={extra}"
        )
    missing = [name for name in space.names if name not in best]
    if missing:
        raise ValueError(f"Best-theta file is missing active parameters: {missing}")
    if not 0 < fraction <= 0.5:
        raise ValueError("fraction must be greater than zero and at most 0.5")

    center = np.array([float(best[name]) for name in space.names], dtype=float)
    vectors = [center.copy()]
    metadata = [{
        "sample_id": 0,
        "candidate": "baseline",
        "parameter": "",
        "direction": "baseline",
        "requested_value": np.nan,
        "requested_delta": 0.0,
        "normalized_delta": 0.0,
    }]
    sample_id = 1
    for i, name in enumerate(space.names):
        ref_lo = float(reference_space.low[i])
        ref_hi = float(reference_space.high[i])
        width = ref_hi - ref_lo
        delta = fraction * width
        low = max(float(space.low[i]), ref_lo, center[i] - delta)
        high = min(float(space.high[i]), ref_hi, center[i] + delta)
        step = (
            float(space.step[i])
            if space.step is not None and np.isfinite(space.step[i]) and space.step[i] > 0
            else np.nan
        )
        if np.isfinite(step):
            effective_center = float(space.to_theta(center)[name])
            low = max(float(space.low[i]), effective_center - step)
            high = min(float(space.high[i]), effective_center + step)
        if not low < center[i]:
            low = max(float(space.low[i]), center[i] - fraction * float(space.high[i] - space.low[i]))
        if not high > center[i]:
            high = min(float(space.high[i]), center[i] + fraction * float(space.high[i] - space.low[i]))
        for direction, value in (("lower", low), ("upper", high)):
            vector = center.copy()
            vector[i] = value
            vectors.append(vector)
            metadata.append({
                "sample_id": sample_id,
                "candidate": f"{name}:{direction}",
                "parameter": name,
                "direction": direction,
                "requested_value": value,
                "requested_delta": value - center[i],
                "normalized_delta": (value - center[i]) / width if width > 0 else np.nan,
            })
            sample_id += 1
    effective_vectors = [
        [theta[name] for name in space.names]
        for theta in (space.to_theta(vector) for vector in vectors)
    ]
    return (
        pd.DataFrame(effective_vectors, columns=space.names),
        pd.DataFrame(metadata).set_index("sample_id", drop=False),
    )


def _objective_components(result, cfg: dict) -> dict[str, float]:
    resid = result.residuals.copy()
    if resid.empty:
        return {}
    resid["_loss"] = _weighted_loss(resid, cfg)
    weights = (cfg.get("objective", {}) or {}).get("weights", {}) or {}
    return {
        str(variable): float(weights.get(variable, 1.0)) * _group_loss(group, cfg)
        for variable, group in resid.groupby("user_var")
    }


def _aligned_max_sim_change(baseline: pd.DataFrame, candidate: pd.DataFrame) -> float:
    keys = ["exp_id", "treatment", "user_var", "dssat", "kind", "date"]
    left = baseline[keys + ["sim"]].rename(columns={"sim": "sim_base"})
    right = candidate[keys + ["sim"]].rename(columns={"sim": "sim_candidate"})
    joined = left.merge(right, on=keys, how="outer")
    if joined[["sim_base", "sim_candidate"]].isna().any(axis=None):
        return float("inf")
    return float(np.max(np.abs(joined["sim_candidate"] - joined["sim_base"])))


def _species_value(path: Path, key: str, index: int) -> float:
    line = next(
        line for line in path.read_text(errors="replace").splitlines()
        if key.lower() in line.lower()
    )
    numbers = [
        float(value)
        for value in re.findall(
            r"(?<![A-Za-z])[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?",
            line,
        )
    ]
    return numbers[index]


def _writer_effective_value(
    cfg: dict,
    spec: dict,
    requested: float,
    temp_dir: Path,
) -> float:
    crop = cfg["crops"][0]
    genotype_dir = resolve_dssat_paths(cfg)["genotype"]
    stem = crop["genotype_stem"]
    group = spec["group"]
    base_name = spec.get("base_name", spec["name"])
    cultivar = spec.get("cultivar")

    if group == "genetic_cultivar":
        path = temp_dir / f"{stem}.CUL"
        shutil.copy(genotype_dir / path.name, path)
        edit_cultivar(path, cultivar, {base_name: requested})
        return float(read_cultivar_values(path, cultivar)[base_name])
    if group == "genetic_ecotype":
        path = temp_dir / f"{stem}.ECO"
        shutil.copy(genotype_dir / path.name, path)
        ecotype = str(crop["cultivar_ecotypes"][cultivar])
        edit_ecotype(path, ecotype, {base_name: requested})
        return float(read_ecotype_values(path, ecotype)[base_name])
    if group == "genetic_species":
        path = temp_dir / f"{stem}.SPE"
        shutil.copy(genotype_dir / path.name, path)
        key = str(spec.get("spe_key", base_name))
        index = int(spec.get("spe_index", spec.get("token_index", 0)))
        edit_species(path, {key: {"value": requested, "index": index}})
        return _species_value(path, key, index)
    return float(requested)


def writer_resolution_table(
    cfg: dict,
    space: ParameterSpace,
    best: dict[str, float],
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Report requested versus file-effective values for each OAT triplet."""
    rows = []
    with tempfile.TemporaryDirectory(prefix="dssatcal_writer_audit_") as tmp:
        temp_dir = Path(tmp)
        for spec in space.specs:
            name = spec["name"]
            subset = metadata[metadata["parameter"] == name]
            requested = {
                "baseline": float(best[name]),
                "lower": float(subset.loc[subset["direction"] == "lower", "requested_value"].iloc[0]),
                "upper": float(subset.loc[subset["direction"] == "upper", "requested_value"].iloc[0]),
            }
            effective = {
                direction: _writer_effective_value(cfg, spec, value, temp_dir)
                for direction, value in requested.items()
            }
            levels = len(set(effective.values()))
            rows.append({
                "parameter": name,
                "group": spec["group"],
                "scope": spec.get("scope", "global"),
                "cultivar": spec.get("cultivar", ""),
                "requested_baseline": requested["baseline"],
                "effective_baseline": effective["baseline"],
                "baseline_write_drift": effective["baseline"] - requested["baseline"],
                "requested_lower": requested["lower"],
                "effective_lower": effective["lower"],
                "requested_upper": requested["upper"],
                "effective_upper": effective["upper"],
                "effective_levels": levels,
                "writer_stuck": levels == 1,
            })
    return pd.DataFrame(rows)


def summarize(
    cfg: dict,
    space: ParameterSpace,
    metadata: pd.DataFrame,
    design: pd.DataFrame,
    objective_results: dict,
    writer_table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    baseline = objective_results[0]
    baseline_components = _objective_components(baseline, cfg)
    missing_penalty = float((cfg.get("objective", {}) or {}).get("missing_simulation_penalty", 1000.0))
    design_by_id = design.set_index("sample_id")
    rows = []
    component_rows = []
    anthesis_rows = []

    for i, spec in enumerate(space.specs):
        name = spec["name"]
        lower_id = int(metadata.index[(metadata["parameter"] == name) & (metadata["direction"] == "lower")][0])
        upper_id = int(metadata.index[(metadata["parameter"] == name) & (metadata["direction"] == "upper")][0])
        lower = objective_results[lower_id]
        upper = objective_results[upper_id]
        requested_low = float(metadata.loc[lower_id, "requested_value"])
        requested_high = float(metadata.loc[upper_id, "requested_value"])
        low_value = float(design_by_id.loc[lower_id, name])
        high_value = float(design_by_id.loc[upper_id, name])
        width = high_value - low_value
        low_change = float(lower.score - baseline.score)
        high_change = float(upper.score - baseline.score)
        max_sim_change = max(
            _aligned_max_sim_change(baseline.residuals, lower.residuals),
            _aligned_max_sim_change(baseline.residuals, upper.residuals),
        )
        writer = writer_table.set_index("parameter").loc[name]
        rows.append({
            "parameter": name,
            "group": spec["group"],
            "scope": spec.get("scope", "global"),
            "cultivar": spec.get("cultivar", ""),
            "baseline_value": float(metadata.iloc[0].get(name, np.nan))
                if name in metadata.columns else float(design_by_id.loc[0, name]),
            "requested_lower_value": requested_low,
            "requested_upper_value": requested_high,
            "lower_value": low_value,
            "upper_value": high_value,
            "baseline_score": float(baseline.score),
            "lower_score": float(lower.score),
            "upper_score": float(upper.score),
            "lower_feasible": bool(np.isfinite(lower.score)),
            "upper_feasible": bool(np.isfinite(upper.score)),
            "feasibility_boundary": bool(
                not (np.isfinite(lower.score) and np.isfinite(upper.score))
            ),
            "lower_score_change": low_change,
            "upper_score_change": high_change,
            "max_abs_score_change": max(abs(low_change), abs(high_change)),
            "central_score_slope": (float(upper.score) - float(lower.score)) / width if width > 0 else np.nan,
            "max_abs_sim_change": max_sim_change,
            "locally_inactive": bool(np.isclose(max_sim_change, 0.0, atol=1e-12, rtol=0.0)),
            "writer_stuck": bool(writer["writer_stuck"]),
            "effective_levels": int(writer["effective_levels"]),
            "lower_n_obs": len(lower.residuals),
            "upper_n_obs": len(upper.residuals),
            "lower_missing_penalties": int(np.isclose(np.abs(lower.residuals["resid"]), missing_penalty).sum()),
            "upper_missing_penalties": int(np.isclose(np.abs(upper.residuals["resid"]), missing_penalty).sum()),
        })

        for variable in sorted(
            set(baseline_components)
            | set(_objective_components(lower, cfg))
            | set(_objective_components(upper, cfg))
        ):
            base_component = baseline_components.get(variable, np.nan)
            low_component = _objective_components(lower, cfg).get(variable, np.nan)
            high_component = _objective_components(upper, cfg).get(variable, np.nan)
            component_rows.append({
                "parameter": name,
                "variable": variable,
                "baseline_component": base_component,
                "lower_component": low_component,
                "upper_component": high_component,
                "max_abs_component_change": np.nanmax(
                    np.abs([low_component - base_component, high_component - base_component])
                ),
            })

        for direction, result in (("lower", lower), ("upper", upper)):
            anthesis = result.residuals[result.residuals["user_var"] == "anthesis"]
            for _, row in anthesis.iterrows():
                anthesis_rows.append({
                    "parameter": name,
                    "direction": direction,
                    "exp_id": row["exp_id"],
                    "treatment": int(row["treatment"]),
                    "observed_anthesis_dap": float(row["obs"]),
                    "simulated_anthesis_dap": float(row["sim"]),
                    "bias": float(row["resid"]),
                })

    ranking = pd.DataFrame(rows).sort_values(
        ["feasibility_boundary", "max_abs_score_change", "max_abs_sim_change"],
        ascending=False,
    ).reset_index(drop=True)
    components = pd.DataFrame(component_rows).sort_values(
        ["variable", "max_abs_component_change"], ascending=[True, False]
    ).reset_index(drop=True)
    anthesis = pd.DataFrame(anthesis_rows)
    return ranking, components, anthesis


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("best_theta")
    parser.add_argument("reference_config")
    parser.add_argument("output_dir")
    parser.add_argument("--fraction", type=float, default=0.05)
    parser.add_argument("--cores", type=int)
    args = parser.parse_args()

    cfg = load_config(args.config)
    reference_cfg = load_config(args.reference_config)
    space = ParameterSpace.from_config(cfg)
    reference_space = ParameterSpace.from_config(reference_cfg)
    best = {
        str(key): float(value)
        for key, value in json.loads(Path(args.best_theta).read_text(encoding="utf-8")).items()
    }
    samples, metadata = build_oat_design(space, best, reference_space, args.fraction)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_cfg = deepcopy(cfg)
    run_cfg["calibrator"]["name"] = output_dir.name
    run_cfg["calibrator"]["workdir"] = str(output_dir.parent / "_workdir_local_sensitivity")
    run_cfg["calibrator"]["cache_spawns"] = False
    run_cfg["calibrator"]["cache_evaluations"] = False
    run_cfg["calibrator"]["keep_run_dirs"] = False
    if args.cores is not None:
        run_cfg["calibrator"]["num_cores"] = int(args.cores)
    run_cfg.setdefault("diagnostics", {})["keep_all_objectives"] = True

    writer_table = writer_resolution_table(run_cfg, space, best, metadata)
    design, objective_results, _ = evaluate_design(run_cfg, samples, progress=True)
    manifest = design.attrs.get("spawn_manifest", pd.DataFrame())
    metadata.reset_index(drop=True).to_csv(output_dir / "candidate_metadata.csv", index=False)
    design.to_csv(output_dir / "design.csv", index=False)
    manifest.to_csv(output_dir / "spawn_manifest.csv", index=False)
    writer_table.to_csv(output_dir / "writer_resolution.csv", index=False)
    failures = manifest[~manifest["status"].isin(["success", "cached", "evaluation_cached"])]
    if not failures.empty:
        failures.to_csv(output_dir / "failed_spawns.csv", index=False)
        raise RuntimeError(
            f"Local sensitivity aborted: {len(failures)} DSSAT spawns failed. "
            f"See {output_dir / 'failed_spawns.csv'}."
        )

    ranking, components, anthesis = summarize(
        run_cfg, space, metadata, design, objective_results, writer_table
    )
    ranking.to_csv(output_dir / "local_sensitivity.csv", index=False)
    components.to_csv(output_dir / "local_sensitivity_by_variable.csv", index=False)
    anthesis.to_csv(output_dir / "anthesis_perturbations.csv", index=False)

    print(f"Evaluated {len(samples)} candidates across {len(run_cfg['experiments'])} FileX runs.")
    print(f"Spawn failures: {len(failures)}")
    print(f"Locally inactive parameters: {int(ranking['locally_inactive'].sum())}")
    print(f"Writer-stuck perturbation triplets: {int(writer_table['writer_stuck'].sum())}")
    print(ranking.head(25).to_string(index=False))


if __name__ == "__main__":
    main()

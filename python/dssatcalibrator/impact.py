"""One-at-a-time DSSAT parameter impact atlas.

The atlas is a real-DSSAT exercise: for each candidate parameter it runs a
baseline plus low/high perturbations, collects broad ``*.OUT`` tables, scores
against observations, and writes human/agent-readable artifacts.
"""

from __future__ import annotations

import json
import math
import re
import shutil
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import dssat_io
from . import objective as obj
from .config import all_parameters, crop_for, resolve_dssat_paths, resolve_exe
from .observations import Observations
from .runner import resolve_cores, run_many
from .spawn import parse_treatments
from .writers import read_cul_calibration_bounds, read_cultivar_values, read_ecotype_values


DEFAULT_GROUPS = [
    "genetic_cultivar",
    "genetic_ecotype",
    "genetic_species",
    "management",
    "initial_conditions",
    "soil",
    "weather",
]


@dataclass
class AtlasResult:
    output_dir: Path
    run_manifest: pd.DataFrame
    file_manifest: pd.DataFrame
    output_effects: pd.DataFrame
    output_long: pd.DataFrame
    capability_map: pd.DataFrame
    parameter_catalog: pd.DataFrame
    score_effects: pd.DataFrame
    output_impact_summary: pd.DataFrame
    parameter_impact_summary: pd.DataFrame


def _safe_id(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")


def _is_num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(float(x))


def _parameter_record(group: str, name: str, spec: dict) -> dict:
    rec = {"group": group, "name": name}
    rec.update(spec)
    return rec


def _relative_bounds(value: float, *, rel_width: float = 0.10) -> tuple[float, float]:
    width = max(abs(float(value)) * float(rel_width), 1e-6 if abs(float(value)) < 0.01 else float(rel_width))
    lo, hi = float(value) - width, float(value) + width
    if lo == hi:
        hi = lo + 1e-6
    return lo, hi


def _genotype_paths(cfg: dict) -> tuple[dict, Path, str]:
    crop = crop_for(cfg, (cfg.get("crops") or [{}])[0].get("code", "HM"))
    stem = crop.get("genotype_stem") or crop.get("model")
    return crop, resolve_dssat_paths(cfg)["genotype"], stem


def _candidate_parameters(cfg: dict, groups: list[str], *, active_only: bool) -> list[dict]:
    specs = []
    for spec in all_parameters(cfg):
        if spec.get("group") not in groups:
            continue
        if active_only and not spec.get("active", False):
            continue
        if not all(_is_num(spec.get(k)) for k in ("min", "max", "start")):
            continue
        if float(spec["min"]) >= float(spec["max"]):
            continue
        specs.append(dict(spec))
    return specs


_NUM_RE = re.compile(r"-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def discover_species_parameters(cfg: dict, *, rel_width: float = 0.10) -> list[dict]:
    """Discover editable species tokens from the crop ``.SPE`` file.

    DSSAT ``.SPE`` rows commonly contain several numeric values followed by a
    comma-separated label list. The returned specs use ``spe_key`` plus
    ``spe_index`` so the writer edits exactly one token.
    """
    _crop, genotype_dir, stem = _genotype_paths(cfg)
    if not stem:
        return []
    spe = genotype_dir / f"{stem}.SPE"
    if not spe.exists():
        return []

    out: list[dict] = []
    seen: set[str] = set()
    for line_no, line in enumerate(spe.read_text(errors="replace").splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped[0] in "*@!":
            continue
        nums = list(_NUM_RE.finditer(line))
        if not nums:
            continue
        label_block = line[nums[-1].end():].split("!", 1)[0].strip()
        if not label_block:
            continue
        labels = [x.strip() for x in label_block.split(",") if x.strip()]
        for idx, m in enumerate(nums):
            if idx < len(labels):
                label = labels[idx]
            else:
                label = f"{_safe_id(label_block)}_{idx + 1}"
            try:
                value = float(m.group(0))
            except ValueError:
                continue
            lo, hi = _relative_bounds(value, rel_width=rel_width)
            name = _safe_id(label)
            if name in seen:
                name = f"{name}_{line_no}_{idx}"
            seen.add(name)
            out.append({
                "group": "genetic_species",
                "name": name,
                "active": True,
                "min": lo,
                "max": hi,
                "start": value,
                "spe_key": label_block,
                "spe_index": idx,
                "source": "discovered_spe",
                "line_no": line_no,
            })
    return out


def discover_cultivar_parameters(cfg: dict, *, rel_width: float = 0.10) -> list[dict]:
    """Discover editable cultivar coefficients from the crop ``.CUL`` row."""
    crop, genotype_dir, stem = _genotype_paths(cfg)
    anchor = crop.get("cultivar_anchor")
    if not stem or not anchor:
        return []
    cul = genotype_dir / f"{stem}.CUL"
    if not cul.exists():
        return []
    values = read_cultivar_values(cul, anchor)
    bounds = read_cul_calibration_bounds(cul)
    out = []
    for name, value in values.items():
        if not _is_num(value):
            continue
        b = bounds.get(name, {})
        lo, hi = b.get("min"), b.get("max")
        if not (_is_num(lo) and _is_num(hi) and float(lo) < float(hi) and float(lo) <= float(value) <= float(hi)):
            lo, hi = _relative_bounds(float(value), rel_width=rel_width)
        out.append({
            "group": "genetic_cultivar",
            "name": name,
            "active": True,
            "min": float(lo),
            "max": float(hi),
            "start": float(value),
            "source": "discovered_cul",
        })
    return out


def discover_ecotype_parameters(cfg: dict, *, rel_width: float = 0.10) -> list[dict]:
    """Discover editable ecotype coefficients from the crop ``.ECO`` row."""
    crop, genotype_dir, stem = _genotype_paths(cfg)
    anchor = crop.get("ecotype")
    if not stem or not anchor:
        return []
    eco = genotype_dir / f"{stem}.ECO"
    if not eco.exists():
        return []
    values = read_ecotype_values(eco, anchor)
    out = []
    for name, value in values.items():
        if not _is_num(value):
            continue
        lo, hi = _relative_bounds(float(value), rel_width=rel_width)
        out.append({
            "group": "genetic_ecotype",
            "name": name,
            "active": True,
            "min": float(lo),
            "max": float(hi),
            "start": float(value),
            "source": "discovered_eco",
        })
    return out


def _variant_specs(spec: dict, levels: list[str]) -> list[dict]:
    values = {
        "start": float(spec["start"]),
        "low": float(spec["min"]),
        "high": float(spec["max"]),
    }
    variants = []
    for level in levels:
        if level not in values:
            continue
        value = values[level]
        rec = dict(spec)
        rec["active"] = True
        rec["scope"] = "global"
        variants.append({
            "variant_kind": "parameter",
            "parameter": spec["name"],
            "group": spec["group"],
            "level": level,
            "value": value,
            "theta": {spec["name"]: value},
            "param_specs": [rec],
        })
    return variants


def _load_observations(cfg: dict, experiments: list[str], crop: dict) -> Observations:
    src_block = cfg.get("source", {}) or {}
    hemp_dir = Path(src_block.get("hemp_dir", ""))
    if "table" in src_block:
        return Observations(src_block["table"])
    if cfg.get("observation_sources"):
        return Observations.from_sources(cfg, experiments)
    src = src_block.get("observations", "dssat")
    if src == "dssat":
        return Observations.from_dssat(hemp_dir, experiments, crop_ext=crop["code"])
    return Observations.from_csv(src)


def _output_stats(long: pd.DataFrame) -> pd.DataFrame:
    if long.empty:
        return pd.DataFrame(columns=["source_file", "variable", "treatment", "stat", "value"])
    df = long.copy()
    if "treatment" not in df.columns:
        df["treatment"] = "__all__"
    else:
        df["treatment"] = df["treatment"].astype("object").where(df["treatment"].notna(), "__all__")
    sort_cols = [c for c in ("date", "dap", "das", "row_index") if c in df.columns]
    group_cols = ["source_file", "variable", "treatment"]
    rows = []
    for keys, g in df.groupby(group_cols, dropna=False):
        vals = pd.to_numeric(g["value"], errors="coerce").dropna()
        if vals.empty:
            continue
        rows.extend([
            (*keys, "min", float(vals.min())),
            (*keys, "max", float(vals.max())),
            (*keys, "mean", float(vals.mean())),
        ])
        if sort_cols:
            last = g.sort_values(sort_cols).iloc[-1]
        else:
            last = g.iloc[-1]
        rows.append((*keys, "final", float(last["value"])))
    return pd.DataFrame(rows, columns=["source_file", "variable", "treatment", "stat", "value"])


def _effects_for_run(stats: pd.DataFrame, baseline: pd.DataFrame) -> pd.DataFrame:
    keys = ["source_file", "variable", "treatment", "stat"]
    if stats.empty or baseline.empty:
        return pd.DataFrame(columns=keys + ["baseline_value", "value", "delta", "relative_delta"])
    merged = stats.merge(baseline, on=keys, how="inner", suffixes=("", "_baseline"))
    merged = merged.rename(columns={"value_baseline": "baseline_value"})
    merged["delta"] = merged["value"] - merged["baseline_value"]
    denom = merged["baseline_value"].replace(0, np.nan).abs()
    merged["relative_delta"] = merged["delta"] / denom
    return merged[keys + ["baseline_value", "value", "delta", "relative_delta"]]


def summarize_score_effects(run_manifest: pd.DataFrame) -> pd.DataFrame:
    """Compare each parameter variant's objective score with its experiment baseline."""
    cols = [
        "variant_id", "exp_id", "group", "parameter", "level", "parameter_value",
        "baseline_score", "score", "delta_score", "abs_delta_score", "improved",
        "status",
    ]
    if run_manifest.empty or "score" not in run_manifest:
        return pd.DataFrame(columns=cols)
    df = run_manifest.copy()
    if "parameter_value" not in df.columns and "value" in df.columns:
        df = df.rename(columns={"value": "parameter_value"})
    base = df[df["variant_kind"] == "baseline"][["exp_id", "score"]].rename(
        columns={"score": "baseline_score"}
    )
    variants = df[df["variant_kind"] == "parameter"].merge(base, on="exp_id", how="left")
    if variants.empty:
        return pd.DataFrame(columns=cols)
    variants["delta_score"] = variants["score"] - variants["baseline_score"]
    variants["abs_delta_score"] = variants["delta_score"].abs()
    variants["improved"] = variants["delta_score"] < 0
    return variants[[c for c in cols if c in variants.columns]].sort_values(
        ["abs_delta_score", "group", "parameter"], ascending=[False, True, True]
    ).reset_index(drop=True)


def summarize_output_effects(output_effects: pd.DataFrame, *, tolerance: float = 1e-9) -> pd.DataFrame:
    """Aggregate raw DSSAT output deltas by parameter and output variable."""
    cols = [
        "group", "parameter", "source_file", "variable", "stat", "n_comparisons",
        "n_changed", "mean_abs_delta", "max_abs_delta", "mean_abs_relative_delta",
        "max_abs_relative_delta",
    ]
    if output_effects.empty:
        return pd.DataFrame(columns=cols)
    df = output_effects.copy()
    df["abs_delta"] = pd.to_numeric(df["delta"], errors="coerce").abs()
    df["abs_relative_delta"] = pd.to_numeric(df["relative_delta"], errors="coerce").abs()
    df["changed"] = df["abs_delta"] > float(tolerance)
    group_cols = ["group", "parameter", "source_file", "variable", "stat"]
    summary = (
        df.groupby(group_cols, dropna=False)
        .agg(
            n_comparisons=("delta", "size"),
            n_changed=("changed", "sum"),
            mean_abs_delta=("abs_delta", "mean"),
            max_abs_delta=("abs_delta", "max"),
            mean_abs_relative_delta=("abs_relative_delta", "mean"),
            max_abs_relative_delta=("abs_relative_delta", "max"),
        )
        .reset_index()
    )
    summary = summary.sort_values(
        ["max_abs_delta", "max_abs_relative_delta", "n_changed"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    return summary[cols]


def summarize_parameter_impacts(
    score_effects: pd.DataFrame,
    output_impact_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse objective and DSSAT-output impacts to one row per parameter."""
    keys = ["group", "parameter"]
    score_cols = [
        "n_score_variants", "max_abs_score_delta", "best_score_delta",
        "worst_score_delta", "n_improved",
    ]
    output_cols = [
        "n_output_summaries", "n_changed_output_summaries", "max_abs_output_delta",
        "max_abs_relative_output_delta",
    ]

    if score_effects.empty:
        score = pd.DataFrame(columns=keys + score_cols)
    else:
        score = (
            score_effects.groupby(keys, dropna=False)
            .agg(
                n_score_variants=("delta_score", "size"),
                max_abs_score_delta=("abs_delta_score", "max"),
                best_score_delta=("delta_score", "min"),
                worst_score_delta=("delta_score", "max"),
                n_improved=("improved", "sum"),
            )
            .reset_index()
        )

    if output_impact_summary.empty:
        output = pd.DataFrame(columns=keys + output_cols)
    else:
        output = (
            output_impact_summary.groupby(keys, dropna=False)
            .agg(
                n_output_summaries=("variable", "size"),
                n_changed_output_summaries=("n_changed", lambda s: int((pd.to_numeric(s, errors="coerce") > 0).sum())),
                max_abs_output_delta=("max_abs_delta", "max"),
                max_abs_relative_output_delta=("max_abs_relative_delta", "max"),
            )
            .reset_index()
        )

    if score.empty and output.empty:
        return pd.DataFrame(columns=keys + score_cols + output_cols)
    if score.empty:
        out = output
        for c in score_cols:
            out[c] = np.nan
    elif output.empty:
        out = score
        for c in output_cols:
            out[c] = np.nan
    else:
        out = score.merge(output, on=keys, how="outer")
    return out.sort_values(
        ["max_abs_score_delta", "max_abs_output_delta", "max_abs_relative_output_delta"],
        ascending=[False, False, False],
        na_position="last",
    ).reset_index(drop=True)


def write_impact_summary_report(
    path: Path,
    run_manifest: pd.DataFrame,
    parameter_impact_summary: pd.DataFrame,
    output_impact_summary: pd.DataFrame,
    score_effects: pd.DataFrame,
) -> None:
    """Write a compact human/agent-readable Markdown summary."""
    def fmt_num(value, digits: str = ".6g") -> str:
        try:
            if pd.isna(value):
                return ""
            return format(float(value), digits)
        except (TypeError, ValueError):
            return ""

    def fmt_int(value) -> str:
        try:
            if pd.isna(value):
                return "0"
            return str(int(value))
        except (TypeError, ValueError):
            return "0"

    total = len(run_manifest)
    ok = int(run_manifest["status"].isin(["success", "cached"]).sum()) if total and "status" in run_manifest else 0
    lines = [
        "# Impact Atlas Summary",
        "",
        f"Spawns succeeded: {ok}/{total}",
        "",
        "## Top Parameters",
        "",
        "| group | parameter | max abs score delta | max abs output delta | changed output summaries |",
        "|---|---|---:|---:|---:|",
    ]
    if parameter_impact_summary.empty:
        lines.append("| | | | | |")
    else:
        for _, r in parameter_impact_summary.head(20).iterrows():
            lines.append(
                f"| {r.get('group', '')} | {r.get('parameter', '')} | "
                f"{fmt_num(r.get('max_abs_score_delta', np.nan))} | "
                f"{fmt_num(r.get('max_abs_output_delta', np.nan))} | "
                f"{fmt_int(r.get('n_changed_output_summaries', 0))} |"
            )

    lines.extend([
        "",
        "## Top Objective Score Changes",
        "",
        "| experiment | group | parameter | level | delta score | score |",
        "|---|---|---|---|---:|---:|",
    ])
    if score_effects.empty:
        lines.append("| | | | | | |")
    else:
        for _, r in score_effects.head(20).iterrows():
            lines.append(
                f"| {r.get('exp_id', '')} | {r.get('group', '')} | {r.get('parameter', '')} | "
                f"{r.get('level', '')} | {fmt_num(r.get('delta_score', np.nan))} | "
                f"{fmt_num(r.get('score', np.nan))} |"
            )

    lines.extend([
        "",
        "## Top DSSAT Output Changes",
        "",
        "| group | parameter | file | variable | stat | max abs delta | max abs relative delta |",
        "|---|---|---|---|---|---:|---:|",
    ])
    if output_impact_summary.empty:
        lines.append("| | | | | | | |")
    else:
        for _, r in output_impact_summary.head(30).iterrows():
            lines.append(
                f"| {r.get('group', '')} | {r.get('parameter', '')} | {r.get('source_file', '')} | "
                f"{r.get('variable', '')} | {r.get('stat', '')} | "
                f"{fmt_num(r.get('max_abs_delta', np.nan))} | "
                f"{fmt_num(r.get('max_abs_relative_delta', np.nan))} |"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def capability_map(parameters: list[dict]) -> pd.DataFrame:
    rows = []
    for spec in parameters:
        group = spec.get("group", "")
        if group in {"genetic_cultivar", "genetic_ecotype", "genetic_species"}:
            owner = "dssatengine"
            note = "Generic genotype file editing and DSSAT run-folder materialization should migrate to dssatengine; calibrator keeps experiment design and scoring."
        elif group in {"management", "initial_conditions"}:
            owner = "dssatengine"
            note = "FileX section editing is reusable run-folder functionality and belongs in dssatengine once stabilized."
        elif group == "soil":
            owner = "dssatengine + dssatutils"
            note = "Local .SOL editing belongs in dssatengine; acquiring/building soil profiles from external sources belongs in dssatutils."
        elif group == "weather":
            owner = "dssatengine + dssatutils"
            note = "Local .WTH perturbation belongs in dssatengine; weather download/conversion/gap-fill belongs in dssatutils."
        else:
            owner = "dssatcalibrator"
            note = "Calibration-specific analysis remains local."
        rows.append({
            "group": group,
            "parameter": spec.get("name"),
            "current_support": "implemented" if group in DEFAULT_GROUPS else "unknown",
            "recommended_owner": owner,
            "notes": note,
        })
    return pd.DataFrame(rows)


def write_capability_report(path: Path, cap: pd.DataFrame) -> None:
    lines = [
        "# Impact Atlas Capability Map",
        "",
        "This file maps the tested parameter groups to the package that should own the reusable capability.",
        "",
        "| group | parameter | recommended owner | notes |",
        "|---|---|---|---|",
    ]
    for _, r in cap.iterrows():
        lines.append(f"| {r['group']} | {r['parameter']} | {r['recommended_owner']} | {r['notes']} |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_impact_atlas(
    cfg: dict,
    *,
    output_dir: str | Path,
    groups: list[str] | None = None,
    experiments: list[str] | None = None,
    levels: list[str] | None = None,
    active_only: bool = False,
    discover_cultivar: bool = False,
    discover_ecotype: bool = False,
    discover_species: bool = False,
    discover_genotype: bool = False,
    allow_species: bool = False,
    max_parameters: int | None = None,
    max_per_group: int | None = None,
    output_files: list[str] | None = None,
    num_cores: int | None = None,
    write_long: bool = True,
    compress_long: bool = False,
    effect_tolerance: float = 1e-9,
    clean: bool = True,
    progress: bool = True,
) -> AtlasResult:
    """Run a one-at-a-time real-DSSAT parameter impact atlas."""
    cfg = deepcopy(cfg)
    groups = groups or list(DEFAULT_GROUPS)
    levels = levels or ["low", "high"]
    if discover_genotype:
        discover_cultivar = discover_ecotype = discover_species = True
    if experiments is not None:
        cfg["experiments"] = list(experiments)
    if num_cores is not None:
        cfg.setdefault("calibrator", {})["num_cores"] = int(num_cores)
    cfg.setdefault("calibrator", {})["keep_run_dirs"] = True
    cfg.setdefault("calibrator", {})["cache_spawns"] = False

    output_dir = Path(output_dir)
    run_root = output_dir / "run_dirs"
    if clean and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_root.mkdir(parents=True, exist_ok=True)

    crop = crop_for(cfg, (cfg.get("crops") or [{}])[0].get("code", "HM"))
    exe = resolve_exe(cfg)
    hemp_dir = Path(cfg["source"]["hemp_dir"])
    obs = _load_observations(cfg, cfg.get("experiments", []), crop)
    obs_exps = set(obs.experiments())
    exps = [e for e in cfg.get("experiments", []) if not obs_exps or e in obs_exps]
    if not exps:
        raise ValueError(
            "Impact atlas has no experiments to run after observation filtering. "
            "Check config.experiments and the observation source."
        )
    treatments = {e: parse_treatments(hemp_dir / f"{e}.{crop['filex_ext']}") for e in exps}

    params = _candidate_parameters(cfg, groups, active_only=active_only)
    if discover_cultivar and "genetic_cultivar" in groups:
        params.extend(discover_cultivar_parameters(cfg))
    if discover_ecotype and "genetic_ecotype" in groups:
        params.extend(discover_ecotype_parameters(cfg))
    if discover_species and "genetic_species" in groups:
        params.extend(discover_species_parameters(cfg))
    # Remove duplicate group/name pairs while preserving order.
    dedup = []
    seen = set()
    for p in params:
        key = (p.get("group"), p.get("name"), p.get("spe_key"), p.get("spe_index"))
        if key not in seen:
            seen.add(key)
            dedup.append(p)
    params = dedup
    if max_per_group is not None:
        limited = []
        counts: dict[str, int] = {}
        for p in params:
            group = str(p.get("group"))
            counts[group] = counts.get(group, 0)
            if counts[group] < max_per_group:
                limited.append(p)
                counts[group] += 1
        params = limited
    params = params[:max_parameters] if max_parameters else params
    parameter_catalog = pd.DataFrame(params)
    species_sweep_requested = any(p.get("group") == "genetic_species" for p in params)
    species_gate_opened = False
    if species_sweep_requested:
        current_gate = str((cfg.get("gating") or {}).get("species", "blocked")).lower()
        if current_gate != "free" and not allow_species:
            raise ValueError(
                "Atlas includes genetic_species parameters but gating.species is "
                f"'{current_gate}'. Set gating.species: free in the config or pass "
                "allow_species=True / --allow-species to explicitly permit .SPE edits."
            )
        if allow_species:
            cfg.setdefault("gating", {})["species"] = "free"
            species_gate_opened = current_gate != "free"

    jobs = []
    meta = []
    for exp_id in exps:
        variant_id = f"{exp_id}__baseline"
        jobs.append({
            "theta": {},
            "exp_id": exp_id,
            "cfg": cfg,
            "crop": crop,
            "param_specs": [],
            "run_root": run_root,
            "treatments": treatments[exp_id],
            "exe": exe,
        })
        meta.append({
            "variant_id": variant_id, "exp_id": exp_id, "variant_kind": "baseline",
            "group": "baseline", "parameter": "baseline", "level": "baseline", "value": np.nan,
        })
        for spec in params:
            for var in _variant_specs(spec, levels):
                variant_id = f"{exp_id}__{_safe_id(spec['group'])}__{_safe_id(spec['name'])}__{var['level']}"
                jobs.append({
                    "theta": var["theta"],
                    "exp_id": exp_id,
                    "cfg": cfg,
                    "crop": crop,
                    "param_specs": var["param_specs"],
                    "run_root": run_root,
                    "treatments": treatments[exp_id],
                    "exe": exe,
                })
                meta.append({
                    "variant_id": variant_id, "exp_id": exp_id, "variant_kind": "parameter",
                    "group": spec["group"], "parameter": spec["name"],
                    "level": var["level"], "value": var["value"],
                })

    n_workers = resolve_cores(cfg["calibrator"].get("num_cores", 0))
    done = {"n": 0}

    def _cb(_res):
        done["n"] += 1
        if progress and (done["n"] % max(1, len(jobs) // 20) == 0 or done["n"] == len(jobs)):
            print(f"  atlas spawns {done['n']}/{len(jobs)}", flush=True)

    results = run_many(jobs, n_workers=n_workers, on_done=_cb if progress else None)

    manifest_rows = []
    file_manifest_frames = []
    long_frames = []
    effects_frames = []
    baselines: dict[str, pd.DataFrame] = {}

    for m, res in zip(meta, results):
        score = np.inf
        n_obs = 0
        if res.status in {"success", "cached"}:
            scored = obj.score({m["exp_id"]: res}, obs.table, cfg)
            score = scored.score
            n_obs = len(scored.residuals)
        collected = dssat_io.collect_run_outputs(res.run_dir, output_files=output_files)
        long = collected["long"].copy()
        if not long.empty:
            for key, value in m.items():
                long["parameter_value" if key == "value" else key] = value
            if write_long:
                long_frames.append(long)
        stats = _output_stats(long)
        if m["variant_kind"] == "baseline":
            baselines[m["exp_id"]] = stats
        else:
            eff = _effects_for_run(stats, baselines.get(m["exp_id"], pd.DataFrame()))
            if not eff.empty:
                for key, value in m.items():
                    eff["parameter_value" if key == "value" else key] = value
                effects_frames.append(eff)

        out_manifest = collected["manifest"]
        if not out_manifest.empty:
            fm = out_manifest.copy()
            for key, value in m.items():
                fm["parameter_value" if key == "value" else key] = value
            fm["run_dir"] = str(res.run_dir)
            file_manifest_frames.append(fm)
        files_present = int(out_manifest["exists"].sum()) if not out_manifest.empty else 0
        manifest_rows.append({
            **{"parameter_value" if key == "value" else key: value for key, value in m.items()},
            "status": res.status,
            "message": res.message,
            "score": score,
            "n_obs": n_obs,
            "run_dir": str(res.run_dir),
            "files_present": files_present,
        })

    run_manifest = pd.DataFrame(manifest_rows)
    file_manifest = pd.concat(file_manifest_frames, ignore_index=True, sort=False) if file_manifest_frames else pd.DataFrame()
    output_long = pd.concat(long_frames, ignore_index=True, sort=False) if long_frames else pd.DataFrame()
    output_effects = pd.concat(effects_frames, ignore_index=True, sort=False) if effects_frames else pd.DataFrame()
    score_effects = summarize_score_effects(run_manifest)
    output_impact_summary = summarize_output_effects(output_effects, tolerance=effect_tolerance)
    parameter_impact_summary = summarize_parameter_impacts(score_effects, output_impact_summary)
    cap = capability_map(params)

    run_manifest.to_csv(output_dir / "run_manifest.csv", index=False)
    if not file_manifest.empty:
        file_manifest.to_csv(output_dir / "file_manifest.csv", index=False)
    parameter_catalog.to_csv(output_dir / "parameter_catalog.csv", index=False)
    output_effects.to_csv(output_dir / "parameter_output_effects.csv", index=False)
    score_effects.to_csv(output_dir / "score_effects.csv", index=False)
    output_impact_summary.to_csv(output_dir / "output_impact_summary.csv", index=False)
    parameter_impact_summary.to_csv(output_dir / "parameter_impact_summary.csv", index=False)
    if write_long and not output_long.empty:
        long_path = output_dir / ("outputs_long.csv.gz" if compress_long else "outputs_long.csv")
        output_long.to_csv(long_path, index=False, compression="gzip" if compress_long else None)
    cap.to_csv(output_dir / "capability_map.csv", index=False)
    write_capability_report(output_dir / "capability_map.md", cap)
    write_impact_summary_report(
        output_dir / "impact_summary.md",
        run_manifest,
        parameter_impact_summary,
        output_impact_summary,
        score_effects,
    )
    (output_dir / "atlas_config.json").write_text(json.dumps({
        "experiments": exps,
        "groups": groups,
        "levels": levels,
        "num_parameters": len(params),
        "max_per_group": max_per_group,
        "num_jobs": len(jobs),
        "output_files": output_files or dssat_io.DEFAULT_OUTPUT_FILES,
        "write_long": bool(write_long),
        "compress_long": bool(compress_long),
        "effect_tolerance": float(effect_tolerance),
        "discover_cultivar": bool(discover_cultivar),
        "discover_ecotype": bool(discover_ecotype),
        "discover_species": bool(discover_species),
        "discover_genotype": bool(discover_genotype),
        "allow_species": bool(allow_species),
        "species_sweep_requested": bool(species_sweep_requested),
        "species_gate_opened": bool(species_gate_opened),
    }, indent=2), encoding="utf-8")

    return AtlasResult(
        output_dir,
        run_manifest,
        file_manifest,
        output_effects,
        output_long,
        cap,
        parameter_catalog,
        score_effects,
        output_impact_summary,
        parameter_impact_summary,
    )

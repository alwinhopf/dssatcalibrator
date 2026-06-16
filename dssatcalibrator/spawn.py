"""Spawn = deterministic function of (theta, experiment) -> a DSSAT run + outputs.

A spawn copies the base FileX and genotype files into an isolated run directory,
writes the perturbed cultivar coefficients (and, when activated, management /
initial conditions), runs ``dscsm048``, and parses PlantGro.OUT + Evaluate.OUT.

Weather and soil are resolved centrally by DSSAT via ``DSSATPRO.V48`` (DSSAT
reads genotype/weather/soil from the run dir first, else the install), so a
spawn only needs the FileX + edited genotype files — no per-spawn weather/soil
copy. Observed FileA/FileT are copied in so DSSAT writes Evaluate.OUT.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import dssat_io
from .writers import edit_cultivar, edit_ecotype

GENETIC_GROUPS = {"genetic_cultivar"}


def theta_hash(theta: dict[str, float]) -> str:
    blob = json.dumps({k: round(float(v), 6) for k, v in sorted(theta.items())})
    return hashlib.sha1(blob.encode()).hexdigest()[:10]


def parse_treatments(filex_path: str | Path) -> list[int]:
    """Return the treatment numbers from a FileX ``*TREATMENTS`` section."""
    lines = Path(filex_path).read_text(errors="replace").splitlines()
    try:
        start = next(i for i, ln in enumerate(lines) if ln.startswith("*TREATMENTS"))
    except StopIteration:
        return [1]
    trts = []
    for ln in lines[start + 1:]:
        if ln.startswith("*"):
            break
        if re.match(r"\s*\d", ln) and not ln.lstrip().startswith("@"):
            trts.append(int(ln.split()[0]))
    return trts or [1]


def write_dssbatch(run_dir: Path, filex_name: str, treatments: list[int]) -> Path:
    """Write a DSSAT ``B``-mode batch file listing the requested treatments."""
    header = ("$BATCH(CALIB)\n!\n@FILEX" + " " * 86 +
              "TRTNO     RP     SQ     OP     CO\n")
    rows = "".join(f"{filex_name:<93}{t:6d}      1      0      1      0\n" for t in treatments)
    p = run_dir / "DSSBatch.v48"
    p.write_text(header + rows)
    return p


@dataclass
class SpawnResult:
    status: str                      # "success" | "error" | "cached"
    run_dir: Path
    theta: dict
    plantgro: pd.DataFrame = field(default_factory=pd.DataFrame)
    evaluate: pd.DataFrame = field(default_factory=pd.DataFrame)
    message: str = ""


def _partition_theta(theta: dict, param_specs: list[dict]) -> dict[str, dict]:
    """Split a flat theta into per-group update dicts using the param specs."""
    group_of = {p["name"]: p["group"] for p in param_specs}
    groups: dict[str, dict] = {}
    for name, val in theta.items():
        g = group_of.get(name, "genetic_cultivar")
        groups.setdefault(g, {})[name] = val
    return groups


def spawn_and_run(
    theta: dict,
    *,
    exp_id: str,
    cfg: dict,
    crop: dict,
    param_specs: list[dict],
    run_root: Path,
    treatments: list[int] | None = None,
    exe: Path,
    timeout: int = 600,
) -> SpawnResult:
    """Materialize and run one spawn; return parsed PlantGro + Evaluate tables."""
    from .config import crop_for  # local import to avoid cycle at module load

    dssat_dir = Path(cfg["calibrator"]["dssat_dir"])
    hemp_dir = Path(cfg["source"]["hemp_dir"])
    geno_dir = dssat_dir / "Genotype"
    stem = crop["genotype_stem"]
    ext = crop["filex_ext"]          # FileX extension, e.g. "HMX"
    code = crop["code"]              # crop code for observed files, e.g. "HM" -> .HMA/.HMT

    run_dir = run_root / exp_id / f"s_{theta_hash(theta)}"
    pg_path = run_dir / "PlantGro.OUT"

    if cfg["calibrator"].get("cache_spawns", True) and pg_path.exists() and pg_path.stat().st_size > 0:
        return SpawnResult(
            status="cached", run_dir=run_dir, theta=theta,
            plantgro=dssat_io.parse_plantgro(pg_path),
            evaluate=dssat_io.parse_evaluate(run_dir / "Evaluate.OUT"),
        )

    run_dir.mkdir(parents=True, exist_ok=True)

    # genotype files (edit a local copy)
    for e in ("CUL", "ECO", "SPE"):
        src = geno_dir / f"{stem}.{e}"
        if src.exists():
            shutil.copy(src, run_dir / f"{stem}.{e}")

    # base FileX + observed files (FileA/FileT for Evaluate.OUT)
    filex_name = f"{exp_id}.{ext}"
    shutil.copy(hemp_dir / filex_name, run_dir / filex_name)
    for obs_ext in (f"{code}A", f"{code}T"):   # observed files use the crop CODE (e.g. .HMA/.HMT)
        src = hemp_dir / f"{exp_id}.{obs_ext}"
        if src.exists():
            shutil.copy(src, run_dir / f"{exp_id}.{obs_ext}")

    # apply parameter perturbations
    groups = _partition_theta(theta, param_specs)
    cul_updates = {}
    for g in GENETIC_GROUPS:
        cul_updates.update(groups.get(g, {}))
    if cul_updates:
        edit_cultivar(run_dir / f"{stem}.CUL", crop["cultivar_anchor"], cul_updates)

    eco_updates = groups.get("genetic_ecotype", {})
    if eco_updates:
        edit_ecotype(run_dir / f"{stem}.ECO", crop["ecotype"], eco_updates)


    # edit management and initial conditions in FileX
    mgt_updates = groups.get("management", {})
    init_updates = groups.get("initial_conditions", {})
    if mgt_updates or init_updates:
        mgt_fields = {}
        for name, val in mgt_updates.items():
            spec = next((s for s in param_specs if s["name"] == name), None)
            if spec and "dssat" in spec:
                mgt_fields[spec["dssat"]] = val
        
        from .writers import edit_filex
        edit_filex(run_dir / filex_name, mgt_fields, init_updates)

    # edit soil (.SOL) and/or weather (.WTH) — DSSAT reads the run dir first, so a
    # local single-profile SOIL.SOL / station .WTH overrides the central copy.
    soil_updates = groups.get("soil", {})
    weather_updates = groups.get("weather", {})
    if soil_updates or weather_updates:
        from .writers import parse_fields
        fields = parse_fields(run_dir / filex_name)

        if soil_updates and fields.get("id_soil"):
            from .writers import extract_soil_profile, edit_soil
            src_sol = dssat_dir / "Soil" / "SOIL.SOL"
            pid = fields["id_soil"]
            if src_sol.exists():
                layer_mults, profile_sets = {}, {}
                for name, val in soil_updates.items():
                    spec = next((s for s in param_specs if s["name"] == name), None)
                    if not spec or "dssat" not in spec:
                        continue
                    (profile_sets if spec.get("op") == "set" else layer_mults)[spec["dssat"]] = val
                (run_dir / "SOIL.SOL").write_text(extract_soil_profile(src_sol, pid))
                edit_soil(run_dir / "SOIL.SOL", pid, layer_mults=layer_mults, profile_sets=profile_sets)

        if weather_updates and fields.get("wsta"):
            from .writers import edit_weather
            wsta = fields["wsta"]
            src_wth = dssat_dir / "Weather" / f"{wsta}.WTH"
            if src_wth.exists():
                ops = {}
                for name, val in weather_updates.items():
                    spec = next((s for s in param_specs if s["name"] == name), None)
                    if not spec or "dssat" not in spec:
                        continue
                    ops[spec["dssat"]] = (spec.get("op", "mult"), val)
                shutil.copy(src_wth, run_dir / f"{wsta}.WTH")
                edit_weather(run_dir / f"{wsta}.WTH", ops)

    if treatments is None:
        treatments = parse_treatments(run_dir / filex_name)
    write_dssbatch(run_dir, filex_name, treatments)

    try:
        subprocess.run(
            [str(exe), crop["model"], "B", "DSSBatch.v48"],
            cwd=str(run_dir), timeout=timeout,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
        )
    except subprocess.TimeoutExpired:
        return SpawnResult(status="error", run_dir=run_dir, theta=theta, message="timeout")

    if not pg_path.exists() or pg_path.stat().st_size == 0:
        return SpawnResult(status="error", run_dir=run_dir, theta=theta,
                           message="no PlantGro.OUT produced")

    pg = dssat_io.parse_plantgro(pg_path)
    ev = dssat_io.parse_evaluate(run_dir / "Evaluate.OUT")

    if not cfg["calibrator"].get("keep_run_dirs", False):
        if not cfg["calibrator"].get("cache_spawns", True):
            shutil.rmtree(run_dir, ignore_errors=True)
        else:
            # keep the parsed outputs only; drop the bulky per-run artifacts
            for f in run_dir.glob("*.OUT"):
                if f.name not in ("PlantGro.OUT", "Evaluate.OUT", "Summary.OUT"):
                    f.unlink(missing_ok=True)

    return SpawnResult(status="success", run_dir=run_dir, theta=theta, plantgro=pg, evaluate=ev)

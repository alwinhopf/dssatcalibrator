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
from .config import resolve_dssat_paths
from .writers import edit_cultivar, edit_ecotype

GENETIC_GROUPS = {"genetic_cultivar"}
BATCH_FILE = "DSSBatch.V48"


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
    p = run_dir / BATCH_FILE
    p.write_text(header + rows, encoding="utf-8")
    return p


def _execution_backend(cfg: dict) -> str:
    backend = str(cfg.get("execution", {}).get("backend", "native")).lower()
    if backend not in {"native", "dssatengine"}:
        raise ValueError("execution.backend must be 'native' or 'dssatengine'.")
    return backend


def _dssatengine_api():
    try:
        from dssatengine import normalize_treatment_list, run_dssat, write_dssbatch as engine_write_dssbatch
    except ImportError as exc:
        raise ImportError(
            "execution.backend: dssatengine requires 'dssatcalibrator[shared]' "
            "or an installed dssatengine package."
        ) from exc
    return normalize_treatment_list, run_dssat, engine_write_dssbatch


def _normalize_treatments(treatments: list[int], backend: str) -> list[int]:
    if backend == "dssatengine":
        normalize_treatment_list, _, _ = _dssatengine_api()
        return normalize_treatment_list(1, 1, treatment_list=treatments)

    seen = set()
    out = []
    for value in treatments:
        trt = int(value)
        if trt < 1:
            raise ValueError("Treatment IDs must be positive integers.")
        if trt not in seen:
            seen.add(trt)
            out.append(trt)
    if not out:
        raise ValueError("No valid treatments selected.")
    return out


def _write_batch(run_dir: Path, filex_name: str, treatments: list[int], backend: str) -> Path:
    if backend == "dssatengine":
        _, _, engine_write_dssbatch = _dssatengine_api()
        batch = run_dir / BATCH_FILE
        engine_write_dssbatch(filex_name, treatments, str(batch), run_mode="experiment")
        return batch
    return write_dssbatch(run_dir, filex_name, treatments)


def _run_native_dssat(run_dir: Path, exe: Path, model: str, timeout: int) -> str:
    try:
        result = subprocess.run(
            [str(exe), model, "B", BATCH_FILE],
            cwd=str(run_dir),
            timeout=timeout,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
    except subprocess.TimeoutExpired:
        return "timeout"

    output = "\n".join(part for part in (result.stdout, result.stderr) if part)
    if output:
        (run_dir / "dssat_B_stdout_stderr.log").write_text(output.rstrip() + "\n",
                                                           encoding="utf-8")
    if result.returncode != 0:
        tail = " | ".join(output.splitlines()[-12:]) if output else "<no stdout/stderr captured>"
        return f"DSSAT exited with status {result.returncode}. Tail: {tail}"
    return ""


def _run_backend_dssat(run_dir: Path, exe: Path, crop: dict,
                       backend: str, timeout: int) -> str:
    if backend == "native":
        return _run_native_dssat(run_dir, exe, crop["model"], timeout)
    try:
        _, run_dssat, _ = _dssatengine_api()
        run_dssat(str(run_dir), str(exe), "B", model=crop.get("model"), timeout=timeout)
    except subprocess.TimeoutExpired:
        return "timeout"
    except Exception as exc:
        return str(exc)
    return ""


def _weather_window(cfg: dict) -> tuple[pd.Timestamp, pd.Timestamp]:
    wcfg = cfg.get("weather", {}) or {}
    start = wcfg.get("start") or wcfg.get("start_date")
    end = wcfg.get("end") or wcfg.get("end_date")
    if start is None and wcfg.get("start_year") is not None:
        start = f"{int(wcfg['start_year'])}-01-01"
    if end is None and wcfg.get("end_year") is not None:
        end = f"{int(wcfg['end_year'])}-12-31"
    if start is None or end is None:
        raise ValueError(
            "weather.provider acquisition requires weather.start/end "
            "or weather.start_year/end_year."
        )
    return pd.Timestamp(start), pd.Timestamp(end)


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
    backend = _execution_backend(cfg)
    dssat_paths = resolve_dssat_paths(cfg)
    hemp_dir = Path(cfg["source"]["hemp_dir"])
    geno_dir = dssat_paths["genotype"]
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

    # DSSAT profile (DSSATPRO) — CSM reads it from the current directory first,
    # then its compiled default path. Copying it from the install root into the
    # run dir lets a non-standard install (one whose root differs from the
    # binary's compiled default, e.g. a relocated or per-user DSSAT48) resolve
    # its Genotype/Weather/Soil paths. Harmless when CSM would find it anyway.
    for pro in ("DSSATPRO.L48", "DSSATPRO.V48", "DSSATPRO.v48", "DSCSM048.CTR"):
        src = dssat_paths["root"] / pro
        if src.exists():
            shutil.copy(src, run_dir / pro)

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

    # species (.SPE) coefficients — gated: only written when gating.species == "free"
    # (physiology-defining; for new-species adaptation from an analog template).
    spe_updates = groups.get("genetic_species", {})
    if spe_updates and str(cfg.get("gating", {}).get("species", "blocked")).lower() == "free":
        from .writers import edit_species
        spe_file = run_dir / f"{stem}.SPE"
        if spe_file.exists():
            updates = {}
            for name, val in spe_updates.items():
                spec = next((s for s in param_specs if s["name"] == name), None)
                updates[(spec or {}).get("spe_key", name)] = val
            edit_species(spe_file, updates)

    # edit management and initial conditions in FileX
    mgt_updates = groups.get("management", {})
    init_updates = groups.get("initial_conditions", {})
    mgt_fields = {}
    for name, val in mgt_updates.items():
        spec = next((s for s in param_specs if s["name"] == name), None)
        if spec and "dssat" in spec:
            mgt_fields[spec["dssat"]] = val
    # per-experiment planting date override (e.g. from farm-management software):
    # set PDATE directly rather than calibrating it. cfg["_planting_dates"] maps
    # exp_id -> a date; written as the DSSAT YYDDD code.
    pdate = (cfg.get("_planting_dates") or {}).get(exp_id)
    if pdate is not None:
        ts = pd.Timestamp(pdate)
        mgt_fields["PDATE"] = int(f"{ts.year % 100:02d}{ts.dayofyear:03d}")
    if mgt_fields or init_updates:
        from .writers import edit_filex
        edit_filex(run_dir / filex_name, mgt_fields, init_updates)

    soil_updates = groups.get("soil", {})
    weather_updates = groups.get("weather", {})
    soil_provider = str(cfg.get("soil", {}).get("provider", "file")).lower()
    weather_provider = str(cfg.get("weather", {}).get("provider", "file")).lower()
    needs_fields = (
        soil_updates or weather_updates
        or soil_provider not in ("", "file", "none")
        or weather_provider not in ("", "file", "none")
    )
    fields = {}
    if needs_fields:
        from .writers import parse_fields
        fields = parse_fields(run_dir / filex_name)

    if soil_provider not in ("", "file", "none"):
        try:
            from .acquisition import acquire_soil_profile

            site_id = fields.get("id_soil") or fields.get("id_field") or exp_id
            lat = fields.get("lat", cfg.get("soil", {}).get("lat"))
            lon = fields.get("lon", cfg.get("soil", {}).get("lon"))
            acquire_soil_profile(
                cfg,
                site_id=str(site_id),
                lat=lat,
                lon=lon,
                out_path=run_dir / "SOIL.SOL",
            )
        except Exception as exc:
            return SpawnResult(status="error", run_dir=run_dir, theta=theta,
                               message=f"soil acquisition failed: {exc}")

    if weather_provider not in ("", "file", "none"):
        try:
            from .weather import acquire_wth

            station = fields.get("wsta") or exp_id
            lat = fields.get("lat", cfg.get("weather", {}).get("lat"))
            lon = fields.get("lon", cfg.get("weather", {}).get("lon"))
            start, end = _weather_window(cfg)
            acquire_wth(
                cfg,
                station=str(station),
                lat=lat,
                lon=lon,
                start=start,
                end=end,
                out_path=run_dir / f"{station}.WTH",
            )
        except Exception as exc:
            return SpawnResult(status="error", run_dir=run_dir, theta=theta,
                               message=f"weather acquisition failed: {exc}")

    # edit soil (.SOL) and/or weather (.WTH) — DSSAT reads the run dir first, so a
    # local single-profile SOIL.SOL / station .WTH overrides the central copy.
    if soil_updates or weather_updates:
        if soil_updates and fields.get("id_soil"):
            from .writers import extract_soil_profile, edit_soil
            local_sol = run_dir / "SOIL.SOL"
            src_sol = local_sol if local_sol.exists() else dssat_paths["soil"] / "SOIL.SOL"
            pid = fields["id_soil"]
            if src_sol.exists():
                layer_mults, profile_sets = {}, {}
                for name, val in soil_updates.items():
                    spec = next((s for s in param_specs if s["name"] == name), None)
                    if not spec or "dssat" not in spec:
                        continue
                    (profile_sets if spec.get("op") == "set" else layer_mults)[spec["dssat"]] = val
                if src_sol != local_sol:
                    local_sol.write_text(extract_soil_profile(src_sol, pid), encoding="utf-8")
                edit_soil(local_sol, pid, layer_mults=layer_mults, profile_sets=profile_sets)

        if weather_updates and fields.get("wsta"):
            from .writers import edit_weather
            wsta = fields["wsta"]
            local_wth = run_dir / f"{wsta}.WTH"
            src_wth = local_wth if local_wth.exists() else dssat_paths["weather"] / f"{wsta}.WTH"
            if src_wth.exists():
                ops = {}
                for name, val in weather_updates.items():
                    spec = next((s for s in param_specs if s["name"] == name), None)
                    if not spec or "dssat" not in spec:
                        continue
                    ops[spec["dssat"]] = (spec.get("op", "mult"), val)
                if src_wth != local_wth:
                    shutil.copy(src_wth, local_wth)
                edit_weather(local_wth, ops)

    if treatments is None:
        treatments = parse_treatments(run_dir / filex_name)
    try:
        treatments = _normalize_treatments(treatments, backend)
        _write_batch(run_dir, filex_name, treatments, backend)
    except Exception as exc:
        return SpawnResult(status="error", run_dir=run_dir, theta=theta,
                           message=f"batch setup failed: {exc}")

    run_error = _run_backend_dssat(run_dir, exe, crop, backend, timeout)
    if run_error:
        return SpawnResult(status="error", run_dir=run_dir, theta=theta, message=run_error)

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

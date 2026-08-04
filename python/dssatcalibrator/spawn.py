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
    def normalize(value):
        if isinstance(value, bool):
            return value
        try:
            return float(value)
        except (TypeError, ValueError):
            return str(value)

    blob = json.dumps({k: normalize(v) for k, v in sorted(theta.items())})
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def _file_digest(path: Path) -> str | None:
    if not path.exists() or not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _spawn_provenance(cfg, crop, param_specs, source_filex, geno_dir,
                      dssat_paths, exe, treatments, effective_theta) -> dict:
    try:
        from .writers import parse_fields
        fields = parse_fields(source_filex)
    except Exception:
        fields = {}
    wsta = fields.get("wsta")
    payload = {
        "schema": 2,
        "theta": {k: v for k, v in sorted(effective_theta.items())},
        "crop": crop,
        "specs": param_specs,
        "treatments": treatments,
        "filex_sha256": _file_digest(source_filex),
        "genotype_sha256": {ext: _file_digest(geno_dir / f"{crop['genotype_stem']}.{ext}")
                             for ext in ("CUL", "ECO", "SPE")},
        "weather_sha256": _file_digest(dssat_paths["weather"] / f"{wsta}.WTH") if wsta else None,
        "soil_sha256": _file_digest(dssat_paths["soil"] / "SOIL.SOL"),
        "exe_sha256": _file_digest(Path(exe)),
        "execution": cfg.get("execution", {}),
        "gating": cfg.get("gating", {}),
        "weather": cfg.get("weather", {}),
        "soil": cfg.get("soil", {}),
        "filex_overrides": cfg.get("filex_overrides", {}),
    }
    return json.loads(json.dumps(payload, sort_keys=True, default=str))


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


def parse_cultivars(filex_path: str | Path) -> list[str]:
    """Return cultivar codes listed in a FileX ``*CULTIVARS`` section."""
    lines = Path(filex_path).read_text(errors="replace").splitlines()
    try:
        start = next(i for i, ln in enumerate(lines) if ln.startswith("*CULTIVARS"))
    except StopIteration:
        return []
    header = None
    out = []
    for ln in lines[start + 1:]:
        if ln.startswith("*"):
            break
        if ln.lstrip().startswith("@"):
            header = ln.lstrip().lstrip("@").split()
            continue
        if header and re.match(r"\s*\d", ln):
            values = ln.split()
            row = dict(zip(header, values))
            code = row.get("INGENO") or row.get("CULTIVAR") or row.get("VAR#")
            if code and code not in out:
                out.append(code)
    return out


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


def _treatment_run_key(treatments: list[int] | None) -> str | None:
    if treatments is None:
        return None
    vals = sorted({int(t) for t in treatments})
    if not vals:
        return None
    return "T" + "-".join(str(t) for t in vals)


def _stamp_single_treatment(outputs: dict[str, pd.DataFrame], treatments: list[int] | None) -> dict[str, pd.DataFrame]:
    vals = sorted({int(t) for t in treatments or []})
    if len(vals) != 1:
        return outputs
    trt = vals[0]
    stamped = {}
    for key, df in (outputs or {}).items():
        if isinstance(df, pd.DataFrame) and not df.empty:
            out = df.copy()
            if "treatment" in out.columns:
                out["treatment"] = pd.to_numeric(out["treatment"], errors="coerce").fillna(trt).astype("Int64")
            else:
                out["treatment"] = trt
            stamped[key] = out
        else:
            stamped[key] = df
    return stamped


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


def _rewrite_dssat_profile_paths(run_dir: Path, root: Path) -> None:
    """Point copied DSSAT profile files at the configured DSSAT root."""
    drive = root.drive.upper()
    tail = str(root.resolve())[len(root.drive):].replace("/", "\\")
    if not drive or not tail:
        return
    dssat_profile_root = f"{drive} {tail}"
    windows_root = f"{drive}{tail}"
    replacements = {
        "C: \\DSSAT48": dssat_profile_root,
        "C:\\DSSAT48": windows_root,
        "c: \\DSSAT48": dssat_profile_root,
        "c:\\DSSAT48": windows_root,
    }
    for name in ("DSSATPRO.L48", "DSSATPRO.V48", "DSSATPRO.v48"):
        path = run_dir / name
        if not path.exists():
            continue
        text = path.read_text(errors="replace")
        for old, new in replacements.items():
            text = text.replace(old, new)
        path.write_text(text, encoding="utf-8")


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
    outputs: dict[str, pd.DataFrame] = field(default_factory=dict)
    message: str = ""
    effective_theta: dict = field(default_factory=dict)


def _spec_applies(spec: dict, exp_id: str | None, cultivars: list[str] | None = None) -> bool:
    scope = spec.get("scope")
    if scope == "experiment":
        return exp_id is None or str(spec.get("exp_id")) == str(exp_id)
    if scope == "cultivar":
        return cultivars is None or str(spec.get("cultivar")) in set(map(str, cultivars))
    return True


def _effective_theta(
    theta: dict,
    param_specs: list[dict],
    exp_id: str | None = None,
    cultivars: list[str] | None = None,
) -> dict[str, float]:
    """Return the coefficient/value pairs that actually affect one spawn."""
    out: dict[str, float] = {}
    matched = set()
    for spec in param_specs:
        name = spec["name"]
        if not _spec_applies(spec, exp_id, cultivars):
            continue
        if name in theta:
            value = theta[name]
        elif spec.get("fixed", False):
            value = spec.get("start")
        else:
            continue
        if value is None:
            continue
        base = spec.get("base_name", name)
        if spec.get("scope") == "cultivar":
            base = f"{base}__{spec.get('cultivar')}"
        out[base] = value
        matched.add(name)
    if not matched:
        return {} if param_specs else dict(theta)
    return out


def _partition_theta(
    theta: dict,
    param_specs: list[dict],
    exp_id: str | None = None,
    cultivars: list[str] | None = None,
) -> dict[str, dict]:
    """Split a flat theta into per-group update dicts using the param specs."""
    groups: dict[str, dict] = {}
    matched = set()
    for spec in param_specs:
        name = spec["name"]
        if not _spec_applies(spec, exp_id, cultivars):
            continue
        if name in theta:
            value = theta[name]
        elif spec.get("fixed", False):
            value = spec.get("start")
        else:
            continue
        if value is None:
            continue
        g = spec.get("group", "genetic_cultivar")
        base = spec.get("base_name", name)
        if spec.get("scope") == "cultivar":
            anchor = str(spec.get("cultivar"))
            groups.setdefault(f"{g}_by_cultivar", {}).setdefault(anchor, {})[base] = value
        else:
            groups.setdefault(g, {})[base] = value
        matched.add(name)
    if not matched and not param_specs:
        for name, val in theta.items():
            groups.setdefault("genetic_cultivar", {})[name] = val
    return groups


def _cultivar_ecotype_map(crop: dict) -> dict[str, str]:
    mapping = {str(k): str(v) for k, v in (crop.get("cultivar_ecotypes") or {}).items()}
    if crop.get("cultivar_anchor") and crop.get("ecotype"):
        mapping.setdefault(str(crop["cultivar_anchor"]), str(crop["ecotype"]))
    return mapping


def _genotype_gate_allows(cfg: dict, level: str) -> bool:
    """Return whether an explicitly configured genotype update may be written."""
    gate = str((cfg.get("gating", {}) or {}).get(level, "free")).lower()
    if level == "species":
        return gate == "free"
    return gate != "blocked"


def _missing_requested_treatments(pg: pd.DataFrame, treatments: list[int]) -> list[int]:
    """Identify requested treatments absent from a parsed PlantGro output."""
    if pg is None or pg.empty or "treatment" not in pg:
        return [int(value) for value in treatments]
    present = {
        int(value) for value in pd.to_numeric(pg["treatment"], errors="coerce").dropna()
    }
    return [int(value) for value in treatments if int(value) not in present]


def _run_identity(theta: dict, effective_theta: dict, cache_spawns: bool) -> dict:
    """Return a collision-free run identity.

    The full candidate is required even with caching: two candidates can differ
    only in a cultivar irrelevant to this experiment, observe the same empty
    cache entry, and otherwise race in one effective-theta directory.
    """
    return theta


def _filex_overrides_for(cfg: dict, exp_id: str) -> list[dict]:
    """Return generic FileX section updates configured for one experiment."""
    block = cfg.get("filex_overrides", {}) or {}
    updates = []
    for rec in block.get("all", []) or []:
        if isinstance(rec, dict):
            updates.append(dict(rec))
    for rec in block.get(exp_id, []) or []:
        if isinstance(rec, dict):
            updates.append(dict(rec))
    return updates


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
    filex_name = f"{exp_id}.{ext}"
    source_filex = hemp_dir / filex_name
    exp_cultivars = parse_cultivars(source_filex)

    effective_theta = _effective_theta(theta, param_specs, exp_id, exp_cultivars)
    provenance = _spawn_provenance(
        cfg, crop, param_specs, source_filex, geno_dir, dssat_paths, exe,
        treatments, effective_theta,
    )
    treatment_key = _treatment_run_key(treatments)
    provenance_blob = json.dumps(provenance, sort_keys=True, separators=(",", ":"), default=str)
    provenance_hash = hashlib.sha256(provenance_blob.encode("utf-8")).hexdigest()[:12]
    run_dir = run_root / exp_id
    if treatment_key:
        run_dir = run_dir / treatment_key
<<<<<<< Updated upstream
    # Isolate runs when any template, forcing, executable, parser, gate, or
    # configuration input changes, even if the coefficient vector is identical.
    run_dir = run_dir / f"s_{theta_hash(effective_theta)}_{provenance_hash}"
=======
    # Always use the full candidate. Effective-theta hashes remain in the
    # manifest for traceability, but are unsafe as writable directory identities
    # when multiple candidates are evaluated concurrently.
    run_identity = _run_identity(
        theta,
        effective_theta,
        bool(cfg["calibrator"].get("cache_spawns", True)),
    )
    run_dir = run_dir / f"s_{theta_hash(run_identity)}"
>>>>>>> Stashed changes
    pg_path = run_dir / "PlantGro.OUT"
    manifest_path = run_dir / "spawn_manifest.json"

    recorded = None
    try:
        recorded = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        pass
    if (cfg["calibrator"].get("cache_spawns", True)
            and recorded == provenance and pg_path.exists() and pg_path.stat().st_size > 0):
        outputs = dssat_io.collect_run_outputs(run_dir)
        outputs = _stamp_single_treatment(outputs, treatments)
<<<<<<< Updated upstream
        pg_cached = outputs.get("plantgro", dssat_io.parse_plantgro(pg_path))
        if pg_cached is None or pg_cached.empty:
            recorded = None
        else:
            if treatments and "treatment" in pg_cached.columns:
                found = set(pd.to_numeric(pg_cached["treatment"], errors="coerce").dropna().astype(int))
                if not found.issubset(set(map(int, treatments))):
                    recorded = None
        if recorded == provenance:
            return SpawnResult(
                status="cached", run_dir=run_dir, theta=theta,
                plantgro=pg_cached,
                evaluate=outputs.get("evaluate", dssat_io.parse_evaluate(run_dir / "Evaluate.OUT")),
                outputs=outputs,
            )
=======
        return SpawnResult(
            status="cached", run_dir=run_dir, theta=theta,
            plantgro=outputs.get("plantgro", dssat_io.parse_plantgro(pg_path)),
            evaluate=outputs.get("evaluate", dssat_io.parse_evaluate(run_dir / "Evaluate.OUT")),
            outputs=outputs,
            effective_theta=effective_theta,
        )
>>>>>>> Stashed changes

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
    _rewrite_dssat_profile_paths(run_dir, dssat_paths["root"])

    # genotype files (edit a local copy)
    for e in ("CUL", "ECO", "SPE"):
        src = geno_dir / f"{stem}.{e}"
        if src.exists():
            shutil.copy(src, run_dir / f"{stem}.{e}")

    # base FileX + observed files (FileA/FileT for Evaluate.OUT)
    shutil.copy(source_filex, run_dir / filex_name)
    filex_overrides = _filex_overrides_for(cfg, exp_id)
    if filex_overrides:
        from .writers import edit_filex
        edit_filex(run_dir / filex_name, {}, {}, section_updates=filex_overrides)
    for obs_ext in (f"{code}A", f"{code}T"):   # observed files use the crop CODE (e.g. .HMA/.HMT)
        src = hemp_dir / f"{exp_id}.{obs_ext}"
        if src.exists():
            shutil.copy(src, run_dir / f"{exp_id}.{obs_ext}")

    # apply parameter perturbations
    groups = _partition_theta(theta, param_specs, exp_id, exp_cultivars)
    cul_updates = {}
    for g in GENETIC_GROUPS:
        cul_updates.update(groups.get(g, {}))
<<<<<<< Updated upstream
    if cul_updates and str(cfg.get("gating", {}).get("cultivar", "free")).lower() != "blocked":
        anchors = crop.get("cultivar_anchors") or [crop["cultivar_anchor"]]
        for anchor in anchors:
            edit_cultivar(run_dir / f"{stem}.CUL", anchor, cul_updates)
    if str(cfg.get("gating", {}).get("cultivar", "free")).lower() != "blocked":
=======
    if cul_updates and _genotype_gate_allows(cfg, "cultivar"):
        anchors = crop.get("cultivar_anchors") or [crop["cultivar_anchor"]]
        for anchor in anchors:
            edit_cultivar(run_dir / f"{stem}.CUL", anchor, cul_updates)
    if _genotype_gate_allows(cfg, "cultivar"):
>>>>>>> Stashed changes
        for anchor, updates in groups.get("genetic_cultivar_by_cultivar", {}).items():
            edit_cultivar(run_dir / f"{stem}.CUL", anchor, updates)

    eco_updates = groups.get("genetic_ecotype", {})
<<<<<<< Updated upstream
    if eco_updates and str(cfg.get("gating", {}).get("ecotype", "free")).lower() != "blocked":
        edit_ecotype(run_dir / f"{stem}.ECO", crop["ecotype"], eco_updates)
    cultivar_ecotypes = _cultivar_ecotype_map(crop)
    for anchor, updates in (groups.get("genetic_ecotype_by_cultivar", {}).items()
                            if str(cfg.get("gating", {}).get("ecotype", "free")).lower() != "blocked" else []):
        eco_anchor = cultivar_ecotypes.get(anchor)
        if eco_anchor is None:
            raise ValueError(
                f"No ecotype mapping for cultivar '{anchor}'. Add crops[].cultivar_ecotypes."
            )
        edit_ecotype(run_dir / f"{stem}.ECO", eco_anchor, updates)
=======
    if eco_updates and _genotype_gate_allows(cfg, "ecotype"):
        edit_ecotype(run_dir / f"{stem}.ECO", crop["ecotype"], eco_updates)
    cultivar_ecotypes = _cultivar_ecotype_map(crop)
    if _genotype_gate_allows(cfg, "ecotype"):
        for anchor, updates in groups.get("genetic_ecotype_by_cultivar", {}).items():
            eco_anchor = cultivar_ecotypes.get(anchor)
            if eco_anchor is None:
                raise ValueError(
                    f"No ecotype mapping for cultivar '{anchor}'. Add crops[].cultivar_ecotypes."
                )
            edit_ecotype(run_dir / f"{stem}.ECO", eco_anchor, updates)
>>>>>>> Stashed changes

    # species (.SPE) coefficients — gated: only written when gating.species == "free"
    # (physiology-defining; for new-species adaptation from an analog template).
    spe_updates = groups.get("genetic_species", {})
<<<<<<< Updated upstream
    if spe_updates and str(cfg.get("gating", {}).get("species", "blocked")).lower() != "blocked":
=======
    if spe_updates and _genotype_gate_allows(cfg, "species"):
>>>>>>> Stashed changes
        from .writers import edit_species
        spe_file = run_dir / f"{stem}.SPE"
        if spe_file.exists():
            for name, val in spe_updates.items():
                spec = next((s for s in param_specs if s.get("base_name", s["name"]) == name), None)
                key = (spec or {}).get("spe_key", name)
                if spec and ("spe_index" in spec or "token_index" in spec):
                    update = {
                        "value": val,
                        "index": int(spec.get("spe_index", spec.get("token_index", 0))),
                    }
                else:
                    update = val
                # Apply one species edit at a time so several calibrated
                # parameters can target different numeric tokens on the same
                # .SPE line via the same spe_key.
                edit_species(spe_file, {key: update})

    # edit management and initial conditions in FileX
    mgt_updates = groups.get("management", {})
    init_updates = groups.get("initial_conditions", {})
    def filex_update_from_spec(name, val, spec, default_section):
        if spec is None:
            return val
        section = spec.get("section", spec.get("filex_section"))
        field = spec.get("field", spec.get("filex_field", spec.get("dssat")))
        is_soil_water_mult = default_section == "INITIAL CONDITIONS" and name == "initial_soil_water_mult"
        if is_soil_water_mult and not field:
            field = "SH2O"
        generic_keys = {
            "header_prefix", "row", "treatment", "trt", "trtno",
            "clip_01", "required", "type", "format",
        }
        if not section and field:
            uses_generic = bool(generic_keys.intersection(spec)) or str(spec.get("op", "set")).lower() != "set"
            if default_section == "PLANTING DETAILS" and not uses_generic:
                return val
            out = {
                "section": default_section,
                "field": field,
                "value": val,
                "op": spec.get("op", "mult" if is_soil_water_mult else "set"),
            }
            if is_soil_water_mult and "clip_01" not in spec:
                out["clip_01"] = True
            for key in generic_keys:
                if key in spec:
                    out[key] = spec[key]
            return out
        out = {
            "section": section or default_section,
            "field": field or name,
            "value": val,
            "op": spec.get("op", "mult" if is_soil_water_mult else "set"),
        }
        if is_soil_water_mult and "clip_01" not in spec:
            out["clip_01"] = True
        for key in ("header_prefix", "row", "treatment", "trt", "trtno", "clip_01", "required", "type", "format"):
            if key in spec:
                out[key] = spec[key]
        return out

    mgt_fields = {}
    for name, val in mgt_updates.items():
        spec = next((s for s in param_specs if s.get("base_name", s["name"]) == name), None)
        if spec:
            key = spec.get("dssat", spec.get("field", spec.get("filex_field", name)))
            mgt_fields[key] = filex_update_from_spec(name, val, spec, "PLANTING DETAILS")
    # per-experiment planting date override (e.g. from farm-management software):
    # set PDATE directly rather than calibrating it. cfg["_planting_dates"] maps
    # exp_id -> a date; written as the DSSAT YYDDD code.
    pdate = (cfg.get("_planting_dates") or {}).get(exp_id)
    if pdate is not None:
        ts = pd.Timestamp(pdate)
        mgt_fields["PDATE"] = int(f"{ts.year % 100:02d}{ts.dayofyear:03d}")
    if mgt_fields or init_updates:
        from .writers import edit_filex
        init_fields = {}
        for name, val in init_updates.items():
            spec = next((s for s in param_specs if s.get("base_name", s["name"]) == name), None)
            if spec:
                key = spec.get("dssat", spec.get("field", spec.get("filex_field", name)))
                init_fields[key] = filex_update_from_spec(name, val, spec, "INITIAL CONDITIONS")
            else:
                init_fields[name] = val
        edit_filex(run_dir / filex_name, mgt_fields, init_fields)

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
                               message=f"soil acquisition failed: {exc}",
                               effective_theta=effective_theta)

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
                               message=f"weather acquisition failed: {exc}",
                               effective_theta=effective_theta)

    # edit soil (.SOL) and/or weather (.WTH) — DSSAT reads the run dir first, so a
    # local single-profile SOIL.SOL / station .WTH overrides the central copy.
    if soil_updates or weather_updates:
        if soil_updates and fields.get("id_soil"):
            from .writers import extract_soil_profile, edit_soil
            local_sol = run_dir / "SOIL.SOL"
            pid = fields["id_soil"]
            candidates = []
            if local_sol.exists():
                candidates.append(local_sol)
            candidates.append(dssat_paths["soil"] / "SOIL.SOL")
            candidates.extend(sorted(dssat_paths["soil"].glob("*.SOL")))
            profile_text = None
            for src_sol in candidates:
                if not src_sol.exists():
                    continue
                try:
                    profile_text = extract_soil_profile(src_sol, pid)
                    break
                except ValueError:
                    continue
            if profile_text is None:
                raise ValueError(f"Soil profile '{pid}' not found in DSSAT soil directory {dssat_paths['soil']}")
            layer_mults, profile_sets = {}, {}
            for name, val in soil_updates.items():
                spec = next((s for s in param_specs if s.get("base_name", s["name"]) == name), None)
                if not spec or "dssat" not in spec:
                    continue
                (profile_sets if spec.get("op") == "set" else layer_mults)[spec["dssat"]] = val
            local_sol.write_text(profile_text, encoding="utf-8")
            edit_soil(local_sol, pid, layer_mults=layer_mults, profile_sets=profile_sets)

        if weather_updates and fields.get("wsta"):
            from .writers import edit_weather
            wsta = fields["wsta"]
            local_wth = run_dir / f"{wsta}.WTH"
            src_wth = local_wth if local_wth.exists() else dssat_paths["weather"] / f"{wsta}.WTH"
            if src_wth.exists():
                ops = {}
                for name, val in weather_updates.items():
                    spec = next((s for s in param_specs if s.get("base_name", s["name"]) == name), None)
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
                           message=f"batch setup failed: {exc}",
                           effective_theta=effective_theta)

    run_error = _run_backend_dssat(run_dir, exe, crop, backend, timeout)
    if run_error:
        return SpawnResult(status="error", run_dir=run_dir, theta=theta,
                           message=run_error, effective_theta=effective_theta)

    if not pg_path.exists() or pg_path.stat().st_size == 0:
        return SpawnResult(status="error", run_dir=run_dir, theta=theta,
                           message="no PlantGro.OUT produced",
                           effective_theta=effective_theta)

    pg = dssat_io.parse_plantgro(pg_path)
    missing_treatments = _missing_requested_treatments(pg, treatments)
    if missing_treatments:
        return SpawnResult(
            status="error",
            run_dir=run_dir,
            theta=theta,
            plantgro=pg,
            message=f"PlantGro.OUT is missing requested treatment(s): {missing_treatments}",
            effective_theta=effective_theta,
        )
    ev = dssat_io.parse_evaluate(run_dir / "Evaluate.OUT")
    outputs = dssat_io.collect_run_outputs(run_dir)
    outputs = _stamp_single_treatment(outputs, treatments)
    manifest_path.write_text(json.dumps(provenance, indent=2, sort_keys=True, default=str) + "\n",
                             encoding="utf-8")

    if not cfg["calibrator"].get("keep_run_dirs", False):
        if not cfg["calibrator"].get("cache_spawns", True):
            shutil.rmtree(run_dir, ignore_errors=True)
        else:
            # keep the parsed outputs only; drop the bulky per-run artifacts
            for f in run_dir.glob("*.OUT"):
                if f.name not in ("PlantGro.OUT", "Evaluate.OUT", "Summary.OUT"):
                    f.unlink(missing_ok=True)

    return SpawnResult(
        status="success", run_dir=run_dir, theta=theta, plantgro=pg,
        evaluate=ev, outputs=outputs, effective_theta=effective_theta,
    )

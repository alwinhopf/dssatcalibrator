"""Persistent objective-result cache for expensive DSSAT evaluations.

The spawn cache avoids rerunning one ``(theta, experiment)`` folder when DSSAT
outputs already exist. This higher-level cache stores the scored objective for a
whole ``theta`` across the requested experiments, so optimizers and samplers can
skip spawning/parsing/scoring duplicate candidates across batches or reruns.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

import pandas as pd

from . import objective as obj
from .config import resolve_dssat_paths

CACHE_SCHEMA_VERSION = 2


def _normalise(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _normalise(v) for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))}
    if isinstance(value, (list, tuple)):
        return [_normalise(v) for v in value]
    if isinstance(value, (str, bool)) or value is None:
        return value
    if isinstance(value, (int,)):
        return int(value)
    if isinstance(value, float):
        if math.isfinite(value):
            return round(float(value), 12)
        return {"__float__": str(value)}
    if hasattr(value, "item"):
        return _normalise(value.item())
    return str(value)


def _file_fingerprint(path: str | Path) -> dict[str, Any] | None:
    p = Path(path)
    if not p.exists() or not p.is_file():
        return None
    st = p.stat()
    h = hashlib.sha1()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return {
        "path": str(p.resolve()),
        "size": int(st.st_size),
        "sha1": h.hexdigest(),
    }


def _frame_digest(df: pd.DataFrame | None) -> str:
    if df is None or df.empty:
        payload = "[]"
    else:
        work = df.copy()
        for col in work.columns:
            if pd.api.types.is_datetime64_any_dtype(work[col]):
                work[col] = work[col].dt.strftime("%Y-%m-%dT%H:%M:%S")
        payload = work.to_json(orient="records", date_format="iso", default_handler=str)
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()


def _encode_scalar(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, pd.Timestamp):
        return None if pd.isna(value) else value.isoformat()
    if hasattr(value, "isoformat") and not isinstance(value, str):
        try:
            return value.isoformat()
        except Exception:
            pass
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    if hasattr(value, "item"):
        return _encode_scalar(value.item())
    if isinstance(value, float):
        if math.isfinite(value):
            return float(value)
        return {"__float__": "inf" if value > 0 else "-inf" if value < 0 else "nan"}
    if isinstance(value, int):
        return int(value)
    return value


def _decode_scalar(value: Any) -> Any:
    if isinstance(value, dict) and "__float__" in value:
        tag = value["__float__"]
        return float("inf") if tag == "inf" else float("-inf") if tag == "-inf" else float("nan")
    return value


def _frame_to_records(df: pd.DataFrame | None) -> list[dict[str, Any]]:
    if df is None or df.empty:
        return []
    return [
        {str(k): _encode_scalar(v) for k, v in row.items()}
        for row in df.to_dict(orient="records")
    ]


def _records_to_frame(records: list[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame([
        {k: _decode_scalar(v) for k, v in row.items()}
        for row in (records or [])
    ])
    if "date" in df.columns:
        df["date"] = pd.to_datetime(df["date"], errors="coerce")
    return df


def _encode_float(value: float) -> float | str:
    value = float(value)
    if math.isfinite(value):
        return value
    return "Inf" if value > 0 else "-Inf" if value < 0 else "NaN"


def _decode_float(value: float | str) -> float:
    if value == "Inf":
        return float("inf")
    if value == "-Inf":
        return float("-inf")
    if value == "NaN":
        return float("nan")
    return float(value)


def _cfg_fingerprint(cfg: dict) -> str:
    keep = {
        "parameters": cfg.get("parameters", {}),
        "crops": cfg.get("crops", []),
        "source": cfg.get("source", {}),
        "engine": cfg.get("engine", {}),
        "objective": cfg.get("objective", {}),
        "execution": cfg.get("execution", {}),
        "templates": cfg.get("templates", {}),
        "gating": cfg.get("gating", {}),
        "management_options": cfg.get("management_options", {}),
        "filex_overrides": cfg.get("filex_overrides", {}),
        "calibration_treatments": cfg.get("calibration_treatments", {}),
        "calibration_treatments_by_experiment": cfg.get(
            "calibration_treatments_by_experiment", {}
        ),
        "weather": cfg.get("weather", {}),
        "soil": cfg.get("soil", {}),
        "observation_sources": cfg.get("observation_sources", {}),
        "fusion": cfg.get("fusion", {}),
        "experiments": cfg.get("experiments", []),
        "cache_salt": (cfg.get("calibrator", {}) or {}).get("evaluation_cache_salt", ""),
    }
    payload = json.dumps(_normalise(keep), sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()


class EvaluationCache:
    """On-disk cache for ``ObjectiveResult`` objects."""

    def __init__(self, root: Path | None, context: dict[str, Any] | None = None):
        self.root = root
        self.context = context or {}
        self.hits = 0
        self.misses = 0
        self.writes = 0

    @property
    def enabled(self) -> bool:
        return self.root is not None

    @classmethod
    def from_setup(
        cls,
        cfg: dict,
        *,
        crop: dict,
        specs: list[dict],
        experiments: list[str],
        treatments: dict,
        obs_table: pd.DataFrame,
        exe: str | Path,
    ) -> "EvaluationCache":
        ccfg = cfg.get("calibrator", {}) or {}
        if not ccfg.get("cache_evaluations", True):
            return cls(None)

        run_root = Path(ccfg.get("workdir", "results/_workdir")) / str(ccfg.get("name", "run"))
        raw_dir = ccfg.get("evaluation_cache_dir") or ccfg.get("eval_cache_dir")
        if raw_dir:
            root = Path(raw_dir)
        else:
            root = run_root / "evaluation_cache"
        if raw_dir and not root.is_absolute():
            root = run_root / root

        hemp_dir = Path((cfg.get("source", {}) or {}).get("hemp_dir", ""))
        input_files: dict[str, Any] = {"exe": _file_fingerprint(exe)}
        try:
            dssat_paths = resolve_dssat_paths(cfg)
            geno_dir = dssat_paths["genotype"]
            stem = crop.get("genotype_stem", "")
            input_files["genotype"] = {
                ext: _file_fingerprint(geno_dir / f"{stem}.{ext}")
                for ext in ("CUL", "ECO", "SPE")
            }
            input_files["soil"] = {
                path.name: _file_fingerprint(path)
                for path in sorted(dssat_paths["soil"].glob("*.SOL"))
            }
        except Exception as exc:  # keep cache usable even when setup is mocked
            input_files["genotype_error"] = str(exc)

        filex_ext = crop.get("filex_ext", "")
        code = crop.get("code", "")
        experiment_inputs = {}
        for exp in experiments:
            filex = hemp_dir / f"{exp}.{filex_ext}"
            rec = {
                "filex": _file_fingerprint(hemp_dir / f"{exp}.{filex_ext}"),
                "filea": _file_fingerprint(hemp_dir / f"{exp}.{code}A"),
                "filet": _file_fingerprint(hemp_dir / f"{exp}.{code}T"),
            }
            try:
                from .writers import parse_fields

                fields = parse_fields(filex)
                station = fields.get("wsta")
                if station and "dssat_paths" in locals():
                    rec["weather"] = _file_fingerprint(
                        dssat_paths["weather"] / f"{station}.WTH"
                    )
            except Exception as exc:
                rec["field_parse_error"] = str(exc)
            experiment_inputs[exp] = rec
        input_files["experiments"] = experiment_inputs
        try:
            dssat_paths = resolve_dssat_paths(cfg)
            input_files["soil_library"] = _file_fingerprint(
                dssat_paths["soil"] / "SOIL.SOL"
            )
            input_files["dssat_profile"] = next(
                (
                    _file_fingerprint(dssat_paths["root"] / name)
                    for name in ("DSSATPRO.V48", "DSSATPRO.L48", "DSCSM048.CTR")
                    if (dssat_paths["root"] / name).exists()
                ),
                None,
            )
        except Exception as exc:
            input_files["support_error"] = str(exc)
        input_files["parser"] = _file_fingerprint(Path(__file__).with_name("dssat_io.py"))
        input_files["objective"] = _file_fingerprint(Path(__file__).with_name("objective.py"))

        context = {
            "schema": CACHE_SCHEMA_VERSION,
            "cfg": _cfg_fingerprint(cfg),
            "crop": _normalise(crop),
            "specs": _normalise(specs),
            "treatments": _normalise(treatments),
            "obs": _frame_digest(obs_table),
            "inputs": _normalise(input_files),
        }
        return cls(root, context)

    def key(self, theta: dict, experiments: list[str]) -> str:
        payload = {
            "context": self.context,
            "experiments": list(experiments),
            "theta": _normalise(theta),
        }
        blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha1(blob.encode("utf-8")).hexdigest()

    def path(self, key: str) -> Path:
        assert self.root is not None
        return self.root / key[:2] / f"{key}.json"

    def get(self, key: str) -> obj.ObjectiveResult | None:
        if not self.enabled:
            return None
        path = self.path(key)
        if not path.exists():
            self.misses += 1
            return None
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("schema") != CACHE_SCHEMA_VERSION:
                self.misses += 1
                return None
            self.hits += 1
            return obj.ObjectiveResult(
                score=_decode_float(payload["score"]),
                loglik=_decode_float(payload["loglik"]),
                residuals=_records_to_frame(payload.get("residuals", [])),
                per_var=payload.get("per_var", {}) or {},
                per_exp_var=_records_to_frame(payload.get("per_exp_var", [])),
            )
        except Exception:
            self.misses += 1
            return None

    def put(self, key: str, result: obj.ObjectiveResult) -> None:
        if not self.enabled:
            return
        path = self.path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema": CACHE_SCHEMA_VERSION,
            "score": _encode_float(result.score),
            "loglik": _encode_float(result.loglik),
            "per_var": _normalise(result.per_var),
            "residuals": _frame_to_records(result.residuals),
            "per_exp_var": _frame_to_records(result.per_exp_var),
        }
        tmp = path.with_suffix(f".tmp-{os.getpid()}")
        tmp.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")), encoding="utf-8")
        tmp.replace(path)
        self.writes += 1

    def stats(self) -> dict[str, int]:
        return {"hits": self.hits, "misses": self.misses, "writes": self.writes}

"""Cloud entrypoint helpers for DSSAT calibration jobs.

The cloud runner stages inputs, rewrites only path-like config fields, then
delegates to ``run_calibration.py``.  That keeps SageMaker/EKS plumbing outside
the scientific calibration code.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml


SAGEMAKER_INPUT = Path("/opt/ml/processing/input")
SAGEMAKER_OUTPUT = Path("/opt/ml/processing/output")
EKS_INPUT = Path("/mnt/dssat/input")
EKS_OUTPUT = Path("/mnt/dssat/output")


def _deep_update(base: dict[str, Any], patch: dict[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for key, value in patch.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_update(out[key], value)
        else:
            out[key] = value
    return out


def _job_index() -> str | None:
    for name in ("DSSATCAL_SHARD_INDEX", "JOB_COMPLETION_INDEX",
                 "AWS_BATCH_JOB_ARRAY_INDEX", "SM_CURRENT_INSTANCE_GROUP"):
        val = os.environ.get(name)
        if val not in (None, ""):
            return str(val)
    return None


def format_indexed(value: str | None, index: str | None = None) -> str | None:
    """Format a shard-aware string containing ``{index}``.

    Strings without the token are returned unchanged.  Missing index leaves the
    value untouched so local single-job runs stay simple.
    """
    if value is None:
        return None
    idx = _job_index() if index is None else index
    if idx is None or "{index}" not in value:
        return value
    return value.format(index=idx)


def default_mounts() -> tuple[Path, Path]:
    """Return the input/output mount pair for SageMaker, EKS, or local Docker."""
    if SAGEMAKER_INPUT.exists():
        return SAGEMAKER_INPUT, SAGEMAKER_OUTPUT
    if EKS_INPUT.exists():
        return EKS_INPUT, EKS_OUTPUT
    return Path(os.environ.get("DSSATCAL_INPUT_DIR", "/work/input")), Path(
        os.environ.get("DSSATCAL_OUTPUT_DIR", "/work/output")
    )


def sync_s3(source: str | None, dest: Path, *, direction: str) -> None:
    """Run ``aws s3 sync`` when an S3 URI is provided."""
    if not source:
        return
    dest.mkdir(parents=True, exist_ok=True)
    if direction == "down":
        cmd = ["aws", "s3", "sync", source, str(dest)]
    elif direction == "up":
        cmd = ["aws", "s3", "sync", str(dest), source]
    else:
        raise ValueError(f"unknown sync direction: {direction}")
    subprocess.run(cmd, check=True)


def find_config(input_dir: Path, explicit: str | None = None,
                pattern: str | None = None, index: str | None = None) -> Path:
    """Find the calibration YAML inside the staged input directory."""
    chosen = format_indexed(pattern, index) or format_indexed(explicit, index)
    if chosen:
        path = Path(chosen)
        return path if path.is_absolute() else input_dir / path

    candidates = sorted(input_dir.glob("*.yaml")) + sorted(input_dir.glob("*.yml"))
    if len(candidates) == 1:
        return candidates[0]
    preferred = input_dir / "config.yaml"
    if preferred.exists():
        return preferred
    names = ", ".join(p.name for p in candidates[:8]) or "none"
    raise FileNotFoundError(
        "Could not choose a calibration config. Set DSSATCAL_CONFIG or "
        f"DSSATCAL_CONFIG_PATTERN. YAML files found in {input_dir}: {names}"
    )


def find_dssat_dir(input_dir: Path, explicit: str | None = None) -> Path:
    """Resolve the mounted DSSAT install directory."""
    if explicit:
        return Path(explicit)
    for candidate in (
        input_dir / "DSSAT48",
        input_dir / "dssat" / "DSSAT48",
        Path("/opt/dssat/DSSAT48"),
    ):
        if candidate.exists():
            return candidate
    return input_dir / "DSSAT48"


def find_dssat_exe(dssat_dir: Path, explicit: str | None = None) -> Path:
    """Resolve the DSSAT executable, preferring the Linux binary when present."""
    if explicit:
        return Path(explicit)
    for name in ("dscsm048", "DSCSM048.EXE", "dscsm048.exe"):
        candidate = dssat_dir / name
        if candidate.exists():
            return candidate
    return dssat_dir / "dscsm048"


def find_source_dir(input_dir: Path, dssat_dir: Path, explicit: str | None = None) -> Path | None:
    """Resolve the experiment/source directory used by the calibrator."""
    if explicit:
        return Path(explicit)
    for candidate in (
        input_dir / "source",
        input_dir / "Hemp",
        dssat_dir / "Hemp",
    ):
        if candidate.exists():
            return candidate
    return None


def find_template_dir(input_dir: Path, explicit: str | None = None) -> Path | None:
    """Resolve optional shared DSSAT templates mounted with the job."""
    if explicit:
        return Path(explicit)
    for candidate in (
        input_dir / "dssat_templates",
        input_dir / "DSSAT_Gridded_Run_Tutorial" / "dssat_templates",
    ):
        if candidate.exists():
            return candidate
    return None


def prepare_config(config_path: Path, output_dir: Path, input_dir: Path,
                   *, dssat_dir: Path | None = None, dssat_exe: Path | None = None,
                   source_dir: Path | None = None, template_dir: Path | None = None,
                   num_cores: int | None = None) -> tuple[Path, dict[str, Any]]:
    """Write a cloud-local config and return ``(path, cfg)``."""
    with open(config_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}

    dssat_dir = dssat_dir or find_dssat_dir(input_dir)
    dssat_exe = dssat_exe or find_dssat_exe(dssat_dir)
    source_dir = source_dir if source_dir is not None else find_source_dir(input_dir, dssat_dir)
    template_dir = template_dir if template_dir is not None else find_template_dir(input_dir)

    patch: dict[str, Any] = {
        "calibrator": {
            "dssat_dir": str(dssat_dir),
            "dssat_exe": str(dssat_exe),
            "workdir": str(output_dir / "_workdir"),
            "results_dir": str(output_dir / "results"),
            "figures_dir": str(output_dir / "figures"),
        }
    }
    if num_cores is not None:
        patch["calibrator"]["num_cores"] = int(num_cores)
    if source_dir is not None:
        patch["source"] = {"hemp_dir": str(source_dir)}
    if template_dir is not None:
        patch["templates"] = {"template_dir": str(template_dir)}

    cloud_cfg = _deep_update(cfg, patch)
    run_dir = output_dir / "_cloud"
    run_dir.mkdir(parents=True, exist_ok=True)
    out_path = run_dir / "config.cloud.yaml"
    with open(out_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cloud_cfg, fh, sort_keys=False)
    return out_path, cloud_cfg


def write_manifest(output_dir: Path, *, input_dir: Path, config_path: Path,
                   cloud_config_path: Path, cfg: dict[str, Any],
                   extra_args: list[str], language: str) -> None:
    """Write a compact launch manifest for reproducibility."""
    manifest = {
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "config": str(config_path),
        "cloud_config": str(cloud_config_path),
        "shard_index": _job_index(),
        "language": language,
        "extra_args": extra_args,
        "calibrator": cfg.get("calibrator", {}),
        "experiments": cfg.get("experiments", []),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "cloud_manifest.json", "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)


def run_calibration(repo_root: Path, cloud_config: Path, extra_args: list[str],
                    *, language: str = "python") -> None:
    """Invoke the existing calibration CLI."""
    language = language.lower()
    if language == "python":
        script = Path(os.environ.get("DSSATCAL_RUNNER", repo_root / "run_calibration.py"))
        cmd = [sys.executable, str(script), str(cloud_config), *extra_args]
    elif language == "r":
        script = Path(os.environ.get("DSSATCAL_RUNNER", repo_root / "run_calibration.R"))
        rscript = os.environ.get("DSSATCAL_RSCRIPT", "Rscript")
        cmd = [rscript, str(script), str(cloud_config), *extra_args]
    else:
        raise ValueError("language must be 'python' or 'r'")
    env = os.environ.copy()
    py_path = str(repo_root / "python")
    env["PYTHONPATH"] = py_path + os.pathsep + env.get("PYTHONPATH", "")
    subprocess.run(cmd, cwd=str(repo_root), env=env, check=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run dssatcalibrator in SageMaker/EKS/local Docker.")
    parser.add_argument("--input-dir", default=os.environ.get("DSSATCAL_INPUT_DIR"))
    parser.add_argument("--output-dir", default=os.environ.get("DSSATCAL_OUTPUT_DIR"))
    parser.add_argument("--config", default=os.environ.get("DSSATCAL_CONFIG"))
    parser.add_argument("--config-pattern", default=os.environ.get("DSSATCAL_CONFIG_PATTERN"))
    parser.add_argument("--dssat-dir", default=os.environ.get("DSSATCAL_DSSAT_DIR"))
    parser.add_argument("--dssat-exe", default=os.environ.get("DSSATCAL_DSSAT_EXE"))
    parser.add_argument("--source-dir", default=os.environ.get("DSSATCAL_SOURCE_DIR"))
    parser.add_argument("--template-dir", default=os.environ.get("DSSAT_TEMPLATE_DIR"))
    parser.add_argument("--language", choices=("python", "r"),
                        default=os.environ.get("DSSATCAL_LANGUAGE", "python").lower())
    parser.add_argument("--num-cores", type=int, default=(
        int(os.environ["DSSATCAL_NUM_CORES"]) if os.environ.get("DSSATCAL_NUM_CORES") else None
    ))
    parser.add_argument("--repo-root", default=os.environ.get("DSSATCAL_REPO_ROOT", "/opt/dssatcalibrator"))
    parser.add_argument("calibration_args", nargs=argparse.REMAINDER,
                        help="Arguments forwarded to run_calibration.py; prefix with --.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    default_input, default_output = default_mounts()
    input_dir = Path(args.input_dir or default_input)
    output_dir = Path(args.output_dir or default_output)
    shard = _job_index()
    input_s3 = format_indexed(os.environ.get("DSSATCAL_INPUT_S3"), shard)
    output_s3 = format_indexed(os.environ.get("DSSATCAL_OUTPUT_S3"), shard)

    sync_s3(input_s3, input_dir, direction="down")
    output_dir.mkdir(parents=True, exist_ok=True)

    config_path = find_config(input_dir, explicit=args.config,
                              pattern=args.config_pattern, index=shard)
    cloud_config, cfg = prepare_config(
        config_path,
        output_dir,
        input_dir,
        dssat_dir=Path(args.dssat_dir) if args.dssat_dir else None,
        dssat_exe=Path(args.dssat_exe) if args.dssat_exe else None,
        source_dir=Path(args.source_dir) if args.source_dir else None,
        template_dir=Path(args.template_dir) if args.template_dir else None,
        num_cores=args.num_cores,
    )
    extra_args = shlex.split(os.environ.get("DSSATCAL_ARGS", "")) + [
        a for a in args.calibration_args if a != "--"
    ]
    write_manifest(output_dir, input_dir=input_dir, config_path=config_path,
                   cloud_config_path=cloud_config, cfg=cfg, extra_args=extra_args,
                   language=args.language)
    run_calibration(Path(args.repo_root), cloud_config, extra_args, language=args.language)
    if not cfg.get("calibrator", {}).get("keep_run_dirs", False):
        shutil.rmtree(output_dir / "_workdir", ignore_errors=True)
    sync_s3(output_s3, output_dir, direction="up")


if __name__ == "__main__":
    main()

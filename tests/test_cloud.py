from pathlib import Path
from types import SimpleNamespace
import importlib.util

import yaml

from dssatcalibrator import cloud


WORKSPACE = Path(__file__).resolve().parents[2]


def _load_deploy_module(name: str):
    path = WORKSPACE / "dssatdocker" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


one_click_sagemaker = _load_deploy_module("one_click_sagemaker")
render_eks_job = _load_deploy_module("render_eks_job")


def _write_min_config(path: Path) -> None:
    path.write_text(
        """
calibrator:
  name: cloud_test
  dssat_dir: C:/DSSAT48
  dssat_exe: C:/DSSAT48/DSCSM048.EXE
  workdir: results/_workdir
  results_dir: results
  figures_dir: figures
  num_cores: 0
source:
  hemp_dir: C:/DSSAT48/Hemp
templates:
  template_dir: ""
parameters:
  genetic_cultivar:
    P1: {active: true, min: 1, max: 2, start: 1.5}
""",
        encoding="utf-8",
    )


def test_prepare_config_patches_cloud_paths(tmp_path):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    dssat_dir = input_dir / "DSSAT48"
    hemp_dir = dssat_dir / "Hemp"
    templates = input_dir / "dssat_templates"
    hemp_dir.mkdir(parents=True)
    templates.mkdir(parents=True)
    exe = dssat_dir / "dscsm048"
    exe.write_text("#!/bin/sh\n", encoding="utf-8")
    cfg_path = input_dir / "config.yaml"
    _write_min_config(cfg_path)

    cloud_path, cfg = cloud.prepare_config(cfg_path, output_dir, input_dir, num_cores=8)

    assert cloud_path == output_dir / "_cloud" / "config.cloud.yaml"
    assert cfg["calibrator"]["dssat_dir"] == str(dssat_dir)
    assert cfg["calibrator"]["dssat_exe"] == str(exe)
    assert cfg["calibrator"]["workdir"] == str(output_dir / "_workdir")
    assert cfg["calibrator"]["results_dir"] == str(output_dir / "results")
    assert cfg["calibrator"]["figures_dir"] == str(output_dir / "figures")
    assert cfg["calibrator"]["num_cores"] == 8
    assert cfg["source"]["hemp_dir"] == str(hemp_dir)
    assert cfg["templates"]["template_dir"] == str(templates)

    persisted = yaml.safe_load(cloud_path.read_text(encoding="utf-8"))
    assert persisted == cfg


def test_find_config_supports_indexed_patterns(tmp_path):
    input_dir = tmp_path / "input"
    (input_dir / "configs").mkdir(parents=True)
    shard = input_dir / "configs" / "config_003.yaml"
    shard.write_text("calibrator: {name: shard}\n", encoding="utf-8")

    found = cloud.find_config(input_dir, pattern="configs/config_{index}.yaml", index="003")

    assert found == shard


def test_format_indexed_leaves_plain_values_unchanged():
    assert cloud.format_indexed("s3://bucket/prefix", index="7") == "s3://bucket/prefix"
    assert cloud.format_indexed("s3://bucket/shard-{index}", index="7") == "s3://bucket/shard-7"


def test_render_eks_job_includes_indexed_shard_settings():
    args = SimpleNamespace(
        template=str(WORKSPACE / "dssatdocker" / "k8s" / "job.yaml.tpl"),
        job_name="dssatcal-test",
        namespace="batch",
        completions=4,
        parallelism=2,
        backoff_limit=1,
        ttl_seconds=3600,
        service_account="dssatcalibrator",
        image_uri="123.dkr.ecr.us-east-1.amazonaws.com/dssatcalibrator:latest",
        input_s3="s3://bucket/input/",
        output_s3="s3://bucket/output/shard-{index}/",
        config="config.yaml",
        config_pattern="configs/config_{index}.yaml",
        language="r",
        num_cores=8,
        calibration_args=["--", "--preset", "C", "--n", "5000"],
        cpu_request="4",
        cpu_limit="8",
        memory_request="16Gi",
        memory_limit="32Gi",
        ephemeral_size="500Gi",
    )

    manifest = render_eks_job.render(args)

    assert "completionMode: Indexed" in manifest
    assert "completions: 4" in manifest
    assert 'value: "configs/config_{index}.yaml"' in manifest
    assert 'value: "s3://bucket/output/shard-{index}/"' in manifest
    assert 'name: DSSATCAL_LANGUAGE' in manifest
    assert 'value: "r"' in manifest
    assert 'value: "--preset C --n 5000"' in manifest


def test_write_manifest_records_language(tmp_path):
    cloud.write_manifest(
        tmp_path,
        input_dir=tmp_path / "input",
        config_path=tmp_path / "input" / "config.yaml",
        cloud_config_path=tmp_path / "output" / "config.cloud.yaml",
        cfg={"calibrator": {"name": "x"}, "experiments": ["A"]},
        extra_args=["--preset", "C"],
        language="r",
    )

    manifest = yaml.safe_load((tmp_path / "cloud_manifest.json").read_text(encoding="utf-8"))
    assert manifest["language"] == "r"
    assert manifest["extra_args"] == ["--preset", "C"]


def test_one_click_s3_join_normalizes_slashes():
    assert one_click_sagemaker.s3_join("s3://bucket/prefix/", "/DSSAT48/") == "s3://bucket/prefix/DSSAT48"

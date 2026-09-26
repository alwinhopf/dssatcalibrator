"""Simulation inputs, including secondary stations, invalidate spawn reuse."""
from pathlib import Path
import pandas as pd
import pytest
from dssatcalibrator.spawn import _spawn_provenance, _spawn_outputs_complete, spawn_and_run


@pytest.mark.parametrize("changed", ["E1.MZX", "fake.exe", "DSSATPRO.V48", "E1.MZA",
                                     "Genotype/MZCER048.CUL", "Genotype/MZCER048.ECO",
                                     "Genotype/MZCER048.SPE", "Soil/CUSTOM.SOL",
                                     "Weather/SECOND.WTH", "Weather/nested/lower.wth"])
def test_spawn_provenance_tracks_inputs(tmp_path, changed):
    paths = {"root": tmp_path, "soil": tmp_path / "Soil", "weather": tmp_path / "Weather",
             "genotype": tmp_path / "Genotype"}
    for path in paths.values():
        path.mkdir(exist_ok=True)
    source = tmp_path / "E1.MZX"
    source.write_text("*EXP.DETAILS\n", encoding="utf-8")
    target = tmp_path / changed
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("original", encoding="utf-8")
    crop = {"code": "MZ", "genotype_stem": "MZCER048"}
    def provenance(cfg=None):
        return _spawn_provenance(cfg or {}, crop, [], source, paths["genotype"], paths,
                                 tmp_path / "fake.exe", [1], {})
    before = provenance()
    target.write_text("modified", encoding="utf-8")
    assert provenance() != before


def test_spawn_provenance_tracks_config_and_metadata(tmp_path):
    paths = {"root": tmp_path, "soil": tmp_path / "Soil", "weather": tmp_path / "Weather",
             "genotype": tmp_path / "Genotype"}
    for path in paths.values():
        path.mkdir(exist_ok=True)
    source = tmp_path / "E1.MZX"
    source.write_text("*EXP.DETAILS\n", encoding="utf-8")
    crop = {"code": "MZ", "genotype_stem": "MZCER048", "model": "MZCER048"}
    base_cfg = {}
    baseline = _spawn_provenance(base_cfg, crop, [], source, paths["genotype"], paths,
                                tmp_path / "fake.exe", [1], {})

    # planting dates
    pdate_cfg = {"_planting_dates": {"E1": "2021-04-15"}}
    assert _spawn_provenance(pdate_cfg, crop, [], source, paths["genotype"], paths,
                             tmp_path / "fake.exe", [1], {}) != baseline

    # filex overrides
    override_cfg = {"filex_overrides": {"all": [{"section": "SIMULATION CONTROLS", "field": "WATER", "value": "Y"}]}}
    assert _spawn_provenance(override_cfg, crop, [], source, paths["genotype"], paths,
                             tmp_path / "fake.exe", [1], {}) != baseline

    # crop model change
    other_crop = dict(crop, model="OTHER048")
    assert _spawn_provenance(base_cfg, other_crop, [], source, paths["genotype"], paths,
                             tmp_path / "fake.exe", [1], {}) != baseline


def test_spawn_outputs_complete():
    assert not _spawn_outputs_complete(pd.DataFrame(), [1])
    assert not _spawn_outputs_complete(pd.DataFrame({"LAID": [1.0]}), [1])
    assert not _spawn_outputs_complete(pd.DataFrame({"treatment": [1]}), [1, 2])
    assert not _spawn_outputs_complete(pd.DataFrame({"treatment": [1, 2]}), [1])
    assert _spawn_outputs_complete(pd.DataFrame({"treatment": [2, 1]}), [1, 2])
    assert _spawn_outputs_complete(pd.DataFrame({"treatment": [1, 1]}), [1])


def test_spawn_cache_requires_complete_output_and_manifest(tmp_path, monkeypatch):
    paths = {"root": tmp_path, "soil": tmp_path / "Soil", "weather": tmp_path / "Weather",
             "genotype": tmp_path / "Genotype"}
    for p in paths.values():
        p.mkdir(exist_ok=True)
    filex = tmp_path / "E1.MZX"
    filex.write_text("*TREATMENTS\n@N R O C TNAME\n 1 1 1 0 first\n*CULTIVARS\n", encoding="utf-8")
    exe = tmp_path / "fake.exe"
    exe.write_text("fake", encoding="utf-8")
    cfg = {"source": {"hemp_dir": str(tmp_path)}, "calibrator": {"cache_spawns": True}}
    crop = {"code": "MZ", "genotype_stem": "MZCER048", "filex_ext": "MZX", "model": "MZCER048"}

    import dssatcalibrator.spawn as sp_mod
    monkeypatch.setattr(sp_mod, "resolve_dssat_paths", lambda c: paths)
    monkeypatch.setattr(sp_mod, "_execution_backend", lambda c: "native")

    calls = {"count": 0}
    complete = {"val": True}
    def fake_run_native(run_dir, *args, **kwargs):
        calls["count"] += 1
        (run_dir / "PlantGro.OUT").write_text("new output", encoding="utf-8")
        return ""

    monkeypatch.setattr(sp_mod, "_run_native_dssat", fake_run_native)
    monkeypatch.setattr(
        sp_mod.dssat_io,
        "collect_run_outputs",
        lambda run_dir: {
            "plantgro": pd.DataFrame({"treatment": [1], "LAID": [2.0]}) if complete["val"] else pd.DataFrame(),
            "evaluate": pd.DataFrame(),
        }
    )

    def run():
        return spawn_and_run({"x": 1.0}, exp_id="E1", cfg=cfg, crop=crop,
                             param_specs=[{"name": "x", "group": "unused"}],
                             run_root=tmp_path / "runs", exe=exe)

    first = run()
    assert first.status == "success"
    assert calls["count"] == 1
    # Cached run
    second = run()
    assert second.status == "cached"
    assert calls["count"] == 1

    # Weather change invalidates cache
    (paths["weather"] / "NEW.WTH").write_text("weather change", encoding="utf-8")
    third = run()
    assert third.status == "success"
    assert third.run_dir != first.run_dir
    assert calls["count"] == 2

    # Incomplete output causes error and does not cache
    complete["val"] = False
    fourth = run()
    assert fourth.status == "error"
    assert not (third.run_dir / "spawn_manifest.json").exists()

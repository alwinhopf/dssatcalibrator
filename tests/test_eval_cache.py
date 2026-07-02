from types import SimpleNamespace

import pandas as pd

from dssatcalibrator import orchestrator
from dssatcalibrator.spaces import ParameterSpace
from dssatcalibrator.spawn import SpawnResult


def _cfg(tmp_path):
    exe = tmp_path / "dscsm048"
    exe.write_text("fake exe", encoding="utf-8")
    return {
        "calibrator": {
            "name": "cache_test",
            "workdir": str(tmp_path / "work"),
            "dssat_dir": str(tmp_path),
            "dssat_exe": str(exe),
            "num_cores": 1,
            "cache_evaluations": True,
        },
        "source": {"hemp_dir": str(tmp_path)},
        "experiments": ["E1"],
        "parameters": {
            "genetic_cultivar": {
                "P1": {"active": True, "min": 0, "max": 10, "start": 5}
            }
        },
        "crops": [{"code": "HM"}],
        "engine": {"scalar_outputs": {"HWAM": "HWAM"}},
        "objective": {"weighting": "unified", "weights": {}, "error_model": {}},
    }


def _setup(cfg, tmp_path):
    space = ParameterSpace.from_config(cfg)
    crop = {"code": "HM", "filex_ext": "HMX", "genotype_stem": "HMGRO048"}
    obs = SimpleNamespace(table=pd.DataFrame())
    return (
        space,
        crop,
        tmp_path / "dscsm048",
        space.specs,
        tmp_path / "runs",
        obs,
        ["E1"],
        {"E1": [1]},
    )


def test_evaluate_thetas_persistent_cache_dedupes_and_reuses(tmp_path, monkeypatch):
    cfg = _cfg(tmp_path)
    setup = _setup(cfg, tmp_path)
    calls = {"jobs": 0}

    def fake_run_many(jobs, n_workers, on_done=None, warmup=0):
        calls["jobs"] += len(jobs)
        out = []
        for job in jobs:
            theta = dict(job["theta"])
            ev = pd.DataFrame({
                "treatment": [1],
                "variable": ["HWAM"],
                "sim": [100.0 + float(theta["P1"])],
                "meas": [100.0],
            })
            res = SpawnResult("success", tmp_path / "run", theta, evaluate=ev)
            out.append(res)
            if on_done:
                on_done(res)
        return out

    monkeypatch.setattr(orchestrator, "run_many", fake_run_many)

    first, setup = orchestrator.evaluate_thetas(
        cfg, [{"P1": 1.0}, {"P1": 1.0}, {"P1": 2.0}], setup=setup, n_workers=1
    )
    assert calls["jobs"] == 2
    assert [r.score for r in first] == [first[0].score, first[0].score, first[2].score]

    second, _ = orchestrator.evaluate_thetas(
        cfg, [{"P1": 1.0}, {"P1": 2.0}], setup=setup, n_workers=1
    )
    assert calls["jobs"] == 2
    assert [r.score for r in second] == [first[0].score, first[2].score]
    assert list((tmp_path / "work" / "cache_test" / "evaluation_cache").glob("*/*.json"))


def test_evaluation_cache_salt_invalidates_persisted_result(tmp_path, monkeypatch):
    cfg = _cfg(tmp_path)
    setup = _setup(cfg, tmp_path)
    calls = {"jobs": 0}

    def fake_run_many(jobs, n_workers, on_done=None, warmup=0):
        calls["jobs"] += len(jobs)
        return [
            SpawnResult(
                "success",
                tmp_path / "run",
                dict(job["theta"]),
                evaluate=pd.DataFrame({
                    "treatment": [1],
                    "variable": ["HWAM"],
                    "sim": [101.0],
                    "meas": [100.0],
                }),
            )
            for job in jobs
        ]

    monkeypatch.setattr(orchestrator, "run_many", fake_run_many)
    orchestrator.evaluate_thetas(cfg, [{"P1": 1.0}], setup=setup, n_workers=1)
    orchestrator.evaluate_thetas(cfg, [{"P1": 1.0}], setup=setup, n_workers=1)
    assert calls["jobs"] == 1

    cfg2 = _cfg(tmp_path)
    cfg2["calibrator"]["evaluation_cache_salt"] = "fresh"
    setup2 = _setup(cfg2, tmp_path)
    orchestrator.evaluate_thetas(cfg2, [{"P1": 1.0}], setup=setup2, n_workers=1)
    assert calls["jobs"] == 2

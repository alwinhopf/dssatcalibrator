"""Offline regressions for the report-oriented cross-validation layer."""

import json
from types import SimpleNamespace

import pandas as pd
import pytest

from dssatcalibrator import cv
from dssatcalibrator.objective import ObjectiveResult


def _objective(score, rmse):
    return ObjectiveResult(
        score=float(score),
        loglik=-0.5 * float(score),
        residuals=pd.DataFrame(),
        per_var={"yield": {"RMSE": float(rmse), "d": 0.9}},
        per_exp_var=pd.DataFrame(),
    )


def test_cross_validation_uses_top_level_experiment_subsets(monkeypatch, tmp_path):
    cfg = {
        "experiments": ["E1", "E2", "E3"],
        "calibrator": {"seed": 7, "num_cores": 1},
        "cross_validation": {
            "enabled": True,
            "strategy": "k_fold",
            "k": 3,
            "final_theta": "best_fold",
        },
    }
    seen_train, seen_test = [], []

    monkeypatch.setattr(
        cv, "_setup",
        lambda work_cfg: (None, None, None, None, None, None,
                          list(work_cfg["experiments"]), None),
    )

    def fake_calibrate(work_cfg, progress=True):
        exps = list(work_cfg["experiments"])
        seen_train.append(exps)
        return SimpleNamespace(
            best_theta={"P1": float(len(exps))},
            best=_objective(score=len(exps), rmse=len(exps)),
        )

    def fake_evaluate(work_cfg, thetas):
        exps = list(work_cfg["experiments"])
        seen_test.append(exps)
        return [_objective(score=10 + len(exps), rmse=2 * len(exps))], None

    monkeypatch.setattr(cv.orchestrator, "calibrate", fake_calibrate)
    monkeypatch.setattr(cv.orchestrator, "evaluate_thetas", fake_evaluate)

    result = cv.run_cross_validation(cfg, progress=False)

    assert len(result.folds) == 3
    assert sorted(e for fold in seen_test for e in fold) == ["E1", "E2", "E3"]
    assert all(len(fold) == 2 for fold in seen_train)
    assert all(set(train).isdisjoint(test)
               for train, test in zip(seen_train, seen_test))
    assert all(set(train) | set(test) == set(cfg["experiments"])
               for train, test in zip(seen_train, seen_test))
    assert "experiments" not in cfg["calibrator"]
    assert result.final_theta_method == "best_fold"
    assert list(result.report_df["fold"]) == [1, 2, 3]
    paths = cv.write_cv_report(result, tmp_path)
    assert [path.name for path in paths] == ["cv_report.csv", "cv_summary.json"]
    assert pd.read_csv(paths[0])["fold"].tolist() == [1, 2, 3]
    summary = json.loads(paths[1].read_text(encoding="utf-8"))
    assert summary["recommendation"] == result.summary["recommendation"]

    seen_train.clear()
    seen_test.clear()
    temporal_cfg = {
        **cfg,
        "experiments": ["E2001", "E2101", "E2201"],
        "cross_validation": {
            **cfg["cross_validation"],
            "strategy": "temporal_forward",
        },
    }

    temporal = cv.run_cross_validation(temporal_cfg, progress=False)

    assert len(temporal.folds) == 2
    assert seen_test == [["E2101"], ["E2201"]]
    assert seen_train == [["E2001"], ["E2001", "E2101"]]


def test_cross_validation_config_is_validated():
    with pytest.raises(ValueError, match="at least 2"):
        cv.parse_cv_config({"cross_validation": {"strategy": "k_fold", "k": 1}})
    with pytest.raises(ValueError, match="Unknown cross-validation strategy"):
        cv.parse_cv_config({"cross_validation": {"strategy": "not-a-strategy"}})
    with pytest.raises(ValueError, match="not enabled"):
        cv.run_cross_validation({"cross_validation": {"enabled": False}}, progress=False)

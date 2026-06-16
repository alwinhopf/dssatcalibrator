"""Tests for combining completed calibration runs."""
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from dssatcalibrator import orchestrator
from dssatcalibrator.config import load_config
from dssatcalibrator.objective import ObjectiveResult
from dssatcalibrator.spaces import ParameterSpace

REPO = Path(__file__).resolve().parents[1]
HEMP_CFG = REPO / "config_hemp.yaml"


def test_combine_runs_synthetic(tmp_path):
    dir1 = tmp_path / "run1"
    dir2 = tmp_path / "run2"
    dir1.mkdir()
    dir2.mkdir()

    # design1
    design1 = pd.DataFrame({
        "sample_id": [0, 1],
        "CSDL": [13.0, 14.0],
        "PPSEN": [0.5, 0.6],
        "score": [1.5, 2.5],
        "loglik": [-1.5, -2.5],
        "n_obs": [10, 10],
        "weight": [0.73, 0.27]
    })
    design1.to_csv(dir1 / "design.csv", index=False)

    # design2 (one duplicate, one new best sample)
    design2 = pd.DataFrame({
        "sample_id": [0, 1],
        "CSDL": [14.0, 15.0],
        "PPSEN": [0.6, 0.7],
        "score": [2.5, 0.5],
        "loglik": [-2.5, -0.5],
        "n_obs": [10, 10],
        "weight": [0.12, 0.88]
    })
    design2.to_csv(dir2 / "design.csv", index=False)

    cfg = load_config(HEMP_CFG)
    # Simplify parameter configuration for testing
    cfg["parameters"] = {
        "genetic_cultivar": {
            "CSDL":  { "active": True,  "role": "obligatory", "min": 12.0, "max": 16.0, "start": 12.80 },
            "PPSEN": { "active": True,  "role": "obligatory", "min": 0.20, "max": 1.00, "start": 0.90 }
        }
    }
    cfg["experiments"] = ["YUKU2101"]

    mock_best = ObjectiveResult(score=0.5, loglik=-0.5, residuals=pd.DataFrame(), per_var={}, per_exp_var=pd.DataFrame())

    with patch("dssatcalibrator.orchestrator._score_theta", return_value=mock_best):
        with patch("dssatcalibrator.orchestrator._setup") as mock_setup:
            space = ParameterSpace.from_config(cfg)
            mock_setup.return_value = (space, {}, Path("."), [], Path("."), None, ["YUKU2101"], {"YUKU2101": []})

            res = orchestrator.combine_runs(cfg, [dir1, dir2])

            assert res.best_theta == {"CSDL": 15.0, "PPSEN": 0.7}
            # Combined length should be 3 (1 duplicate removed)
            assert len(res.design) == 3
            assert list(res.design["sample_id"]) == [0, 1, 2]
            assert res.best.score == 0.5

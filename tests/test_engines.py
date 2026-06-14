"""Tests for calibration engines: GLUE (fast/synthetic) + NSGA-II & validation (tiny, real)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from dssatcalibrator.config import load_config
from dssatcalibrator.engines import run_glue

REPO = Path(__file__).resolve().parents[1]
HEMP_CFG = REPO / "config_hemp.yaml"


def test_run_glue_synthetic():
    # 5 parameter sets, scores increasing -> best is sample 0
    design = pd.DataFrame({
        "sample_id": range(5),
        "P1": [1, 2, 3, 4, 5], "P2": [5, 4, 3, 2, 1],
        "score": [0.1, 0.5, 1.0, 2.0, 5.0],
        "loglik": [-1.0, -3.0, -6.0, -10.0, -20.0],
    })
    cfg = {"method": {"bayesian": {"behavioural_quantile": 0.4}}}
    g = run_glue(design, ["P1", "P2"], cfg)
    assert g.best_sample_id == 0
    assert g.best_theta == {"P1": 1.0, "P2": 5.0}
    assert abs(g.design["weight"].sum() - 1.0) < 1e-9
    assert g.design["weight"].iloc[0] == g.design["weight"].max()  # best has most weight
    assert len(g.behavioural) >= 2 and g.ess > 0


@pytest.mark.slow
def test_nsga2_smoke(hemp_dir):
    cfg = load_config(HEMP_CFG)
    cfg["experiments"] = ["YUKU2101"]
    cfg["calibrator"]["num_cores"] = 4
    cfg["method"]["sample"]["n"] = 2
    cfg["method"]["bayesian"] = {"engine": "glue", "behavioural_quantile": 0.1}
    cfg["method"]["multiobjective"] = {"engine": "nsga2", "variables": ["biomass", "height"],
                                       "pop_size": 4, "n_gen": 1}
    from dssatcalibrator import orchestrator
    res = orchestrator.calibrate(cfg, progress=False)
    assert res.nsga2 is not None
    front = res.nsga2.front()
    assert len(front) >= 1
    assert {"nRMSE_biomass", "nRMSE_height"}.issubset(front.columns)


@pytest.mark.slow
def test_validation_loeo_smoke(hemp_dir):
    cfg = load_config(HEMP_CFG)
    cfg["experiments"] = ["YUKU2101", "YUFE2201"]
    cfg["calibrator"]["num_cores"] = 4
    cfg["method"]["sample"]["n"] = 2
    cfg["method"]["bayesian"] = {"engine": "glue", "behavioural_quantile": 0.1}
    from dssatcalibrator import orchestrator
    df = orchestrator.validate_loeo(cfg, progress=False)
    assert not df.empty
    assert set(df["split"]) == {"calibration", "evaluation"}
    assert set(df["held_out"]) == {"YUKU2101", "YUFE2201"}

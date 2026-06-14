"""Tests for the SMC Particle Filter engine (run_smc_pf) using synthetic mocks."""
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from dssatcalibrator.engines.smc_pf import run_smc_pf
from dssatcalibrator.spawn import SpawnResult
from dssatcalibrator.spaces import ParameterSpace


def mock_setup(cfg):
    specs = [
        {"group": "genetic_cultivar", "name": "P1", "min": 1.0, "max": 10.0, "start": 5.0, "active": True},
        {"group": "genetic_cultivar", "name": "P2", "min": 10.0, "max": 100.0, "start": 50.0, "active": True}
    ]
    space = ParameterSpace(
        names=["P1", "P2"],
        low=np.array([1.0, 10.0]),
        high=np.array([10.0, 100.0]),
        start=np.array([5.0, 50.0]),
        specs=specs
    )
    
    crop = {"genotype_stem": "HMGRO048", "cultivar_anchor": "IB0008", "filex_ext": "HMX", "code": "HM", "model": "CRGRO"}
    exe = Path("fake_exe")
    run_root = Path("fake_root")
    
    from dssatcalibrator.observations import Observations
    obs_df = pd.DataFrame([
        # exp_id | treatment | variable | kind | date | value | sigma | weight
        ("EXP1", 1, "biomass", "timeseries", pd.Timestamp("2021-06-20"), 100.0, 10.0, 1.0),
        ("EXP1", 1, "biomass", "timeseries", pd.Timestamp("2021-06-25"), 200.0, 20.0, 1.0),
        ("EXP1", 1, "grain_yield", "scalar", pd.NaT, 1500.0, 150.0, 1.0),
    ], columns=["exp_id", "treatment", "variable", "kind", "date", "value", "sigma", "weight"])
    obs = Observations(obs_df)
    
    experiments = ["EXP1"]
    treatments = {"EXP1": [1]}
    return space, crop, exe, specs, run_root, obs, experiments, treatments


def mock_spawn_and_run(theta, exp_id, cfg, crop, param_specs, run_root, treatments, exe):
    p1 = theta.get("P1", 5.0)
    p2 = theta.get("P2", 50.0)
    
    pg = pd.DataFrame([
        {"treatment": 1, "date": pd.Timestamp("2021-06-20"), "CWAD": p1 * 20.0},
        {"treatment": 1, "date": pd.Timestamp("2021-06-25"), "CWAD": p1 * 40.0},
    ])
    
    ev = pd.DataFrame([
        {"treatment": 1, "variable": "HWAM", "sim": p2 * 30.0, "meas": 1500.0}
    ])
    
    return SpawnResult(status="success", run_dir=Path("fake_dir"), theta=theta, plantgro=pg, evaluate=ev)


@patch("dssatcalibrator.orchestrator._setup", side_effect=mock_setup)
@patch("dssatcalibrator.runner.spawn_and_run", side_effect=mock_spawn_and_run)
def test_run_smc_pf_synthetic(mock_setup_fn, mock_spawn_fn):
    cfg = {
        "calibrator": {
            "name": "test_smc",
            "seed": 42,
            "workdir": "fake_workdir",
            "dssat_dir": "fake_dssat",
            "num_cores": 1
        },
        "method": {
            "sample": {"engine": "lhs", "n": 8},
            "bayesian": {"engine": "smc_pf", "n_particles": 8, "ess_frac": 0.8, "mutation_scale": 0.05, "behavioural_quantile": 0.5}
        },
        "engine": {
            "timeseries_outputs": {"biomass": "CWAD"},
            "scalar_outputs": {"grain_yield": "HWAM"}
        },
        "objective": {
            "weighting": "unified",
            "weights": {"biomass": 1.0, "grain_yield": 1.0},
            "error_model": {
                "biomass": {"type": "absolute", "value": 10.0},
                "grain_yield": {"type": "absolute", "value": 150.0}
            }
        }
    }
    
    result = run_smc_pf(cfg, progress=False)
    
    assert result.best_theta is not None
    assert "P1" in result.best_theta
    assert "P2" in result.best_theta
    
    # The optimal values should yield simulated values close to observed.
    # EXP1 has observed biomass 100 on 2021-06-20 (sim = P1 * 20.0 => P1 = 5.0)
    # and grain_yield 1500 (sim = P2 * 30.0 => P2 = 50.0)
    assert abs(result.best_theta["P1"] - 5.0) < 1.0
    assert abs(result.best_theta["P2"] - 50.0) < 5.0
    
    assert result.ess > 0
    assert len(result.design) == 9
    assert "weight" in result.design.columns
    assert "score" in result.design.columns
    
    # Verify that the best particle indeed has the minimum score
    best_idx = result.best_sample_id
    assert result.design.loc[best_idx, "score"] == result.design["score"].min()

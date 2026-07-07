"""End-to-end pipeline wiring tests (no DSSAT binary).

These mock the DSSAT spawn + experiment setup exactly like ``test_smc_pf.py`` but
drive the *full* ``orchestrator.calibrate()`` dispatch, so they verify the preset
routing and every estimator branch (glue / smc_pf / mcmc / optimizer) hang together.
"""
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from dssatcalibrator import orchestrator
from dssatcalibrator.observations import Observations
from dssatcalibrator.spaces import ParameterSpace
from dssatcalibrator.spawn import SpawnResult


def mock_setup(cfg):
    specs = [
        {"group": "genetic_cultivar", "name": "P1", "min": 1.0, "max": 10.0, "start": 5.0, "active": True},
        {"group": "genetic_cultivar", "name": "P2", "min": 10.0, "max": 100.0, "start": 50.0, "active": True},
    ]
    space = ParameterSpace(names=["P1", "P2"], low=np.array([1.0, 10.0]),
                           high=np.array([10.0, 100.0]), start=np.array([5.0, 50.0]), specs=specs)
    crop = {"genotype_stem": "HMGRO048", "cultivar_anchor": "IB0008", "filex_ext": "HMX",
            "code": "HM", "model": "CRGRO"}
    obs_df = pd.DataFrame([
        ("EXP1", 1, "biomass", "timeseries", pd.Timestamp("2021-06-20"), 100.0, 10.0, 1.0),
        ("EXP1", 1, "biomass", "timeseries", pd.Timestamp("2021-06-25"), 200.0, 20.0, 1.0),
        ("EXP1", 1, "grain_yield", "scalar", pd.NaT, 1500.0, 150.0, 1.0),
    ], columns=["exp_id", "treatment", "variable", "kind", "date", "value", "sigma", "weight"])
    return space, crop, Path("fake"), specs, Path("fake"), Observations(obs_df), ["EXP1"], {"EXP1": [1]}


def mock_spawn(theta, exp_id, cfg, crop, param_specs, run_root, treatments, exe):
    p1, p2 = theta.get("P1", 5.0), theta.get("P2", 50.0)
    pg = pd.DataFrame([{"treatment": 1, "date": pd.Timestamp("2021-06-20"), "CWAD": p1 * 20.0},
                       {"treatment": 1, "date": pd.Timestamp("2021-06-25"), "CWAD": p1 * 40.0}])
    ev = pd.DataFrame([{"treatment": 1, "variable": "HWAM", "sim": p2 * 30.0, "meas": 1500.0}])
    return SpawnResult(status="success", run_dir=Path("fake"), theta=theta, plantgro=pg, evaluate=ev)


BASE_CFG = {
    "calibrator": {
        "name": "t",
        "seed": 1,
        "workdir": "fake",
        "dssat_dir": "fake",
        "num_cores": 1,
        "cache_evaluations": False,
    },
    "engine": {"timeseries_outputs": {"biomass": "CWAD"}, "scalar_outputs": {"grain_yield": "HWAM"}},
    "objective": {"weighting": "unified", "weights": {"biomass": 1.0, "grain_yield": 1.0},
                  "error_model": {"biomass": {"type": "absolute", "value": 10.0},
                                  "grain_yield": {"type": "absolute", "value": 150.0}}},
}

METHODS = {
    "glue": {"preset": "C", "sample": {"engine": "lhs", "n": 6},
             "bayesian": {"engine": "glue", "behavioural_quantile": 0.5}},
    "smc_pf": {"preset": "A", "sample": {"engine": "lhs", "n": 6},
               "bayesian": {"engine": "smc_pf", "n_particles": 6, "ess_frac": 0.8,
                            "behavioural_quantile": 0.5}},
    "mcmc": {"bayesian": {"engine": "mcmc", "n_walkers": 4, "n_steps": 20, "burn_in": 10,
                          "proposal_scale": 0.2, "behavioural_quantile": 0.5}},
    "optimizer": {"bayesian": {"engine": "none"},
                  "optimizer": {"engine": "nelder_mead", "restarts": 1, "maxiter": 15}},
}


def _calibrate(method_block):
    cfg = {**BASE_CFG, "method": method_block}
    with patch("dssatcalibrator.orchestrator._setup", side_effect=mock_setup), \
         patch("dssatcalibrator.runner.spawn_and_run", side_effect=mock_spawn):
        return orchestrator.calibrate(cfg, progress=False)


def _check(res):
    assert set(res.best_theta) == {"P1", "P2"}
    assert res.design is not None and len(res.design) >= 1
    assert np.isfinite(res.best.score)


def test_pipeline_glue():
    _check(_calibrate(METHODS["glue"]))


def test_pipeline_smc_pf():
    _check(_calibrate(METHODS["smc_pf"]))


def test_pipeline_mcmc():
    res = _calibrate(METHODS["mcmc"])
    _check(res)
    assert (res.extras or {}).get("engine") == "mcmc"


def test_pipeline_optimizer():
    res = _calibrate(METHODS["optimizer"])
    _check(res)
    assert (res.extras or {}).get("engine") == "optimizer"

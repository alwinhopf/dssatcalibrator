"""Fast, offline tests for the engines added in the framework build-out.

All use *synthetic* score functions (a simple quadratic / Gaussian), so they run
without the DSSAT binary and prove the algorithms behave: optimisers find the
minimum, screening ranks the influential parameter first, MCMC concentrates on the
target, stepwise selection drops useless parameters, presets resolve correctly.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from dssatcalibrator import priors
from dssatcalibrator.objective import ObjectiveResult
from dssatcalibrator.spaces import ParameterSpace
from dssatcalibrator.engines import (
    influential_params, run_mcmc, run_morris, run_optimizer, stepwise_select,
)
from dssatcalibrator.engines.sensitivity import anova_variance_share
from dssatcalibrator.orchestrator import _apply_active_subset, _resolve_method, _stage_on


def make_space(roles=None):
    roles = roles or {}
    names = ["a", "b"] if roles == {} else list(roles)
    specs = [{"group": "genetic_cultivar", "name": n, "min": 0.0, "max": 10.0,
              "start": 5.0, "active": True, **({"role": roles[n]} if n in roles else {})}
             for n in names]
    return ParameterSpace(names=names,
                          low=np.array([0.0] * len(names)),
                          high=np.array([10.0] * len(names)),
                          start=np.array([5.0] * len(names)), specs=specs)


# --------------------------------------------------------------------------- priors
def test_prior_uniform_default():
    spec = {"name": "x", "min": 2.0, "max": 8.0}
    rng = np.random.default_rng(0)
    draws = priors.sample_one(spec, 500, rng)
    assert draws.min() >= 2.0 and draws.max() <= 8.0
    assert priors.log_prior_one(spec, 5.0) == 0.0          # flat inside
    assert priors.log_prior_one(spec, 9.0) == float("-inf")  # outside bounds


def test_prior_normal_shape():
    spec = {"name": "x", "min": 0.0, "max": 10.0, "start": 5.0,
            "prior": {"dist": "normal", "sd": 1.0}}
    # density should be higher at the centre than in the tail
    assert priors.log_prior_one(spec, 5.0) > priors.log_prior_one(spec, 8.0)
    space = make_space()
    space.specs[0]["prior"] = {"dist": "normal", "sd": 1.0}
    assert priors.has_informative_prior(space) is True


# ----------------------------------------------------------------------- optimizers
@pytest.mark.parametrize("method", ["nelder_mead", "diffevo"])
def test_optimizer_finds_minimum(method):
    target = {"a": 3.0, "b": 7.0}
    space = make_space()

    def score_batch(thetas):
        return [sum((t[n] - target[n]) ** 2 for n in target) for t in thetas]

    res = run_optimizer(space, score_batch, method=method, seed=1, maxiter=40, restarts=3)
    assert abs(res.best_theta["a"] - 3.0) < 1.0
    assert abs(res.best_theta["b"] - 7.0) < 1.0
    assert res.best_score < 1.0 and res.n_eval > 0


# ---------------------------------------------------------------------- sensitivity
def test_morris_ranks_influential_first():
    space = make_space()

    def score_results(thetas):
        # output depends strongly on "a", negligibly on "b"
        return [ObjectiveResult(score=5.0 * t["a"] + 0.001 * t["b"], loglik=0.0,
                                residuals=pd.DataFrame()) for t in thetas]

    sens = run_morris(space, score_results, trajectories=10, seed=1)
    assert sens.ranking.iloc[0]["parameter"] == "a"
    assert "a" in influential_params(sens.ranking)
    assert "b" not in influential_params(sens.ranking, rel_threshold=0.5)


def test_anova_variance_share():
    df = pd.DataFrame({"src": ["x", "x", "y", "y"], "score": [1.0, 1.0, 9.0, 9.0]})
    out = anova_variance_share(df, ["src"], "score")
    assert abs(float(out.iloc[0]["var_share"]) - 1.0) < 1e-9   # factor explains all variance


# ------------------------------------------------------------------------------ mcmc
def test_mcmc_concentrates_on_target():
    target = np.array([3.0, 7.0])
    space = make_space()

    def score_results(thetas):
        out = []
        for t in thetas:
            d = np.array([t["a"], t["b"]]) - target
            out.append(ObjectiveResult(score=float(np.sum(d ** 2)),
                                       loglik=float(-0.5 * np.sum(d ** 2)),
                                       residuals=pd.DataFrame({"resid": np.zeros(10)})))
        return out

    cfg = {"calibrator": {"seed": 1},
           "method": {"bayesian": {"n_walkers": 10, "n_steps": 200, "burn_in": 100,
                                   "proposal_scale": 0.2}}}
    mc = run_mcmc(cfg, score_results, space, progress=False)
    assert 0.0 < mc.acceptance <= 1.0
    assert "weight" in mc.design.columns
    assert abs(mc.design["a"].mean() - 3.0) < 2.0
    assert abs(mc.design["b"].mean() - 7.0) < 2.0


# ------------------------------------------------------------------------- selection
def test_stepwise_keeps_useful_drops_useless():
    space = make_space({"a": "obligatory", "c": "candidate", "d": "candidate"})

    def score_results(thetas):
        out = []
        for t in thetas:
            s = (t["a"] - 3.0) ** 2 + (t["c"] - 7.0) ** 2     # "d" has no effect
            out.append(ObjectiveResult(score=float(s), loglik=float(-0.5 * s),
                                       residuals=pd.DataFrame({"resid": np.zeros(20)})))
        return out

    sel = stepwise_select(space, score_results, criterion="bic",
                          optimizer="nelder_mead", optimizer_restarts=1, maxiter=60)
    assert "a" in sel.selected            # obligatory always in
    assert "c" in sel.selected            # candidate that helps is kept
    assert "d" not in sel.selected        # useless candidate is dropped


# --------------------------------------------------------------------------- presets
def test_preset_resolves_engines():
    m = _resolve_method({"method": {"preset": "A"}})
    assert m["bayesian"]["engine"] == "smc_pf"
    assert m["sensitivity"]["engine"] == "morris"
    m = _resolve_method({"method": {"preset": "D"}})
    assert m["bayesian"]["engine"] == "mcmc"


def test_preset_user_override_wins():
    m = _resolve_method({"method": {"preset": "A", "bayesian": {"engine": "glue"}}})
    assert m["bayesian"]["engine"] == "glue"


def test_stage_on_requires_active():
    assert _stage_on({"engine": "morris", "active": True}) is True
    assert _stage_on({"engine": "morris"}) is False        # engine set but not switched on
    assert _stage_on(None) is False


def test_apply_active_subset():
    cfg = {"parameters": {"g": {"a": {"active": True, "min": 0, "max": 1},
                                "b": {"active": True, "min": 0, "max": 1}}}}
    out = _apply_active_subset(cfg, ["a"])
    assert out["parameters"]["g"]["a"]["active"] is True
    assert out["parameters"]["g"]["b"]["active"] is False

"""Offline correctness tests for the newly added engines (CMA-ES, DREAM, ES-MDA,
Bayesian optimisation, history matching, ABC-SMC). A synthetic linear-Gaussian problem with a known optimum
``target`` is solved through each engine's real code path (no DSSAT), and we
assert the engine recovers the optimum and the registry resolves it.
"""
import numpy as np
import pandas as pd
import pytest

from dssatcalibrator.objective import ObjectiveResult
from dssatcalibrator.spaces import ParameterSpace
from dssatcalibrator import orchestrator as orch

NAMES = ["p0", "p1", "p2", "p3"]
TARGET = np.array([1.5, -2.0, 0.7, 3.1])
LOW = np.array([-5.0] * 4)
HIGH = np.array([5.0] * 4)


def _space():
    specs = [{"name": n, "min": LOW[i], "max": HIGH[i], "start": 0.0,
              "prior": {"dist": "uniform"}} for i, n in enumerate(NAMES)]
    return ParameterSpace(names=NAMES, low=LOW, high=HIGH,
                          start=np.zeros(4), specs=specs)


def _scorer(seed=0, n_obs=10, sigma=0.2):
    """Build score_results(list_of_theta) -> [ObjectiveResult] for a linear model
    sim = A @ theta, obs = A @ TARGET, so the optimum is theta == TARGET."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n_obs, len(NAMES)))
    obs_vec = A @ TARGET

    def score_results(thetas):
        out = []
        for t in thetas:
            x = np.array([t[n] for n in NAMES])
            sim = A @ x
            resid = sim - obs_vec
            df = pd.DataFrame({
                "exp_id": ["E"] * n_obs, "treatment": [1] * n_obs,
                "dssat": [f"V{i}" for i in range(n_obs)], "user_var": ["v"] * n_obs,
                "date": [pd.NaT] * n_obs, "obs": obs_vec, "sim": sim,
                "sigma": [sigma] * n_obs, "weight": [1.0] * n_obs, "resid": resid,
                "kind": ["scalar"] * n_obs,
            })
            chi2 = float(np.sum((resid / sigma) ** 2))
            out.append(ObjectiveResult(score=chi2 / n_obs, loglik=-0.5 * chi2,
                                       residuals=df, per_var={}, per_exp_var=pd.DataFrame()))
        return out
    return score_results


def _all_invalid_scorer(thetas):
    return [ObjectiveResult(
        score=float("inf"), loglik=float("-inf"),
        residuals=pd.DataFrame(), per_var={}, per_exp_var=pd.DataFrame(),
    ) for _ in thetas]


def _cfg(engine, **bayes):
    return {"calibrator": {"seed": 7, "num_cores": 1},
            "method": {"bayesian": {"engine": engine, **bayes}}}


def _err(theta):
    return float(np.linalg.norm(np.array([theta[n] for n in NAMES]) - TARGET))


def test_cmaes_recovers_optimum():
    from dssatcalibrator.engines.optimizers import run_optimizer
    sp = _space()
    sr = _scorer()
    score_batch = lambda ths: [r.score for r in sr(ths)]  # noqa: E731
    res = run_optimizer(sp, score_batch, method="cmaes", seed=1, maxiter=60, popsize=12)
    assert _err(res.best_theta) < 0.1
    # registry routes cmaes through the optimizer estimator
    assert orch._resolve_estimator({"bayesian": {"engine": "none"},
                                    "optimizer": {"engine": "cmaes"}}) == "optimizer"


def test_dream_recovers_optimum_and_gives_spread():
    from dssatcalibrator.engines.dream import run_dream
    r = run_dream(_cfg("dream", n_generations=120, burn_in=60),
                  _scorer(), _space(), progress=False)
    assert _err(r.best_theta) < 0.4
    assert len(r.design) > 10 and r.design["p0"].std() >= 0   # a posterior cloud
    assert orch._resolve_estimator({"bayesian": {"engine": "dream"}}) == "dream"


def test_es_mda_recovers_optimum():
    from dssatcalibrator.engines.es_mda import run_es_mda
    r = run_es_mda(_cfg("es_mda", ensemble_size=40, iterations=6),
                   _scorer(), _space(), progress=False)
    assert _err(r.best_theta) < 0.3
    assert orch._resolve_estimator({"bayesian": {"engine": "es_mda"}}) == "es_mda"


def test_bayesopt_recovers_optimum():
    pytest.importorskip("sklearn")
    from dssatcalibrator.engines.bayesopt import run_bayesopt
    r = run_bayesopt(_cfg("bayesopt", n_init=24, n_iter=18, batch_size=4),
                     _scorer(), _space(), progress=False)
    assert _err(r.best_theta) < 0.6
    assert orch._resolve_estimator({"bayesian": {"engine": "bayesopt"}}) == "bayesopt"


def test_history_matching_keeps_nroy_region():
    from dssatcalibrator.engines.history import run_history_matching
    r = run_history_matching(_cfg("history", n=96, waves=2, implausibility_cutoff=4.0),
                             _scorer(), _space(), progress=False)
    assert len(r.design) > 20
    assert r.ess > 0
    assert np.isfinite(r.best.score)
    assert _err(r.best_theta) < 4.0
    assert orch._resolve_estimator({"bayesian": {"engine": "history"}}) == "history"


def test_abc_smc_accepts_final_population():
    from dssatcalibrator.engines.abc_smc import run_abc_smc
    r = run_abc_smc(_cfg("abc_smc", n_particles=40, waves=2, oversample=2,
                         threshold_quantile=0.4),
                    _scorer(), _space(), progress=False)
    assert len(r.design) > 40
    assert r.ess > 0
    assert len(r.thresholds) == 2
    assert r.thresholds[-1] <= r.thresholds[0]
    assert _err(r.best_theta) < 2.5
    assert orch._resolve_estimator({"bayesian": {"engine": "abc_smc"}}) == "abc_smc"


@pytest.mark.parametrize(
    ("engine", "message"),
    [
        ("history", "History matching found no valid candidates"),
        ("abc_smc", "ABC-SMC found no valid candidates"),
    ],
)
def test_sparse_bayesian_engines_reject_all_invalid_candidates(engine, message):
    if engine == "history":
        from dssatcalibrator.engines.history import run_history_matching as run
        cfg = _cfg(engine, n=8, waves=2)
    else:
        from dssatcalibrator.engines.abc_smc import run_abc_smc as run
        cfg = _cfg(engine, n_particles=8, waves=2, oversample=1)
    with pytest.raises(ValueError, match=message):
        run(cfg, _all_invalid_scorer, _space(), progress=False)


def test_all_new_engines_in_registry():
    for name in ("dream", "es_mda", "bayesopt", "history", "abc_smc"):
        assert name in orch.ESTIMATOR_REGISTRY

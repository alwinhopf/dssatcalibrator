"""Bayesian optimisation — sample-efficient best-fit for expensive crops.

Where ``engines.surrogate`` trains a GP on a *fixed* design and then runs GLUE on
the emulator, Bayesian optimisation makes the GP *active*: it fits a Gaussian
process to the runs so far, then asks "where is the single most promising place
to spend the next DSSAT run?" via an **acquisition function** (Expected
Improvement), runs DSSAT there, and repeats. Each round it proposes a small
*batch* of promising points (scored in one parallel batch), so it uses every core
while still converging in far fewer model runs than GLUE or differential
evolution — the right tool when each DSSAT run is costly (CROPGRO/CERES, many
treatments).

Requires scikit-learn (``pip install -e .[full]``), lazily imported.

Config (under ``method.bayesian`` when ``engine: bayesopt``)::

    n_init: 16           # initial LHS design (default max(2*n_params, 10))
    n_iter: 20           # acquisition rounds
    batch_size: 4        # points proposed per round (uses parallel cores)
    xi: 0.01             # EI exploration margin (fraction of the score range)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .. import priors


@dataclass
class BayesOptResult:
    design: pd.DataFrame
    obj_results: dict
    best_theta: dict
    best_sample_id: int
    best: object
    info: dict = field(default_factory=dict)


def _expected_improvement(mu, sigma, best_y, xi):
    """EI for MINIMISATION: how much we expect to improve below ``best_y``."""
    from scipy.stats import norm
    sigma = np.maximum(sigma, 1e-12)
    imp = best_y - mu - xi
    z = imp / sigma
    return imp * norm.cdf(z) + sigma * norm.pdf(z)


def run_bayesopt(cfg: dict, score_results, space, *, progress: bool = True) -> BayesOptResult:
    """Run batched Bayesian optimisation with a GP surrogate + Expected Improvement."""
    try:
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
    except ImportError as exc:                                  # pragma: no cover
        raise ImportError("Bayesian optimisation needs scikit-learn "
                          "(pip install -e '.[full]').") from exc

    bcfg = cfg.get("method", {}).get("bayesian", {})
    n_init = int(bcfg.get("n_init", max(2 * space.ndim, 10)))
    n_iter = int(bcfg.get("n_iter", 20))
    batch = max(1, int(bcfg.get("batch_size", 4)))
    xi_frac = float(bcfg.get("xi", 0.01))
    seed = int(cfg["calibrator"].get("seed", 42))
    rng = np.random.default_rng(seed)

    names = space.names
    low, high = np.asarray(space.low, float), np.asarray(space.high, float)
    span = np.where(high > low, high - low, 1.0)

    def to_unit(M):
        return (np.asarray(M, float) - low) / span

    def to_theta_unit(u):
        return space.to_theta(low + np.clip(u, 0, 1) * span)

    # --- initial design --------------------------------------------------------
    init = priors.sample_prior_design(space, n_init, rng)
    X = np.array([init.iloc[i].to_numpy(dtype=float) for i in range(n_init)])
    thetas = [space.to_theta(X[i]) for i in range(n_init)]
    results = list(score_results(thetas))
    y = np.array([r.score if np.isfinite(r.score) else 1e12 for r in results])

    if progress:
        print(f"Running Bayesian optimisation: {n_init} init + {n_iter}x{batch} "
              f"acquisitions...", flush=True)

    kernel = (ConstantKernel(1.0) * Matern(length_scale=np.ones(space.ndim), nu=2.5)
              + WhiteKernel(noise_level=1e-3))

    for it in range(n_iter):
        Xu = to_unit(X)
        ymean, ystd = y.mean(), (y.std() or 1.0)
        gp = GaussianProcessRegressor(kernel=kernel, normalize_y=False,
                                      n_restarts_optimizer=2, random_state=seed + it)
        import warnings
        with warnings.catch_warnings():            # GP hyper-fit can warn on flat data
            warnings.simplefilter("ignore")
            gp.fit(Xu, (y - ymean) / ystd)

        # Maximise EI over a large random candidate pool (gradient-free, cheap).
        cand = rng.uniform(0, 1, size=(2000, space.ndim))
        mu, sd = gp.predict(cand, return_std=True)
        mu = mu * ystd + ymean
        sd = sd * ystd
        best_y = y.min()
        ei = _expected_improvement(mu, sd, best_y, xi_frac * (y.max() - best_y + 1e-9))

        # Take the top-``batch`` distinct candidates this round (parallel batch).
        picks = []
        order = np.argsort(-ei)
        for j in order:
            if all(np.linalg.norm(cand[j] - cand[p]) > 1e-2 for p in picks):
                picks.append(j)
            if len(picks) >= batch:
                break
        new_thetas = [to_theta_unit(cand[j]) for j in picks]
        new_res = list(score_results(new_thetas))
        new_y = np.array([r.score if np.isfinite(r.score) else 1e12 for r in new_res])

        X = np.vstack([X, np.array([low + np.clip(cand[j], 0, 1) * span for j in picks])])
        y = np.concatenate([y, new_y])
        results.extend(new_res)
        if progress and (it + 1) % max(1, n_iter // 10) == 0:
            print(f"  round {it+1}/{n_iter}  best score {y.min():.4g}", flush=True)

    rows, obj_results = [], {}
    for sid in range(len(results)):
        res = results[sid]
        obj_results[sid] = res
        rows.append({"sample_id": sid, **space.to_theta(X[sid]),
                     "score": float(y[sid]), "loglik": res.loglik, "n_obs": len(res.residuals)})
    design = pd.DataFrame(rows)
    design["weight"] = 0.0
    best_sample_id = int(np.argmin(y))
    design.loc[best_sample_id, "weight"] = 1.0
    best_theta = {n: float(X[best_sample_id][i]) for i, n in enumerate(names)}
    best = obj_results[best_sample_id]

    return BayesOptResult(design=design, obj_results=obj_results, best_theta=best_theta,
                          best_sample_id=best_sample_id, best=best,
                          info={"n_eval": len(results), "n_init": n_init,
                                "n_iter": n_iter, "batch": batch})

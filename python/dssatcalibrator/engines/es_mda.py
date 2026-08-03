"""ES-MDA — Ensemble Smoother with Multiple Data Assimilation.

A derivative-free ensemble method (Emerick & Reynolds 2013) that estimates many
parameters *with uncertainty* in only a handful of iterations — and every
iteration scores its whole ensemble in **one parallel batch**, so it maps
perfectly onto the framework's parallel evaluator.

How it works
------------
Start from an ensemble drawn from the prior. Each iteration runs DSSAT for every
member, then nudges every member's parameters toward the observations using the
ensemble's own cross-covariance between parameters and simulated outputs
(a Kalman-style update). The observation error is *inflated* by a factor
``alpha`` and the assimilation is repeated ``Na`` times (with ``sum 1/alpha = 1``),
which makes the update gentle and far more robust to nonlinearity than a single
ensemble-smoother step. The spread of the final ensemble is the posterior
uncertainty.

Cost is fixed and predictable: ``Na x ensemble_size`` DSSAT evaluations total
(default ``4 x 32``), regardless of dimension — attractive for moderately
expensive crops where a long MCMC is impractical.

Config (under ``method.bayesian`` when ``engine: es_mda``)::

    ensemble_size: 32     # members (default max(4*n_params, 24))
    iterations: 4         # Na assimilation steps (inflation alpha = Na each)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import priors
from .mcmc import McmcResult


def _obs_vectors(results, names_key="user_var"):
    """Align each member's residual table to a common (exp,trt,var,date) key set.

    Returns ``(keys, d_obs, sigma, d_sim)`` where ``d_sim`` is ``(n_members x n_obs)``
    of simulated values, restricted to keys present in EVERY member so the matrices
    line up. Members with an empty residual table contribute no keys (handled by the
    caller)."""
    per_member = []
    for r in results:
        rd = getattr(r, "residuals", None)
        if rd is None or rd.empty:
            per_member.append({})
            continue
        m = {}
        for _, row in rd.iterrows():
            key = (row["exp_id"], int(row["treatment"]), row["dssat"],
                   str(row["date"]) if not pd.isna(row["date"]) else "NA")
            weight = max(float(row.get("weight", 1.0)), 1e-12)
            # ES-MDA uses an observation-error covariance; objective weights are
            # equivalent to shrinking sigma by sqrt(weight).
            m[key] = (float(row["sim"]), float(row["obs"]),
                      float(row["sigma"]) / np.sqrt(weight))
        per_member.append(m)

    common = None
    for m in per_member:
        if not m:
            continue
        ks = set(m)
        common = ks if common is None else (common & ks)
    common = sorted(common) if common else []
    if not common:
        return [], np.array([]), np.array([]), np.zeros((len(results), 0))

    d_obs = np.array([per_member[next(i for i, m in enumerate(per_member) if common[0] in m)][k][1]
                      for k in common])
    sigma = np.array([per_member[next(i for i, m in enumerate(per_member) if common[0] in m)][k][2]
                      for k in common])
    d_sim = np.full((len(results), len(common)), np.nan)
    for i, m in enumerate(per_member):
        if not m:
            continue
        for j, k in enumerate(common):
            if k in m:
                d_sim[i, j] = m[k][0]
    return common, d_obs, sigma, d_sim


def run_es_mda(cfg: dict, score_results, space, *, progress: bool = True) -> McmcResult:
    """Run ES-MDA. ``score_results(list_of_theta) -> list[ObjectiveResult]``."""
    bcfg = cfg.get("method", {}).get("bayesian", {})
    ne = int(bcfg.get("ensemble_size", max(4 * space.ndim, 24)))
    na = int(bcfg.get("iterations", 4))
    seed = int(cfg["calibrator"].get("seed", 42))
    rng = np.random.default_rng(seed)
    names = space.names

    # Initial ensemble from the prior (truncated to bounds).
    init = priors.sample_prior_design(space, ne, rng)
    ens = np.array([init.iloc[i].to_numpy(dtype=float) for i in range(ne)])
    alpha = float(na)                                   # equal inflation; sum(1/alpha)=1

    if progress:
        print(f"Running ES-MDA: {ne} members x {na} iterations "
              f"(alpha={alpha:g})...", flush=True)

    results = score_results([space.to_theta(ens[i]) for i in range(ne)])
    for it in range(na):
        keys, d_obs, sigma, d_sim = _obs_vectors(results)
        nd = len(keys)
        if nd == 0:
            if progress:
                print("  no common observations across members; stopping ES-MDA.", flush=True)
            break
        # Failed members must be penalised, never made identical to observations.
        bad = ~np.isfinite(d_sim)
        if bad.any():
            cols = np.where(bad)[1]
            d_sim[bad] = np.take(d_obs + 10.0 * np.maximum(sigma, 1e-6), cols)

        theta_mean = ens.mean(axis=0)
        d_mean = d_sim.mean(axis=0)
        Ta = ens - theta_mean                           # (ne x np)
        Da = d_sim - d_mean                             # (ne x nd)
        C_td = (Ta.T @ Da) / (ne - 1)                   # (np x nd)
        C_dd = (Da.T @ Da) / (ne - 1)                   # (nd x nd)
        R = alpha * np.diag(sigma ** 2)
        # Solve (C_dd + R) X = (d_uc - d_sim)^T robustly via pseudo-inverse.
        Cinv = np.linalg.pinv(C_dd + R)
        K = C_td @ Cinv                                 # (np x nd) Kalman gain

        # Perturbed observations per member.
        pert = d_obs[None, :] + np.sqrt(alpha) * sigma[None, :] * rng.standard_normal((ne, nd))
        innov = pert - d_sim                            # (ne x nd)
        ens = ens + innov @ K.T                         # (ne x np)
        ens = np.clip(ens, space.low, space.high)

        results = score_results([space.to_theta(ens[i]) for i in range(ne)])
        if progress:
            sc = np.array([r.score if np.isfinite(r.score) else np.nan for r in results])
            print(f"  iter {it+1}/{na}  median score {np.nanmedian(sc):.4g}", flush=True)

    rows, obj_results = [], {}
    for sid in range(ne):
        res = results[sid]
        obj_results[sid] = res
        rows.append({"sample_id": sid, **space.to_theta(ens[sid]),
                     "score": res.score, "loglik": res.loglik, "n_obs": len(res.residuals)})
    design = pd.DataFrame(rows)
    design["weight"] = 1.0 / len(design)                # final ensemble ~ posterior

    valid = design[np.isfinite(design["score"])]
    best_sample_id = int(valid["score"].idxmin()) if not valid.empty else 0
    best_theta = {n: float(design.loc[best_sample_id, n]) for n in names}
    best = obj_results[best_sample_id]

    q = float(bcfg.get("behavioural_quantile", 0.1))
    valid = design[np.isfinite(design["score"])]
    threshold = float(valid["score"].quantile(q)) if not valid.empty else float("inf")
    behavioural = design[design["score"] <= threshold].copy()

    return McmcResult(design=design, behavioural=behavioural, best_theta=best_theta,
                      best_sample_id=best_sample_id, threshold=threshold,
                      ess=float(len(design)), obj_results=obj_results, best=best,
                      acceptance=float("nan"), chain=pd.DataFrame(),
                      initial_design=pd.DataFrame(
                          [{"sample_id": i, **space.to_theta(init.iloc[i].to_numpy())}
                           for i in range(ne)]))

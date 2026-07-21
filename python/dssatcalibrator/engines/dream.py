"""DREAM (DE-MC) — robust Bayesian posterior for multimodal/correlated parameters.

Why this over plain Metropolis (``engines.mcmc``)
-------------------------------------------------
Random-walk Metropolis proposes jumps from a *fixed* Gaussian whose shape you
must tune; on a curved or correlated posterior (very common for crop genetic
coefficients — e.g. a phenology coefficient trades off against a growth one) it
mixes slowly and can get stuck in one mode. DREAM builds each proposal from the
*differences between other chains* (differential evolution), so the proposal
automatically adopts the scale and orientation of the posterior and several
chains can discover several modes. It is the de-facto standard for hydrologic
and crop-model Bayesian calibration (Vrugt 2016).

This is the DE-MC core of DREAM (ter Braak & Vrugt 2008): each chain's proposal
is built from the *difference between two other chains*, scaled by
``gamma = 2.38 / sqrt(2 d)``, with periodic ``gamma = 1`` jumps for mode-hopping.
Because the difference vectors come from the live population, the proposal
self-adapts to the posterior's scale and orientation. NumPy only; uses the
framework's ``prior`` declarations directly.

Config (under ``method.bayesian`` when ``engine: dream``)::

    n_chains: 8           # parallel chains (default max(2*n_params, 6); min 4)
    n_generations: 400    # generations
    burn_in: 200          # discarded warm-up (default half)
    thin: 1               # keep every k-th sample
    snooker: 0.1          # fraction of jumps that are gamma=1 mode-hops
    eps: 1e-4             # jitter added to each proposal (as a fraction of range)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import priors
from .mcmc import McmcResult


def run_dream(cfg: dict, score_results, space, *, progress: bool = True) -> McmcResult:
    """Run DREAM(ZS). ``score_results(list_of_theta) -> list[ObjectiveResult]`` is
    the framework's parallel evaluator; one generation is one parallel batch."""
    bcfg = cfg.get("method", {}).get("bayesian", {})
    n_chains = int(bcfg.get("n_chains", max(2 * space.ndim, 6)))
    n_gen = int(bcfg.get("n_generations", bcfg.get("n_steps", 400)))
    burn = int(bcfg.get("burn_in", n_gen // 2))
    thin = max(1, int(bcfg.get("thin", 1)))
    snooker = float(bcfg.get("snooker", 0.1))
    eps = float(bcfg.get("eps", 1e-4))
    seed = int(cfg["calibrator"].get("seed", 42))
    rng = np.random.default_rng(seed)

    names = space.names
    d = space.ndim
    ranges = space.high - space.low
    gamma_default = 2.38 / np.sqrt(2 * d)

    def vec(theta):
        return np.array([theta[n] for n in names], dtype=float)

    if n_chains < 4:
        n_chains = 4                                      # DE-MC needs >=2 other chains

    # Initialise chains from the prior.
    init = priors.sample_prior_design(space, n_chains, rng)
    cur_theta = [space.to_theta(init.iloc[c].to_numpy()) for c in range(n_chains)]
    cur_res = score_results(cur_theta)
    cur_lp = np.array([priors.log_prior_vec(space, t) for t in cur_theta])
    cur_logpost = np.array([lp + (r.loglik if np.isfinite(r.loglik) else -1e300)
                            for lp, r in zip(cur_lp, cur_res)])
    initial_design = pd.DataFrame([{"sample_id": c, **cur_theta[c]} for c in range(n_chains)])

    if progress:
        print(f"Running DREAM(ZS): {n_chains} chains x {n_gen} generations "
              f"(burn-in {burn}, thin {thin})...", flush=True)

    chain_rows, samples = [], []
    accepts = proposals = 0

    for gen in range(n_gen):
        # Snapshot the current chain states; difference vectors come from OTHER
        # chains (classic DE-MC), so the proposal self-scales as chains contract.
        cur_vecs = [vec(t) for t in cur_theta]
        prop = []
        for c in range(n_chains):
            others = [k for k in range(n_chains) if k != c]
            a, b = rng.choice(others, size=2, replace=False)
            gamma = 1.0 if rng.uniform() < snooker else gamma_default
            jump = gamma * (cur_vecs[a] - cur_vecs[b]) + eps * ranges * rng.standard_normal(d)
            prop.append(space.to_theta(np.clip(cur_vecs[c] + jump, space.low, space.high)))

        lp_prop = np.array([priors.log_prior_vec(space, t) for t in prop])
        idx_in = [c for c in range(n_chains) if np.isfinite(lp_prop[c])]
        res_in = score_results([prop[c] for c in idx_in]) if idx_in else []
        logpost_prop = np.full(n_chains, -np.inf)
        res_prop = [None] * n_chains
        for c, r in zip(idx_in, res_in):
            ll = r.loglik if np.isfinite(r.loglik) else -1e300
            logpost_prop[c] = lp_prop[c] + ll
            res_prop[c] = r

        for c in range(n_chains):
            proposals += 1
            if np.log(rng.uniform()) < (logpost_prop[c] - cur_logpost[c]):
                cur_theta[c], cur_logpost[c], cur_res[c] = prop[c], logpost_prop[c], res_prop[c]
                accepts += 1
            chain_rows.append({"step": gen, "walker": c, "logpost": float(cur_logpost[c]),
                               **cur_theta[c]})
            if gen >= burn and ((gen - burn) % thin == 0):
                samples.append((cur_theta[c], cur_res[c]))

        if progress and (gen + 1) % max(1, n_gen // 10) == 0:
            print(f"  gen {gen+1}/{n_gen}  acceptance {accepts/max(proposals,1):.2f}", flush=True)

    if not samples:
        samples = [(cur_theta[c], cur_res[c]) for c in range(n_chains)]

    rows, obj_results = [], {}
    for sid, (theta, res) in enumerate(samples):
        obj_results[sid] = res
        rows.append({"sample_id": sid, **theta, "score": res.score, "loglik": res.loglik,
                     "n_obs": len(res.residuals)})
    design = pd.DataFrame(rows)
    design["weight"] = 1.0 / len(design)

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
                      acceptance=accepts / max(proposals, 1),
                      chain=pd.DataFrame(chain_rows), initial_design=initial_design)

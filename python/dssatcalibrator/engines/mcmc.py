"""MCMC — full Bayesian posterior by Markov-chain Monte-Carlo (preset D).

What you get
------------
Where GLUE and SMC-PF give a weighted cloud of parameter sets, MCMC builds a chain
of samples whose density *is* the posterior — giving clean credible intervals and
parameter correlations. We use an **adaptive random-walk Metropolis** ensemble:

* several chains ("walkers") explore in parallel;
* each step proposes a small Gaussian jump, then accepts it with probability
  ``min(1, posterior(new) / posterior(current))`` where
  ``posterior = likelihood x prior``;
* the proposal size self-tunes during *burn-in* toward a healthy ~23% acceptance
  rate, then freezes (so the kept samples are a valid posterior).

It is dependency-free (NumPy only) and uses the framework's ``prior`` declarations
directly, so ``prior: {dist: normal, sd: ...}`` genuinely informs the posterior.

Cost & tips
-----------
Every step evaluates ``n_walkers`` parameter sets, each requiring DSSAT runs across
all experiments — MCMC is the most expensive engine. All walkers of a step are
scored in **one parallel batch**. For expensive crops, screen first
(``method.sensitivity``) to cut the parameter count, or train the surrogate
(:mod:`engines.surrogate`) and run MCMC on that.

Config (under ``method.bayesian`` when ``engine: mcmc``)::

    n_walkers: 16        # chains (default max(2*n_params, 8))
    n_steps: 400         # iterations per walker
    burn_in: 200         # discarded warm-up (default half of n_steps)
    thin: 1              # keep every k-th sample
    proposal_scale: 0.1  # initial jump size as a fraction of each range
    target_accept: 0.234 # acceptance the adaptation aims for
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .. import priors


@dataclass
class McmcResult:
    design: pd.DataFrame          # posterior samples (+ "weight", "score", "loglik")
    behavioural: pd.DataFrame
    best_theta: dict
    best_sample_id: int
    threshold: float
    ess: float
    obj_results: dict
    best: object
    acceptance: float
    chain: pd.DataFrame = field(default_factory=pd.DataFrame)
    initial_design: pd.DataFrame = None


def run_mcmc(cfg: dict, score_results, space, *, progress: bool = True) -> McmcResult:
    """Run adaptive random-walk Metropolis. ``score_results(list_of_theta) -> list
    of ObjectiveResult`` is the framework's parallel evaluator."""
    bcfg = cfg.get("method", {}).get("bayesian", {})
    n_walkers = int(bcfg.get("n_walkers", max(2 * space.ndim, 8)))
    n_steps = int(bcfg.get("n_steps", 400))
    burn = int(bcfg.get("burn_in", n_steps // 2))
    thin = max(1, int(bcfg.get("thin", 1)))
    target = float(bcfg.get("target_accept", 0.234))
    scale = float(bcfg.get("proposal_scale", 0.1))
    adapt_every = int(bcfg.get("adapt_interval", 20))
    seed = int(cfg["calibrator"].get("seed", 42))
    rng = np.random.default_rng(seed)

    ranges = space.high - space.low
    names = space.names

    def vec(theta):
        return np.array([theta[n] for n in names], dtype=float)

    # --- initialise walkers from the prior (uniform prior -> uniform draws) -------
    init = priors.sample_prior_design(space, n_walkers, rng)
    cur_theta = [space.to_theta(init.iloc[w].to_numpy()) for w in range(n_walkers)]
    init_res = score_results(cur_theta)
    cur_lp = np.array([priors.log_prior_vec(space, t) for t in cur_theta])
    cur_logpost = np.array([
        lp + (r.loglik if np.isfinite(r.loglik) else -1e300)
        for lp, r in zip(cur_lp, init_res)
    ])
    cur_res = list(init_res)
    initial_design = pd.DataFrame([{"sample_id": w, **cur_theta[w]} for w in range(n_walkers)])

    if progress:
        print(f"Running MCMC: {n_walkers} walkers x {n_steps} steps "
              f"(burn-in {burn}, thin {thin})...", flush=True)

    chain_rows, samples = [], []          # samples: list of (theta, ObjectiveResult)
    accepts = proposals = 0

    for step in range(n_steps):
        sd = scale * ranges
        prop = [space.to_theta(vec(cur_theta[w]) + rng.normal(0.0, sd)) for w in range(n_walkers)]
        lp_prop = np.array([priors.log_prior_vec(space, t) for t in prop])

        # Only run DSSAT for in-bounds proposals; out-of-bounds auto-reject.
        idx_in = [w for w in range(n_walkers) if np.isfinite(lp_prop[w])]
        res_in = score_results([prop[w] for w in idx_in]) if idx_in else []
        logpost_prop = np.full(n_walkers, -np.inf)
        res_prop = [None] * n_walkers
        for w, r in zip(idx_in, res_in):
            ll = r.loglik if np.isfinite(r.loglik) else -1e300
            logpost_prop[w] = lp_prop[w] + ll
            res_prop[w] = r

        for w in range(n_walkers):
            proposals += 1
            if np.log(rng.uniform()) < (logpost_prop[w] - cur_logpost[w]):
                cur_theta[w], cur_logpost[w], cur_res[w] = prop[w], logpost_prop[w], res_prop[w]
                accepts += 1
            chain_rows.append({"step": step, "walker": w, "logpost": float(cur_logpost[w]),
                               **cur_theta[w]})
            if step >= burn and ((step - burn) % thin == 0):
                samples.append((cur_theta[w], cur_res[w]))

        # Adapt the jump size during burn-in only (keeps the kept samples valid).
        if step < burn and (step + 1) % adapt_every == 0 and proposals > 0:
            ar = accepts / proposals
            scale = float(np.clip(scale * np.exp(ar - target), 1e-3, 1.0))
        if progress and (step + 1) % max(1, n_steps // 10) == 0:
            print(f"  step {step+1}/{n_steps}  acceptance {accepts/max(proposals,1):.2f}  "
                  f"scale {scale:.3g}", flush=True)

    if not samples:                          # degenerate (burn >= n_steps): use final states
        samples = [(cur_theta[w], cur_res[w]) for w in range(n_walkers)]

    rows, obj_results = [], {}
    for sid, (theta, res) in enumerate(samples):
        obj_results[sid] = res
        rows.append({"sample_id": sid, **theta, "score": res.score, "loglik": res.loglik,
                     "n_obs": len(res.residuals)})
    design = pd.DataFrame(rows)
    design["weight"] = 1.0 / len(design)     # posterior samples are already posterior-distributed

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

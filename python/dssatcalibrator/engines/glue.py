"""GLUE / Monte-Carlo engine (preset C, the v1 default).

Post-processes an evaluated design (one row per sampled parameter set, with a
``score`` and ``loglik``) into:

* posterior-like **likelihood weights** ``w = softmax(loglik)`` (GLUE), and
* a **behavioural set** — the best ``behavioural_quantile`` fraction by score,

plus the single best (MAP-ish) parameter set. This mirrors DSSAT's own GLUE tool
(Monte-Carlo sampling + Gaussian likelihood) but over the framework's joint
multi-variable / multi-experiment objective.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class GlueResult:
    design: pd.DataFrame          # design + a "weight" column (posterior weights)
    behavioural: pd.DataFrame     # the behavioural subset
    best_theta: dict
    best_sample_id: int
    threshold: float
    ess: float                    # effective sample size of the posterior weights


def run_glue(design: pd.DataFrame, param_names: list[str], cfg: dict,
             space=None) -> GlueResult:
    """Turn an evaluated design into GLUE weights + a behavioural set.

    ``space`` is optional: when given and any parameter declares a non-uniform
    prior, each run's weight becomes ``likelihood x prior`` (proper Bayesian
    importance weighting). With uniform priors (or ``space=None``) the prior is a
    constant and the weights are pure likelihood — the original behaviour.
    """
    d = design.copy().reset_index(drop=True)
    ll = d["loglik"].to_numpy(dtype=float)

    # Fold the log-prior into the log-likelihood when an informative prior exists.
    if space is not None:
        from .. import priors
        if priors.has_informative_prior(space):
            lp = np.array([
                priors.log_prior_vec(space, {n: float(d.loc[i, n]) for n in param_names})
                for i in range(len(d))
            ])
            ll = ll + lp

    finite = np.isfinite(ll)
    w = np.zeros(len(d))
    if finite.any():
        m = np.max(ll[finite])
        w[finite] = np.exp(ll[finite] - m)
        total = w.sum()
        w = w / total if total > 0 else w
    d["weight"] = w
    ess = float(1.0 / np.sum(w ** 2)) if w.sum() > 0 else 0.0

    q = float(cfg.get("method", {}).get("bayesian", {}).get("behavioural_quantile", 0.1))
    valid = d[np.isfinite(d["score"])]
    threshold = float(valid["score"].quantile(q)) if not valid.empty else float("inf")
    behavioural = d[d["score"] <= threshold].copy()

    best_id = int(valid["score"].idxmin()) if not valid.empty else 0
    best_theta = {n: float(d.loc[best_id, n]) for n in param_names}

    return GlueResult(design=d, behavioural=behavioural, best_theta=best_theta,
                      best_sample_id=best_id, threshold=threshold, ess=ess)


def posterior_summary(glue: GlueResult, param_names: list[str]) -> pd.DataFrame:
    """Weighted posterior mean / sd / quantiles per parameter (for reporting)."""
    d = glue.design
    w = d["weight"].to_numpy()

    def wquantile(x, probs):
        if w.sum() <= 0:
            return np.quantile(x, probs)
        order = np.argsort(x)
        xs = x[order]
        ws = w[order] / w[order].sum()
        cdf = np.cumsum(ws)
        return np.interp(probs, cdf, xs)

    rows = []
    for n in param_names:
        x = d[n].to_numpy(dtype=float)
        if w.sum() > 0:
            mean = float(np.sum(w * x))
            var = float(np.sum(w * (x - mean) ** 2))
            sd = float(np.sqrt(max(var, 0.0)))
        else:
            mean, sd = float(np.mean(x)), float(np.std(x))
        q05, q95 = wquantile(x, [0.05, 0.95])
        rows.append({"parameter": n, "best": glue.best_theta[n],
                     "post_mean": mean, "post_sd": sd,
                     "p05": float(q05), "p95": float(q95)})
    return pd.DataFrame(rows)

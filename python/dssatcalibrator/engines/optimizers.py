"""Classical optimisers — find a single best-fit parameter set (preset B).

When you do not need full uncertainty and just want *good calibrated numbers*,
these search for the parameter vector that minimises the objective:

============  ================================================================
``method``    what it does
============  ================================================================
``nelder_mead``  Downhill-simplex search from a start point. Fast and local; set
                 ``restarts`` > 1 to launch several from random points and keep
                 the best (a cheap way to escape poor local minima).
``diffevo``      Differential evolution — a global, population-based search.
                 Robust but uses more model runs. The whole population of each
                 generation is scored in **one parallel batch**.
``cmaes``        Covariance-Matrix-Adaptation Evolution Strategy — the strongest
                 gradient-free optimiser for 5-30 continuous coefficients. It
                 learns the local shape of the objective (an evolving covariance)
                 so it follows curved valleys that trip up Nelder-Mead, and
                 typically reaches a good fit in fewer DSSAT runs than diffevo.
                 Each generation's population is scored in **one parallel batch**.
============  ================================================================

Both call a single ``score_batch(list_of_theta) -> list_of_float`` callback, so
every DSSAT run is dispatched through the framework's parallel evaluator and uses
all cores. (The AgMIP stepwise protocol lives in :mod:`engines.selection`; it
repeatedly calls ``nelder_mead`` while growing the parameter set.)

Only SciPy is required (already a core dependency).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

_FAIL = 1e12  # finite stand-in for a failed/infinite score, so the search continues


@dataclass
class OptimizerResult:
    best_theta: dict
    best_score: float
    method: str
    n_eval: int
    history: list = field(default_factory=list)   # [{"iter": int, "score": float}]


def run_optimizer(space, score_batch: Callable[[list[dict]], list[float]], *,
                  method: str = "diffevo", seed: int = 42, maxiter: int | None = None,
                  popsize: int = 15, restarts: int = 1, tol: float = 1e-4,
                  progress: bool = False) -> OptimizerResult:
    """Minimise ``score_batch`` over ``space`` and return the best parameter set."""
    from scipy.optimize import differential_evolution, minimize

    bounds = list(zip(space.low, space.high))
    state = {"n": 0, "best_theta": None, "best_score": np.inf, "hist": []}

    def _note(theta: dict, score: float):
        if score < state["best_score"]:
            state["best_score"] = score
            state["best_theta"] = theta
        state["hist"].append({"iter": state["n"], "score": state["best_score"]})
        if progress and state["n"] % 25 == 0:
            print(f"  opt eval {state['n']}: best score so far {state['best_score']:.4g}", flush=True)

    def cost_one(x):
        theta = space.to_theta(np.clip(x, space.low, space.high))
        s = score_batch([theta])[0]
        s = float(s) if np.isfinite(s) else _FAIL
        state["n"] += 1
        _note(theta, s)
        return s

    method = method.lower()
    if method in ("nelder_mead", "neldermead", "nm"):
        rng = np.random.default_rng(seed)
        starts = [np.asarray(space.start, float)]
        starts += [rng.uniform(space.low, space.high) for _ in range(max(0, restarts - 1))]
        for x0 in starts:
            minimize(cost_one, x0, method="Nelder-Mead", bounds=bounds,
                     options={"maxiter": maxiter or 200 * space.ndim,
                              "xatol": tol, "fatol": tol})

    elif method in ("diffevo", "differential_evolution", "de"):
        def cost_vec(X):
            # SciPy hands us X with shape (n_params, n_population): one column per
            # candidate. Score the whole population in a single parallel batch.
            thetas = [space.to_theta(np.clip(X[:, j], space.low, space.high))
                      for j in range(X.shape[1])]
            scores = score_batch(thetas)
            arr = np.array([float(s) if np.isfinite(s) else _FAIL for s in scores])
            state["n"] += len(thetas)
            j = int(np.argmin(arr))
            _note(thetas[j], float(arr[j]))
            return arr

        differential_evolution(cost_vec, bounds, seed=seed, vectorized=True,
                               updating="deferred", maxiter=maxiter or 30,
                               popsize=popsize, tol=tol, polish=False, init="sobol")

    elif method in ("cmaes", "cma_es", "cma"):
        def score_pop(thetas):
            scores = score_batch(thetas)
            arr = np.array([float(s) if np.isfinite(s) else _FAIL for s in scores])
            state["n"] += len(thetas)
            j = int(np.argmin(arr))
            _note(thetas[j], float(arr[j]))
            return arr
        _cma_es(space, score_pop, seed=seed, maxiter=maxiter,
                popsize=popsize if popsize and popsize > 4 else None, progress=progress)

    else:
        raise ValueError(f"unknown optimizer method '{method}' (use nelder_mead | diffevo | cmaes)")

    if state["best_theta"] is None:                # nothing scored finitely
        state["best_theta"] = space.to_theta(space.start)
        state["best_score"] = float("inf")
    return OptimizerResult(best_theta=state["best_theta"], best_score=state["best_score"],
                           method=method, n_eval=state["n"], history=state["hist"])


def _cma_es(space, score_pop, *, seed=42, maxiter=None, popsize=None, progress=False):
    """Minimal CMA-ES (Hansen) on the unit cube, mapped to the parameter bounds.

    Dependency-free reference implementation of the (mu/mu_w, lambda)-CMA-ES with
    the standard rank-mu + rank-one covariance update and cumulative step-size
    adaptation. Candidates are sampled in a normalised ``[0, 1]`` box (clipped to
    stay feasible) and mapped to real bounds, which keeps the step size meaningful
    across coefficients with very different units. ``score_pop(list_of_theta) ->
    np.ndarray`` scores a whole generation in one parallel batch.
    """
    rng = np.random.default_rng(seed)
    n = space.ndim
    low = np.asarray(space.low, float)
    high = np.asarray(space.high, float)
    span = np.where(high > low, high - low, 1.0)

    # Start at the normalised start point (fall back to the box centre).
    start = np.asarray(space.start, float)
    mean = np.clip((start - low) / span, 0.0, 1.0)
    sigma = 0.3

    lam = int(popsize) if popsize else 4 + int(3 * np.log(n))
    mu = lam // 2
    w = np.log(mu + 0.5) - np.log(np.arange(1, mu + 1))
    w /= w.sum()
    mueff = 1.0 / np.sum(w ** 2)

    # Strategy parameters (Hansen 2016 defaults).
    cc = (4 + mueff / n) / (n + 4 + 2 * mueff / n)
    cs = (mueff + 2) / (n + mueff + 5)
    c1 = 2 / ((n + 1.3) ** 2 + mueff)
    cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((n + 2) ** 2 + mueff))
    damps = 1 + 2 * max(0, np.sqrt((mueff - 1) / (n + 1)) - 1) + cs
    chiN = np.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n ** 2))

    pc = np.zeros(n)
    ps = np.zeros(n)
    B = np.eye(n)
    D = np.ones(n)
    C = np.eye(n)
    n_gen = int(maxiter) if maxiter else 100 + 50 * n

    def to_theta(x_unit):
        return space.to_theta(low + np.clip(x_unit, 0.0, 1.0) * span)

    for gen in range(n_gen):
        # Sample lambda candidates: y ~ N(0, C), x = mean + sigma*y (in unit box).
        z = rng.standard_normal((lam, n))
        y = z @ (B * D).T
        x = mean + sigma * y
        x = np.clip(x, 0.0, 1.0)
        scores = score_pop([to_theta(x[k]) for k in range(lam)])

        order = np.argsort(scores)
        x_sorted = x[order][:mu]
        old_mean = mean.copy()
        mean = w @ x_sorted

        # Step-size control + covariance adaptation.
        y_w = (mean - old_mean) / sigma
        C_invsqrt = B @ np.diag(1.0 / D) @ B.T
        ps = (1 - cs) * ps + np.sqrt(cs * (2 - cs) * mueff) * (C_invsqrt @ y_w)
        hsig = (np.linalg.norm(ps) / np.sqrt(1 - (1 - cs) ** (2 * (gen + 1))) / chiN
                < 1.4 + 2 / (n + 1))
        pc = (1 - cc) * pc + hsig * np.sqrt(cc * (2 - cc) * mueff) * y_w

        y_k = (x_sorted - old_mean) / sigma
        C = ((1 - c1 - cmu) * C
             + c1 * (np.outer(pc, pc) + (1 - hsig) * cc * (2 - cc) * C)
             + cmu * (y_k.T * w) @ y_k)
        sigma *= np.exp((cs / damps) * (np.linalg.norm(ps) / chiN - 1))
        sigma = float(np.clip(sigma, 1e-4, 1.0))

        # Eigen-decompose C for the next sampling (symmetrise first).
        C = np.triu(C) + np.triu(C, 1).T
        try:
            evals, B = np.linalg.eigh(C)
            evals = np.clip(evals, 1e-14, None)
            D = np.sqrt(evals)
        except np.linalg.LinAlgError:
            B, D, C = np.eye(n), np.ones(n), np.eye(n)

        if progress and (gen + 1) % max(1, n_gen // 10) == 0:
            print(f"  cma-es gen {gen+1}/{n_gen}  sigma {sigma:.3g}  "
                  f"best {scores[order][0]:.4g}", flush=True)
        if sigma < 1e-3 and np.max(D) * sigma < 1e-3:
            break

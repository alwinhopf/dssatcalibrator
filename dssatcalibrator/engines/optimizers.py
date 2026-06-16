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
    else:
        raise ValueError(f"unknown optimizer method '{method}' (use nelder_mead | diffevo)")

    if state["best_theta"] is None:                # nothing scored finitely
        state["best_theta"] = space.to_theta(space.start)
        state["best_score"] = float("inf")
    return OptimizerResult(best_theta=state["best_theta"], best_score=state["best_score"],
                           method=method, n_eval=state["n"], history=state["hist"])

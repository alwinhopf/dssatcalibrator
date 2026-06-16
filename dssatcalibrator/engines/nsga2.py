"""NSGA-II multi-objective engine (per-variable Pareto front).

What it does
------------
Instead of collapsing LAI / biomass / yield / phenology into one weighted score,
NSGA-II treats each chosen variable's misfit (nRMSE%) as a **separate** objective
and returns the *trade-off front* — the parameter sets where you cannot improve
one variable's fit without worsening another. Useful when targets genuinely
conflict (e.g. the parameters that nail yield slightly spoil phenology).

Parallelism
-----------
The whole population of a generation is evaluated in **one** parallel batch via the
``evaluate_batch`` callback (which fans every ``candidate x experiment`` DSSAT run
across all cores). This is the fix for the old version, which evaluated population
members one at a time and left most cores idle.

Cost = ``pop_size x n_gen x n_experiments`` DSSAT runs. Keep ``pop_size``/``n_gen``
small unless you have the compute (or pair it with the surrogate engine).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass
class Nsga2Result:
    objective_vars: list[str]
    X: np.ndarray                 # Pareto-optimal parameter vectors
    F: np.ndarray                 # their objective values (nRMSE% per variable)
    param_names: list[str]

    def front(self):
        import pandas as pd
        cols = {n: self.X[:, i] for i, n in enumerate(self.param_names)}
        for j, v in enumerate(self.objective_vars):
            cols[f"nRMSE_{v}"] = self.F[:, j]
        return pd.DataFrame(cols)


def run_nsga2(evaluate_batch: Callable[[list[dict]], list[dict]], space,
              objective_vars: list[str], pop_size: int = 16, n_gen: int = 5,
              seed: int = 42) -> Nsga2Result:
    """Run NSGA-II.

    Parameters
    ----------
    evaluate_batch
        ``evaluate_batch(list_of_theta) -> list_of_{user_var: nRMSE%}`` — scores a
        whole list of parameter sets at once (so the population runs in parallel).
        Missing variables for a candidate should be reported as a large value.
    space
        the :class:`spaces.ParameterSpace` (gives bounds + name order).
    objective_vars
        the variables to trade off, one objective each.
    """
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.core.problem import Problem
    from pymoo.optimize import minimize

    names = space.names
    big = 1e6  # penalty for a candidate that fails to produce a comparison

    class _Problem(Problem):
        def __init__(self):
            super().__init__(n_var=space.ndim, n_obj=len(objective_vars),
                             xl=space.low, xu=space.high)

        def _evaluate(self, X, out, *args, **kwargs):
            # pymoo hands us the entire population at once -> evaluate it in one
            # parallel batch rather than looping candidate-by-candidate.
            thetas = [space.to_theta(x) for x in np.atleast_2d(X)]
            per_var = evaluate_batch(thetas)
            out["F"] = np.array(
                [[float(pv.get(v, big)) for v in objective_vars] for pv in per_var]
            )

    res = minimize(_Problem(), NSGA2(pop_size=pop_size),
                   ("n_gen", n_gen), seed=seed, verbose=False)
    X = np.atleast_2d(res.X)
    F = np.atleast_2d(res.F)
    return Nsga2Result(objective_vars=objective_vars, X=X, F=F, param_names=names)

"""NSGA-II multi-objective engine (per-variable Pareto front).

Instead of collapsing LAI / biomass / yield / phenology into one weighted score,
NSGA-II treats each chosen variable's misfit (nRMSE%) as a separate objective and
returns the **trade-off front** — the parameter sets where you cannot improve one
variable's fit without worsening another. Useful when targets genuinely conflict.

DSSAT runs in the optimisation loop, so cost = pop_size x n_gen x n_experiments.
Keep pop/gen small unless you have the compute.
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


def run_nsga2(evaluate_theta: Callable[[dict], dict], space, objective_vars: list[str],
              pop_size: int = 16, n_gen: int = 5, seed: int = 42) -> Nsga2Result:
    """Run NSGA-II. ``evaluate_theta(theta) -> {user_var: nRMSE%}``.

    Missing variables for a candidate are penalised with a large value so the
    optimiser avoids parameter sets that fail to produce a comparison.
    """
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.core.problem import Problem
    from pymoo.optimize import minimize

    names = space.names

    class _Problem(Problem):
        def __init__(self):
            super().__init__(n_var=space.ndim, n_obj=len(objective_vars),
                             xl=space.low, xu=space.high)

        def _evaluate(self, X, out, *args, **kwargs):
            F = []
            for x in np.atleast_2d(X):
                pv = evaluate_theta(space.to_theta(x))
                F.append([float(pv.get(v, 1e6)) for v in objective_vars])
            out["F"] = np.array(F)

    res = minimize(_Problem(), NSGA2(pop_size=pop_size),
                   ("n_gen", n_gen), seed=seed, verbose=False)
    X = np.atleast_2d(res.X)
    F = np.atleast_2d(res.F)
    return Nsga2Result(objective_vars=objective_vars, X=X, F=F, param_names=names)

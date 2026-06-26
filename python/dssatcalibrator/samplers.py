"""Design-of-experiment samplers over a :class:`ParameterSpace`.

Backed by ``scipy.stats.qmc`` (no extra dependency): Latin-Hypercube and Sobol
low-discrepancy sequences, plain Monte-Carlo, and a coarse grid. Each returns a
``pandas.DataFrame`` of samples in native parameter units, one column per param.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy.stats import qmc

from .spaces import ParameterSpace


def sample(space: ParameterSpace, n: int, engine: str = "lhs",
           seed: int = 42, include_start: bool = True) -> pd.DataFrame:
    """Draw ``n`` parameter sets from ``space`` using the chosen design engine.

    engine: ``lhs`` | ``sobol`` | ``montecarlo`` | ``grid``. When
    ``include_start`` the configured start point is prepended as the first row
    (so a calibration always evaluates the prior/default).
    """
    d = space.ndim
    rng = np.random.default_rng(seed)
    engine = engine.lower()

    if engine == "lhs":
        unit = qmc.LatinHypercube(d=d, seed=seed).random(n)
    elif engine == "sobol":
        unit = qmc.Sobol(d=d, scramble=True, seed=seed).random(n)
    elif engine == "montecarlo":
        unit = rng.random((n, d))
    elif engine == "grid":
        per = max(2, int(round(n ** (1.0 / d))))
        axes = [np.linspace(0, 1, per) for _ in range(d)]
        unit = np.array(list(itertools.product(*axes)))[:n]
    else:
        raise ValueError(f"unknown sampler engine: {engine}")

    scaled = qmc.scale(unit, space.low, space.high)
    df = pd.DataFrame(scaled, columns=space.names)

    if include_start:
        start_row = pd.DataFrame([space.start], columns=space.names)
        df = pd.concat([start_row, df], ignore_index=True)
    return df

"""Design-of-experiment samplers over a :class:`ParameterSpace`.

Backed by ``scipy.stats.qmc`` (no extra dependency): Latin-Hypercube and Sobol
low-discrepancy sequences, plain Monte-Carlo, and a coarse grid. Each returns a
``pandas.DataFrame`` of samples in native parameter units, one column per param.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
try:
    from scipy.stats import qmc
except ImportError:
    qmc = None

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
        if qmc is not None:
            unit = qmc.LatinHypercube(d=d, seed=seed).random(n)
        else:
            unit = np.zeros((n, d))
            for i in range(d):
                unit[:, i] = (rng.permutation(n) + rng.random(n)) / n
    elif engine == "sobol":
        if qmc is not None:
            unit = qmc.Sobol(d=d, scramble=True, seed=seed).random(n)
        else:
            unit = rng.random((n, d))
    elif engine == "montecarlo":
        unit = rng.random((n, d))
    elif engine == "grid":
        per = max(2, int(np.ceil(n ** (1.0 / d))))
        axes = [np.linspace(0, 1, per) for _ in range(d)]
        rows = itertools.islice(itertools.product(*axes), n)
        unit = np.fromiter((value for row in rows for value in row),
                           dtype=float, count=n * d).reshape(-1, d)
    else:
        raise ValueError(f"unknown sampler engine: {engine}")

    scaled = space.low + unit * (space.high - space.low)
    df = pd.DataFrame(scaled, columns=space.names)

    if include_start:
        start_row = pd.DataFrame([space.start], columns=space.names)
        df = pd.concat([start_row, df], ignore_index=True)
    return df

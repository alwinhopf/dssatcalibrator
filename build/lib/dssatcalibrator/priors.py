"""Prior distributions for parameters — sampling and log-density.

Why this module exists
----------------------
The config lets every parameter declare a ``prior`` (see ``CONCEPT.md`` §4)::

    P5: { active: true, min: 300, max: 700, start: 505, prior: {dist: normal, sd: 50} }

A *prior* is your belief about a parameter **before** looking at the data: a
plausible value (``start``), a spread (``sd``), and hard bounds (``min``/``max``).
The Bayesian engines (GLUE, SMC-PF, MCMC) combine this prior with how well each
parameter set fits the observations to produce a *posterior* (your belief
**after** seeing the data).

Earlier versions of the framework silently ignored the ``prior`` field and
always assumed a flat (uniform) prior between ``min`` and ``max``. This module
makes the declared prior actually count, while keeping ``uniform`` as the safe
default when no ``prior`` is given.

Supported distributions (all truncated to ``[min, max]``):

======================  ==========================================================
``dist``                meaning
======================  ==========================================================
``uniform`` (default)   every value in ``[min, max]`` equally likely
``normal``              bell curve centred on ``start`` (or ``prior.mean``),
                        width ``prior.sd``; values outside the bounds are excluded
``lognormal``           skewed, positive-only; ``prior.sd`` is the sd of log(x)
``triangular``          peak at ``start`` (or ``prior.mode``), zero at the bounds
======================  ==========================================================

Everything is deliberately small and dependency-light: it uses only SciPy
(already a core dependency). To add a new distribution, add one ``elif`` branch
in both :func:`sample_one` and :func:`log_prior_one` — nothing else changes.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def _bounds(spec: dict) -> tuple[float, float]:
    return float(spec["min"]), float(spec["max"])


def _center(spec: dict) -> float:
    """The prior's central value: explicit ``prior.mean``/``mode`` else ``start``
    else the midpoint of the range."""
    lo, hi = _bounds(spec)
    prior = spec.get("prior") or {}
    if "mean" in prior:
        return float(prior["mean"])
    if "mode" in prior:
        return float(prior["mode"])
    if "start" in spec and spec["start"] is not None:
        return float(spec["start"])
    return 0.5 * (lo + hi)


def _dist_name(spec: dict) -> str:
    prior = spec.get("prior") or {}
    return str(prior.get("dist", "uniform")).lower()


def sample_one(spec: dict, n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw ``n`` samples for a single parameter from its (truncated) prior."""
    lo, hi = _bounds(spec)
    dist = _dist_name(spec)
    prior = spec.get("prior") or {}

    if dist == "uniform":
        return rng.uniform(lo, hi, size=n)

    if dist == "normal":
        mu = _center(spec)
        sd = float(prior.get("sd", 0.25 * (hi - lo)))
        # Truncated normal: a, b are the bounds expressed in standard deviations.
        a, b = (lo - mu) / sd, (hi - mu) / sd
        return stats.truncnorm.rvs(a, b, loc=mu, scale=sd, size=n, random_state=rng)

    if dist == "lognormal":
        # prior.sd is the sd of log(x); centre is the median (exp of log-mean).
        center = max(_center(spec), 1e-9)
        sigma = float(prior.get("sd", 0.5))
        draws = rng.lognormal(mean=np.log(center), sigma=sigma, size=n)
        return np.clip(draws, lo, hi)

    if dist == "triangular":
        mode = _center(spec)
        return rng.triangular(lo, np.clip(mode, lo, hi), hi, size=n)

    raise ValueError(f"unknown prior dist '{dist}' for parameter '{spec.get('name')}'")


def log_prior_one(spec: dict, value: float) -> float:
    """Log prior density of one parameter value (``-inf`` outside the bounds).

    Normalising constants are dropped where they cancel in Metropolis-Hastings
    ratios; what matters is the *shape* and that out-of-bounds values are
    impossible.
    """
    lo, hi = _bounds(spec)
    if not (lo <= value <= hi):
        return float("-inf")
    dist = _dist_name(spec)
    prior = spec.get("prior") or {}

    if dist == "uniform":
        return 0.0  # constant inside the box

    if dist == "normal":
        mu = _center(spec)
        sd = float(prior.get("sd", 0.25 * (hi - lo)))
        a, b = (lo - mu) / sd, (hi - mu) / sd
        return float(stats.truncnorm.logpdf(value, a, b, loc=mu, scale=sd))

    if dist == "lognormal":
        center = max(_center(spec), 1e-9)
        sigma = float(prior.get("sd", 0.5))
        return float(stats.lognorm.logpdf(value, s=sigma, scale=center))

    if dist == "triangular":
        mode = float(np.clip(_center(spec), lo, hi))
        c = (mode - lo) / (hi - lo) if hi > lo else 0.5
        return float(stats.triang.logpdf(value, c, loc=lo, scale=(hi - lo)))

    raise ValueError(f"unknown prior dist '{dist}' for parameter '{spec.get('name')}'")


def sample_prior_design(space, n: int, rng: np.random.Generator):
    """Draw ``n`` parameter sets, one column per active parameter.

    Returns a ``pandas.DataFrame`` in native parameter units (same shape the
    samplers in :mod:`samplers` produce), so engines can use prior draws or
    space-filling draws interchangeably.
    """
    import pandas as pd

    cols = {s["name"]: sample_one(s, n, rng) for s in space.specs}
    return pd.DataFrame(cols, columns=space.names)


def log_prior_vec(space, theta: dict) -> float:
    """Total log prior of a parameter vector = sum over independent parameters."""
    total = 0.0
    for s in space.specs:
        total += log_prior_one(s, float(theta[s["name"]]))
        if total == float("-inf"):
            break
    return total


def has_informative_prior(space) -> bool:
    """True if *any* active parameter declares a non-uniform prior.

    Lets engines skip the prior term entirely (it is a no-op constant) when every
    parameter is uniform — keeping the old, simpler behaviour bit-for-bit.
    """
    return any(_dist_name(s) != "uniform" for s in space.specs)

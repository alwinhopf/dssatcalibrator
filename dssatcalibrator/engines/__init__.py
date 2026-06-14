"""Pluggable calibration engines.

Implemented now:
    glue   — Monte-Carlo / GLUE over a sampled design (preset C, the v1 default)
    smc_pf — sequential Monte-Carlo particle filter + Metropolis-Hastings move
             (preset A); assimilates time-ordered observations, resamples on ESS
    nsga2  — NSGA-II multi-objective Pareto front (per-variable trade-offs)

Documented as future work in CONCEPT.md (not yet implemented):
    mcmc, optimizers (nelder_mead/diffevo/agmip_stepwise), surrogate.
"""

from .glue import GlueResult, run_glue  # noqa: F401
from .nsga2 import Nsga2Result, run_nsga2  # noqa: F401
from .smc_pf import SmcResult, run_smc_pf  # noqa: F401


"""Pluggable calibration engines.

Each engine is a self-contained module with a ``run_*`` function and extensive,
plain-language docstrings. They share the framework's parallel evaluator
(``orchestrator.evaluate_thetas``) so every engine uses all cores.

Stages (the order they typically run in a pipeline — see ``method.preset``):

    sensitivity  — morris / sobol / anova : which parameters matter?      (screen)
    selection    — stepwise BIC/AICc (AgMIP)  : how many to calibrate?    (guard)
    sampling     — lhs / sobol / montecarlo / grid  (in ``samplers.py``)  (map)
    glue         — Monte-Carlo / GLUE pseudo-posterior                    (estimate)
    smc_pf       — sequential Monte-Carlo particle filter + MH move       (estimate)
    mcmc         — adaptive random-walk Metropolis posterior              (estimate)
    optimizers   — nelder_mead / diffevo : single best-fit point          (estimate)
    nsga2        — NSGA-II multi-objective Pareto front                   (estimate)
    surrogate    — GP/RF emulator acceleration (wraps any of the above)   (accelerate)

Optional dependencies (lazy-imported, only when that engine is used):
    sobol sensitivity -> SALib;  surrogate -> scikit-learn.
Everything else needs only NumPy/SciPy (already core dependencies).
"""

from .bayesopt import BayesOptResult, run_bayesopt  # noqa: F401
from .dream import run_dream  # noqa: F401
from .es_mda import run_es_mda  # noqa: F401
from .glue import GlueResult, posterior_summary, run_glue  # noqa: F401
from .mcmc import McmcResult, run_mcmc  # noqa: F401
from .nsga2 import Nsga2Result, run_nsga2  # noqa: F401
from .optimizers import OptimizerResult, run_optimizer  # noqa: F401
from .selection import SelectionResult, stepwise_select  # noqa: F401
from .sensitivity import (  # noqa: F401
    SensitivityResult,
    anova_variance_share,
    influential_params,
    run_morris,
    run_sensitivity,
    run_sobol,
)
from .smc_pf import SmcResult, run_smc_pf  # noqa: F401
from .surrogate import SurrogateResult, run_surrogate  # noqa: F401
from .recalibration import InSeasonRecalibrator  # noqa: F401
from .forcing import ForcingAssimilator  # noqa: F401
from .enkf import EnsembleKalmanFilter  # noqa: F401


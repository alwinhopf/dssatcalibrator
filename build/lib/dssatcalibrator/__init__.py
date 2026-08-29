"""dssatcalibrator — Monte-Carlo & Bayesian calibration of DSSAT-CSM.

New users: read ``WALKTHROUGH.md`` (plain-language guide). Architecture/design:
``CONCEPT.md`` (its §0 is the feature-status table).

Module map:
    config        — load the YAML config; enumerate active parameters
    observations  — read observed data (FileA/FileT, long-format CSV)
    spaces        — the parameter search space (bounds + start)
    samplers      — design-of-experiment samplers (lhs / sobol / montecarlo / grid)
    priors        — prior distributions (uniform / normal / lognormal / triangular)
    spawn         — run one DSSAT simulation for a parameter set
    runner        — run many spawns in parallel (thread-pool over subprocesses)
    objective     — align sim-vs-observed, metrics, and the weighted score / likelihood
    impact        — one-at-a-time real-DSSAT parameter impact atlas
    orchestrator  — the pipeline driver: presets -> screen -> select -> estimate
    engines/      — calibration engines (glue, smc_pf, mcmc, optimizers, sensitivity,
                    selection, surrogate, nsga2)
    viz           — figures + CSV report
"""

__version__ = "0.0.1"

from .config import active_parameters, load_config, validate_config

__all__ = ["__version__", "active_parameters", "load_config", "validate_config",
           "calibrate", "calibrate_fixed_design", "validate_loeo"]


def __getattr__(name):
    if name in {"calibrate", "calibrate_fixed_design", "validate_loeo"}:
        from . import orchestrator
        return getattr(orchestrator, name)
    raise AttributeError(name)

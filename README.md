# dssatcalibrator

> **AI agents & maintainers:** Read [`../AGENTS.md`](../AGENTS.md) before editing this repo.

`dssatcalibrator` is a Monte-Carlo and Bayesian calibration framework for the DSSAT-CSM crop model. It calibrates cultivar (and optionally management / soil / weather) parameters against observed crop data — fitting **LAI, biomass, grain yield, and phenology** jointly, across one or many experiments.

> **New here? Read [`WALKTHROUGH.md`](WALKTHROUGH.md)** — a step-by-step, plain-language guide that takes you from "what is calibration?" to a finished run, with no statistics background assumed. [`CONCEPT.md`](CONCEPT.md) is the deeper architecture/design document (its §0 is the feature-status table).

## Key Features

- **One config, pluggable engines.** Pick a *preset* and the framework runs the right pipeline:
  - **A** `morris → lhs → smc_pf` — screen, map, then a Bayesian **particle filter** (posterior + credible intervals).
  - **B** `morris → diffevo` — screen, then a global **optimiser** for a single best-fit point.
  - **C** `lhs → glue` *(default)* — **GLUE** pseudo-posterior from one big parallel batch; simplest.
  - **D** `morris → sobol → mcmc` — full **MCMC** posterior (pair with the surrogate for expensive crops).
- **Engines**: sampling (LHS / Sobol / Monte-Carlo / grid), GLUE, SMC-PF, MCMC, Nelder-Mead & differential-evolution optimisers, NSGA-II multi-objective, Morris/Sobol sensitivity screening, AgMIP stepwise BIC/AICc selection, and GP/RF surrogate acceleration.
- **Priors that count**: declare `uniform` / `normal` / `lognormal` / `triangular` priors per parameter; the Bayesian engines use them.
- **Honest objective**: RMSE/nRMSE/MBE/Willmott-d/EF/R² metrics, four weighting modes (`unified`, `sigma`, `count_scale`, `user`), `agmip_wls` reweighting, and optional `obs_autocorr` down-weighting of dense time-series.
- **Validation**: leave-one-environment-out cross-validation.
- **Multi-source & in-season**: pluggable observation adapters (satellite, UAV, IoT, farm software, field) fused by inverse-variance/priority, plus an **in-season recalibration** mode that re-estimates parameters as data arrives (`--assimilate` / `--combined`). See [`WALKTHROUGH.md`](WALKTHROUGH.md) §14 and [`CONCEPT.md`](CONCEPT.md) §17. *(EnKF/forcing state-assimilation modes are uncoupled prototypes, gated behind `allow_uncoupled`.)*
- **In-season LAI nowcast**: forecast LAI forward with an ensemble uncertainty band and last-observation anchoring (`--nowcast DATE --forecast`); optional NASA POWER weather driver with latency gap-fill. See [`WALKTHROUGH.md`](WALKTHROUGH.md) §15.
- **New crop / cultivar / species**: scaffold from an analog DSSAT module (`scaffold_crop.py`) with a gated `.SPE` writer, parameter **staging** (freeze what the data can't constrain), **identifiability/structural-adequacy** diagnostics (`--diagnostics`), and `year`/`site`/`random` cross-validation (`--cv-scheme`). See [`WALKTHROUGH.md`](WALKTHROUGH.md) §16.
- **Shared-stack plumbing**: optional `execution.backend: dssatengine` delegates DSSAT spawning and `DSSBatch.V48` writing to the shared engine; optional `weather.provider: dssatutils` / `soil.provider: dssatutils` acquire new-site inputs through the shared download layer. The calibration-specific writers, PlantGro/Evaluate parsers, objective, and engines stay local.
- **Parallel by default**: every engine fans its DSSAT runs across all cores (`num_cores`).
- **Visualization**: posterior distributions, observed-vs-simulated fits, sensitivity tornado, MCMC traces, ESS trajectory, Pareto front.

## Installation

Requires Python 3.10+.

```bash
pip install -e .          # core (numpy/scipy/pandas/matplotlib/pymoo) — enough for
                          # sampling, GLUE, SMC-PF, MCMC, optimisers, Morris screening
pip install -e .[shared]  # + pinned dssatengine@v0.3.0 execution backend
pip install -e .[acquire] # + pinned dssatutils@v0.2.0 weather/soil acquisition
pip install -e .[full]    # + SALib (Sobol sensitivity) and scikit-learn (surrogate)
pip install -e .[dev]     # + pytest
```

## Usage

Run calibration using the command-line entry point:

```bash
# Quick GLUE run with 300 samples (preset C)
python run_calibration.py config_hemp.yaml --preset C --n 300

# Bayesian particle filter with 250 particles (preset A)
python run_calibration.py config_hemp.yaml --preset A --n-particles 250

# Single best-fit via differential evolution (preset B)
python run_calibration.py config_hemp.yaml --optimizer diffevo

# Full MCMC posterior (preset D)
python run_calibration.py config_hemp.yaml --bayesian-engine mcmc

# Add a screening stage first (keep only influential parameters)
python run_calibration.py config_hemp.yaml --sensitivity morris

# AgMIP stepwise selection / surrogate acceleration
python run_calibration.py config_hemp.yaml --select bic
python run_calibration.py config_hemp.yaml --surrogate gp

# Subset of experiments / leave-one-environment-out cross-validation
python run_calibration.py config_hemp.yaml --n 50 --experiments YUKU2101 YUFE2201
python run_calibration.py config_hemp.yaml --validate
```

Calibration outputs (figures, design matrices, and best parameters) are written to `results/<calibrator_name>/` or custom directory paths.

## Testing

Run unit tests via `pytest`:

```bash
# Run fast offline tests
python -m pytest -m "not slow"

# Run all tests (including slow E2E tests that run DSSAT CSM)
python -m pytest
```

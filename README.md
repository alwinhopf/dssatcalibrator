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
- **Engines**: sampling (LHS / Sobol / Monte-Carlo / grid), GLUE, SMC-PF, MCMC, **DREAM (DE-MC)** and **ES-MDA** ensemble posteriors, Nelder-Mead / differential-evolution / **CMA-ES** optimisers, **Bayesian optimisation** (GP + Expected Improvement), NSGA-II multi-objective, Morris/Sobol sensitivity screening, AgMIP stepwise BIC/AICc selection, and GP/RF surrogate acceleration.

  Pick the estimator with `method.bayesian.engine` (`glue` | `smc_pf` | `mcmc` | `dream` | `es_mda` | `bayesopt`) or `method.optimizer.engine` (`nelder_mead` | `diffevo` | `cmaes`). Each engine is one entry in the estimator registry, so adding another is a function + one registration. Rules of thumb: **CMA-ES** for the fastest single best-fit, **DREAM** for an honest posterior on correlated/multimodal coefficients, **ES-MDA** for many parameters with uncertainty in a few iterations, **Bayesian optimisation** when each DSSAT run is expensive.
- **Priors that count**: declare `uniform` / `normal` / `lognormal` / `triangular` priors per parameter; the Bayesian engines use them.
- **Pooled calibrations**: keep species parameters shared across all experiments by default, or set `scope: experiment` / `pooling: per_experiment` on cultivar/ecotype parameters to give each experiment its own coefficient while fitting them in one run.
- **Parameter impact atlas**: run real-DSSAT one-at-a-time sweeps across discovered cultivar/ecotype/species coefficients plus management, soil, and weather parameters, collecting broad `*.OUT` tables, compact impact summaries, and a support-package capability map.
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

# One-at-a-time real-DSSAT impact atlas
python run_impact_atlas.py config_hemp.yaml --experiments UFCI2101 --discover-genotype --allow-species --max-per-group 1 --cores 2 --no-long

# R front end to the same atlas runner
Rscript -e "library(dssatcalibrator); run_impact_atlas('config_hemp.yaml', experiments='UFCI2101', discover_genotype=TRUE, allow_species=TRUE, max_per_group=1, num_cores=2, write_long=FALSE)"
```

Calibration outputs (`design.csv`, `best_theta.json`, `objective_breakdown.csv`, `manifest.csv`/`.json`, and summary tables) are written to `results/<calibrator_name>/` or custom directory paths; figures go under `figures/<calibrator_name>/` by default.

## Testing

Run unit tests via `pytest`:

```bash
# Run fast offline tests
python -m pytest -m "not slow"

# Run all tests (including slow E2E tests that run DSSAT CSM)
python -m pytest
```

## R interface (full parity)

`dssatcalibrator` ships as a **dual-language package**: every public function has an
R twin with the same name, and both languages read the same `config.yaml`. The
layout mirrors the workspace's other shared packages (`dssatengine`, `dssatutils`):

```
python/dssatcalibrator/   # the Python package
R/                        # the R twin (DESCRIPTION / NAMESPACE / R/*.R)
config_*.yaml             # one config, read by BOTH languages
run_calibration.{py,R}    # mirrored CLIs
scaffold_crop.{py,R}      # mirrored CLIs
```

```r
# install (from the repo root)
# install.packages("remotes"); remotes::install_local(".")
library(dssatcalibrator)
cfg <- load_config("config_hemp.yaml")
result <- calibrate(cfg)                 # preset C / GLUE by default
result$best_theta
```

```bash
Rscript run_calibration.R config_hemp.yaml --preset A --n-particles 250
Rscript run_calibration.R config_hemp.yaml --validate --cv-scheme year
Rscript run_calibration.R config_hemp.yaml --nowcast 2021-07-15 --forecast
```

**Engines.** GLUE, AgMIP stepwise selection, Nelder-Mead/`DEoptim`, CMA-ES,
Morris screening, the SMC particle filter, adaptive-Metropolis MCMC, DREAM (DE-MC),
ES-MDA, Bayesian optimisation (`DiceKriging`), NSGA-II (`mco`), and GP/RF
surrogates (`DiceKriging`/`ranger`) are all available in R. This R path **supersedes
the former `DSSAT_Calibration` repo** (its CroptimizR / AgMIP-stepwise capability now
lives here with full Python parity).

**LAI assimilation mode.** In-season LAI data assimilation (formerly the planned
`DSSAT_LAI_Assimilation` repo) is a **mode of this package**, not a separate project:
a satellite LAI observation source (`sentinel2_lai` / `modis_lai`) feeds the coupled
recalibration path (`assimilate()` with `mode: recalibration`) or the `nowcast()` +
`forecast_lai()` forward-LAI product. The uncoupled EnKF/forcing prototypes remain
gated behind `assimilation.allow_uncoupled: true`.

## Parity testing

R↔Python parity is a contract. Python is the source of truth; golden fixtures are
generated from it and checked by the R suite:

```bash
PYTHONPATH=python python tests/generate_parity_fixtures.py   # regenerate goldens
Rscript -e 'devtools::test()'                                # run R parity tests
python -m pytest -m "not slow"                               # run Python tests
```

Deterministic surfaces (config, priors' log-density, objective metrics, the
fixed-width DSSAT file writers, output parsers, `theta_hash`, GLUE post-processing)
are checked to machine precision. Stochastic engines (MCMC, SMC-PF, NSGA-II,
surrogate) are validated statistically, since RNG streams cannot match bit-for-bit
across languages.

## Dependencies (R)

R packages: `yaml` (core); plus, per engine, `lhs` + `randtoolbox` (sampling),
`DEoptim` (differential evolution), `sensitivity` (Sobol), `BayesianTools` (MCMC
helpers), `mco` (NSGA-II), `DiceKriging`/`ranger` (surrogate), `ggplot2` (figures),
`digest` (theta hashing), `jsonlite` (state/IO), and the workspace packages
`dssatengine` (execution backend) and `dssatutils` (weather/soil acquisition).
Pins follow the workspace `DEPENDENCIES.md`: `dssatengine@v0.3.0`,
`dssatutils@v0.2.0`.


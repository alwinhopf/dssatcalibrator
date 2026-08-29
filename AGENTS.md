# AGENTS.md — dssatcalibrator

> **Workspace context:** Read the root [`../AGENTS.md`](../AGENTS.md) first. This document
> holds rules and guidance specific to the `dssatcalibrator` repository.

## 1. Role in the Workspace

`dssatcalibrator` is the **canonical Monte-Carlo and Bayesian calibration framework** for DSSAT-CSM:
- It calibrates cultivar, ecotype, and species coefficients (plus optional management, soil, and weather parameters) against observed crop data (LAI, biomass, yield, phenology) across one or many experiments.
- It provides in-season observation assimilation and LAI nowcast forecasting.
- It supersedes the former `DSSAT_Calibration` and `DSSAT_LAI_Assimilation` repositories.

## 2. 1:1 R ↔ Python Parity Contract

`dssatcalibrator` maintains strict 1:1 parity between its R and Python implementations:
- R implementation: `R/` (22 mirrored files: `config.R`, `spaces.R`, `priors.R`, `dssat_io.R`, `writers.R`, `spawn.R`, `runner.R`, `observations.R`, `objective.R`, `fusion.R`, `acquisition.R`, `weather.R`, `forecast.R`, `diagnostics.R`, `scaffold.R`, `viz.R`, `orchestrator.R`, `cv.R`, `sparse.R`, `eval_cache.R`, `engines.R`, `sources.R`).
- Python implementation: `python/dssatcalibrator/` (matching modules and engine subpackages).
- Both languages read the **same `config.yaml`** and observation datasets.
- Mirrored CLIs: `run_calibration.py` / `run_calibration.R` and `scaffold_crop.py` / `scaffold_crop.R`.

## 3. Critical Implementation Rules & Common Pitfalls

### A. Parameter Bounds & Physiological Plausibility
- Parameter search spaces (`spaces.py` / `spaces.R`) must enforce strict physical and biological bounds.
- Never allow calibration engines to sample unphysical values (e.g. negative thermal times, harvest index > 1.0, extinction coefficients outside [0.1, 1.2]).
- Never tune cultivar coefficients to compensate for uncalibrated soil, bad weather, or erroneous management inputs.

### B. Genotype File Formatting (.CUL, .ECO, .SPE)
- Cultivar and ecotype files use strict fixed-width Fortran tables.
- Parameter substitution via `writers.py` / `writers.R` must preserve column positions, header comments, and trailing field definitions.
- Always isolate parameter mutations in temporary run directories; **never overwrite base files in `DSSAT48/`**.

### C. Validation & Overfitting Prevention
- Calibrate in stages: phenology first (P1, P2, P5), then vegetative growth/LAI (G1, G2), then grain filling and yield.
- Never report calibration performance based only on training fit. Always run leave-one-environment-out cross-validation (`validate_loeo` / `--validate`) or k-fold spatial/temporal cross-validation.

### D. Shared-Stack Integration
- Optional `execution.backend: dssatengine` delegates DSSBatch writing and DSSAT execution to the shared `dssatengine` library.
- Optional `weather.provider: dssatutils` / `soil.provider: dssatutils` acquire inputs for new calibration sites through `dssatutils`.
- Keep these integrations clean and non-breaking when dependencies are loaded or absent.

### E. Evaluation Caching
- `eval_cache.py` / `eval_cache.R` cache simulation objective evaluations indexed by `hash(parameter_vector, experiment)`.
- Invalidate or clear the cache when modifying experiment files, observation weights, or objective functions.

## 4. Verification & Testing

```bash
# Python tests
pytest tests/test_config.py tests/test_dependency_pins.py
pytest tests/test_new_features.py tests/test_impact_atlas.py

# R tests
Rscript -e "testthat::test_dir('tests/testthat')"
```

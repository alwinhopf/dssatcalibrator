# dssatcalibrator — Concept & Architecture

**Status:** implemented prototype plus architecture notes
**Author:** drafted 2026-06-13
**Audience:** maintainers of the DSSAT crop-modeling workspace

A general-purpose framework for calibrating **DSSAT-CSM** to a new *environment*,
*crop*, or *scenario* using **Monte-Carlo sampling** and **Bayesian** techniques. It
spawns large numbers of DSSAT runs with perturbed settings (genetic coefficients,
management, initial soil/water conditions), runs them in parallel, and identifies the
parameter combination(s) that best reproduce a central observations file containing
**LAI, biomass, grain yield, and phenology** — for one experiment or many across sites,
years and crops.

This document started as the design brief and now records both the target
architecture and the implementation choices in the current prototype. The Python
package exists in `dssatcalibrator/`; see `README.md` for the runnable status and
**`WALKTHROUGH.md` for a step-by-step, plain-language guide** aimed at first-time
users.

---

## 0. Implementation status (updated 2026-06-18)

Earlier drafts of this document over-claimed ("all four presets implemented") while
only GLUE + SMC-PF existed. The table below is the **honest, current** map of
feature → module → how you switch it on. Everything marked *implemented* has fast
offline tests under `tests/`.

| Capability | Status | Module | Turn it on with |
|---|---|---|---|
| Sampling: LHS / Sobol / Monte-Carlo / grid | implemented | `samplers.py` | `method.sample.engine` |
| Priors (uniform / normal / lognormal / triangular) | implemented | `priors.py` | per-parameter `prior:` block |
| Objective + metrics (RMSE/nRMSE/MBE/d/EF/R²) | implemented | `objective.py` | always on |
| Weighting: unified / sigma / count_scale / user | implemented | `objective.py` | `objective.weighting` |
| Weighting: `agmip_wls` (iterative reweighted LS) | implemented | `orchestrator.py` | `objective.weighting: agmip_wls` |
| `obs_autocorr` time-series down-weighting | implemented | `objective.py` | `objective.obs_autocorr: true` |
| GLUE (preset C) | implemented | `engines/glue.py` | `method.bayesian.engine: glue` |
| SMC particle filter + MH (preset A) | implemented | `engines/smc_pf.py` | `method.bayesian.engine: smc_pf` |
| MCMC posterior (preset D) | implemented | `engines/mcmc.py` | `method.bayesian.engine: mcmc` |
| Optimisers: Nelder-Mead, differential evolution (preset B) | implemented | `engines/optimizers.py` | `method.optimizer.engine` |
| Sensitivity: Morris (NumPy), Sobol (SALib) | implemented | `engines/sensitivity.py` | `method.sensitivity.active: true` |
| ANOVA factor variance-share | implemented (helper) | `engines/sensitivity.py` | `anova_variance_share()` |
| AgMIP stepwise BIC/AICc selection | implemented | `engines/selection.py` | `method.select.active: true` |
| Surrogate (GP / Random-Forest) acceleration | implemented | `engines/surrogate.py` | `method.surrogate.active: true` |
| NSGA-II multi-objective Pareto | implemented | `engines/nsga2.py` | `method.multiobjective.engine: nsga2` |
| Preset pipelines A/B/C/D + custom | implemented | `orchestrator.py` | `method.preset` |
| Leave-one-environment-out validation | implemented | `orchestrator.py` | `--validate` |
| Parallel execution (thread-pool over DSSAT subprocesses) | implemented | `runner.py` | `calibrator.num_cores` |
| Serial warm-up → parallel schedule | hook present (no-op here) | `runner.py` | `run_many(..., warmup=k)` |
| New-site weather/soil acquisition via `dssatutils` | implemented (optional extra; real experiments unchanged) | `weather.py`, `acquisition.py` | `weather.provider: dssatutils`, `soil.provider: dssatutils` |
| R parity layer (wrap CroptimizR) | **future** (Python-only today) | — | — |
| Multi-source observation adapters (satellite/UAV/IoT/farm/field) | implemented | `sources/` | `observation_sources:` block |
| Multi-source fusion (keep_all / inverse-variance / priority) | implemented | `fusion.py` | `fusion.conflict_resolution` |
| In-season **recalibration** (coupled: re-estimate params as data arrives) | implemented | `engines/recalibration.py` | `assimilation.mode: recalibration` + `--assimilate` |
| In-season **EnKF / forcing** state assimilation | **prototype — UNCOUPLED** (no DSSAT state re-injection; gated by `allow_uncoupled`) | `engines/enkf.py`, `engines/forcing.py` | `assimilation.mode: enkf\|forcing`, `allow_uncoupled: true` |
| In-season LAI **forecast/nowcast** (ensemble band + anchor continuity) | implemented | `forecast.py` | `forecast.active: true` / `--forecast` |
| Operational **nowcast** (as-of date, persist, warm-start, forecast) | implemented | `orchestrator.nowcast` | `--nowcast YYYY-MM-DD` |
| **Weather driver** layer (file / NASA POWER / dssatutils) + gap-fill | implemented (file default; dssatutils is optional) | `weather.py` | `weather.provider` |
| **Planting date** ingestion → FileX PDATE (input, not fit) | implemented | `spawn.py` + `observations.planting_dates` | `management_options.use_source_planting_date` |
| Satellite **cloud masking** + LAI observation operator | implemented | `sources/satellite.py` | `max_cloud_fraction`, `obs_operator` |
| New-crop **scaffolding** from an analog module (+ starter config) | implemented | `scaffold.py`, `scaffold_crop.py` | `python scaffold_crop.py …` |
| **`.SPE` writer** for new-species adaptation (gated) | implemented | `writers.edit_species` | `group: genetic_species` + `gating.species: free` |
| **Identifiability + structural-adequacy** diagnostics | implemented | `diagnostics.py` | `diagnostics.active: true` / `--diagnostics` |
| **Parameter staging** (freeze groups/params) | implemented | `orchestrator._apply_staging` | `method.staging.freeze_groups` |
| Cross-validation schemes (loeo / **year / site / random**) | implemented | `orchestrator.validate_cv` | `--validate --cv-scheme …` |
| In-season yield/maturity **forecasting** (true NWP-driven look-ahead) | **future** (LAI nowcast done; needs a forecast-weather provider) | — | — |

**Optional dependencies** (lazy-imported, only when that engine runs): Sobol needs
`SALib`; the surrogate needs `scikit-learn`. Install both with
`pip install -e .[full]`. Everything else uses only NumPy/SciPy/pymoo.

---

## 1. Motivation — what's new vs. the existing repos

`DSSAT_ML_Phenology_Prediction` already contains the *seed* of this framework: its
`01_particle_filter.R` spawns a 100-particle ensemble where every particle perturbs
genetic coefficients (P1V, P1D, P5, G1–G3, PHINT, P1–P4), management (planting-date
offset, planting depth), initial water (`sh2o_mult`), soil (SALB, SLDR) and weather
(`srad_mult`), runs them in parallel via `run_parallel()`, and applies a sequential
**SMC particle filter + Metropolis-Hastings (SMC-MCMC)** likelihood update against
observations. That is a working Bayesian calibrator — but it is **narrow** in five ways
this framework removes:

| Limitation today (phenology repo) | dssatcalibrator generalizes to |
|---|---|
| Single target: **phenology (BBCH/GSTD)** only | Multi-variable: **LAI, biomass, grain yield, phenology** (and any DSSAT output) |
| Single crop: **CERES-Wheat**, hard-coded coefficient block | **Any DSSAT crop/model**; parameters declared in config, not code |
| Parameter set & ranges **hard-coded** in the particle loop | **Central config** lists every varyable parameter with range + start + prior; activation is a flag |
| One method (**SMC-MCMC**) baked into the loop | **Pluggable engines**: sensitivity screening → MC/LHS/GLUE → Bayesian (SMC/MCMC) → optimizers |
| Observation = one BBCH series per field-season | **Long-format observation store**: time-series *and* end-of-season scalars, multi-experiment, per-point weights & error model |

The goal is to lift that proven spawning/likelihood machinery out of the wheat-specific
script and turn it into a configurable, crop-agnostic, multi-objective calibration
engine that sits beside `dssatengine` and `dssatutils` as a reusable package.

---

## 2. Position in the workspace stack

```
                 ┌─────────────────────────────────────────────┐
                 │              dssatcalibrator                 │  ← NEW
                 │  config → sample → spawn → run → score → fit │
                 └───────┬───────────────┬───────────────┬─────┘
                         │ spawns/parses │ weather/soil  │ likelihood/sampling
                         ▼               ▼               ▼
            ┌────────────────┐  ┌────────────────┐  ┌──────────────────────┐
            │  dssatengine   │  │   dssatutils   │  │  SMC-PF / MCMC / GLUE │
            │ run + parse    │  │ .WTH / .SOL    │  │ Morris / Sobol / LHS  │
            │ FileX/DSSBatch │  │ download layer │  │ (calibration backends)│
            └───────┬────────┘  └────────────────┘  └──────────────────────┘
                    ▼
            ┌────────────────┐
            │ DSSAT48/dscsm  │  compiled DSSAT-CSM v4.8 binary (read-only)
            └────────────────┘
```

**Reused, borrowed, or targeted for integration:**

- **`dssatengine`** (`engine.py` / `engine.R`) — the workspace run/parse layer.
  `execution.backend: dssatengine` delegates `DSSBatch.V48` writing and DSSAT
  spawning to the shared public API (`run_dssat`, `write_dssbatch`,
  `normalize_treatment_list`). The calibrator still parses `PlantGro.OUT`,
  `Evaluate.OUT`, and `Summary.OUT` through `dssatcalibrator.dssat_io`, because
  the shared engine's summary-CSV contract does not cover daily calibration outputs.
- **`dssatutils`** — the workspace weather/soil acquisition layer
  (`process_weather_*`, `process_soils_*`). Real hemp experiments still use files
  already present in the DSSAT install / experiment folder. New-site/synthesized
  experiments can opt into `weather.provider: dssatutils` and
  `soil.provider: dssatutils` when the optional `acquire` extra is installed.
- **`DSSAT_Gridded_Run_Tutorial`** — donor of the **spawn-across-settings orchestration**,
  which is the calibrator's outer loop almost verbatim:
  - `run_experiment.R` / `config_loader.{R,py}` — **non-invasive config injection**: write
    a private *merged* config per scenario to a temp file, run the pipeline as a fresh
    subprocess with `DSSAT_CONFIG_FILE` pointing at it, never touching the user's
    `config.yml`. Layered merge `config < base < factor values`, unique `project_name`
    per scenario so parallel workers never collide, a **serial warm-up pass then parallel
    remainder** (populate shared download caches once, then workers only read them),
    pre-flight coverage validation, `reuse_existing` resumability, `dry_run`. The
    calibrator reuses this exact orchestration for its θ-spawns.
  - `validate_against_observed.R` — already the **sim-vs-observed layer**: `parse_filea`
    (read observed `HWAM` from a DSSAT FileA), `parse_filex` (treatments / planting year /
    WSTA / soil id / coordinates), `adapt_template` (insert the pipeline's
    `ID_FIELD/WSTA/SOIL_ID/LATITUDE/LONGITUDE` placeholders into a real FileX), feasible
    weather/soil `select_sources` + `source_status` pre-flight, and a `metrics()` function
    (RMSE, nRMSE, MBE/bias, Willmott *d*, modelling efficiency EF, R²) with 1:1 scatter and
    RMSE-bar figures. This is the calibrator's observation reader, objective metrics, and
    report — already prototyped, for yield, against the bundled DSSAT example experiments.
  - `dssat_main_pipeline.{R,py}` — the 3-step points→weather+soil→FileX→run→parse engine
    (now `dssatengine`) used to set up each experiment's inputs.
- **`DSSAT_ML_Phenology_Prediction`** — donor of concrete, working code patterns:
  `run_parallel()` (cross-platform Windows-socket/Unix-fork pool), `modify_cultivar_file()`,
  `modify_ecotype_file()`, `modify_soil_sol()`, `parse_soil_scom()`, `compile_weather_file()`,
  `cache_soil_profiles()`, `date_to_yydoy()`, and the whole SMC-PF + Metropolis-Hastings
  update in `01_particle_filter.R` (becomes the `smc_pf` engine).
- **`DSSAT_Calibration`** (`dssatcal`, in the workspace index but **not checked out
  locally**) — conceptual sibling that already does the full protocol: **Morris → Sobol
  screening → stepwise AICc/BIC selection → optional DREAM-zs Bayesian UQ**, via a
  CroptimizR-compatible `wrapper.R`, with **parameter-file gating** (`.CUL` free, `.ECO`
  gated, `.SPE` blocked by default). Two ideas adopted directly: this engine pipeline
  ordering, and the **gating safety model** (a config-level allow/deny per genotype file
  so a user can't silently calibrate species-wide constants). Its routines become an
  alternative optimizer engine and a results cross-check where present.
- **`pythia`** — optional large-scale parallel executor; relevant only if calibration is
  scaled to HPC/grid. Not a core dependency.

---

## 3. Design principles

1. **Config is the contract.** Everything a user might vary lives in one file with a
   `min`/`max`/`start`/`prior` and an `active` flag. Switching crop, environment or
   scenario, and turning parameters on/off, is editing config — never code.
2. **Spawn = pure function of a parameter vector.** Given a vector θ and a base
   experiment, the spawner deterministically materializes a run directory. Same θ →
   byte-identical inputs → cacheable.
3. **Engines are pluggable and stageable.** Sensitivity screening, Monte-Carlo/LHS,
   Bayesian (SMC/MCMC), and classical optimizers share one interface
   (`propose θ → run → score`) so they compose into a pipeline and can be swapped.
4. **Objective is multi-variable and multi-experiment by construction.** A single
   scalar/log-likelihood aggregates LAI + biomass + yield + phenology across all
   experiments with explicit per-variable weights and an error model.
5. **Crop-agnostic core, crop-specific only at the file writers.** The one place that
   knows CUL/ECO column layouts is a small, table-driven writer module.
6. **Reproducible & resumable.** Seeded sampling, hashed run cache (mirrors the
   phenology repo's `plantgro.csv` cache check), and a manifest of every spawn.

---

## 4. The central config file (heart of the request)

A single `config.yaml`. Precedence mirrors the phenology repo:
`environment variable > config.yaml > built-in default`. It declares **every** parameter
a user could ever vary; only those with `active: true` enter the calibration. Inactive
parameters are held at `start` (their value still documents the run).

Parameters are grouped by **where they are written in DSSAT**, because that determines
the file writer used:

```yaml
# ============================================================================
# dssatcalibrator config — declares EVERY varyable parameter. `active` decides
# which actually enter calibration; the rest are fixed at `start`.
# ============================================================================

calibrator:
  name: "carinata_DE_2024"
  seed: 42
  workdir: "runs/"            # per-spawn run directories live here
  results_dir: "results/"
  dssat_exe: ""               # blank => platform auto-detect (as phenology repo)
  dssat_dir: ""
  num_cores: 0                # 0 => physical cores - 2
  cache_spawns: true          # skip a spawn whose output already exists (hash-guarded)
  keep_run_dirs: false        # true = keep run folders for audit

# --- which DSSAT outputs are scored, and from which file --------------------
engine:
  run_mode: "experiment"      # experiment | sequence  (dssatengine modes)
  timeseries_outputs:         # parsed from PlantGro.OUT/plantgro.csv (daily)
    LAI: "LAID"               # leaf area index
    biomass: "CWAD"           # tops dry weight kg/ha
    phenology: "GSTD"         # growth stage (BBCH/Zadoks or crop scale)
  summary_outputs:            # parsed from Summary.OUT/summary.csv (end of season)
    grain_yield: "HWAM"
    biomass_maturity: "CWAM"
    anthesis_date: "ADAT"
    maturity_date: "MDAT"

# --- calibration method ------------------------------------------------------
# All four named pipelines (§14a) are implemented; pick one with `preset`, or set
# preset: custom and compose the stages by hand below. preset just pre-fills the
# stage blocks; any explicit stage block you set overrides the preset.
method:
  preset: "A"        # A = morris->lhs->smc_pf (default) | B = morris->diffevo (fast)
                     # C = lhs->glue (quick, fully parallel) | D = morris->sobol->mcmc (full UQ)
                     # custom = use exactly the stage blocks below
  # Stage blocks (shown with preset A's values). Each stage is optional.
  sensitivity: { engine: "morris", active: true,  trajectories: 20 }   # morris | sobol | anova | none
  select:      { engine: "none" }   # stepwise_bic | stepwise_aicc | none  (AgMIP overfit guard, §16)
  sample:      { engine: "lhs",    active: true,  n: 500 }             # lhs | sobol | montecarlo | grid | none
  surrogate:   { engine: "none" }   # gp | rf | none — emulator-accelerated calibration (§16)
  bayesian:    { engine: "smc_pf", active: true,  n_particles: 200, ess_frac: 0.5,
                 mutation: "metropolis" }                              # smc_pf | mcmc | glue | none
  optimizer:   { engine: "none" }   # none | nelder_mead | diffevo | agmip_stepwise | gencalc
  validation:  { scheme: "none" }   # none | loeo (leave-one-environment-out) | year | random
  # `sensitivity`/`select` may auto-set `active` flags below (keep impactful params).

# --- objective / likelihood (weighting modes implemented; §14b) ---------------
objective:
  weighting: "unified"   # unified (default) | sigma | count_scale | user | agmip_wls
  # agmip_wls = iterative reweighted LS: run OLS once, set each variable group's
  # weight = 1/var(residuals), refit (AgMIP protocol). See §16.
  obs_autocorr: false    # true = down-weight dense time-series for serial correlation
  weights:   { LAI: 1.0, biomass: 1.0, grain_yield: 1.0, phenology: 1.0 }  # priority multipliers
  # default per-variable measurement error / scale (used as sigma by Bayesian engines,
  # as the normalizing scale otherwise); overridable per-observation via the CSV `sigma` col.
  error_model:
    LAI:         { type: relative, value: 0.15 }
    biomass:     { type: relative, value: 0.15 }
    grain_yield: { type: relative, value: 0.10 }
    phenology:   { type: absolute, value: 3 }   # days (or BBCH units)

# ============================================================================
# PARAMETERS — grouped by destination file. Each: active, min, max, start,
# prior (for Bayesian), optional transform/type, optional per-crop override.
# ============================================================================
parameters:

  # -- Genetic: cultivar coefficients -> <MODEL>.CUL (modify_cultivar_file) --
  genetic_cultivar:
    P1V:   { active: true,  min: 0,    max: 60,  start: 48.45, prior: {dist: normal, sd: 10} }
    P1D:   { active: true,  min: 0,    max: 150, start: 73.50, prior: {dist: normal, sd: 15} }
    P5:    { active: true,  min: 300,  max: 700, start: 505.0, prior: {dist: normal, sd: 50} }
    G1:    { active: true,  min: 10,   max: 50,  start: 35.42, prior: {dist: normal, sd: 5} }
    G2:    { active: true,  min: 10,   max: 80,  start: 22.60, prior: {dist: normal, sd: 4} }
    G3:    { active: false, min: 0.5,  max: 8.0, start: 0.78,  prior: {dist: normal, sd: 0.15} }
    PHINT: { active: false, min: 50,   max: 150, start: 95.0,  prior: {dist: normal, sd: 10} }

  # -- Genetic: ecotype coefficients -> <MODEL>.ECO (modify_ecotype_file) ----
  genetic_ecotype:
    P1: { active: false, min: 200, max: 600, start: 400, prior: {dist: normal, sd: 40} }
    P2: { active: false, min: 150, max: 450, start: 285, prior: {dist: normal, sd: 30} }
    P3: { active: false, min: 100, max: 300, start: 190, prior: {dist: normal, sd: 20} }
    P4: { active: false, min: 100, max: 300, start: 200, prior: {dist: normal, sd: 20} }

  # -- Management -> FileX *PLANTING DETAILS / *SIMULATION CONTROLS ----------
  management:
    planting_date_offset: { active: true,  min: -10, max: 10, start: 0, type: int, unit: days }   # shifts PDATE
    plant_population:     { active: false, min: 50,  max: 400, start: 200, dssat: PPOP, unit: "plants/m2" }
    plant_emergence_pop:  { active: false, min: 50,  max: 400, start: 200, dssat: PPOE }
    row_spacing:          { active: false, min: 15,  max: 75,  start: 17,  dssat: PLRS, unit: cm }
    planting_depth:       { active: false, min: 2.0, max: 7.0, start: 4.0, dssat: PLDP, unit: cm }

  # -- Initial conditions -> FileX *INITIAL CONDITIONS -----------------------
  initial_conditions:
    initial_soil_water_mult: { active: true,  min: 0.6, max: 1.4, start: 1.0 }  # scales ICBL SH2O (cf. sh2o_mult)
    initial_no3_ppm:         { active: false, min: 0,   max: 30,  start: 5 }    # SNO3 per layer
    initial_nh4_ppm:         { active: false, min: 0,   max: 10,  start: 1 }    # SNH4 per layer
    initial_residue_kg_ha:   { active: false, min: 0,   max: 8000,start: 1000 } # ICRES

  # -- Soil profile -> .SOL (modify_soil_sol + layer multipliers) -----------
  soil:
    SALB: { active: false, min: 0.05, max: 0.30, start: 0.13 }   # albedo
    SLDR: { active: false, min: 0.10, max: 0.90, start: 0.50 }   # drainage rate
    SLLL_mult: { active: false, min: 0.8, max: 1.2, start: 1.0 } # lower limit scaler (all layers)
    SDUL_mult: { active: false, min: 0.8, max: 1.2, start: 1.0 } # drained upper limit scaler
    runoff_CN: { active: false, min: 50, max: 95, start: 75 }    # SLRO curve number

  # -- Weather perturbation -> .WTH (add_weather_noise) ----------------------
  weather:
    srad_mult:   { active: false, min: 0.8, max: 1.2, start: 1.0 }
    tmax_offset: { active: false, min: -2,  max: 2,   start: 0 }
    tmin_offset: { active: false, min: -2,  max: 2,   start: 0 }
    rain_mult:   { active: false, min: 0.7, max: 1.3, start: 1.0 }

# --- Discrete factors -> swept as a full factorial (run_experiment.R style) -
# Not continuous params: enumerated choices crossed with the θ samples. Lets the
# user ALSO calibrate/quantify which input dataset best matches observations
# (exactly what validate_against_observed.R does for yield).
discrete_factors:
  weather_source: { active: false, levels: ["NASA_POWER", "DAYMET", "OPEN_METEO"] }
  soil_source:    { active: false, levels: ["SSURGO", "SOILGRIDS_ONLINE"] }
  exclude: []     # drop specific (weather_source, soil_source) pairs

# --- Parameter-file gating (safety; from dssatcal: CUL free, ECO gated, SPE blocked)
gating:
  cultivar: "free"     # free | gated | blocked  — may write <MODEL>.CUL
  ecotype:  "gated"    # gated => allowed only if explicitly opted-in per run
  species:  "blocked"  # blocked => never written, even if a param targets it

# ============================================================================
# Crop blocks — one per crop in the observations. Selects model files + scale.
# (Generalizes the phenology repo's single `crop:` block to a list.)
# ============================================================================
crops:
  - code: "WH"                    # DSSAT crop code
    model: "WHCER048"             # CUL/ECO/SPE stem
    ecotype: "USWH01"
    cultivar_anchor: "IB0488"     # row edited in the .CUL (cf. modify_cultivar_file)
    filex_ext: "WHX"
    stage_scale: "bbch"           # bbch | ceres_int | ...  (phenology comparison scale)
```

**Notes on the schema**

- A parameter's **group** chooses its writer; its **`dssat`/`unit`** keys carry the exact
  DSSAT field name where ambiguous. This keeps the spawner table-driven.
- **Ranges + start + prior** satisfy the explicit request: a plausible range *and* a
  starting point, plus a prior distribution for the Bayesian stage (default
  `uniform(min,max)` if `prior` omitted; `start` is the prior mean for `normal`).
- **`active: false`** parameters are *declared but dormant* — the "all variables a user
  might vary even if only a few are activated" requirement.
- The **sensitivity stage** can flip `active` automatically (screen all, keep the most
  influential), so a user can start with everything declared and let the framework prune.
- An optional **`role: obligatory | candidate`** per parameter supports the AgMIP protocol
  (§16): *obligatory* params are estimated at every step; *candidate* params are added only
  if they lower the BIC/AICc — a principled overfit guard that supersedes pure sensitivity
  pruning when `method.select` is enabled.

---

## 5. The central observations file

One **long-format** CSV/Parquet handles single-experiment and multi-experiment cases
uniformly, mixing daily time-series (LAI, biomass) with end-of-season scalars (yield)
and phenological events (a stage reached on a date):

| column | meaning |
|---|---|
| `exp_id` | experiment key (site × year × treatment); ties an obs to a base FileX |
| `site`, `year`, `crop` | grouping / crop-block selector |
| `treatment` | DSSAT treatment number within the FileX (optional) |
| `variable` | one of `LAI`, `biomass`, `grain_yield`, `phenology` (extensible) |
| `date` *or* `dap` | observation date (or days-after-planting); for `phenology` this is the date the stage was observed |
| `value` | measured value (units per `variable`) |
| `stage` | for `phenology`: the observed BBCH/Zadoks/CERES stage |
| `sigma` | measurement uncertainty (defaults to config error model if blank) |
| `weight` | optional per-observation weight (defaults to 1) |

This directly subsumes the phenology repo's input (its `variety/lat/lon/planting_date`
become an `exp_id`; its `Obs_date`/`GS` become `variable=phenology` rows). Each `exp_id`
maps to a **base experiment**: either a real DSSAT FileX shipped with the dataset, or a
synthetic one the calibrator builds from site coordinates + planting date using
`dssatutils` (weather/soil) — exactly the `compile_weather_file` + `cache_soil_profiles`
path the phenology repo already uses.

**Native DSSAT observation ingestion (optional).** DSSAT experiments already carry their
own measurements: a **FileA** (`.??A`, end-of-season scalars — `HWAM`, `CWAM`, `ADAT`,
`MDAT`, …) and a **FileT** (`.??T`, in-season time-series — `LAID`, `CWAD`, growth stage
by date). When an `exp_id` points at a real DSSAT experiment, the calibrator can read its
observations straight from FileA/FileT instead of (or in addition to) the central CSV —
reusing `validate_against_observed.R::parse_filea` for the scalar side and a parallel
FileT reader for the time-series side. FileA → `grain_yield`/`biomass`/phenology-date
rows; FileT → `LAI`/`biomass` time-series rows. So the four target variables the user
named map cleanly onto DSSAT's two observation files, and a user calibrating against the
bundled example experiments needs no separate observations file at all. The central CSV
remains the path for field data that doesn't live in DSSAT files.

---

## 6. Parameter space, sampling, and the spawn model

### 6.1 Building the design

Active parameters define the search space Θ. Samplers (each an engine):

- **`lhs`** — Latin-Hypercube over the active ranges (space-filling MC).
- **`sobol`** — low-discrepancy quasi-MC (good for Sobol sensitivity indices).
- **`montecarlo`** — plain random draws from priors (GLUE-style).
- **`grid`** — full factorial (small/debug spaces).
- Bayesian engines draw their own proposals (priors → SMC particles / MCMC chain).

### 6.2 Two-level execution: per-experiment setup, then per-spawn runs

Calibration has a natural two-level structure that maps onto the two donor patterns:

**Level A — per-experiment setup (serial warm-up).** For each `exp_id`, *once*: resolve
the base FileX (read a real one via `parse_filex`, or build one; insert the per-point
placeholders via `adapt_template`), acquire its weather + soil — either from the
experiment's shipped files or, for a bare site, by picking a feasible source
(`select_sources` + `source_status` pre-flight) and downloading via `dssatutils`
(`process_weather_*` / `process_soils_*`), then caching. This mirrors `run_experiment.R`'s
**serial warm-up pass that populates the shared download caches before any parallel work**,
so the many θ-spawns never trigger cold-download races on the same site.

**Level B — per-spawn run (deterministic, parallel).** A spawn is a pure function of
(θ, exp_id). For each one:

1. Create `runs/<exp_id>/sample_<hash(θ)>/` with a unique name (hash + experiment), so
   parallel workers never collide — the `project_name` namespacing trick from
   `run_experiment.R`. The hash also drives the **spawn cache** (same θ skips, mirroring
   the phenology repo's `plantgro.csv` existence/size check and `reuse_existing`).
2. Materialize inputs from the **cached** Level-A weather/soil + base FileX (copy, don't
   re-download), plus DSSAT support files (`.CDE/.ERR/.L48/.CTR/.TXT`).
3. Write the perturbed inputs via the table-driven writers (subject to `gating`):
   - `genetic_cultivar` → `<MODEL>.CUL` (`modify_cultivar_file` generalized to a CUL column map)
   - `genetic_ecotype` → `<MODEL>.ECO` (`modify_ecotype_file`; only if `gating.ecotype != blocked`)
   - `management` + `initial_conditions` → FileX section edits (PDATE, PPOP/PPOE, PLRS,
     PLDP in *PLANTING; SH2O/SNO3/SNH4/ICRES in *INITIAL CONDITIONS)
   - `soil` → `.SOL` (`modify_soil_sol` + layer multipliers)
   - `weather` → `.WTH` (`add_weather_noise` generalized to T/rain offsets)
4. Write `DSSBatch.V48` and invoke `dscsm048`. `execution.backend: native` uses
   the local fallback; `execution.backend: dssatengine` delegates batch writing
   and spawning to the shared engine executor.
5. Parse outputs: `Summary.OUT` / `PlantGro.OUT` / `Evaluate.OUT` through
   `dssat_io.py`; this remains calibration-specific until the shared engine grows
   daily PlantGro and simulated-vs-measured Evaluate readers.

θ-perturbed genotype/management/initial-condition edits are **in-process file writes inside
a parallel worker** (the phenology-repo pattern) — far cheaper than launching a fresh
subprocess per spawn. The heavier `DSSAT_CONFIG_FILE` subprocess-injection pattern is kept
for the *outer* discrete-factor sweep (Level A scenarios), where it earns its cost.

### 6.3 Parallel execution & orchestration controls

The unit of parallelism is one Level-B spawn. The framework flattens (experiments ×
samples × discrete factors) into one task list and maps it across `num_cores` with the
proven cross-platform pool — `run_parallel()` (R, PSOCK on Windows / fork on Unix) or
`multiprocessing`/`joblib` over `_run_one_point` (Python). For HPC the same task list goes
to `pythia` / the `hpc/dssat_mpi_runner.py` work-stealing MPI runner.

Orchestration options are lifted straight from `experiment.yml` and exposed in the
`calibrator:`/`method:` blocks: **`max_parallel`** (with the serial-warm-up-then-parallel
schedule), **`reuse_existing`** (resume — skip cached spawns), **`stop_on_error`** (log and
continue vs. abort), **`dry_run`** (print the spawn plan and exit), and **`validate`**
(pre-flight coverage checks that drop infeasible discrete-factor combos, e.g. NASA_POWER
before 1984 or AGERA5 without a CDS key).

---

## 7. Calibration engines (pluggable backends)

All implement one interface: `propose() -> [θ]`, then receive back `score(θ)` (or the full
simulated-vs-observed table) and `update()`. They chain into a pipeline:

1. **Sensitivity screening** — `morris` (elementary effects, cheap) or `sobol` (variance
   decomposition) for **continuous** parameters; plus **`anova`** variance decomposition for
   **discrete** factors (management levels, input-source choices), reusing
   `run_experiment.R`'s ANOVA "share of variance per factor" routine almost verbatim.
   Output: parameter/factor influence ranking → optionally auto-set `active` flags so only
   impactful parameters are calibrated. (This is the recommended *first* stage for a new
   crop/environment, and the Morris→Sobol ordering `DSSAT_Calibration`/AgMIP use up front.)
2. **Monte-Carlo / GLUE** — run the LHS/Sobol/MC design, keep behavioural runs (those
   under an objective threshold), report posterior-like parameter distributions. Cheap,
   embarrassingly parallel, good first map of the space.
3. **Bayesian — `smc_pf`** — the **port of the phenology repo's SMC particle filter +
   Metropolis-Hastings mutation**, generalized from BBCH-only to the multi-variable
   likelihood (§8). Sequential assimilation of time-ordered observations, ESS-triggered
   resampling, parameter mutation with acceptance ratio — the algorithm already exists in
   `01_particle_filter.R`; this lifts it out and feeds it the general objective.
4. **Bayesian — `mcmc`** — full posterior via DREAM/`emcee`/affine-invariant samplers for
   richer uncertainty (credible intervals on coefficients). Heavier; surrogate-assisted
   (GP/RF emulator trained on the MC design) variant for expensive crops.
5. **Optimizers** — `nelder_mead`, `diffevo` (SciPy/`DEoptim`), or `agmip_stepwise`
   (the AgMIP phenology-then-growth ordering used by `dssatcal`/CroptimizR) for a single
   best-fit point estimate when full UQ isn't needed.

A typical run: **morris → lhs → smc_pf** (screen, map, then Bayesian refine).

---

## 8. Objective / likelihood (multi-variable, multi-experiment)

For a parameter vector θ, simulate every experiment, align simulated to observed
(by date/DAP for time-series; at maturity for scalars; by stage-hit date for phenology —
the phenology repo's "DAP when each particle first hits the target BBCH" logic), and
aggregate:

- **Per-observation residual** normalized by `sigma` (config error model: absolute or
  relative). Gaussian log-likelihood `-0.5·(sim−obs)²/σ²` — identical to the phenology
  filter's term, now summed over **all variables and experiments**.
- **Per-variable weights** (`weights: {LAI, biomass, grain_yield, phenology}`) balance
  units/counts so yield (1 point) isn't swamped by LAI (many points). Default weighting
  normalizes by observation count and variable variance.
- **Aggregate score** = weighted sum of normalized RMSE (for optimizers/GLUE) **or** total
  log-likelihood (for SMC/MCMC). Both come from the same residual table; the four
  `objective.weighting` modes (§14b — `unified`/`sigma`/`count_scale`/`user`) are all
  implemented and just change how that table is reduced to one number.
- Standard reported metrics per variable reuse `validate_against_observed.R::metrics()`
  verbatim: **RMSE, nRMSE%, MBE (bias), Willmott *d*, modelling efficiency EF, R²** — plus
  MAE and, for phenology, the ±n-day accuracy the phenology repo already reports. The same
  function backs both the optimizer objective and the final report, so "what's optimized"
  and "what's reported" never drift.

---

## 9. Outputs

- `posterior_parameters.csv` — calibrated values: MAP/mean + credible intervals (Bayesian)
  or best point (optimizer); GLUE behavioural distributions.
- **Updated genotype/soil files** — a ready-to-use `<MODEL>.CUL`/`.ECO` (and `.SOL`) with
  the calibrated coefficients written in, so the result drops straight into other pipelines.
- `objective_breakdown.csv` — fit per variable × experiment (RMSE/nRMSE/bias/d-index).
- `sim_vs_obs.csv` + plots — time-series overlays (LAI, biomass), 1:1 yield/phenology
  scatter, parameter posteriors / pair-plots, sensitivity tornado, convergence/ESS traces.
- `manifest.json` — every spawn (θ, paths, status, score) for full reproducibility/resume.

---

## 10. Workflow / pipeline steps

```
0  load config + observations (CSV and/or FileA/FileT); resolve crop blocks; auto-detect DSSAT
1  LEVEL A (serial warm-up): per exp_id resolve/adapt base FileX, acquire+cache weather+soil
                             (dssatutils), pre-flight validate discrete-factor feasibility
2  [sensitivity]  Morris/Sobol (continuous) + ANOVA (discrete) -> ranking -> (opt) set active
3  [sample]       LHS/Sobol/MC design over active params (× discrete factors)
4  LEVEL B (parallel): spawn + run all (experiments × samples × factors)  [spawn.py + runner.py]
                       resumable (skip cached), namespaced run dirs, stop_on_error policy
5  parse outputs (summary + PlantGro) -> sim table
6  score against observations (multi-variable likelihood / metrics())
7  [bayesian/optimizer] SMC-PF / MCMC / optimizer refine -> posterior / best-fit
8  write calibrated CUL/ECO/SOL, diagnostics, plots (1:1, RMSE, posteriors), manifest
```

Driven like the phenology repo's `run_pipeline.R --steps=… --skip=…`, so any stage can run
standalone or be skipped (e.g. skip sensitivity, go straight to SMC).

---

## 11. Proposed module layout

```
dssatcalibrator/
  config.yaml                 # the central config (§4) — the only file most users touch
  observations.csv            # the central observations store (§5)
  dssatcalibrator/            # package (Python primary; R parity per workspace convention)
    config_loader.{py,R}      # env > yaml > default loader (DSSAT_CONFIG_FILE-aware)
    observations.{py,R}       # central CSV + native FileA/FileT readers (parse_filea/parse_filet)
    setup.{py,R}              # Level A: resolve base FileX, acquire+cache weather/soil
    spaces.{py,R}             # active-param -> sampling space; transforms/bounds
    samplers.{py,R}           # lhs / sobol / montecarlo / grid
    spawn.{py,R}              # θ + cached base exp -> run dir (writers below, gating-aware)
    writers.{py,R}            # CUL/ECO/FileX/SOL/WTH editors (table-driven, crop-agnostic)
    orchestrate.{py,R}        # warm-up→parallel schedule, resume, dry-run, validate, namespacing
    runner.{py,R}             # parallel map over spawns -> native / optional dssatengine backend
    parse.{py,R}              # summary.csv + PlantGro -> tidy sim table
    objective.{py,R}          # align sim vs obs; metrics() (RMSE/nRMSE/MBE/d/EF/R²) + likelihood
    engines/
      sensitivity.{py,R}      # morris / sobol (continuous) + anova (discrete factors)
      montecarlo.{py,R}       # MC / GLUE
      smc_pf.{py,R}           # ← port of 01_particle_filter.R, generalized
      mcmc.{py,R}             # DREAM / emcee (+ surrogate option)
      optimize.{py,R}         # nelder_mead / diffevo / agmip_stepwise
    report.{py,R}             # metrics, 1:1 + RMSE plots, calibrated-file writeout
  run_calibration.{py,R}      # orchestrator (--steps / --skip)
  templates/                  # FileX templates + stock genotype files per crop
```

Follows workspace `CONVENTIONS.md`: snake_case repo name, R↔Python parity with mirrored
names, pins `dssatengine`/`dssatutils` to tags, ignores `runs/` + `results/` scratch.

---

## 12. Generalizing to a new crop / environment / scenario

- **New environment (same crop):** point `observations.csv` at the new site-years; the
  calibrator pulls weather/soil via `dssatutils`. No code change.
- **New crop:** add a `crops:` block (code, model stem, ecotype, cultivar anchor, FileX
  ext, stage scale) and drop `<MODEL>.CUL/.ECO/.SPE` + a FileX template into `templates/`.
  The CUL/ECO column maps live in `writers` as data tables, so only a table entry is
  needed for a new CERES/CROPGRO genotype layout — the one genuinely crop-specific piece,
  matching the phenology repo's "Adapting to a new crop" checklist.
- **New scenario:** flip `active` flags / edit ranges (e.g. calibrate management instead
  of genetics, or add soil hydraulics) — pure config.

---

## 13. Reuse map (what comes from where)

| Need | Source | Mechanism |
|---|---|---|
| Spawn + run + parse outputs | `spawn.py` / `dssat_io.py`, with optional `dssatengine` executor | `DSSBatch` subprocess, `dssatengine.run_dssat`, `parse_summary`, `parse_plantgro`, `parse_evaluate` |
| Parallel pool | phenology `utils.R` | `run_parallel()` (Win socket / Unix fork) |
| Edit CUL / ECO / SOL | phenology `utils.R` | `modify_cultivar_file`, `modify_ecotype_file`, `modify_soil_sol` |
| Weather/soil for a site | `dssatutils` | `process_weather_*`, `process_soils_*` |
| Per-experiment input setup | gridded `dssat_main_pipeline.{R,py}` | points→weather+soil→FileX (now `dssatengine`) |
| Spawn orchestration (warm-up→parallel, resume, dry-run, validate, no-clobber config) | gridded `run_experiment.R` + `config_loader.{R,py}` | `DSSAT_CONFIG_FILE` injection, `project_name` namespacing |
| Read base FileX / observed FileA / adapt template | gridded `validate_against_observed.R` | `parse_filex`, `parse_filea`, `adapt_template`, `select_sources`/`source_status` |
| Sim-vs-obs metrics (+ANOVA factor variance) | gridded `validate_against_observed.R` + `run_experiment.R` | `metrics()` (RMSE/nRMSE/MBE/d/EF/R²); ANOVA variance share |
| SMC particle filter + MCMC | phenology `01_particle_filter.R` | generalized into `engines/smc_pf` |
| Config (env>yaml>default) | phenology `config.R` + `config.yaml` | extended schema |
| AgMIP Morris→Sobol→stepwise→DREAM-zs + file gating | `DSSAT_Calibration` (`dssatcal`) | optional optimizer engine / cross-check; gating safety model |
| HPC scale-out (optional) | `pythia` / `hpc/dssat_mpi_runner.py` | hand the spawn task list to it |

---

## 14. Decisions

**Decided**

1. **Primary language: Python core + R parity.** `dssatengine`'s run/parse layer is most
   complete in Python and the MC/Bayesian/sensitivity ecosystem (SALib, emcee/PyMC, SciPy,
   Optuna) is richer there; an R parity layer follows the workspace convention. The donor
   SMC-PF code (R) is ported to Python as the reference implementation.
2. **Source of truth for an experiment:** the **real DSSAT FileX** when one exists for an
   `exp_id` (read via `parse_filex`/`adapt_template`, with observations taken natively from
   its FileA/FileT) — **otherwise a user-supplied observations CSV** (the long-format store
   of §5). For a CSV experiment with no shipped weather/soil, inputs are acquired from a
   feasible source via `dssatutils`.
3. **Spawn orchestration:** `run_experiment.R`'s warm-up→parallel schedule +
   `DSSAT_CONFIG_FILE` config injection + `reuse_existing` resume + `project_name`
   namespacing.
4. **Engine pipeline: all four (A–D) are implemented and user-selectable** via
   `method.preset` (plus `custom` for hand-composed stages). The default if you set
   nothing is **C** (`lhs → glue`, fully parallel, simplest); **A**
   (`morris → lhs → smc_pf`) is the recommended choice when you want uncertainty.
   A preset selects the *main estimator* automatically; the optional screening /
   selection / surrogate stages each turn on with their own `active: true`. See §14a.
5. **Weighting: implement all four modes as user-selectable** via `objective.weighting`.
   **Default = `unified`** (count+scale for optimizers/GLUE/report; σ for Bayesian engines;
   user `weights` always applied). See §14b.

### 14a. Engine pipeline presets (all implemented; `method.preset` selects)

The default for a user who sets nothing is **A**; every preset is built and selectable.

| Option | Pipeline | Gives you | Cost | Best when |
|---|---|---|---|---|
| **A — Screen → Map → Bayesian** *(recommended)* | `morris → lhs → smc_pf` | Posterior + credible intervals; native time-series assimilation | High (3 stages) | You want uncertainty and are fitting LAI/biomass/phenology over the season |
| **B — Screen → Optimize** | `morris → diffevo` (or `agmip_stepwise`) | Single best-fit coefficient set | Low–Med | "Just give me calibrated numbers"; AgMIP-style workflows |
| **C — Map → GLUE** | `lhs → glue` | Pseudo-posterior from one big batch | Low, fully parallel | Quick first answer; simplest to ship/debug; high-throughput |
| **D — Full UQ** | `morris → sobol → mcmc(DREAM)` | Rigorous posterior + Sobol indices | Very high | Publication-grade UQ (usually needs a surrogate) |

All four are built and selectable today. **C** (`lhs → glue`) is the zero-config
default because it is the simplest and embarrassingly parallel; **A** is the
recommended upgrade when you want a posterior with credible intervals; **B** is the
`fast` best-fit preset; **D** is the full-UQ preset (pair it with the surrogate for
expensive crops). Switching is one line: `method.preset`.

### 14b. Weighting modes (all implemented; `objective.weighting` selects)

How LAI + biomass + yield + phenology, with very different units and observation counts,
combine into one objective. The default for a user who sets nothing is **`unified`**; all
four are built and selectable.

| `weighting` | Rule | Pros | Cons |
|---|---|---|---|
| **`unified`** *(default)* | count+scale for optimizers/GLUE/report; Bayesian engines reinterpret each variable's `scale` as σ — same residual table both ways | Best of both; no extra user input needed; consistent | A composite (two interpretations of one table) |
| **`sigma`** | `Σ −0.5·((sim−obs)/σ)²`, σ per variable from `error_model` | Statistically correct; **what smc_pf/mcmc need**; encodes trust per variable | User must supply/assume σ |
| **`count_scale`** | per-variable mean of `((sim−obs)/scale)²`, averaged so point-rich vars don't dominate, then × `weights` | Robust; every target gets a fair say; units cancel | Not a strict likelihood |
| **`user`** | raw normalized RMSE × explicit `weights`, no auto-balancing | Full control, transparent | User must hand-manage unit scale + count imbalance |

In every mode the per-variable `weights` apply as priority multipliers, and the per-row
`sigma` column (or the `error_model` defaults — relative ~15% for LAI/biomass, ~10% for
yield, absolute days/BBCH for phenology) is honoured.

---

## 15. Phased roadmap

All four engine presets (A–D) and all four weighting modes are first-class deliverables;
the roadmap orders *when each engine lands*, not whether. Each phase ships its engines fully
selectable via `method.preset` / `objective.weighting`.

- **P0 — Skeleton & config (Python).** Loader, full schema (incl. `method.preset` +
  `objective.weighting` switches), observations reader (CSV + FileA/FileT), base-experiment
  resolver. Validate against the phenology dataset.
- **P1 — Spawn + run + score → preset C.** Table-driven writers, native parallel runner,
  PlantGro parser, all four weighting modes, `lhs → glue` end-to-end. Ships
  the v1 default.
- **P2 — Bayesian → preset A.** Port SMC-PF (generalized) + Morris/Sobol/ANOVA screening
  (SALib); auto-activate; **identifiability diagnostics** + **leave-one-environment-out
  validation**. Makes A (the headline default) available.
- **P3 — Optimizers, AgMIP & full UQ → presets B, D.** Nelder-Mead/diffevo/AgMIP stepwise
  with **BIC/AICc parameter selection** + **agmip_wls** weighting (B); MCMC/DREAM +
  **GP/RF surrogate** acceleration (D); optional **pymoo NSGA-II** multi-objective.
- **P4 — Polish.** R parity layer (**wrap CroptimizR** rather than reimplement),
  reporting/plots, HPC hook, tag-pinned deps, docs.
- **P5 — Multi-source in-season assimilation (partly landed; see §17).** Pluggable
  observation adapters + fusion and a **coupled recalibration** mode are implemented;
  EnKF/forcing remain **uncoupled prototypes** and in-season forecasting is still future.

A natural v1 milestone: reproduce the phenology repo's wheat phenology result through the
generalized engine (regression check), then extend the same run to also fit LAI + biomass
+ yield on a multi-crop dataset — proving the generalization.

---

## 16. Literature review — positioning, corrections, and additions

An online review (June 2026) of DSSAT calibration tooling and the crop-model calibration
literature. Sources listed at the end.

### 16.1 Prior art — and why a new framework is still justified

The concept must be positioned against tools that already exist, or it risks reinventing them:

- **DSSAT ships its own calibrators.** **GENCALC** (sequential, one-parameter-at-a-time
  fine-tuning) and **GLUE** (Generalized Likelihood Uncertainty Estimation — Bayesian
  Monte-Carlo: sample genetic coefficients from priors defined by the `.CUL/.ECO`
  MINIMA/MAXIMA, score with a Gaussian likelihood) are bundled with DSSAT. GLUE is pure R,
  run via `Rscript GLUE.r`, and a 2024 **GLUEP** variant adds chunked multi-core
  parallelism (87–95% faster). **Our preset C (`lhs → glue`) is essentially DSSAT-GLUE** —
  good, but the concept previously failed to acknowledge the built-in tool.
- **CroptimizR** (SticsRPacks) already calibrates DSSAT (4.7 CERES/CropSim) via a model
  wrapper, with frequentist (Nelder-Mead simplex multistart) and **Bayesian (DREAM-zs via
  BayesianTools)** algorithms, the **AgMIP protocol**, and user-defined multi-step
  workflows. Our presets B and D substantially overlap with it.

**The defensible niche** (state this explicitly in the README): DSSAT-GLUE and CroptimizR
calibrate **genetic coefficients** against mostly **end-of-season** targets, per crop. This
framework's value-add is the union of axes neither covers in one place: **(1)** jointly
fitting **LAI + biomass + grain yield + phenology**, time-series *and* scalars; **(2)**
**multi-experiment across sites, years, and crops** in one objective; **(3)** calibrating
**management, initial soil/water, soil, and weather** axes alongside genetics; **(4)** a
**pluggable engine** menu (GLUE, SMC-PF, MCMC, optimizers, AgMIP) under one config; **(5)**
the **spawn orchestration** (warm-up→parallel, resume, namespacing). **Corollary
recommendation:** don't reimplement DREAM-zs / the AgMIP protocol in R — the R parity layer
should **wrap CroptimizR** (caveat: its official DSSAT wrapper targets 4.7; we run 4.8, so
the wrapper needs a version bump or our own `dssatengine`-backed wrapper).

### 16.2 Corrections

- **Acknowledge GLUE/GENCALC/CroptimizR** (done above) — the concept now cites them as
  baselines, reuse targets, and cross-checks rather than ignoring them.
- **Calibration ≠ data assimilation.** The concept conflates two tasks. Estimating fixed
  coefficients from observations is **parameter estimation** (what SMC-PF/GLUE/MCMC do
  here). Updating model **state** in-season from streaming remote-sensing LAI is **data
  assimilation** (EnKF / particle filter), a different operation — and exactly the planned
  `DSSAT_LAI_Assimilation` repo. Add a clear scope line, and optionally an **EnKF state
  mode** for in-season LAI/biomass updating (literature shows EnKF LAI assimilation into
  DSSAT improves yield estimates, with caveats: satellite-LAI error, weak LAI→yield link,
  nonlinear-observation linearization error).
- **Time-series residuals are autocorrelated.** A naive per-observation Gaussian likelihood
  over a dense LAI/biomass series over-weights it (each day is not independent). Add an
  `obs_autocorr` option (effective-sample-size down-weighting or an AR(1) error model);
  flagged in §4.

### 16.3 Additions worth implementing

> **Status (2026-06-16):** items 1–7 below are now **implemented** — AgMIP stepwise
> selection (`engines/selection.py`), `agmip_wls` weighting (`orchestrator.py`),
> identifiability via the sensitivity ranking + posterior pair data, leave-one-
> environment-out validation, surrogate acceleration (`engines/surrogate.py`),
> SALib/scikit-learn/pymoo libraries, and the NSGA-II multi-objective front.
> Item 8 (run-budget guidance) is documented in `WALKTHROUGH.md`. The discrete-
> factor sweep and EnKF state-assimilation remain future work (see §0).

1. **AgMIP calibration protocol as a first-class engine** (`method.select`). Its core idea
   beats pure sensitivity pruning: split params into **obligatory** (degree-day / nearly
   additive, always estimated) vs **candidate** (added one at a time, **kept only if BIC
   (or AICc) drops** — a formal overfit guard); estimate with Nelder-Mead multistart; treat
   one variable group at a time (phenology before growth). Shown to **reduce both
   between-modeler variability and prediction error**. Mapped to the new `role:` field and
   `method.select`.
2. **Iterative reweighted (WLS) weighting** (`weighting: agmip_wls`). Run OLS once, set each
   variable group's weight to **1/variance(residuals)**, refit. More principled than a
   hand-set `unified` scale and is the AgMIP default. Added to §4/§14b.
3. **Identifiability & equifinality diagnostics** — the central risk in genetic-coefficient
   calibration (many parameter sets fit equally well). Report a **parameter
   correlation/collinearity matrix**, posterior **pair-plots**, and a practical-
   identifiability index from the sensitivity stage; **warn and optionally auto-deactivate**
   non-identifiable params. Mitigations the literature recommends, all config-expressible:
   **priors as regularization** (use `.CUL` MINIMA/MAXIMA as default bounds/priors — the
   GLUE convention), and **data design** (calibrating across *contrasting* environments
   sharply reduces equifinality vs. similar seasons — surface this as guidance).
4. **Held-out validation** (`method.validation`). Report **calibration vs. evaluation**
   error separately via **leave-one-environment-out** (or year/random) CV — reusing the
   phenology repo's CV-blocking machinery. AgMIP stresses this; without it, multi-axis
   calibration overfits silently.
5. **Surrogate / emulator acceleration** promoted to first-class (`method.surrogate: gp|rf`).
   Train a GP/RF emulator on the LHS design, calibrate on the *emulator* (cheap), validate
   top candidates on the real model; optionally **active-learning / Bayesian-optimization**
   to place new runs where they most reduce posterior uncertainty. Strongly supported for
   expensive simulators and makes preset D (MCMC/DREAM) tractable.
6. **Concrete libraries** (Python core): **SALib** (Morris / Sobol / FAST) for §7
   sensitivity; **SciPy** / **pymoo** (NSGA-II) for optimizers; **emcee** or
   **BayesianTools-equivalent** for MCMC/DREAM; **scikit-learn**/**GPy/BoTorch** for
   surrogates. R parity via **CroptimizR** + **BayesianTools** + **sensitivity**.
7. **Multi-objective (Pareto) option** (`pymoo` NSGA-II) as an alternative to scalarized
   weighting — exposes the **trade-off front** between fitting yield vs LAI vs phenology
   instead of collapsing them, useful when targets conflict.
8. **Run-budget guidance & GLUEP-style chunking.** Calibration cost is dominated by DSSAT
   runs; document rules of thumb (Morris ≈ r·(k+1) runs; GLUE typically 10³–10⁴; MCMC ≫),
   and keep the chunked-parallel + resume design (we already have it) front-and-centre.

### 16.4 Sources

- DSSAT GLUE (tool): https://github.com/DSSAT/glue · DSSAT tools: https://dssat.net/tools/
- GLUE + parallel computing (GLUEP), 2024: https://www.sciencedirect.com/science/article/abs/pii/S0168169924009049
- GenCalc vs GLUE (spring wheat): https://scialert.net/fulltext/?doi=ja.2016.130.135
- AgMIP calibration protocol (soil-crop models): https://www.sciencedirect.com/science/article/pii/S1364815224002081 · phenology protocol: https://link.springer.com/article/10.1007/s13593-023-00900-0 · AgMIP page: https://agmip.org/crop-model-calibration-3/
- CroptimizR: https://sticsrpacks.github.io/CroptimizR/ · DREAM-zs: https://sticsrpacks.github.io/CroptimizR/articles/Parameter_estimation_DREAM.html · AgMIP phenology vignette: https://sticsrpacks.github.io/CroptimizR/articles/AgMIP_Calibration_Phenology_protocol.html
- Practical identifiability for calibration data design: https://www.sciencedirect.com/science/article/abs/pii/S0168169921004749 · data requirements: https://www.sciencedirect.com/science/article/abs/pii/S0168192316307420
- SALib: https://salib.readthedocs.io/ · GSA of a cropping-system model: https://www.sciencedirect.com/science/article/pii/S1364815223003183
- GP emulators for expensive simulators: https://link.springer.com/chapter/10.1007/978-3-031-66085-6_15
- DSSAT + remote sensing via MCMC (rice, 2025): https://www.mdpi.com/2223-7747/14/8/1206 · EnKF LAI assimilation: https://www.sciencedirect.com/science/article/abs/pii/S030438001300416X

---

## 17. Multi-Source In-Season Observation Assimilation Framework

This extends the calibrator from "fit fixed coefficients to a finished season" toward
"ingest whatever observations are streaming in during the season, from whatever sources,
and update the model." It has three layers; **be precise about which are coupled to DSSAT
and which are not**, because §16.2 flagged exactly this calibration-vs-assimilation trap.

### 17.1 Source adapters (`sources/`) — implemented

A pluggable `ObservationSource` base (`sources/base.py`) with a registry
(`sources/registry.py`) and concrete adapters, each emitting the **extended schema**
`exp_id | treatment | variable | kind | date | value | sigma | weight | source |
quality_flag | spatial_res_m` and carrying a **source-specific error model**:

| adapter | source_type | emits (DSSAT var) | error model |
|---|---|---|---|
| `sentinel2_lai`, `modis_lai` | satellite | `LAID` | base RMSE, LAI-saturation & cloud/QC inflation |
| `uav_multispectral` | uav | `LAID`, `canopy_cover`, `canopy_height` | per-variable, flight-quality inflation |
| `field_measurements` | field | FileA/FileT vars | relative/absolute per variable |
| `farm_phenology`, `farm_management` | farm_software | `GSTD`/`ADAT`/`MDAT`, events | date-precision dependent |
| `soil_moisture_iot`, `canopy_temperature` | iot | `SW`, `TMEAN` | sensor-type & calibration dependent |

Activate them in the `observation_sources:` config block (`adapter:` selects the class).

> **Scoring caveat (enforced by a warning).** A variable is only scored if it maps to an
> `engine.timeseries_outputs`/`scalar_outputs` column *and* DSSAT outputs it. `SW`,
> `TMEAN`, `canopy_cover`, `canopy_height` are **not** produced by the PlantGro/Evaluate
> parse today, so those rows are ingested but ignored unless you add the output mapping
> (and a parser for it). `orchestrator` logs the unmatched variables so this is never silent.

### 17.2 Fusion (`fusion.py`) — implemented

`ObservationFuser` gathers every active source for an experiment, applies each source's
quality filter, and reconciles coincident `(exp_id, treatment, variable, date)` records by
one of three strategies (`fusion.conflict_resolution`): **`keep_all`**, **`inverse_variance`**
(Bayesian-optimal 1/σ² combine — the recommended default), or **`priority`** (highest-ranked
source wins, per `fusion.source_priority`). The fused long table flows into the *exact same*
objective/engines as ordinary calibration via `Observations.from_sources`.

### 17.3 In-season assimilation modes (`assimilation.mode`)

| mode | engine | coupled to DSSAT? | what it does |
|---|---|---|---|
| **`recalibration`** *(default)* | `engines/recalibration.py` | **Yes** | At each observation checkpoint, filter obs to that date and re-run the calibration pipeline; `warm_start` seeds each checkpoint with the previous best fit. This is parameter estimation through time — the genuinely working in-season path. |
| `enkf` | `engines/enkf.py` | **No (prototype)** | Stochastic EnKF *update* math, but no DSSAT forecast step and no state re-injection — DSSAT-CSM has no clean state-restart hook. As wired, the ensemble is synthetic; output is illustrative. |
| `forcing` | `engines/forcing.py` | **No (prototype)** | Direct state replacement on a dict that is never fed back to DSSAT. Illustrative. |

`enkf`/`forcing` **refuse to run** unless `assimilation.allow_uncoupled: true` — so they
cannot silently return meaningless numbers. `combined_mode()` runs a base calibration and
then seeds the assimilation step with its result.

**Driving it:** `python run_calibration.py config.yaml --assimilate`
(or `--combined`, or `--assim-mode recalibration|enkf|forcing`). Recalibration writes
`assimilation_trace.csv` (best θ per checkpoint) plus `assimilation.json`.

### 17.4 Honest status & what's still missing

- **Coupled EnKF/forcing** would require a DSSAT forecast→update→re-inject loop (e.g. via
  `.WTH`/initial-condition forcing or a restart build). Until then they stay gated.
- **In-season forecasting** (project yield/maturity forward with uncertainty from the
  current best θ) is not implemented; the old `assimilation.forecast` config block was
  removed to avoid a dead knob.
- The fusion + recalibration spine is tested offline (`tests/test_assimilation.py`); the
  full recalibration loop needs a DSSAT install (the slow E2E tier).

---

## 18. In-season forecasting & new-crop adaptation (implemented layer)

Two further use cases drove the modules below. Everything here is **optional and off by
default** (a bare config behaves exactly as before), and the pure logic is unit-tested in
`tests/test_new_features.py`.

### 18.1 In-season LAI nowcast (`forecast.py`, `orchestrator.nowcast`)

Once calibrated, DSSAT already simulates the whole season, so a forward LAI estimate is the
calibrated run read past the last observation. `forecast.py` turns that into a product:

* **`ensemble_percentiles`** — propagate the behavioural parameter sets (`n_ensemble`) and
  take daily P10/P50/P90, so the forecast carries parameter uncertainty that widens with lead
  time (`lead_time_table` summarises this).
* **`anchor_correction`** — shift the forward curve to start from the last observation and
  decay the correction to zero over `decay_days`. A seam-free nowcast without injecting state
  into DSSAT (which CSM does not support) — the pragmatic alternative to coupled EnKF.
* **`orchestrator.nowcast(cfg, as_of_date)`** — the operational loop: filter obs ≤ date,
  (re)calibrate, **persist** `nowcast_state.json`, **warm-start** the next call from it, and
  forecast. CLI: `--nowcast YYYY-MM-DD [--forecast]`.

### 18.2 Weather driver (`weather.py`)

Current implementation: `provider: file` is a no-op, `provider: nasa_power`
is the lightweight single-point fallback, and `provider: dssatutils` delegates
whole-year acquisition to `dssatutils.process_weather_*` before applying the
calibrator's in-season `fill_gap` / `horizon` layer. The matching new-site soil
path lives in `acquisition.py` behind `soil.provider: dssatutils`.

A pluggable provider layer (default `file` = use DSSAT's own `.WTH`, a no-op). `NasaPowerProvider`
acquires daily weather from the keyless NASA POWER API and `write_wth` renders a `.WTH`;
`fill_gap` extends the record to (and a `horizon` past) today by `persistence`/`climatology`.
NASA POWER's ~1–2 week latency means even reaching "today" needs gap-fill — documented as a
stand-in for a true forecast-weather provider (the remaining future work for genuine
look-ahead forecasting).

### 18.3 Multi-source refinements

* **Cloud masking** — `sources/satellite.py` now *drops* scenes above `max_cloud_fraction`
  (and `quality_filter` drops bad-QC rows) rather than only inflating σ.
* **Observation operator** — an optional linear `obs_operator: {scale, offset}` corrects the
  satellite "effective LAI" vs model-LAI bias so it does not leak into the coefficients.
* **Planting date as input** — `observations.planting_dates()` reads farm-management sowing
  dates; `spawn` writes them as the FileX `PDATE` (opt-in `management_options.use_source_planting_date`)
  rather than calibrating an offset.

### 18.4 New crop / cultivar / species (`scaffold.py`, `diagnostics.py`, staging, CV)

* **New cultivar of an existing species** — the core competency; declare the `.CUL`
  coefficients and calibrate against dates + yield + biomass across site-years.
* **New species** — only via the **analog-template path**. `scaffold_crop.py` clones the most
  similar module's `.CUL/.ECO/.SPE` under a new code and emits a starter `parameters:` block
  (bounds from MINIMA/MAXIMA rows, normal priors, phenology/growth role split). The `.SPE`
  writer (`writers.edit_species`) enables species-coefficient adaptation but stays **gated**
  (`gating.species: free` + `group: genetic_species`). The framework supplies the calibration
  engine and scaffolding — it does **not** choose the analog, build a new module, or detect
  structural inadequacy for you.
* **Sparse-data discipline** — `method.staging` freezes groups/params the data can't
  constrain (e.g. seed/yield params mid-season); `diagnostics.identifiability` reports
  posterior-vs-prior width and collinearity (what is actually pinned); `diagnostics.structural_adequacy`
  flags variables the best fit still misses badly (wrong analog/observations, not a calibration
  problem); and `validate_cv` adds `year`/`site`/`random` folds for an honest transfer test on
  few site-years.

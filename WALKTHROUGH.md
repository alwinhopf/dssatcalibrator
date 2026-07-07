# dssatcalibrator — Walkthrough (start here)

A friendly, step-by-step guide. **No statistics background assumed.** By the end you
will understand what the tool does, how to run it, how to read the results, and how
to adapt it to your own crop or field site.

If you want the deep design rationale instead, read [`CONCEPT.md`](CONCEPT.md). For a
one-screen feature list, read [`README.md`](README.md).

---

## 1. What problem does this solve?

DSSAT is a crop simulation model. To make it match a *real* field, it needs
**cultivar coefficients** — numbers describing how fast your variety develops, how
big its leaves get, how much it yields, and so on. Out of the box those numbers are
generic. **Calibration** is the process of nudging them until the model's output
lines up with what you actually measured in the field.

Think of it like tuning a guitar:

- the **strings** are the parameters (e.g. `CSDL`, `PPSEN`, `LFMAX`);
- the **target pitch** is your field data (measured LAI, biomass, yield, flowering date);
- **calibration** turns the tuning pegs until the model "plays in tune" with reality.

Doing this by hand is slow and subjective. `dssatcalibrator` does it automatically:
it runs DSSAT hundreds or thousands of times with different parameter values, scores
each run against your data, and reports the values that fit best — *plus* how
confident you can be in them.

---

## 2. The mental model (five words)

> **config → sample → run → score → fit**

| Step | What happens | Where it lives |
|---|---|---|
| **config** | You declare which parameters to vary and what data to fit | `config_hemp.yaml` |
| **sample** | The tool picks many candidate parameter sets | `samplers.py` |
| **run** | Each candidate is fed to DSSAT (in parallel) | `spawn.py`, `runner.py` |
| **score** | Each run is compared to your observations | `objective.py` |
| **fit** | An *engine* turns all those scores into an answer | `engines/` |

You mostly touch the **config**. Everything else is automatic.

---

## 3. Quick start

```bash
# 1. install (one time)
pip install -e .          # core
pip install -e .[full]    # optional: Sobol sensitivity + surrogate

# 2. run the bundled hemp example (uses preset C = the simplest engine)
python run_calibration.py config_hemp.yaml --preset C --n 200

# 3. look in results/<name>/ for tables and figures/<name>/ for plots
```

That's a complete calibration. The sections below explain how to steer it.

---

## 4. The config, the only file you usually edit

Open `config_hemp.yaml`. It has a few blocks. Here are the ones you'll actually change.

### 4a. Which parameters to calibrate

```yaml
parameters:
  genetic_cultivar:
    CSDL:  { active: true,  role: obligatory, min: 11.5, max: 18.0, start: 12.8 }
    LFMAX: { active: true,  role: candidate,  min: 0.3,  max: 2.5,  start: 1.4 }
    PHINT: { active: false, min: 50, max: 150, start: 95.0 }
```

For each parameter:

- **`active`** — `true` means "calibrate this"; `false` means "leave it at `start`".
  Declare everything you *might* vary; switch on only the few you want this run.
- **`min` / `max`** — the allowed range to search. Keep it physically plausible.
- **`start`** — the default/best-guess value (also the value used when `active: false`).
- **`prior`** *(optional)* — your belief about the value before seeing the data
  (see §7). Omit it and the tool assumes "anything in `[min, max]` is equally likely".
- **`role`** *(optional)* — `obligatory` or `candidate`, used only by the stepwise
  selection stage (§6d).
- **`scope`** *(optional)* — omit it or use `global` to fit one shared value across
  all experiments. Use `scope: experiment` (or `pooling: per_experiment`) when each
  experiment should get its own value inside one pooled calibration.

For example, this fits one shared species parameter while letting every experiment
have its own cultivar and ecotype offsets:

```yaml
parameters:
  genetic_species:
    LFMAX: { active: true, min: 0.3, max: 2.5, start: 1.4 }      # shared by default
  genetic_cultivar:
    CSDL:  { active: true, min: 11.5, max: 18.0, start: 12.8, scope: experiment }
  genetic_ecotype:
    ECOA:  { active: true, min: 0.8, max: 1.2, start: 1.0, pooling: per_experiment }
```

> **Tip for non-experts:** start with 3–6 active parameters. Calibrating 15 at once is
> slow and tends to *over-fit* (see §9). Use the **sensitivity** stage (§6c) to find
> which ones actually matter.

### 4b. What data to fit

```yaml
engine:
  timeseries_outputs:        # things measured repeatedly through the season
    biomass: "CWAD"          #   your name : the DSSAT output column
    LAI:     "LAID"
  scalar_outputs:            # things measured once (end of season)
    grain_yield: "HWAM"
    anthesis:    "ADAP"
```

The tool reads your measurements from the DSSAT `FileA`/`FileT` files (or a CSV) and
matches them to these outputs automatically.

### 4c. Before calibrating: build a parameter impact atlas

When you need to understand what each DSSAT input actually changes, run one-at-a-time
real-DSSAT sweeps before a full calibration:

```bash
python run_impact_atlas.py config_hemp.yaml \
  --experiments UFCI2101 \
  --groups genetic_cultivar genetic_ecotype genetic_species management initial_conditions soil weather \
  --discover-genotype --allow-species --max-per-group 1 --grid-points 3 --cores 2 --no-long
```

From R, call the same real-runner:

```r
library(dssatcalibrator)
atlas <- run_impact_atlas(
  "config_hemp.yaml",
  experiments = "UFCI2101",
  groups = c("genetic_cultivar", "genetic_ecotype", "genetic_species", "soil", "weather"),
  discover_genotype = TRUE,
  allow_species = TRUE,
  max_per_group = 1,
  grid_points = 3,
  num_cores = 2,
  write_long = FALSE
)
```

The atlas writes `impact_summary.md`, `run_manifest.csv`, `run_manifest.json`,
`failed_runs.csv`, `file_manifest.csv`, `parameter_catalog.csv`,
`score_effects.csv`, `parameter_impact_summary.csv`, `output_impact_summary.csv`,
`parameter_output_effects.csv`, `capability_map.md`, and per-group PNGs under
`plots/`. If you keep the full long table, it also writes `outputs_long.csv` or
`outputs_long.csv.gz`; for broad or longer sweeps, prefer `--no-long` /
`write_long = FALSE` until you know you need the row-level table. Use the
manifests to spot failed input edits, the summary tables to rank which
parameters move the objective or DSSAT outputs, and the capability map to decide
what should be upstreamed into `dssatengine` or `dssatutils`.

FileX parameters can target any section/field, including row-scoped management
and starting-condition edits. For example, a management parameter can set
`section: IRRIGATION`, `field: IRVAL`, `op: mult`, and `row: 2`, while an
initial-condition parameter can set `field: SNH4` or multiply `SH2O` via
`initial_soil_water_mult`. Text/code fields are supported with `type: code`;
use `required: true` when a missing section or column should fail loudly.

### 4d. Which engine / pipeline (the one knob that matters most)

```yaml
method:
  preset: "C"     # A | B | C | D  — see the decision table in §5
```

That's it. A preset wires up the whole pipeline for you.

---

## 5. Which preset should I pick?

| You want… | Use | What it does | Cost |
|---|---|---|---|
| A quick first answer | **C** | LHS sampling → **GLUE** | Low (fully parallel) |
| Just "good numbers", no uncertainty | **B** | screen → **optimiser** (best-fit point) | Low–Medium |
| Best-fit **with** uncertainty (recommended) | **A** | screen → map → **particle filter** | High |
| Publication-grade uncertainty | **D** | screen → **MCMC** posterior | Very high |

Rule of thumb: **start with C** to sanity-check your setup, then move to **A** for a
real result. Use **D** only when you need rigorous credible intervals (and consider
turning on the surrogate, §6f, to make it affordable).

A preset chooses the **main estimator** automatically. The optional helper stages
(sensitivity, selection, surrogate) stay **off** until you switch them on with
`active: true` — so nothing surprising happens behind your back.

---

## 6. The engines, in plain language

### 6a. GLUE (preset C) — "try many, keep the good ones"
Runs a big batch of random-but-spread-out parameter sets, scores them all, and keeps
the best fraction as a cloud of plausible answers. Simple, fast, and every run is
independent so it uses all your cores. Great first look.

### 6b. SMC particle filter (preset A) — "a swarm that learns from the data"
Keeps a *swarm* of parameter guesses and feeds it your observations one date at a
time. Guesses that match get more weight; when the swarm gets too concentrated it
*resamples* (copies the good ones) and *jiggles* them to keep exploring. The end
result is a posterior — best estimate **plus** how uncertain it is.

### 6c. Sensitivity screening (Morris / Sobol) — "which knobs even matter?"
Before calibrating, ask which parameters actually change the output. Turn it on:

```yaml
method:
  sensitivity: { engine: morris, active: true, trajectories: 12, auto_activate: true }
```

With `auto_activate: true` the tool keeps only the influential parameters and freezes
the rest — fewer knobs, faster runs, less over-fitting. Morris is cheap and needs no
extra packages; Sobol is more thorough but needs `SALib` (`pip install -e .[full]`).

### 6d. AgMIP stepwise selection — "don't add a knob unless it earns its place"
Calibrating too many parameters fits the calibration data but predicts new seasons
badly. This stage estimates the `obligatory` parameters first, then adds `candidate`
ones **one at a time, keeping each only if it improves an information criterion**
(BIC or AICc) that penalises complexity. Turn it on:

```yaml
method:
  select: { engine: stepwise_bic, active: true }
```

### 6e. Optimisers (preset B) — "just give me the best single set"
`nelder_mead` (fast, local) or `diffevo` (global, robust). No uncertainty, just one
best-fit parameter vector. Differential evolution scores its whole population in
parallel each generation.

### 6f. MCMC (preset D) — "the full posterior"
A Markov chain that samples parameter sets in proportion to how well they fit
(times your prior). Gives clean credible intervals and shows which parameters trade
off against each other. The most expensive engine — see §6g.

### 6g. Surrogate acceleration — "let a fast stand-in do the heavy lifting"
DSSAT is slow. The surrogate runs a modest design on the *real* model, fits a fast
statistical emulator (Gaussian Process or Random Forest), searches the emulator for
free, and validates only the best handful on the real model. This is what makes MCMC
practical for slow crops. Needs `scikit-learn` (`pip install -e .[full]`).

```yaml
method:
  surrogate: { engine: gp, active: true, n_train: 64, top_k: 10 }
```

### 6h. NSGA-II multi-objective — "show me the trade-offs, don't average them"
Instead of collapsing yield + LAI + phenology into one score, it finds the
*trade-off front*: the parameter sets where you can't improve one target without
hurting another. Useful when targets genuinely conflict.

```yaml
method:
  multiobjective: { engine: nsga2, variables: [grain_yield, LAI], pop_size: 16, n_gen: 5 }
```

---

## 7. Priors (optional, but useful)

A *prior* says what you believe about a parameter **before** looking at the data.

```yaml
P5: { active: true, min: 300, max: 700, start: 505, prior: { dist: normal, sd: 50 } }
```

- `uniform` (default) — "anything in the range is equally likely". Safe and unbiased.
- `normal` — "probably near `start`, within about `sd`". Use when you trust a value.
- `lognormal` / `triangular` — for positive-only or peaked beliefs.

Priors gently pull the answer toward sensible values and stop the search wasting time
in implausible corners. The Bayesian engines (GLUE, SMC-PF, MCMC) use them; if you
leave them out, you simply get the flat default.

---

## 8. How the fit is scored (the objective)

Every run is compared to your data and reduced to one number (lower = better). You
choose **how** the variables are combined:

| `objective.weighting` | Plain meaning |
|---|---|
| `unified` *(default)* | Balanced: each variable's average error, summed. |
| `count_scale` | Like unified but averaged, so a 100-point LAI series can't drown a 1-point yield. |
| `sigma` | Strict statistical misfit; what the Bayesian engines use internally. |
| `user` | You set the weights by hand. |
| `agmip_wls` | Auto-weights each variable by 1/its-noise (the AgMIP default). |

Extra options:

- `weights: { LAI: 1.0, grain_yield: 2.0 }` — make some variables count more.
- `obs_autocorr: true` — down-weights dense daily time-series so they don't unfairly
  dominate one-off measurements like final yield.

The reported metrics per variable are the standard set: **RMSE, nRMSE%, bias (MBE),
Willmott's d, modelling efficiency (EF), R²**.

---

## 9. The big trap: over-fitting & equifinality

The single most important thing to understand:

> With enough free parameters you can fit your calibration data *perfectly* and still
> predict next year badly. And often **many different parameter sets fit equally
> well** (this is called *equifinality*).

How this tool helps you avoid it:

1. **Screen first** (§6c) and calibrate only influential parameters.
2. **Use stepwise selection** (§6d) so extra parameters must earn their place.
3. **Validate on held-out data**: `python run_calibration.py config.yaml --validate`
   runs leave-one-environment-out CV and reports *calibration* vs *evaluation* error
   separately. If evaluation error is much worse, you're over-fitting.
4. **Calibrate across contrasting environments** (different sites/years), not several
   near-identical seasons — this sharply reduces equifinality.
5. **Look at the posterior pair-plots / correlations** (in the figures) — parameters
   that are strongly correlated are not separately identifiable; consider fixing one.

---

## 10. Reading the outputs

After a run, `results/<name>/` holds the **tables** and `figures/<name>/` holds the
**plots**.

**Tables**

| File | What it tells you |
|---|---|
| `best_theta.json` | the best-fit parameter values (drop straight into DSSAT) |
| `summary_fit.csv` | fit quality per variable (RMSE, nRMSE%, d, EF, R²) |
| `objective_breakdown.csv` | objective components per experiment, variable, and observation kind |
| `posterior_summary.csv` | best + mean ± spread + 5–95% range per parameter |
| `design.csv` | every parameter set tried, with its score |
| `manifest.csv` / `manifest.json` | every DSSAT spawn with sample, experiment, theta, run directory, and status |
| `fit_by_experiment.csv` | fit broken down per experiment |
| `sensitivity_ranking.csv` | parameter influence (if screening ran) |
| `validation_loeo.csv` | calibration vs evaluation error (if `--validate`) |

**Plots**

| Figure | What to look for |
|---|---|
| `fig_obs_vs_sim*` | points on the 1:1 line = good fit |
| `fig_timeseries` | model curves should track the measured points |
| `fig_experiment_<EXP>_T<TRT>_3x3` | one diagnostic panel per calibrated treatment: LAI, stages, canopy size, biomass, stress, soil, weather, and tissue N when DSSAT outputs/observations are available |
| `fig_param_posteriors` | narrow posterior (vs flat prior) = well-identified parameter |
| `fig_fit_bars` | low nRMSE%, high Willmott d per variable |
| `fig_sensitivity` | which parameters matter most |
| `fig_ess_trajectory` | SMC health: ESS shouldn't crash to near-zero |
| `fig_mcmc_trace` | MCMC health: walker lines should overlap once "mixed" |
| `fig_pareto` | the trade-off front (NSGA-II) |

---

## 11. How long will it take? (run-budget guidance)

Cost is dominated by DSSAT runs. Rough rules of thumb (× number of experiments):

| Stage | Approx. model runs |
|---|---|
| Morris screening | `trajectories × (n_params + 1)` (e.g. 12 × 14 ≈ 170) |
| GLUE / LHS | the `n` you set (typically 200–10,000) |
| SMC-PF | `n_particles × (1 + #resample steps)` |
| MCMC | `n_walkers × n_steps` (the most expensive) |
| Optimiser (diffevo) | `popsize × maxiter` |

Everything runs in parallel across `calibrator.num_cores` (0 = all cores but two). To
go faster: screen to fewer parameters, lower the sample/particle count for a draft,
or turn on the surrogate for the expensive engines.

---

## 12. Adapting to your own crop, site, or scenario

- **New site, same crop:** point the config's experiments/observations at your
  site-years. No code change.
- **New crop:** add a `crops:` block (DSSAT code, model file stem, ecotype, cultivar
  anchor, FileX extension) and make sure the cultivar coefficients you want appear in
  the `parameters` block. The only crop-specific code is the column layout in
  `writers.py`.
- **New scenario:** flip `active` flags or edit ranges — e.g. calibrate management or
  soil instead of genetics. Pure config.

---

## 13. Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| "No active parameters in config" | No parameter has `active: true`. |
| All scores are `inf` / no fit | DSSAT produced no output — check `dssat_exe`/`dssat_dir` paths and that the experiment runs standalone. |
| "needs SALib" / "needs scikit-learn" | Install the optional engines: `pip install -e .[full]`. |
| Calibration great, validation poor (`--validate`) | Over-fitting — screen + select fewer parameters (§9). |
| Posterior as wide as the prior | That parameter isn't identifiable from your data; fix it or add more contrasting data. |
| MCMC walkers don't overlap in `fig_mcmc_trace` | Not converged — increase `n_steps`/`burn_in`, or reduce parameters. |
| Runs feel slow | Lower `n`/`n_particles` for a draft, screen first, or enable the surrogate. |

---

## 14. In-season calibration & multiple data sources

Everything above fits a *finished* season. You can also calibrate **during** the season,
and pull observations from **more than just field measurements** (satellite, drone, farm
software, sensors).

### 14a. Many sources, one fused dataset

Instead of `source.observations: dssat`, declare an `observation_sources:` block — each
entry is an adapter that reads one kind of data:

```yaml
observation_sources:
  field_measurements: { active: true,  hemp_dir: "C:/Users/alwin/Documents/GitHub/DSSAT48Hemp/Hemp", crop_ext: "HM" }
  sentinel2_lai:      { active: true,  data_path: "obs/sentinel2_lai.csv" }
  uav_multispectral:  { active: false, data_path: "obs/uav.csv" }

fusion:
  conflict_resolution: "inverse_variance"   # combine overlapping obs by 1/σ² (recommended)
  source_priority: ["field_measurements", "uav_multispectral", "sentinel2_lai"]
```

The framework fetches each active source, attaches a **source-specific error bar** (a cloudy
Sentinel LAI pixel is trusted less than a hand measurement), and **fuses** overlapping
observations. Then calibration proceeds exactly as before. Adapters available:
`sentinel2_lai`, `modis_lai`, `uav_multispectral`, `field_measurements`, `farm_phenology`,
`farm_management`, `soil_moisture_iot`, `canopy_temperature`.

> **Heads-up:** a source can only help if its variable is one you score (in
> `engine.timeseries_outputs`/`scalar_outputs`) *and* DSSAT outputs it. Soil-moisture (`SW`)
> and canopy-temperature (`TMEAN`) are ingested but **not scored** out of the box — the tool
> prints a warning listing any such ignored variables.

### 14b. In-season calibration (recalibration)

"What are this field's coefficients **given the data so far**?" — re-estimated at each new
observation date:

```yaml
assimilation:
  active: true
  mode: "recalibration"      # the path that actually drives DSSAT
  recalibration: { recal_sample_size: 100, warm_start: true, update_frequency: "on_observation" }
```

```bash
python run_calibration.py config_hemp.yaml --assimilate     # writes assimilation_trace.csv
python run_calibration.py config_hemp.yaml --combined        # full calibration, THEN in-season
```

You get `assimilation_trace.csv` — the best-fit parameters at each checkpoint — so you can
watch the estimate sharpen as the season unfolds. `warm_start` makes each checkpoint refine
the previous one instead of starting cold.

> **About `enkf` / `forcing` modes:** these are **prototypes that are *not* coupled to
> DSSAT** — they do the data-assimilation math but never feed the updated state back into a
> running simulation, so their numbers are illustrative only. They refuse to run unless you
> set `assimilation.allow_uncoupled: true`. For real in-season work use `recalibration`.

---

## 15. In-season nowcasting: weather, forecast & continuity

Calibrating "to date" is half the job; the other half is **projecting the crop forward**
(e.g. estimating LAI across a stretch of cloudy days when no satellite image is available).
All of the following are **optional and off by default**.

### 15a. The operational nowcast

One command does the whole loop — (re)calibrate on everything observed up to a date, save
the calibration, and forecast forward:

```bash
python run_calibration.py config_hemp.yaml --nowcast 2024-07-15 --forecast
```

It writes `nowcast_state.json` (the latest best parameters) and `forecast_LAID.csv` (the
forward LAI). The saved state **warm-starts the next call**, so between satellite passes you
reuse the calibration, and when a fresh cloud-free image lands you just re-run with a later
date.

### 15b. The forecast (with uncertainty and continuity)

```yaml
forecast:
  active: true
  variables: ["LAID"]
  n_ensemble: 30          # propagate the best 30 parameter sets -> P10/P50/P90 band
  anchor_continuity: true # start the forecast FROM the last observation...
  decay_days: 21          # ...and relax back to the pure model over 3 weeks
```

- **Ensemble band** — running the behavioural parameter sets forward gives a P10/P50/P90
  fan, so you see *how uncertain* the projection is (it widens with lead time).
- **Anchor continuity** — the model's LAI on the last-observation day rarely equals the
  observed value; `anchor_continuity` shifts the forward curve to start exactly at the
  observation and fades that correction out over `decay_days`. This gives a seam-free
  nowcast without needing to inject state into DSSAT.

### 15c. Acquiring weather (and filling the gap to "today")

By default the tool uses the `.WTH` files DSSAT already has (`provider: file`). To pull
weather for a bare site, or to extend the record toward today for a forecast:

```yaml
weather:
  provider: nasa_power      # keyless NASA POWER daily API (needs internet)
  gap_fill: climatology     # none | persistence | climatology
  horizon: 10               # extend this many days past the last record
```

> **Reality check:** NASA POWER is reanalysis with a ~1–2 week lag, so even reaching *today*
> needs `gap_fill`. `climatology`/`persistence` are honest stand-ins, **not** a true weather
> forecast — the forecast is only as good as the weather driving it.

---

## 16. Calibrating a new crop, cultivar, or species

### 16a. A new cultivar of an existing crop — the easy case

Point the config at your site-years (planting/flowering/maturity dates, yield, biomass),
list the cultivar coefficients in `parameters`, and calibrate. With sparse data, **calibrate
few coefficients**: screen first, freeze the rest.

### 16b. A new species — scaffold from an analog

DSSAT can't model a species it doesn't have; you adapt the **most similar existing module**
(the carinata-from-canola pattern). `scaffold_crop.py` clones the analog genotype files and
writes a starter parameter block:

```bash
python scaffold_crop.py --dssat-dir C:/DSSAT48 \
    --analog-stem SBGRO048 --new-stem QUGRO048 --new-code QU \
    --source-anchor IB0001 --out-dir templates/quinoa
```

You get `templates/quinoa/QUGRO048.CUL/.ECO/.SPE` and a `parameters_block.yaml` (bounds from
the file's MINIMA/MAXIMA rows, normal priors, a phenology→`obligatory` / growth→`candidate`
split). **Review every bound against literature** — the scaffold copies physiology, it does
not invent it. Editing the species file (`.SPE`) stays **gated**: set `gating: { species:
free }` and tag the parameters with `group: genetic_species` only when you truly mean to
adapt species physiology.

### 16c. Sparse-data discipline (the three guards)

With only a few site-years you can fit the calibration data and still predict badly. Use:

```yaml
method:
  sensitivity: { engine: morris, active: true, auto_activate: true }   # keep only what matters
  select:      { engine: stepwise_bic, active: true }                  # add a param only if it earns it
  staging:     { freeze_groups: ["genetic_ecotype"], freeze_params: ["XFRT", "WTPSD"] }  # freeze yield params
  validation:  { scheme: site }                                        # honest transfer test
diagnostics: { active: true }                                          # report what's identifiable
```

- **`staging`** freezes coefficients your data can't constrain (e.g. seed/yield params when
  you only have mid-season LAI + phenology).
- **`diagnostics`** writes `identifiability.csv` (posterior-vs-prior width; which coefficients
  are actually pinned) and `structural_adequacy.csv` (a warning when *no* parameter set fits —
  a sign the analog module is wrong, which calibration can't fix).
- **Validation schemes** (`--validate --cv-scheme site|year|loeo|random`) report calibration
  vs evaluation error; with few site-years prefer `site` or `year` folds.

---

## 17. One-line recipes

```bash
# sanity check (fast)
python run_calibration.py config_hemp.yaml --preset C --n 100

# recommended real run (uncertainty)
python run_calibration.py config_hemp.yaml --preset A --n-particles 250

# screen first, then calibrate only the parameters that matter
python run_calibration.py config_hemp.yaml --preset A --sensitivity morris

# just the best numbers, fast
python run_calibration.py config_hemp.yaml --optimizer diffevo

# check you're not over-fitting
python run_calibration.py config_hemp.yaml --validate

# in-season: recalibrate as observations arrive (writes assimilation_trace.csv)
python run_calibration.py config_hemp.yaml --assimilate

# in-season nowcast: calibrate to a date, then forecast LAI forward
python run_calibration.py config_hemp.yaml --nowcast 2024-07-15 --forecast

# honest transfer test with site folds + identifiability/structural diagnostics
python run_calibration.py config_hemp.yaml --validate --cv-scheme site
python run_calibration.py config_hemp.yaml --diagnostics

# scaffold a new crop from an analog DSSAT module
python scaffold_crop.py --dssat-dir C:/DSSAT48 --analog-stem SBGRO048 \
    --new-stem QUGRO048 --new-code QU --source-anchor IB0001 --out-dir templates/quinoa
```

Happy calibrating. When in doubt: **screen, calibrate few parameters, and validate.**

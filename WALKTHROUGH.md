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

### 4c. Which engine / pipeline (the one knob that matters most)

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
| `posterior_summary.csv` | best + mean ± spread + 5–95% range per parameter |
| `design.csv` | every parameter set tried, with its score |
| `fit_by_experiment.csv` | fit broken down per experiment |
| `sensitivity_ranking.csv` | parameter influence (if screening ran) |
| `validation_loeo.csv` | calibration vs evaluation error (if `--validate`) |

**Plots**

| Figure | What to look for |
|---|---|
| `fig_obs_vs_sim*` | points on the 1:1 line = good fit |
| `fig_timeseries` | model curves should track the measured points |
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

## 14. One-line recipes

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
```

Happy calibrating. When in doubt: **screen, calibrate few parameters, and validate.**

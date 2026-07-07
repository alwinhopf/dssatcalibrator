# Europe hemp calibration audit

This note captures the July 2026 staged hemp calibration trial on:

- Italy/Piacenza: `ITPI2001`, `ITPI2101`
- United Kingdom/Aberystwyth: `UKAB2101`

## Approach reviewed

The current China hemp workflow is staged:

1. Stage 1 fits phenology only (`emergence`, `anthesis`) using cultivar/ecotype
   coefficients.
2. Stage 2 fits biomass/growth after freezing the Stage 1 phenology values.
3. Management, initial conditions, soil, and weather calibration remain inactive;
   the spawned DSSAT runs use the source FileX management as encoded.
4. Cultivar/ecotype coefficients are scoped by cultivar when multiple anchors
   must be calibrated in the same run.

After correcting the UK FileX cultivar row, all three Europe experiments now use
`HM IB0007 Futura75`. The Europe configs therefore use one shared cultivar anchor
(`IB0007`) and one shared ecotype anchor (`HM0005`) rather than the earlier
Italy/UK split.

## Configs added

- `stage1_phenology.yaml`: pooled Europe phenology, one shared `IB0007`
  cultivar/ecotype parameter vector across `ITPI2001`, `ITPI2101`, and
  `UKAB2101`.
- `stage2_biomass_after_phenology.yaml`: pooled Europe biomass/growth, fixed
  shared Stage 1 phenology values and one shared `IB0007` growth vector.
- `stage1_anthesis_lhs.yaml`: pooled Europe anthesis-only GLUE run using Latin
  hypercube sampling. Emergence is ignored, and only timing-related
  cultivar/ecotype coefficients expected to influence time to anthesis are
  active.
- `stage1_anthesis_lhs_wide.yaml`: same anthesis-only target as above, but with
  much wider bounds for the pre-flowering timing, daylength, and flowering
  temperature response coefficients.
- `stage1_anthesis_lhs_wide_no_photo.yaml`: same wider anthesis-only search, but
  with cultivar photoperiod response forced off by fixing `PPSEN = 0`, while
  holding `CSDL = 14.3` and `R1PPO = 0`.
- `stage1_phenology_emergence_temp_lhs_1h.yaml`: bounded one-hour test of
  emergence + anthesis phenology. It scores `EDAP` and `ADAP`, varies `PL-EM`,
  and opens gated `.SPE` edits for the vegetative temperature-response row used
  by planting-to-emergence.
- `stage1_phenology_emergence_temp_lhs_anthesis_priority.yaml`: same 512-sample
  emergence + anthesis setup, but with anthesis weighted much more strongly
  (`anthesis = 4.0`, `emergence = 0.5`).

The earlier Italy-only and UK-only Stage 2 configs were removed because they were
only needed when UK used a different cultivar anchor.

## Results

Stage 1 pooled phenology:

- Output: `results/europe_hemp_calibration/europe_hemp_stage1_phenology`
- Search budget: CMA-ES, 24 generations, population 16
- Anthesis: n = 12, RMSE = 16.391 d, nRMSE = 15.710%, MBE = -3.333 d
- Emergence: n = 12, RMSE = 3.162 d, nRMSE = 28.748%, MBE = -2.667 d

The shared `IB0007` phenology vector is a compromise. Italy anthesis is too
early, while UK anthesis is too late:

- `ITPI2001` anthesis: RMSE = 18.000 d, MBE = -18.000 d
- `ITPI2101` anthesis: RMSE = 11.000 d, MBE = -11.000 d
- `UKAB2101` anthesis: RMSE = 19.000 d, MBE = 19.000 d

A 96-point Latin-hypercube probe plus the previous Italy-tuned and UK-tuned
vectors did not find a lower Stage 1 score than the CMA-ES shared-cultivar
result. A later 48-generation attempt reached generation 24/48 and slightly
improved the score, but exited before writing final artifacts; the completed
24-generation run above preserves that larger-search result.

Anthesis-only Latin hypercube:

- Config: `stage1_anthesis_lhs.yaml`
- Output: `results/europe_hemp_calibration/europe_hemp_stage1_anthesis_lhs`
- Search budget: 1024 Latin-hypercube samples, GLUE scoring, no optimizer
- Objective: anthesis (`ADAP`) only; emergence and all growth variables ignored
- Active timing parameters: `CSDL`, `PPSEN`, `EM-FL`, `THVAR`, `PL-EM`,
  `EM-V1`, `V1-JU`, `JU-R0`, `R1PPO`, `OPTBI`, `SLOBI`
- Anthesis: n = 12, RMSE = 16.503 d, nRMSE = 15.817%, MBE = -4.333 d

The best LHS sample was close to, but did not beat, the completed CMA-ES
anthesis RMSE:

- `ITPI2001` anthesis: RMSE = 18.000 d, MBE = -18.000 d
- `ITPI2101` anthesis: RMSE = 13.000 d, MBE = -13.000 d
- `UKAB2101` anthesis: RMSE = 18.000 d, MBE = 18.000 d

Best anthesis-only LHS vector:

- `CSDL = 16.791585817198346`
- `PPSEN = 0.1825770546319404`
- `EM-FL = 67.22159582079411`
- `THVAR = 0.04363804701443986`
- `PL-EM = 9.077477280866864`
- `EM-V1 = 7.89578369740768`
- `V1-JU = 14.292784623408997`
- `JU-R0 = 16.062737852938106`
- `R1PPO = 0.2844931019353082`
- `OPTBI = 20.457855703781078`
- `SLOBI = 0.02673691348312649`

Wider anthesis-only Latin hypercube:

- Config: `stage1_anthesis_lhs_wide.yaml`
- Output: `results/europe_hemp_calibration/europe_hemp_stage1_anthesis_lhs_wide`
- Search budget: 1024 Latin-hypercube samples, GLUE scoring, no optimizer
- Anthesis: n = 12, RMSE = 12.356 d, nRMSE = 11.843%, MBE = -1.333 d

The wider photoperiod-active search improved the pooled fit, mainly by matching
`ITPI2101`, but it still did not make all sites flower at the observed similar
DAP:

- `ITPI2001` anthesis: observed 106 DAP, simulated 89 DAP, bias = -17 d
- `ITPI2101` anthesis: observed 104 DAP, simulated 104 DAP, bias = 0 d
- `UKAB2101` anthesis: observed 103 DAP, simulated 116 DAP, bias = 13 d

Best wide photoperiod-active vector:

- `CSDL = 14.278927815058173`
- `PPSEN = 1.7217836008892762`
- `EM-FL = 19.79444700137449`
- `THVAR = 0.20587871111881126`
- `PL-EM = 18.22514658854301`
- `EM-V1 = 7.988675599629143`
- `V1-JU = 17.366376609689183`
- `JU-R0 = 55.71841587243925`
- `R1PPO = 2.627170064059839`
- `OPTBI = 2.7820496193846607`
- `SLOBI = 0.18156054107154096`

Wider anthesis-only LHS with photoperiod response off:

- Config: `stage1_anthesis_lhs_wide_no_photo.yaml`
- Output:
  `results/europe_hemp_calibration/europe_hemp_stage1_anthesis_lhs_wide_no_photo`
- Fixed values: `CSDL = 14.3`, `PPSEN = 0.0`, `R1PPO = 0.0`
- Anthesis: n = 12, RMSE = 16.381 d, nRMSE = 15.701%, MBE = -2.333 d

Turning cultivar photoperiod response off did not resolve the site conflict:

- `ITPI2001` anthesis: observed 106 DAP, simulated 88 DAP, bias = -18 d
- `ITPI2101` anthesis: observed 104 DAP, simulated 95 DAP, bias = -9 d
- `UKAB2101` anthesis: observed 103 DAP, simulated 123 DAP, bias = 20 d

Temperature forcing check:

- `ITPI2001.WTH` and `ITPI2101.WTH` use local `UNICATT, Piacenza, Italy`
  weather at latitude 45.008 N.
- `UKAB2101.WTH` uses `NASA` weather at latitude 52.434 N, longitude -4.018,
  with NASA/POWER tile elevation 255.11 m.
- From planting to observed anthesis:
  - `ITPI2001`: 106 d, mean Tmean = 19.47 C, base-1 GDD = 1976.3
  - `ITPI2101`: 104 d, mean Tmean = 18.95 C, base-1 GDD = 1885.2
  - `UKAB2101`: 103 d, mean Tmean = 12.46 C, base-1 GDD = 1191.7

The observed flowering DAP is almost identical across sites, but the weather
forcing gives Aberystwyth far less thermal time. This means a normal
temperature-driven phenology response will tend to delay UK flowering relative
to Italy. The wide photoperiod-active run partly compensated by selecting a very
low `OPTBI` (2.78 C), which nearly removes the extra TMIN-based slowing toward
flowering, but the fixed species temperature response still makes UK much cooler
developmentally.

Emergence + anthesis one-hour LHS test:

- Config: `stage1_phenology_emergence_temp_lhs_1h.yaml`
- Output:
  `results/europe_hemp_calibration/europe_hemp_stage1_phenology_emergence_temp_lhs_1h`
- Search budget: 512 Latin-hypercube samples, GLUE scoring, no optimizer
- Runtime: about 40 minutes on this machine
- Active dimensions: 15
  - Cultivar/ecotype timing and flowering response: `CSDL`, `PPSEN`, `EM-FL`,
    `THVAR`, `PL-EM`, `EM-V1`, `V1-JU`, `JU-R0`, `R1PPO`, `OPTBI`, `SLOBI`
  - Species vegetative temperature response: `VEG_TB`, `VEG_TO1`, `VEG_TO2`,
    `VEG_TM`, targeting the `VEGETATIVE DEVELOPMENT` `.SPE` row
- Harness change: species `.SPE` edits now apply one token at a time so multiple
  calibrated parameters can target different numeric tokens on the same `.SPE`
  line without overwriting each other.

Best fit:

- Emergence: n = 12, RMSE = 1.291 d, nRMSE = 11.736%, MBE = -1.000 d
- Anthesis: n = 12, RMSE = 17.616 d, nRMSE = 16.885%, MBE = 0.333 d

Per-site best fit:

- `ITPI2001` emergence: observed 8 DAP, simulated 7 DAP, bias = -1 d
- `ITPI2101` emergence: observed 11 DAP, simulated 9 DAP, bias = -2 d
- `UKAB2101` emergence: observed 14 DAP, simulated 14 DAP, bias = 0 d
- `ITPI2001` anthesis: observed 106 DAP, simulated 91 DAP, bias = -15 d
- `ITPI2101` anthesis: observed 104 DAP, simulated 95 DAP, bias = -9 d
- `UKAB2101` anthesis: observed 103 DAP, simulated 128 DAP, bias = 25 d

Best vector:

- `CSDL = 20.399860163905814`
- `PPSEN = 1.2136871565677554`
- `EM-FL = 69.99076976663743`
- `THVAR = 0.255298013570393`
- `PL-EM = 22.10815576904204`
- `EM-V1 = 10.60557775921988`
- `V1-JU = 21.594532239870716`
- `JU-R0 = 75.52209825588342`
- `R1PPO = 0.14793100080155613`
- `OPTBI = 11.472526488534456`
- `SLOBI = 0.04145413756730815`
- `VEG_TB = 1.860050165816975`
- `VEG_TO1 = 20.832928636779926`
- `VEG_TO2 = 33.226980080833854`
- `VEG_TM = 44.62496797020597`

Anthesis-priority emergence + anthesis test:

- Config: `stage1_phenology_emergence_temp_lhs_anthesis_priority.yaml`
- Output:
  `results/europe_hemp_calibration/europe_hemp_stage1_phenology_emergence_temp_lhs_anthesis_priority`
- Search budget: 512 Latin-hypercube samples, GLUE scoring, no optimizer
- Weights: `anthesis = 4.0`, `emergence = 0.5`
- The run reached the one-hour cap while writing final outputs; result tables
  were written successfully.

Best fit:

- Anthesis: n = 12, RMSE = 17.253 d, nRMSE = 16.536%, MBE = -6.333 d
- Emergence: n = 12, RMSE = 5.066 d, nRMSE = 46.057%, MBE = 4.333 d

Per-site best fit:

- `ITPI2001` anthesis: observed 106 DAP, simulated 85 DAP, bias = -21 d
- `ITPI2101` anthesis: observed 104 DAP, simulated 90 DAP, bias = -14 d
- `UKAB2101` anthesis: observed 103 DAP, simulated 119 DAP, bias = 16 d
- `ITPI2001` emergence: observed 8 DAP, simulated 10 DAP, bias = 2 d
- `ITPI2101` emergence: observed 11 DAP, simulated 14 DAP, bias = 3 d
- `UKAB2101` emergence: observed 14 DAP, simulated 22 DAP, bias = 8 d

Best anthesis-priority vector:

- `CSDL = 19.613104283509323`
- `PPSEN = 0.02193693599013577`
- `R1PPO = 2.925244595027315`
- `EM-FL = 62.187970153001224`
- `THVAR = 0.09722997253552969`
- `PL-EM = 13.539417309646351`
- `EM-V1 = 4.977638775153402`
- `V1-JU = 43.9010716685267`
- `JU-R0 = 26.137616720201883`
- `OPTBI = 14.297618886290165`
- `SLOBI = 0.1584473425309307`
- `VEG_TB = 2.9961675270439057`
- `VEG_TO1 = 20.912945106996222`
- `VEG_TO2 = 31.21242921146042`
- `VEG_TM = 40.593128876047814`

Stage 2 pooled biomass/growth baseline:

- Config: `stage2_biomass_after_phenology.yaml`
- Fixed phenology values were updated from the completed 24-generation Stage 1.
- A one-theta score of the configured Stage 2 starts gave:
  - Biomass: n = 80, RMSE = 2180.611, nRMSE = 26.459%, MBE = -1236.685
  - LAI: n = 47, RMSE = 2.975, nRMSE = 70.604%, MBE = 1.930
  - Height: n = 76, RMSE = 0.614, nRMSE = 35.499%, MBE = 0.260
  - Leaf: n = 80, RMSE = 1146.715, nRMSE = 55.468%, MBE = 620.806
  - Stem: n = 80, RMSE = 2677.505, nRMSE = 44.008%, MBE = -1781.745

The Stage 2 config is prepared for a larger CMA-ES pass (24 generations,
population 16), but this was not launched after the long Stage 1 run because the
completed phenology improvement was modest and the 24-generation Stage 1 run
already exceeded the expected 1-2 hour runtime on this machine.

- `FL-LF = 35.0`
- `LFMAX = 1.4`
- `SLAVR = 300.0`
- `SIZLF = 250.0`
- `FL-VS = 10.0`
- `TRIFL = 0.3`
- `RWDTH = 1.0`
- `RHGHT = 1.5`

## Caveats

- The shared Futura75 setup is now mechanically consistent with the FileX files,
  but the current Stage 1 search exposes a site conflict: the same phenology
  coefficients make Italy flower early and UK flower late.
- The anthesis-only LHS search confirms the same conflict even after removing
  emergence from the objective and varying only anthesis-relevant timing
  parameters.
- The wider photoperiod-active LHS reduced the pooled anthesis RMSE, but the
  best vector is biologically stretched at several bounds/near-bounds
  (`PPSEN`, `THVAR`, `R1PPO`, `SLOBI`, `OPTBI`). This looks more like
  compensation for structural/site forcing issues than a stable cultivar
  estimate.
- The photoperiod-off run performed worse than the wide active run, so simply
  disabling cultivar daylength response is not enough to recover the observed
  similar DAP across Piacenza and Aberystwyth.
- The weather forcing check suggests UK temperature forcing should be audited
  before further cultivar tuning. The UK file uses NASA/POWER tile weather while
  the Italy files use local site weather; the observed equal DAP is difficult to
  reconcile with the much cooler UK thermal-time forcing.
- Adding emergence scoring and varying the vegetative temperature-response row
  can fit emergence well, but in the 512-sample test it worsened anthesis. This
  reinforces that emergence and flowering should probably be staged or weighted
  carefully, rather than letting emergence dominate a shared phenology search.
- Increasing the anthesis weight in the same 512-sample design did not recover
  the previous anthesis-only fit. The best solution drove `PPSEN` nearly to zero
  and used a large `R1PPO`, but anthesis remained too early in Italy and too late
  in the UK while emergence degraded.
- The larger Stage 1 search reduced the anthesis bias for Italy but increased
  the late UK anthesis bias and worsened emergence. This reinforces that a
  shared cultivar-only phenology vector is probably absorbing site, weather, or
  management effects.
- Because Stage 2 freezes the imperfect shared phenology vector, growth fitting
  should be treated as exploratory.
- `GWAD` and `FLWAD` observations are not used in this Stage 2 objective because
  the China biomass stage also focuses on biomass, leaf, stem, height, LAI, and
  phenology. Mapping these would be a Stage 3/all-observed extension.

## Verification

- Confirmed `ITPI2001.HMX`, `ITPI2101.HMX`, and `UKAB2101.HMX` all contain
  `HM IB0007 Futura75`.
- Loaded the Europe configs with Python `load_config`.
- Confirmed dimensions:
  - Stage 1: 14 active, 0 fixed
  - Stage 2: 8 active, 14 fixed
- `python -m pytest tests/test_config.py tests/test_observations.py tests/test_spawn.py -q --basetemp results/_pytest_tmp -o cache_dir=results/_pytest_cache`
  passed with 25 tests.
- After adding same-line `.SPE` token edits and the emergence-temperature config,
  `python -m pytest tests/test_config.py tests/test_observations.py tests/test_spawn.py tests/test_writers.py -q --basetemp results/_pytest_tmp -o cache_dir=results/_pytest_cache`
  passed with 30 tests.

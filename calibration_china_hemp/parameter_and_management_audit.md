# China hemp parameter and management audit

This audit captures the July 2026 calibration setup for the China hemp data in
`C:/Users/alwin/Documents/GitHub/DSSAT48Hemp/Hemp`.

## Active genetic parameters

The base hemp configs now expose the full expanded cultivar/ecotype search
space. Staged calibration narrows that set by data type:

- Stage 1 phenology/growth stages:
  - Cultivar: `CSDL`, `PPSEN`, `EM-FL`, `FL-SH`, `FL-SD`, `SD-PM`
  - Ecotype: `THVAR`, `PL-EM`, `EM-V1`, `V1-JU`, `JU-R0`, `R1PPO`,
    `OPTBI`, `SLOBI`
- Stage 2 biomass/vegetative growth:
  - Cultivar: `FL-LF`, `LFMAX`, `SLAVR`, `SIZLF`
  - Ecotype: `FL-VS`, `TRIFL`, `RWDTH`, `RHGHT`
  - IB0008 phenology and early-stage parameters from Stage 1 are fixed during
    this stage.
- Stage 3 all observed variables:
  - Cultivar: `CSDL`, `PPSEN`, `EM-FL`, `FL-SH`, `FL-SD`, `SD-PM`,
    `FL-LF`, `LFMAX`, `SLAVR`, `SIZLF`, `XFRT`, `WTPSD`, `SFDUR`,
    `SDPDV`, `PODUR`, `THRSH`
  - Ecotype: `THVAR`, `PL-EM`, `EM-V1`, `V1-JU`, `JU-R0`, `PM09`,
    `LNGSH`, `FL-VS`, `TRIFL`, `RWDTH`, `RHGHT`, `R1PPO`, `OPTBI`,
    `SLOBI`

The new stage files are rebuilt from the expanded `config_hemp.yaml` ranges
rather than inheriting the narrower post-hoc ranges from the previous run.
`YUKU2202` and `YUKU2203` are excluded from the stage configs.

## Management inputs

Planting date, emergence date, planting population, emergence population,
row spacing, planting depth, fertilizer rows, and irrigation rows were extracted
from the FileX treatment sections into:

`calibration_china_hemp/management_inputs_from_filex.csv`

The calibration parameter specs do not currently activate any `management` or
`initial_conditions` parameters. During a spawn, `spawn_and_run` copies the
source HMX file into the run directory and only calls `edit_filex` when such
parameters are active, when an explicit `_planting_dates` override is present,
or when `filex_overrides` are configured. For the China stage configs,
management parameters are inactive, so DSSAT uses the planting date, planting
population, row spacing, depth, fertilizer, and irrigation practices already
encoded in each treatment of the experimental HMX files.

The one deliberate FileX environment override is `CNKU2101`: the source HMX
references Florida/Citra weather and soil identifiers, so spawned calibration
runs replace `WSTA` with `CNKU2101` and `ID_SOIL` with `YUKU2101`.

## Final corrected Stage 2 run

After fixing the `CNKU2101` weather/soil override and fixed-column code
alignment, the final Stage 2 biomass/growth run is:

`results/china_hemp_calibration/china_hemp_stage2_biomass_after_phenology_cultivar`

All six best spawns succeeded. The final best-fit summary is:

- `LAI`: n = 16, RMSE = 2.513, nRMSE = 100.757%, MBE = -1.803
- `biomass`: n = 51, RMSE = 5356.000, nRMSE = 64.263%, MBE = -3271.300
- `height`: n = 79, RMSE = 1.702, nRMSE = 72.273%, MBE = -1.422
- `leaf`: n = 51, RMSE = 1118.392, nRMSE = 57.299%, MBE = -226.035
- `stem`: n = 51, RMSE = 5223.795, nRMSE = 82.671%, MBE = -3633.400

The calibrated genotype copies are saved in the result directory as
`HMGRO048.CUL`, `HMGRO048.ECO`, and `calibrated_genotype_parameters.csv` with
new IDs `IBCN02`/`HMCN02` for `IB0002` and `IBCN08`/`HMCN08` for `IB0008`.

A spot check on `YUKU2101` confirmed that the spawned `YUKU2101.HMX` remained
byte-identical to the source file after calibration setup and a DSSAT run.

## Verification

- `python -m pytest tests/test_config.py tests/test_observations.py tests/test_spawn.py -q`
  passed with 16 tests.
- `python -m pytest tests/test_config.py tests/test_spawn.py tests/test_observations.py tests/test_writers.py tests/test_writers_ecotype.py tests/test_writers_filex.py tests/test_writers_soil_weather.py tests/test_new_engines.py -q`
  passed with 47 tests and 4 skipped tests.
- `Rscript -e "pkgload::load_all('.', quiet=TRUE); testthat::test_file('tests/testthat/test-pooled-calibration.R', reporter='summary'); testthat::test_file('tests/testthat/test-parity-writers.R', reporter='summary')"`
  passed.
- R source parsing passed for all files under `R/`.
- Config loading confirms 14 active Stage 1 parameter specs and 8 active
  Stage 2 biomass/growth specs, with the Stage 2 specs expanding to cultivar-
  specific dimensions for `IB0002` and `IB0008`.
- `YUKU2202` and `YUKU2203` are absent from the China calibration stage configs.

# China hemp experimental translation audit

This audit checks that the original DSSAT hemp experimental files are translated
into calibration inputs consistently.

## Findings

- The China calibration configs use six trials: `CNKU2101`, `YUBA2101`,
  `YUFE2101`, `YUFE2201`, `YUKU2101`, and `YUKU2201`.
- `YUKU2202` and `YUKU2203` remain excluded.
- FileX management is copied into each spawn unchanged unless management
  parameters are explicitly active. Current China phenology configs keep
  management inactive.
- `YUBA2101.HMA` encodes `EDATE` and `ADAT` as three-digit day-of-year values
  (`137`, `216`) rather than `YYDDD`. The observation readers now interpret
  these using the experiment year, giving `2021-05-17` and `2021-08-04`.
- For every treatment with anthesis observations, FileA/FileX-derived anthesis
  days after planting match DSSAT `Evaluate.OUT` measured `ADAPM` exactly.
- `CNKU2101` has FileT growth observations but no FileA anthesis observations,
  so it contributes no anthesis residuals in anthesis-only calibration.
- `CNKU2101` references cultivar `IB0002`; the other anthesis trials reference
  `IB0008`. The spawn code edits all configured China cultivar anchors
  (`IB0001`, `IB0002`, and `IB0008`) so the shared China coefficient vector is
  applied consistently when those cultivars occur in FileX files.
- `CNKU2101.HMX` points to weather station `CTRA2101` and soil `IBSB910015`,
  which are Florida/Citra inputs. The China calibration configs now apply a
  spawn-time FileX override for `CNKU2101` so run copies use station
  `CNKU2101` and soil `YUKU2101` (Kunming/Chenggong) without modifying the
  original HMX file.

The row-level audit is saved as:

`calibration_china_hemp/experimental_translation_audit.csv`

## Parser/writer corrections

- `read_filea` in both Python and R now treats `EDATE` as a phenology date
  column.
- FileA date columns now accept either `YYDDD` or bare day-of-year values when
  an experiment year can be inferred from the experiment ID.
- The Python cultivar/ecotype writers now format replacement values like the
  existing field token. This avoids changing DSSAT behavior merely by rewriting
  a coefficient at its existing start value.
- `edit_filex` now supports generic section-row updates for code fields such as
  `WSTA` and `ID_SOIL`; the spawn runner applies these from the optional
  `filex_overrides` config block.
- FileX code replacement preserves the original leading offset inside widened
  cells, so `CNKU2101` spawned rows now carry `ID_SOIL` as
  `YUKU2101` in the same fixed-column position used by the working `YUKU2101`
  experiment.

## Flowering parameters now included

The expanded anthesis calibration includes the direct and plausible flowering
controls:

- Cultivar: `CSDL`, `PPSEN`, `EM-FL`, `FL-SH`, `FL-SD`, `SD-PM`
- Ecotype: `THVAR`, `PL-EM`, `EM-V1`, `V1-JU`, `JU-R0`, `R1PPO`, `OPTBI`,
  `SLOBI`

`THVAR`, `OPTBI`, and `SLOBI` were added to the active base ecotype parameter
space because the ECO header identifies them as reproductive development /
flowering-temperature controls.

## Expanded anthesis result

The expanded anthesis-only CMA-ES run is:

`calibration_china_hemp/stage1_anthesis_flowering_expanded.yaml`

It wrote results to:

`results/china_hemp_calibration/china_hemp_stage1_anthesis_flowering_expanded`

Best anthesis fit:

- `n = 16`
- `RMSE = 8.396 days`
- `bias = -4.0 days`
- `R2 = 0.965`

The calibrated genotype copies are saved in that result directory with new IDs:

- Cultivar: `IBCN01` (`YUNMA8_CN_ANTH`)
- Ecotype: `HMCN01` (`China hemp anth`)

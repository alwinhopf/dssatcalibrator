# Germany Hemp Calibration Run Notes

## Runnable scaffold

- FileX files installed in `C:/Users/alwin/Documents/GitHub/DSSAT48Hemp/Hemp`:
  `GEMQ18S1.HMX` for Santhica 27 and `GEMQ18I1.HMX` for Ivory.
- Provisional cultivar/ecotype rows installed if missing:
  `DE0001/HMDE01` for Santhica 27 and `DE0002/HMDE02` for Ivory.
- Weather installed as `DSSAT48Hemp/Weather/GEMQ2018.WTH`; source is DWD
  Potsdam station 03987, with the DSSAT weather header INSI harmonized to `GEMQ`.
- Soil profile `*GEMQ2018` appended to `DSSAT48Hemp/Soil/SOIL.SOL` from the
  reported low loamy sand profile, cross-checked against SoilGrids.
- Kieserite/MgSO4 is intentionally omitted from FileX.

## Management assumptions currently encoded

- Planting depth: 2 cm.
- Planned population: 200 seeds/plants m-2 in `PPOP`; emerged stand in `PPOE`
  from reported start-density values.
- Row spacing: 12.5 cm.
- Initial soil water: 75 percent relative soil moisture, converted per layer as
  `SLLL + 0.75 * (SDUL - SLLL)`, giving SH2O about 0.219 cm3 cm-3.
- Initial residual mineral N: 100 kg ha-1 total, approximated as 0.5 ppm NH4-N
  plus 2.9 ppm NO3-N in each profile layer through 200 cm.
- Fertilizer: 70 kg N ha-1 as calcium ammonium nitrate, assumed one day before
  each cultivar-specific planting date.
- Irrigation: 10 mm sprinkler event on 2018-05-30.
- No tillage and no residue/organic amendment sections are used.

## Phenology handling

The PDFs do not report direct flowering dates for either cultivar. The only
development information is indirect: Ivory is discussed as reaching maximum
photosynthesis around 52 DAS, possibly end flowering/beginning senescence, and
Santhica 27 is described as staying steady until about 70 DAS, consistent with
the beginning of senescence. Because these are not direct observed stage dates,
the current configs keep phenology coefficients active and let LAI/height/biomass
shape constrain development timing.

## Calibration runs started

- `stage1_phenology_lai_shape_smoke.yaml` completed successfully. It verified
  the Germany FileX/weather/soil/genotype chain and produced PlantGro output for
  both experiments. Best aggregate fit was poor but useful as a smoke test:
  LAI nRMSE 67.6 percent, height nRMSE 66.8 percent, biomass nRMSE 74.6 percent,
  stem nRMSE 81.3 percent.
- `stage2_biomass_lai_followup.yaml` completed successfully with phenology and
  growth/leaf parameters active. It improved biomass/stem but worsened LAI and
  height: biomass nRMSE 55.6 percent, stem nRMSE 70.4 percent, LAI nRMSE 75.7
  percent, height nRMSE 76.3 percent.
- `stage1_phenology_lai_shape_deep.yaml` completed successfully. Adding
  post-flowering duration parameters (`FL-SH`, `FL-SD`, `SD-PM`) and running a
  deeper CMA-ES pass improved canopy/height timing: LAI nRMSE 44.4 percent and
  height nRMSE 34.9 percent. Biomass remained weak: biomass nRMSE 62.4 percent,
  stem nRMSE 67.0 percent.
- `stage2_biomass_lai_expanded_from_stage1.yaml` was generated from the Stage 1
  best phenology starts and completed successfully. It widened growth bounds
  (`LFMAX`, `SLAVR`, `SIZLF`) and added `XFRT`. It improved harvest biomass
  relative to Stage 1 but degraded height: biomass nRMSE 49.8 percent, stem
  nRMSE 62.6 percent, LAI nRMSE 46.3 percent, height nRMSE 56.2 percent.

A compact comparison table is written to
`results/germany_hemp_run_comparison.csv`.

## Water-Stress Diagnostic Split

The diagnostic no-water-stress config
`stage2_biomass_lai_expanded_from_stage1_no_waterstress.yaml` was created from
the same seeded expanded Stage 2 config, with only one FileX override:
`*SIMULATION CONTROLS / @N OPTIONS / WATER = N`.

It completed successfully and changed the fit pattern:

- Current-water expanded Stage 2: biomass nRMSE 49.8 percent, stem nRMSE 62.6
  percent, LAI nRMSE 46.3 percent, height nRMSE 56.2 percent.
- No-water-stress expanded Stage 2: biomass nRMSE 37.8 percent, stem nRMSE 45.2
  percent, height nRMSE 32.4 percent, but LAI worsened to 69.3 percent.

Best-run review tables are written to
`results/germany_hemp_waterstress_diagnostic_best_review.csv`.

For Santhica 27, removing water stress improved harvest CWAD from about
8.2 t ha-1 to 12.0 t ha-1 and harvest SWAD from about 5.6 t ha-1 to 10.7 t ha-1.
This confirms that the Santhica biomass gap is partly water-stress/soil-water
driven. It does not fully close the gap to the reported 17.9 t ha-1 CWAD and
16.0 t ha-1 straw target.

## Flowering-Date Check

A fresh keyword audit is written to `derived/flowering_date_audit.md`, with
contexts in `derived/flowering_keyword_contexts.txt`.

No direct observed flowering dates, anthesis dates, or BBCH-style flowering
stage dates were found for Santhica 27 or Ivory. The only usable development
evidence is indirect:

- Ivory maximum photosynthesis around 52 DAS may indicate end flowering or
  beginning senescence by analogy to Chameleon.
- Santhica 27 is discussed as remaining steady until about 70 DAS, consistent
  with reported beginning of senescence timing.
- LAI decline after the mid-season peak is a growth-shape/senescence constraint,
  not a direct stage observation.

## Late-Flowering Probe

A targeted Santhica 27 probe is written to
`results/germany_hemp_late_flowering_probe_santhica.csv`.

The probe supports the hypothesis that later flowering can retain much more
end-season biomass, but only when water stress is also relaxed:

- Current-water best: Santhica harvest CWAD about 8.2 t ha-1; first flower around
  82 DAP.
- Current-water forced-late-flowering probe: harvest CWAD stayed about
  8.2 t ha-1; water stress still prevented recovery.
- No-water-stress best: Santhica harvest CWAD about 12.0 t ha-1; first flower
  around 53 DAP.
- No-water-stress forced-late-flowering probe: Santhica harvest CWAD reached
  about 16.0 t ha-1 and SWAD about 12.1 t ha-1, but LAI and height became too
  large (harvest LAI about 12.7, height about 3.2 m). This is close to the
  reported biomass target but biologically over-vegetative.

## High Initial Soil-Water Calibration

The next diagnostic increased initial soil water from the 75 percent relative
assumption (`SH2O = .219`) to an approximate field-capacity value
(`SH2O = .240`) while keeping water balance active. The config is
`stage2_biomass_lai_high_initial_water_long.yaml`.

Implementation note: the FileX override must write `.240` as text, not numeric
`0.240`, because the narrow DSSAT `SH2O` column otherwise becomes malformed
(`50.2400` instead of `5  .240`). The fixed writer check produced valid rows:
`5  .240`, `15  .240`, and `30  .240`.

Two runs were completed:

- `germany_hemp_stage2_biomass_lai_high_initial_water_smoke`: fixed high initial
  water at the previous Stage 2 best point, with only a short optimizer pass.
  This improved biomass nRMSE to 40.7 percent and raised Santhica harvest CWAD
  to about 9.9 t ha-1, but LAI remained only moderate.
- `germany_hemp_stage2_biomass_lai_high_initial_water_fc_long`: a longer CMA-ES
  run (`maxiter=56`, `popsize=28`, one restart; about 51 minutes). It produced
  the best LAI fit so far (LAI nRMSE 29.7 percent), but biomass nRMSE was
  49.4 percent because the optimizer favored canopy shape over final Santhica
  biomass.

Best-run review is written to
`results/germany_hemp_high_initial_water_fc_best_review.csv`.

For the long high-water run:

- Santhica 27 peaked at 14.3 t ha-1 CWAD and 9.4 t ha-1 SWAD, but harvested at
  only 8.3 t ha-1 CWAD and 7.1 t ha-1 SWAD. LAI peaked at 6.95 and declined to
  1.39 at harvest.
- Ivory harvested at 8.5 t ha-1 CWAD and 5.9 t ha-1 SWAD, close to its total
  biomass target but still low for straw.
- Santhica modeled first flower was 36 DAP, first seed 73 DAP, physiological
  maturity 100 DAP, end leaf 89 DAP, and harvest 138 DAP. Ivory modeled first
  flower was inferred around 49 DAP, first seed 65 DAP, physiological maturity
  68 DAP, harvest maturity 93 DAP, and harvest 120 DAP.
- Seasonal mean water stress remained substantial for Santhica (`WSPD` about
  0.21, `WSGD` about 0.23) and almost absent for Ivory. Higher initial water
  helped early growth potential but did not remove the late Santhica collapse.

## Current interpretation

The scaffold is runnable, but the first follow-up parameter space was not yet a
good calibration. Several parameters in that earlier Stage 2 run hit bounds
(`CSDL__DE0001`, `SIZLF__DE0001`, `FL-VS__DE0001`), and Santhica 27 remained
under-predicted at harvest. Across that candidate pool, the highest final
simulated CWAD was about 12.8 t ha-1 for Santhica 27 and 11.6 t ha-1 for Ivory,
compared with reported harvest dry mass targets of 17.9 and 10.0 t ha-1.

After the deeper runs, the interpretation is sharper:

- Stage 1 deep gives the best canopy timing and height behavior so far.
- The expanded Stage 2 can fit Ivory reasonably on harvest biomass and LAI, but
  Santhica still loses too much biomass after peak growth. In the reviewed best
  run, Santhica peaks at about 13.8 t ha-1 CWAD and harvests at 8.2 t ha-1,
  while the paper target is 17.9 t ha-1. Ivory peaks at about 11.2 t ha-1 and
  harvests at 8.2 t ha-1, closer to the 10.0 t ha-1 target.
- Santhica's reviewed best run flowers very late (first flower around 82 DAP)
  and then hits high late water stress after first seed. This aligns with the
  observed late-season senescence problem, but the modeled collapse is too large.
- Only `CSDL__DE0002` is still exactly at a Stage 2 expanded bound. The poor fit
  is therefore not just a simple bound problem; it is a structural/management
  tradeoff between late canopy duration, water stress, and vegetative biomass
  accumulation.
- The next calibration should not use a hard observed flowering target, because
  none exists. Instead, use flowering as a latent parameter and constrain it with
  LAI decline, height plateau, and harvest biomass. A promising next run is a
  bounded late-flowering Santhica scenario with relaxed, but not disabled, water
  stress. The high-initial-water run shows that initial soil water alone is not
  enough; the remaining issue is late-season Santhica water stress/senescence.
  The next test should therefore adjust soil-water supply later in the season
  (soil profile/drainage/runoff or a modest supplemental irrigation diagnostic)
  while penalizing excessive LAI/height.

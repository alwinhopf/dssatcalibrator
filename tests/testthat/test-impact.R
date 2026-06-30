test_that("impact atlas CLI args expose R output controls", {
  args <- .impact_atlas_args(
    "run_impact_atlas.py", "config.yaml", output_dir = "out",
    experiments = c("EXP1", "EXP2"),
    groups = c("genetic_cultivar", "soil"),
    levels = c("low", "high"),
    discover_cultivar = TRUE,
    discover_ecotype = TRUE,
    discover_species = TRUE,
    discover_genotype = TRUE,
    allow_species = TRUE,
    max_per_group = 1L,
    output_files = c("PlantGro.OUT", "Summary.OUT"),
    num_cores = 4L,
    write_long = FALSE,
    compress_long = TRUE,
    effect_tolerance = 0.01,
    progress = FALSE
  )

  expect_equal(args[1:2], c("run_impact_atlas.py", "config.yaml"))
  expect_true("--outdir" %in% args)
  expect_true("--discover-cultivar" %in% args)
  expect_true("--discover-ecotype" %in% args)
  expect_true("--discover-species" %in% args)
  expect_true("--discover-genotype" %in% args)
  expect_true("--allow-species" %in% args)
  expect_true("--max-per-group" %in% args)
  expect_true("--outputs" %in% args)
  expect_true("--cores" %in% args)
  expect_true("--no-long" %in% args)
  expect_true("--compress-long" %in% args)
  expect_true("--effect-tolerance" %in% args)
  expect_true("--no-progress" %in% args)
})

test_that("impact atlas default output directory mirrors Python CLI", {
  cfg <- list(calibrator = list(name = "demo", results_dir = "results"))
  expect_equal(.impact_default_outdir(cfg), file.path("results", "demo_impact_atlas"))
})

test_that("impact atlas script path can come from environment", {
  old <- Sys.getenv("DSSATCAL_IMPACT_SCRIPT", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("DSSATCAL_IMPACT_SCRIPT") else Sys.setenv(DSSATCAL_IMPACT_SCRIPT = old)
  }, add = TRUE)
  Sys.setenv(DSSATCAL_IMPACT_SCRIPT = "C:/tools/run_impact_atlas.py")
  expect_equal(.impact_script_path(), normalizePath("C:/tools/run_impact_atlas.py", mustWork = FALSE))
})

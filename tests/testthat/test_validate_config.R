# Behaviour tests for validate_config (R twin of config.py:validate_config).
# Mirrors tests/test_config.py so the two languages reject the same configs.

test_that("the reference config validates", {
  cfg_path <- "../../config_hemp.yaml"
  if (!file.exists(cfg_path)) cfg_path <- "config_hemp.yaml"
  skip_if_not(file.exists(cfg_path), "config_hemp.yaml not found")
  cfg <- load_config(cfg_path)                           # load_config validates
  expect_identical(validate_config(cfg), invisible(cfg))
  expect_gt(length(active_parameters(cfg)), 0L)
})

test_that("validate_config collects every problem", {
  bad <- list(
    method = list(preset = "Z", validation = list(scheme = "weekly")),
    objective = list(weighting = "bogus"),
    execution = list(backend = "spark"),
    assimilation = list(mode = "kalman"),
    gating = list(cultivar = "loose"),
    calibrator = list(num_cores = -2L),
    parameters = list(g = list(
      A = list(active = TRUE, min = 5, max = 1),                 # inverted bounds
      B = list(active = TRUE, min = 0, max = 10, start = 50),    # start out of bounds
      C = list(active = TRUE, min = 0, max = 1, prior = list(dist = "cauchy"))
    ))
  )
  msg <- tryCatch(validate_config(bad), error = function(e) conditionMessage(e))
  for (needle in c("method.preset", "validation.scheme", "objective.weighting",
                   "execution.backend", "assimilation.mode", "gating.cultivar",
                   "num_cores", "min (5) must be < max (1)", "start (50) is outside",
                   "prior.dist 'cauchy'")) {
    expect_true(grepl(needle, msg, fixed = TRUE), info = needle)
  }
})

test_that("validate_config requires at least one active parameter", {
  cfg <- list(parameters = list(g = list(A = list(active = FALSE, min = 0, max = 1))))
  expect_error(validate_config(cfg), "No active parameters")
})

test_that("cultivar-scoped parameters and zero-observation defaults match Python", {
  cfg <- .dssatcal_defaults()
  expect_true("ignore_zero_observations" %in% names(cfg$objective))
  cfg$experiments <- list("E1")
  cfg$parameters <- list(g = list(A = list(active = TRUE, min = 0, max = 1,
                                               scope = "per_cultivar")))
  expect_identical(validate_config(cfg), invisible(cfg))
})

test_that("portable DSSAT path resolution uses a canonical executable name", {
  cfg_path <- "../../config_hemp.yaml"
  if (!file.exists(cfg_path)) cfg_path <- "config_hemp.yaml"
  skip_if_not(file.exists(cfg_path), "config_hemp.yaml not found")
  cfg <- load_config(cfg_path, validate = FALSE)
  exe <- resolve_exe(cfg)
  expect_true(tolower(basename(exe)) %in% c("dscsm048", "dscsm048.exe"))
  if (file.exists(exe)) {
    expect_true(dir.exists(resolve_dssat_paths(cfg)$genotype))
  }
})

test_that("DSSAT path resolution retains custom name when native discovery fails", {
  target_env <- environment(resolve_exe)
  if (identical(target_env, .GlobalEnv)) {
    original_discovery <- .workspace_dssat_root
    on.exit(assign(".workspace_dssat_root", original_discovery, envir = .GlobalEnv),
            add = TRUE)
    assign(".workspace_dssat_root", function() NULL, envir = .GlobalEnv)
  } else {
    local_mocked_bindings(
      .workspace_dssat_root = function() NULL,
      .package = environmentName(target_env)
    )
  }
  custom <- paste0(
    "C:/Users/example/DSSAT48Hemp/",
    "dscsm048_compiled_4.8.2.with_HM_code.exe"
  )
  cfg <- list(calibrator = list(
    dssat_dir = "C:/Users/example/DSSAT48Hemp",
    dssat_exe = custom
  ))

  expect_identical(resolve_exe(cfg), custom)
})

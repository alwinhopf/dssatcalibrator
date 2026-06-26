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

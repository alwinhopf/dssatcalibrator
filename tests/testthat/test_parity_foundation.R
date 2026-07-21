# Cross-language parity tests for the Phase 1 foundation modules.
# These load golden fixtures generated from the Python implementation
# (tests/generate_parity_fixtures.py) and assert the R twin reproduces them.
#
# Run from the repo root (R available):
#   Rscript -e 'testthat::test_dir("tests/testthat")'
# Regenerate fixtures first if Python changed:
#   PYTHONPATH=python python tests/generate_parity_fixtures.py

fixture_dir <- file.path("..", "fixtures")
if (!dir.exists(fixture_dir)) fixture_dir <- file.path("tests", "fixtures")
read_fix <- function(name) jsonlite::fromJSON(file.path(fixture_dir, name),
                                              simplifyVector = FALSE)

TOL <- 1e-9

test_that("priors::log_prior_one matches Python for every distribution", {
  gold <- read_fix("priors_logdensity.json")
  for (dist in names(gold)) {
    spec <- gold[[dist]]$spec
    pts <- unlist(gold[[dist]]$points)
    expected <- vapply(gold[[dist]]$log_prior, function(x) {
      if (is.character(x) && x %in% c("-inf", "-Inf")) -Inf else as.numeric(x)
    }, numeric(1))
    got <- vapply(pts, function(v) log_prior_one(spec, v), numeric(1))
    for (i in seq_along(pts)) {
      if (is.infinite(expected[i])) {
        expect_true(is.infinite(got[i]) && got[i] < 0,
                    info = sprintf("%s @ %s expected -Inf", dist, pts[i]))
      } else {
        expect_equal(got[i], expected[i], tolerance = TOL,
                     info = sprintf("%s @ %s", dist, pts[i]))
      }
    }
  }
})

test_that("objective::metrics matches Python", {
  gold <- read_fix("metrics.json")
  for (case in names(gold)) {
    obs <- unlist(gold[[case]]$input$obs)
    sim <- unlist(gold[[case]]$input$sim)
    exp_m <- gold[[case]]$metrics
    got <- metrics(obs, sim)
    expect_equal(got$n, as.integer(exp_m$n), info = case)
    for (k in c("RMSE", "nRMSE_pct", "MBE", "d", "EF", "R2")) {
      e <- exp_m[[k]]
      if (is.null(e) || (is.character(e) && grepl("nan|NA", e, ignore.case = TRUE))) {
        expect_true(is.na(got[[k]]), info = sprintf("%s$%s expected NA", case, k))
      } else {
        expect_equal(got[[k]], as.numeric(e), tolerance = 1e-9,
                     info = sprintf("%s$%s", case, k))
      }
    }
  }
})

test_that("config::load_config + active_parameters + ParameterSpace match Python", {
  gold <- read_fix("config_space.json")
  cfg <- load_config(file.path(fixture_dir, "_sample_config.yaml"))

  act <- active_parameters(cfg)
  expect_equal(vapply(act, function(a) a$name, character(1)),
               unlist(gold$active_names))
  expect_equal(vapply(act, function(a) a$group, character(1)),
               unlist(gold$active_groups))

  # defaults survive the merge untouched
  expect_equal(cfg$method$preset, gold$merged_preset)
  expect_equal(cfg$calibrator$seed, as.integer(gold$merged_seed))
  expect_equal(cfg$gating$ecotype, gold$default_gating_ecotype)

  space <- parameter_space_from_config(cfg)
  expect_equal(space$names, unlist(gold$space_names))
  expect_equal(space$low, unlist(gold$space_low), tolerance = TOL)
  expect_equal(space$high, unlist(gold$space_high), tolerance = TOL)
  expect_equal(space$start, unlist(gold$space_start), tolerance = TOL)
  expect_equal(ps_ndim(space), gold$ndim)
})

test_that("configured zero observations are filtered like Python", {
  resid <- data.frame(
    user_var = c("biomass", "biomass", "grain"),
    dssat = c("CWAD", "CWAD", "GWAD"),
    obs = c(0, 8000, 0), sim = c(100, 8000, 0)
  )
  cfg <- list(objective = list(ignore_zero_observations = list("biomass")))
  filtered <- .drop_configured_zero_observations(resid, cfg)
  expect_equal(filtered$obs, c(8000, 0))
  expect_equal(filtered$user_var, c("biomass", "grain"))
})

test_that("objective scores only configured PlantGro time-series columns", {
  dates <- as.Date(c("2021-09-17", "2021-10-05"))
  cfg <- list(
    engine = list(
      timeseries_outputs = list(biomass = "CWAD"),
      scalar_outputs = list(anthesis = "ADAP", grain_yield = "HWAM")
    ),
    objective = list(weighting = "unified", weights = list(), error_model = list())
  )
  result <- list(
    evaluate = data.frame(
      treatment = c(1L, 1L), variable = c("ADAP", "HWAM"),
      sim = c(75, 1000), meas = c(75, 1000)
    ),
    plantgro = data.frame(
      treatment = c(1L, 1L), date = dates, DAP = c(75, 93),
      CWAD = c(5000, 8000), RWAD = c(250, 300)
    )
  )
  obs <- data.frame(
    exp_id = "E1", treatment = 1L,
    variable = c("CWAD", "CWAD", "RWAD"), kind = "timeseries",
    date = c(dates, dates[1]), value = c(5000, 8000, 250),
    sigma = NA_real_, weight = 1, stringsAsFactors = FALSE
  )

  resid <- build_residuals(list(E1 = result), obs, cfg)

  expect_false("RWAD" %in% resid$dssat)
  expect_equal(sum(resid$kind == "timeseries"), 2L)
})

test_that("FileA phenology dates map to configured DAP outputs", {
  dates <- as.Date(c("2021-09-17", "2021-10-05"))
  cfg <- list(
    engine = list(timeseries_outputs = list(biomass = "CWAD"),
                  scalar_outputs = list(anthesis = "ADAP")),
    objective = list(weighting = "unified", weights = list(), error_model = list())
  )
  result <- list(
    evaluate = data.frame(treatment = 1L, variable = "ADAP", sim = 75, meas = NA_real_),
    plantgro = data.frame(treatment = c(1L, 1L), date = dates,
                          DAP = c(75, 93), CWAD = c(5000, 8000))
  )
  obs <- data.frame(
    exp_id = "E1", treatment = 1L, variable = "ADAT", kind = "phenology",
    date = dates[1], value = 21260, sigma = NA_real_, weight = 1,
    stringsAsFactors = FALSE
  )

  resid <- build_residuals(list(E1 = result), obs, cfg)

  expect_equal(unmatched_variables(obs, cfg), character(0))
  expect_equal(resid$user_var, "anthesis")
  expect_equal(resid$dssat, "ADAP")
  expect_equal(resid$obs, 75)
  expect_equal(resid$sim, 75)
})

test_that("observation sigma and weight override objective defaults", {
  date <- as.Date("2021-09-17")
  cfg <- list(
    engine = list(timeseries_outputs = list(biomass = "CWAD"),
                  scalar_outputs = list()),
    objective = list(
      weighting = "unified", weights = list(biomass = 0.5),
      error_model = list(biomass = list(type = "relative", value = 0.2))
    )
  )
  result <- list(
    evaluate = data.frame(),
    plantgro = data.frame(treatment = 1L, date = date, DAP = 75, CWAD = 5000)
  )
  obs <- data.frame(
    exp_id = "E1", treatment = 1L, variable = "CWAD", kind = "timeseries",
    date = date, value = 5000, sigma = 50, weight = 0.25,
    stringsAsFactors = FALSE
  )

  resid <- build_residuals(list(E1 = result), obs, cfg)

  expect_equal(resid$sigma, 50)
  expect_equal(resid$weight, 0.25)
})

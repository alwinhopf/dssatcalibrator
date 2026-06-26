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

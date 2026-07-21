test_that("report CV uses disjoint top-level experiment subsets", {
  seen <- new.env(parent = emptyenv())
  seen$train <- list()
  seen$test <- list()
  fake_setup <- function(work_cfg) {
    list(experiments = as.character(work_cfg$experiments))
  }
  fake_calibrate <- function(work_cfg, progress = TRUE) {
    exps <- as.character(work_cfg$experiments)
    seen$train[[length(seen$train) + 1L]] <- exps
    list(
      best_theta = list(P1 = as.numeric(length(exps))),
      best = list(
        score = as.numeric(length(exps)), loglik = -0.5 * length(exps),
        residuals = data.frame(),
        per_var = list(yield = list(RMSE = as.numeric(length(exps)), d = 0.9))
      )
    )
  }
  fake_evaluate <- function(work_cfg, thetas, setup = NULL,
                            n_workers = NULL, progress = FALSE) {
    exps <- as.character(work_cfg$experiments)
    seen$test[[length(seen$test) + 1L]] <- exps
    list(results = list(list(
      score = 10 + length(exps), loglik = -0.5 * (10 + length(exps)),
      residuals = data.frame(),
      per_var = list(yield = list(RMSE = 2 * length(exps), d = 0.8))
    )), setup = setup)
  }
  target_env <- environment(run_cross_validation)
  if (identical(target_env, .GlobalEnv)) {
    original_setup <- .setup
    original_calibrate <- calibrate
    original_evaluate <- evaluate_thetas
    on.exit({
      assign(".setup", original_setup, envir = .GlobalEnv)
      assign("calibrate", original_calibrate, envir = .GlobalEnv)
      assign("evaluate_thetas", original_evaluate, envir = .GlobalEnv)
    }, add = TRUE)
    assign(".setup", fake_setup, envir = .GlobalEnv)
    assign("calibrate", fake_calibrate, envir = .GlobalEnv)
    assign("evaluate_thetas", fake_evaluate, envir = .GlobalEnv)
  } else {
    local_mocked_bindings(
      .setup = fake_setup,
      calibrate = fake_calibrate,
      evaluate_thetas = fake_evaluate,
      .package = environmentName(target_env)
    )
  }

  cfg <- list(
    experiments = c("E1", "E2", "E3"),
    calibrator = list(seed = 7L, num_cores = 1L),
    cross_validation = list(
      enabled = TRUE, strategy = "k_fold", k = 3L, final_theta = "best_fold"
    )
  )
  result <- run_cross_validation(cfg, progress = FALSE)

  expect_s3_class(result, "cv_result")
  expect_length(result$folds, 3L)
  expect_equal(sort(unlist(seen$test)), c("E1", "E2", "E3"))
  expect_true(all(vapply(seen$train, length, integer(1)) == 2L))
  for (i in seq_along(seen$train)) {
    expect_length(intersect(seen$train[[i]], seen$test[[i]]), 0L)
    expect_equal(sort(union(seen$train[[i]], seen$test[[i]])), cfg$experiments)
  }
  expect_null(cfg$calibrator$experiments)
  expect_equal(result$final_theta_method, "best_fold")
  expect_equal(result$report_df$fold, 1:3)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    report_dir <- tempfile("cv-report-")
    paths <- write_cv_report(result, report_dir)
    expect_true(all(file.exists(unlist(paths))))
    expect_equal(utils::read.csv(paths$csv)$fold, 1:3)
    summary <- jsonlite::fromJSON(paths$json)
    expect_equal(summary$recommendation, result$summary$recommendation)
  }

  seen$train <- list()
  seen$test <- list()
  temporal_cfg <- cfg
  temporal_cfg$experiments <- c("E2001", "E2101", "E2201")
  temporal_cfg$cross_validation$strategy <- "temporal_forward"

  temporal <- run_cross_validation(temporal_cfg, progress = FALSE)

  expect_length(temporal$folds, 2L)
  expect_equal(seen$test, list("E2101", "E2201"))
  expect_equal(seen$train, list("E2001", c("E2001", "E2101")))
})

test_that("report CV configuration is validated", {
  expect_error(
    parse_cv_config(list(cross_validation = list(strategy = "k_fold", k = 1L))),
    "at least 2"
  )
  expect_error(
    parse_cv_config(list(cross_validation = list(strategy = "not-a-strategy"))),
    "Unknown cross-validation strategy"
  )
  expect_error(
    run_cross_validation(list(cross_validation = list(enabled = FALSE)), progress = FALSE),
    "not enabled"
  )
})

test_that("persistent evaluation cache round-trips objective_result objects", {
  tmp <- tempfile("eval_cache_roundtrip_")
  dir.create(tmp, recursive = TRUE)
  exe <- file.path(tmp, "dscsm048")
  writeLines("fake exe", exe)
  cfg <- list(
    calibrator = list(name = "cache_test", workdir = file.path(tmp, "work"),
                      dssat_dir = tmp, dssat_exe = exe, cache_evaluations = TRUE),
    source = list(hemp_dir = tmp),
    experiments = list("E1"),
    parameters = list(genetic_cultivar = list(
      P1 = list(active = TRUE, min = 0, max = 10, start = 5)
    )),
    crops = list(list(code = "HM")),
    engine = list(scalar_outputs = list(HWAM = "HWAM")),
    objective = list(weighting = "unified", weights = list(), error_model = list())
  )
  space <- parameter_space_from_config(cfg)
  cache <- .evaluation_cache_from_setup(
    cfg,
    crop = list(code = "HM", filex_ext = "HMX", genotype_stem = "HMGRO048"),
    specs = space$specs,
    experiments = c("E1"),
    treatments = list(E1 = 1L),
    obs_table = data.frame(),
    exe = exe
  )
  key <- .eval_cache_key(cache, list(P1 = 1), c("E1"))
  result <- structure(list(
    score = 0.25,
    loglik = -0.125,
    residuals = data.frame(exp_id = "E1", treatment = 1L, user_var = "HWAM",
                           resid = 1, stringsAsFactors = FALSE),
    per_var = list(HWAM = list(n = 1L, RMSE = 1)),
    per_exp_var = data.frame(exp_id = "E1", user_var = "HWAM", n = 1L,
                             stringsAsFactors = FALSE)
  ), class = "objective_result")

  expect_true(.eval_cache_put(cache, key, result))
  loaded <- .eval_cache_get(cache, key)
  expect_s3_class(loaded, "objective_result")
  expect_equal(loaded$score, result$score)
  expect_equal(loaded$residuals$resid, 1)

  result2 <- result
  result2$score <- 0.125
  expect_true(.eval_cache_put(cache, key, result2))
  loaded2 <- .eval_cache_get(cache, key)
  expect_equal(loaded2$score, result2$score)
})

test_that("persistent evaluation cache dedupes and reuses objective results", {
  runner_env <- environment(evaluate_thetas)
  skip_if(bindingIsLocked("run_many", runner_env),
          "run_many binding is locked; exercised by source-based parity runner")

  tmp <- tempfile("eval_cache_")
  dir.create(tmp, recursive = TRUE)
  exe <- file.path(tmp, "dscsm048")
  writeLines("fake exe", exe)

  cfg <- list(
    calibrator = list(name = "cache_test", workdir = file.path(tmp, "work"),
                      dssat_dir = tmp, dssat_exe = exe, num_cores = 1L,
                      cache_evaluations = TRUE),
    source = list(hemp_dir = tmp),
    experiments = list("E1"),
    parameters = list(genetic_cultivar = list(
      P1 = list(active = TRUE, min = 0, max = 10, start = 5)
    )),
    crops = list(list(code = "HM")),
    engine = list(scalar_outputs = list(HWAM = "HWAM")),
    objective = list(weighting = "unified", weights = list(), error_model = list())
  )
  space <- parameter_space_from_config(cfg)
  setup <- list(
    space = space,
    crop = list(code = "HM", filex_ext = "HMX", genotype_stem = "HMGRO048"),
    exe = exe,
    specs = space$specs,
    run_root = file.path(tmp, "runs"),
    obs = observations(data.frame()),
    experiments = c("E1"),
    treatments = list(E1 = 1L)
  )

  old_run_many <- get("run_many", envir = runner_env)
  calls <- 0L
  fake_run_many <- function(jobs, n_workers, on_done = NULL, warmup = 0) {
    calls <<- calls + length(jobs)
    lapply(jobs, function(job) {
      theta <- job$theta
      spawn_result(
        "success",
        tmp,
        theta,
        evaluate = data.frame(
          treatment = 1L,
          variable = "HWAM",
          sim = 100 + as.numeric(theta$P1),
          meas = 100,
          stringsAsFactors = FALSE
        )
      )
    })
  }
  assign("run_many", fake_run_many, envir = runner_env)
  on.exit(assign("run_many", old_run_many, envir = runner_env), add = TRUE)

  first <- evaluate_thetas(cfg, list(list(P1 = 1), list(P1 = 1), list(P1 = 2)),
                           setup = setup, n_workers = 1L)$results
  expect_equal(calls, 2L)
  expect_equal(first[[1]]$score, first[[2]]$score)

  second <- evaluate_thetas(cfg, list(list(P1 = 1), list(P1 = 2)),
                            setup = setup, n_workers = 1L)$results
  expect_equal(calls, 2L)
  expect_equal(second[[1]]$score, first[[1]]$score)
  expect_equal(second[[2]]$score, first[[3]]$score)

  cache_files <- list.files(file.path(tmp, "work", "cache_test", "evaluation_cache"),
                            pattern = "\\.rds$", recursive = TRUE)
  expect_gt(length(cache_files), 0)
})

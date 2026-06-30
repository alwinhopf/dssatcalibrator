test_that("run_many uses Windows PSOCK workers when requested", {
  skip_if_not(.Platform$OS.type == "windows")
  skip_if_not(requireNamespace("parallel", quietly = TRUE))

  runner_env <- environment(run_many)
  skip_if(bindingIsLocked("spawn_and_run", runner_env),
          "spawn_and_run binding is locked; exercised by source-based parity runner")

  old_spawn <- get("spawn_and_run", envir = runner_env)
  fake_spawn <- function(theta, ..., sleep = 0.25) {
    Sys.sleep(sleep)
    spawn_result("success", ".", theta,
                 plantgro = data.frame(pid = Sys.getpid()),
                 evaluate = data.frame())
  }
  assign("spawn_and_run", fake_spawn, envir = runner_env)
  on.exit(assign("spawn_and_run", old_spawn, envir = runner_env), add = TRUE)

  jobs <- lapply(seq_len(8), function(i) list(theta = list(i = i), sleep = 0.25))

  t1 <- system.time(res1 <- run_many(jobs, n_workers = 1L))[["elapsed"]]
  t4 <- system.time(res4 <- run_many(jobs, n_workers = 4L))[["elapsed"]]

  expect_equal(vapply(res4, function(r) r$theta$i, integer(1)), seq_len(8))
  expect_true(all(vapply(res4, function(r) r$status, character(1)) == "success"))
  expect_gt(length(unique(vapply(res4, function(r) r$plantgro$pid[[1]], integer(1)))), 1)
  expect_lt(t4, t1)
})

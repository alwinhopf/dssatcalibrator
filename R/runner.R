# Parallel execution of spawns. R twin of python/dssatcalibrator/runner.py.
#
# The unit of work is one spawn_and_run() call (an isolated DSSAT subprocess).
# A serial warm-up pass runs the first `warmup` jobs before the remainder are
# fanned across workers, mirroring the Python runner. Results are returned in
# job order. Parallelism uses parallel::mclapply where available (fork), falling
# back to serial lapply on platforms without fork (e.g. Windows) — the result
# ordering and contents are identical either way.

#' How many workers to use. `num_cores > 0` taken as-is; `0` (or less) leaves 2
#' logical cores free. Mirrors runner.py:resolve_cores.
#' @export
resolve_cores <- function(num_cores) {
  if (!is.null(num_cores) && num_cores > 0) return(as.integer(num_cores))
  max(1L, (parallel::detectCores() %||% 2L) - 2L)
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

#' Run a list of spawn jobs and return their results in job order.
#'
#' Each job is a list forwarded as arguments to [spawn_and_run] (must contain
#' `theta` plus the spawn arguments). A spawn that errors is captured as an
#' "error" result so one bad run never kills the batch. `on_done` (if given) is
#' called as each job finishes. `warmup` runs the first N jobs serially first.
#' Mirrors runner.py:run_many.
#' @export
run_many <- function(jobs, n_workers, on_done = NULL, warmup = 0) {
  n <- length(jobs)
  results <- vector("list", n)
  n_workers <- max(1L, as.integer(n_workers))

  run_one <- function(i) {
    job <- jobs[[i]]
    theta <- job$theta
    args <- job[setdiff(names(job), "theta")]
    tryCatch(
      do.call(spawn_and_run, c(list(theta = theta), args)),
      error = function(e) spawn_result("error", ".", theta %||% list(),
                                       message = paste0("<", conditionMessage(e), ">"))
    )
  }

  warm <- min(max(warmup, 0L), n)
  if (warm > 0) {
    for (i in seq_len(warm)) {
      results[[i]] <- run_one(i)
      if (!is.null(on_done)) on_done(results[[i]])
    }
  }

  rest <- if (warm < n) (warm + 1L):n else integer(0)
  if (length(rest) > 0) {
    use_fork <- n_workers > 1L && .Platform$OS.type != "windows" &&
      requireNamespace("parallel", quietly = TRUE)
    rest_results <- if (use_fork) {
      parallel::mclapply(rest, run_one, mc.cores = n_workers)
    } else {
      lapply(rest, run_one)
    }
    for (k in seq_along(rest)) {
      results[[rest[k]]] <- rest_results[[k]]
      if (!is.null(on_done)) on_done(results[[rest[k]]])
    }
  }
  results
}

# Parallel execution of spawns. R twin of python/dssatcalibrator/runner.py.
#
# The unit of work is one spawn_and_run() call (an isolated DSSAT subprocess).
# A serial warm-up pass runs the first `warmup` jobs before the remainder are
# fanned across workers, mirroring the Python runner. Results are returned in
# job order. Parallelism uses parallel::mclapply on fork-capable platforms and a
# PSOCK cluster on Windows, so the R interface can fan out DSSAT subprocesses on
# every supported OS. The result ordering and contents are identical either way.

#' How many workers to use. `num_cores > 0` taken as-is; `0` (or less) leaves 2
#' logical cores free. Mirrors runner.py:resolve_cores.
#' @export
resolve_cores <- function(num_cores) {
  if (!is.null(num_cores) && num_cores > 0) return(as.integer(num_cores))
  max(1L, (parallel::detectCores() %||% 2L) - 2L)
}

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 1L && is.atomic(a) && is.na(a)) return(b)
  a
}

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

  error_result <- function(theta, e) {
    spawn_result("error", ".", theta %||% list(),
                 message = paste0("<", conditionMessage(e), ">"))
  }

  run_one <- function(job) {
    theta <- job$theta
    args <- job[setdiff(names(job), "theta")]
    tryCatch(
      do.call(spawn_and_run, c(list(theta = theta), args)),
      error = function(e) error_result(theta, e)
    )
  }

  warm <- min(max(warmup, 0L), n)
  if (warm > 0) {
    for (i in seq_len(warm)) {
      results[[i]] <- run_one(jobs[[i]])
      if (!is.null(on_done)) on_done(results[[i]])
    }
  }

  rest <- if (warm < n) (warm + 1L):n else integer(0)
  if (length(rest) > 0) {
    worker_count <- min(n_workers, length(rest))
    use_parallel <- worker_count > 1L && requireNamespace("parallel", quietly = TRUE)
    rest_jobs <- jobs[rest]
    rest_results <- if (!use_parallel) {
      lapply(rest_jobs, run_one)
    } else if (.Platform$OS.type != "windows") {
      parallel::mclapply(rest_jobs, run_one, mc.cores = worker_count)
    } else {
      cl <- parallel::makeCluster(worker_count, type = "PSOCK")
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterCall(cl, setwd, getwd())

      # PSOCK workers start clean R sessions. Export the package/source
      # environment so spawn_and_run and its helpers are available whether this
      # code is installed as a package or sourced by run_r_tests.R.
      pkg_env <- parent.env(environment())
      export_names <- ls(envir = pkg_env, all.names = TRUE)
      export_names <- export_names[vapply(export_names, function(nm) {
        obj <- get(nm, envir = pkg_env)
        is.function(obj) || is.atomic(obj) || is.list(obj)
      }, logical(1))]
      parallel::clusterExport(cl, export_names, envir = pkg_env)

      pkg_name <- environmentName(pkg_env)
      if (nzchar(pkg_name) && pkg_name %in% loadedNamespaces()) {
        parallel::clusterEvalQ(cl, suppressPackageStartupMessages(
          require(pkg_name, character.only = TRUE)
        ))
      }

      run_one_psock <- function(job) {
        theta <- job$theta
        args <- job[setdiff(names(job), "theta")]
        tryCatch(
          do.call(get("spawn_and_run", envir = .GlobalEnv),
                  c(list(theta = theta), args)),
          error = function(e) {
            structure(
              list(status = "error", run_dir = ".", theta = theta %||% list(),
                   plantgro = data.frame(), evaluate = data.frame(),
                   message = paste0("<", conditionMessage(e), ">")),
              class = "spawn_result"
            )
          }
        )
      }
      parallel::parLapplyLB(cl, rest_jobs, run_one_psock)
    }
    for (k in seq_along(rest)) {
      results[[rest[k]]] <- rest_results[[k]]
      if (!is.null(on_done)) on_done(results[[rest[k]]])
    }
  }
  results
}

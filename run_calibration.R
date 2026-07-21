#!/usr/bin/env Rscript
# Command-line entry point for a DSSAT calibration run (R twin of run_calibration.py).
#
#   Rscript run_calibration.R config_hemp.yaml --n 300
#   Rscript run_calibration.R config_hemp.yaml --preset A --n-particles 250
#   Rscript run_calibration.R config_hemp.yaml --validate --cv-scheme year
#   Rscript run_calibration.R config_hemp.yaml --nowcast 2021-07-15 --forecast
#
# Writes figures + CSV summaries to results/<name>/ (or --outdir).

.script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  getwd()
}

.load_dssatcalibrator <- function() {
  r_dir <- file.path(.script_dir(), "R")
  if (dir.exists(r_dir)) {
    r_files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
    invisible(lapply(r_files, function(f) sys.source(f, envir = globalenv())))
  } else {
    suppressMessages(library(dssatcalibrator))
  }
}

.load_dssatcalibrator()

# --- tiny flag parser: --key value, or --flag (boolean) ---------------------
.parse_args <- function(argv) {
  out <- list(`_positional` = character(0))
  bool_flags <- c("validate", "assimilate", "combined", "forecast", "diagnostics",
                  "no-progress")
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (startsWith(a, "--")) {
      key <- substring(a, 3)
      if (key %in% bool_flags) { out[[key]] <- TRUE; i <- i + 1L }
      else { out[[key]] <- argv[i + 1L]; i <- i + 2L }
    } else { out[["_positional"]] <- c(out[["_positional"]], a); i <- i + 1L }
  }
  out
}

args <- .parse_args(commandArgs(trailingOnly = TRUE))
if (length(args[["_positional"]]) < 1) stop("usage: run_calibration.R <config.yaml> [options]")
cfg <- load_config(args[["_positional"]][1])

setk <- function(cfg, path, val) { cfg[[path[1]]][[path[2]]] <- val; cfg }
if (!is.null(args$n))               cfg$method$sample$n <- as.integer(args$n)
if (!is.null(args$engine))          cfg$method$sample$engine <- args$engine
if (!is.null(args$preset))          cfg$method$preset <- args$preset
if (!is.null(args[["bayesian-engine"]])) cfg$method$bayesian$engine <- args[["bayesian-engine"]]
if (!is.null(args$optimizer)) { cfg$method$optimizer$engine <- args$optimizer; cfg$method$bayesian$engine <- "none" }
if (!is.null(args$sensitivity)) { cfg$method$sensitivity <- modifyList(cfg$method$sensitivity %||% list(), list(engine = args$sensitivity, active = TRUE, auto_activate = TRUE)) }
if (!is.null(args$select))          cfg$method$select <- modifyList(cfg$method$select %||% list(), list(engine = paste0("stepwise_", args$select), active = TRUE))
if (!is.null(args$surrogate))       cfg$method$surrogate <- modifyList(cfg$method$surrogate %||% list(), list(engine = args$surrogate, active = TRUE))
if (!is.null(args[["n-particles"]])) cfg$method$bayesian$n_particles <- as.integer(args[["n-particles"]])
if (!is.null(args$experiments))     cfg$experiments <- as.list(strsplit(args$experiments, ",")[[1]])
if (!is.null(args[["assim-mode"]])) cfg$assimilation$mode <- args[["assim-mode"]]
if (isTRUE(args$forecast))          cfg$forecast$active <- TRUE

progress <- !isTRUE(args[["no-progress"]])
name <- cfg$calibrator$name
outdir <- args$outdir %||% file.path(.cfg_get(cfg$calibrator, "results_dir", "results"), name)
figdir <- file.path(.cfg_get(cfg$calibrator, "figures_dir", "figures"), name)

if (!is.null(args$nowcast)) {
  res <- nowcast(cfg, args$nowcast, progress = progress)
  cat(sprintf("\nNowcast as of %s -> %s\n", res$as_of, normalizePath(outdir, mustWork = FALSE)))
  for (k in names(res$best_theta)) cat(sprintf("  %-8s %.4f\n", k, res$best_theta[[k]]))
} else if (isTRUE(args$combined)) {
  res <- combined_mode(cfg, progress = progress)
  cat("\n=== Calibrated (base) parameters ===\n")
  for (k in names(res$calibration$best_theta)) cat(sprintf("  %-8s %.4f\n", k, res$calibration$best_theta[[k]]))
} else if (isTRUE(args$assimilate) || isTRUE(.cfg_get(.cfg_get(cfg, "assimilation", list()), "active", FALSE))) {
  res <- assimilate(cfg, progress = progress)
  cat(sprintf("\nAssimilation (%s) complete.\n", res$mode))
} else if (isTRUE(args$validate)) {
  scheme <- args[["cv-scheme"]] %||% .cfg_get(.cfg_get(cfg$method, "validation", list()), "scheme", "loeo")
  df <- validate_cv(cfg, scheme = scheme, progress = progress)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, file.path(outdir, sprintf("validation_%s.csv", scheme)), row.names = FALSE)
  cat(sprintf("\nCross-validation (%s) written.\n", scheme))
} else {
  result <- if (!is.null(args$combine)) combine_runs(cfg, strsplit(args$combine, ",")[[1]]) else calibrate(cfg, progress = progress)
  best_spawns <- spawn_results_for(cfg, result$best_theta, result$experiments)
  paths <- make_report(result, outdir, best_spawns = best_spawns, figdir = figdir)
  cat("\n=== Best-fit parameters ===\n")
  for (k in names(result$best_theta)) cat(sprintf("  %-8s %.4f\n", k, result$best_theta[[k]]))
  cat("\n=== Fit summary (best) ===\n"); print(summary_fit_table(result))
  if (isTRUE(args$diagnostics) || isTRUE(.cfg_get(.cfg_get(cfg, "diagnostics", list()), "active", FALSE))) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(identifiability(result), file.path(outdir, "identifiability.csv"), row.names = FALSE)
    utils::write.csv(structural_adequacy(result), file.path(outdir, "structural_adequacy.csv"), row.names = FALSE)
  }
  if (isTRUE(.cfg_get(.cfg_get(cfg, "forecast", list()), "active", FALSE))) {
    for (v in unlist(.cfg_get(.cfg_get(cfg, "forecast", list()), "variables", list("LAID")))) forecast_lai(cfg, result, variable = v)
  }
  cat(sprintf("\nData tables -> %s\nFigures     -> %s\n", normalizePath(outdir, mustWork = FALSE), normalizePath(figdir, mustWork = FALSE)))
}

#!/usr/bin/env Rscript

# Manual integration test for pooled calibration against a local DSSAT install.
# It uses the existing hemp example experiments in C:/DSSAT48/Hemp and does not
# modify the DSSAT installation; all genotype edits happen in per-spawn run dirs.

source_package <- function() {
  r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
  invisible(lapply(r_files, function(f) sys.source(f, envir = globalenv())))
}

make_pooled_cfg <- function(num_cores = 4L, workdir = "results/_actual_pooled_r_eval",
                            keep_run_dirs = FALSE,
                            experiments = c("UFCI2101", "UFCI2201", "UFJA2101", "UFJA2201")) {
  cfg <- load_config("config_hemp.yaml", validate = FALSE)
  cfg$calibrator$dssat_exe <- "C:/DSSAT48/DSCSM048.EXE"
  cfg$calibrator$dssat_dir <- "C:/DSSAT48"
  cfg$calibrator$workdir <- workdir
  cfg$calibrator$cache_spawns <- FALSE
  cfg$calibrator$keep_run_dirs <- keep_run_dirs
  cfg$calibrator$num_cores <- as.integer(num_cores)
  cfg$source$hemp_dir <- "C:/DSSAT48/Hemp"
  cfg$experiments <- as.list(experiments)
  cfg$gating$species <- "free"

  for (group in names(cfg$parameters)) {
    for (name in names(cfg$parameters[[group]])) {
      spec <- cfg$parameters[[group]][[name]]
      if (is.list(spec)) {
        spec$active <- FALSE
        spec$scope <- NULL
        spec$pooling <- NULL
        cfg$parameters[[group]][[name]] <- spec
      }
    }
  }

  if (is.null(cfg$parameters$genetic_species)) cfg$parameters$genetic_species <- list()
  cfg$parameters$genetic_species$SLWREF <- list(
    active = TRUE, min = 0.0038, max = 0.0052, start = 0.0046,
    scope = "global", spe_key = "SLWREF,SLWSLO,NSLOPE,LNREF,PGREF"
  )
  cfg$parameters$genetic_cultivar$CSDL$active <- TRUE
  cfg$parameters$genetic_cultivar$CSDL$scope <- "global"
  cfg$parameters$genetic_cultivar[["EM-FL"]]$active <- TRUE
  cfg$parameters$genetic_cultivar[["EM-FL"]]$scope <- "experiment"
  cfg$parameters$genetic_ecotype[["PL-EM"]]$active <- TRUE
  cfg$parameters$genetic_ecotype[["PL-EM"]]$scope <- "experiment"
  cfg
}

make_samples <- function(space) {
  samples <- sample_design(space, n = 4L, engine = "montecarlo", seed = 730L, include_start = TRUE)
  samples[2, "EM-FL__UFCI2101"] <- 28.0
  samples[2, "EM-FL__UFCI2201"] <- 32.0
  samples[2, "EM-FL__UFJA2101"] <- 36.0
  samples[2, "EM-FL__UFJA2201"] <- 40.0
  samples[2, "PL-EM__UFCI2101"] <- 3.5
  samples[2, "PL-EM__UFCI2201"] <- 4.5
  samples[2, "PL-EM__UFJA2101"] <- 5.5
  samples[2, "PL-EM__UFJA2201"] <- 6.5
  samples[2, "SLWREF"] <- 0.0048
  samples
}

stop_if_missing_dssat <- function() {
  if (!file.exists("C:/DSSAT48/DSCSM048.EXE")) stop("DSSAT executable not found: C:/DSSAT48/DSCSM048.EXE")
  missing <- character(0)
  for (exp_id in c("UFCI2101", "UFCI2201", "UFJA2101", "UFJA2201")) {
    for (ext in c("HMX", "HMA", "HMT")) {
      p <- file.path("C:/DSSAT48/Hemp", sprintf("%s.%s", exp_id, ext))
      if (!file.exists(p)) missing <- c(missing, p)
    }
  }
  if (length(missing)) stop(sprintf("Missing DSSAT hemp example files:\n%s", paste(missing, collapse = "\n")))
}

source_package()
stop_if_missing_dssat()

debug_dir <- "results/_actual_pooled_r_debug"
if (dir.exists(debug_dir)) unlink(debug_dir, recursive = TRUE, force = TRUE)
debug_cfg <- make_pooled_cfg(1L, debug_dir, keep_run_dirs = TRUE, experiments = "UFCI2101")
debug_space <- parameter_space_from_config(debug_cfg)
debug_theta <- as.list(debug_space$start)
names(debug_theta) <- debug_space$names
debug_theta[["EM-FL__UFCI2101"]] <- 28.0
debug_theta[["PL-EM__UFCI2101"]] <- 3.5
debug_theta[["SLWREF"]] <- 0.0048

debug_res <- spawn_and_run(
  theta = debug_theta, exp_id = "UFCI2101", cfg = debug_cfg, crop = crop_for(debug_cfg, "HM"),
  param_specs = debug_space$specs, run_root = debug_dir, treatments = 1L,
  exe = debug_cfg$calibrator$dssat_exe
)
cat(sprintf("DIRECT_SPAWN status=%s plantgro=%d evaluate=%d run_dir=%s\n",
            debug_res$status, nrow(debug_res$plantgro), nrow(debug_res$evaluate), debug_res$run_dir))
if (!identical(debug_res$status, "success")) stop(debug_res$message)
spe_path <- file.path(debug_res$run_dir, "HMGRO048.SPE")
spe_line <- grep("SLWREF,SLWSLO,NSLOPE,LNREF,PGREF", readLines(spe_path, warn = FALSE), value = TRUE)
spe_line <- spe_line[!startsWith(trimws(spe_line, which = "left"), "!")][1]
cat(sprintf("DIRECT_SPE_LINE %s\n", spe_line))
if (!startsWith(spe_line, " .0048 .0004")) stop("R species writer did not preserve .SPE leading-dot decimal")

eval_dir <- "results/_actual_pooled_r_eval"
if (dir.exists(eval_dir)) unlink(eval_dir, recursive = TRUE, force = TRUE)
cfg <- make_pooled_cfg(4L, eval_dir)
space <- parameter_space_from_config(cfg)
samples <- make_samples(space)
t0 <- proc.time()[["elapsed"]]
ed <- evaluate_design(cfg, samples, progress = TRUE)
elapsed <- proc.time()[["elapsed"]] - t0
design <- ed$design

cat(sprintf("EVAL space_names=%s\n", paste(space$names, collapse = ",")))
cat(sprintf("EVAL design_shape=%dx%d samples_shape=%dx%d elapsed_sec=%.3f\n",
            nrow(design), ncol(design), nrow(samples), ncol(samples), elapsed))
cat(sprintf("EVAL n_obs=%s\n", paste(as.integer(design$n_obs), collapse = ",")))
cat(sprintf("EVAL scores=%s\n", paste(sprintf("%.6f", as.numeric(design$score)), collapse = ",")))
cat(sprintf("EVAL best_sample_id=%d best_score=%.6f\n",
            as.integer(design$sample_id[which.min(design$score)]), as.numeric(min(design$score))))

if (any(!is.finite(design$score))) stop("Non-finite score in R pooled DSSAT evaluation")
if (any(as.integer(design$n_obs) <= 0L)) stop("No matched observations in R pooled DSSAT evaluation")

parallel_samples <- samples[1:3, , drop = FALSE]
parallel_runs <- list()
for (cores in c(1L, 4L)) {
  pdir <- sprintf("results/_actual_pooled_r_parallel_%d", cores)
  if (dir.exists(pdir)) unlink(pdir, recursive = TRUE, force = TRUE)
  pcfg <- make_pooled_cfg(cores, pdir)
  pt0 <- proc.time()[["elapsed"]]
  ped <- evaluate_design(pcfg, parallel_samples, progress = FALSE)
  pelapsed <- proc.time()[["elapsed"]] - pt0
  parallel_runs[[as.character(cores)]] <- list(
    elapsed = pelapsed,
    scores = as.numeric(ped$design$score),
    n_obs = as.integer(ped$design$n_obs)
  )
  cat(sprintf("PARALLEL cores=%d elapsed_sec=%.3f scores=%s n_obs=%s\n",
              cores, pelapsed,
              paste(sprintf("%.6f", parallel_runs[[as.character(cores)]]$scores), collapse = ","),
              paste(parallel_runs[[as.character(cores)]]$n_obs, collapse = ",")))
}
if (!isTRUE(all.equal(parallel_runs[["1"]]$scores, parallel_runs[["4"]]$scores, tolerance = 1e-9))) {
  stop("R pooled DSSAT scores differ between 1 and 4 cores")
}
cat(sprintf("PARALLEL speedup_4_vs_1=%.3f\n", parallel_runs[["1"]]$elapsed / parallel_runs[["4"]]$elapsed))
cat("R_ACTUAL_POOLED_DSSAT_TEST_OK\n")

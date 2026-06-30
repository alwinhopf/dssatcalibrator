# Real-DSSAT parameter impact atlas.
#
# The atlas implementation is Python-first because it already collects broad
# DSSAT *.OUT tables and writes analysis artifacts. This R front end keeps the
# workflow available from R while returning ordinary data frames.

.impact_default_outdir <- function(cfg) {
  cal <- .cfg_get(cfg, "calibrator", list())
  name <- .cfg_get(cal, "name", "run")
  results_dir <- .cfg_get(cal, "results_dir", "results")
  file.path(results_dir, paste0(name, "_impact_atlas"))
}

.impact_script_path <- function(script = NULL, config_path = NULL) {
  if (!is.null(script) && nzchar(script)) {
    return(normalizePath(script, mustWork = FALSE))
  }
  env_script <- Sys.getenv("DSSATCAL_IMPACT_SCRIPT", unset = "")
  if (nzchar(env_script)) {
    return(normalizePath(env_script, mustWork = FALSE))
  }
  candidates <- c(
    file.path(getwd(), "run_impact_atlas.py"),
    if (!is.null(config_path)) file.path(dirname(normalizePath(config_path, mustWork = FALSE)), "run_impact_atlas.py") else character(0)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Could not find run_impact_atlas.py. Pass script = '<path>/run_impact_atlas.py'.")
  }
  normalizePath(hit[1], mustWork = TRUE)
}

.impact_add_values <- function(args, flag, values) {
  if (is.null(values) || length(values) == 0) return(args)
  c(args, flag, as.character(values))
}

.impact_atlas_args <- function(script, config_path, output_dir = NULL,
                               experiments = NULL, groups = NULL, levels = NULL,
                               active_only = FALSE, discover_cultivar = FALSE,
                               discover_ecotype = FALSE, discover_species = FALSE,
                               discover_genotype = FALSE,
                               allow_species = FALSE,
                               max_parameters = NULL, max_per_group = NULL,
                               output_files = NULL, num_cores = NULL,
                               dssat_exe = NULL, dssat_dir = NULL,
                               hemp_dir = NULL, keep_existing = FALSE,
                               write_long = TRUE, compress_long = FALSE,
                               effect_tolerance = NULL, progress = TRUE) {
  args <- c(script, config_path)
  if (!is.null(output_dir)) args <- c(args, "--outdir", output_dir)
  args <- .impact_add_values(args, "--experiments", experiments)
  args <- .impact_add_values(args, "--groups", groups)
  args <- .impact_add_values(args, "--levels", levels)
  if (isTRUE(active_only)) args <- c(args, "--active-only")
  if (isTRUE(discover_cultivar)) args <- c(args, "--discover-cultivar")
  if (isTRUE(discover_ecotype)) args <- c(args, "--discover-ecotype")
  if (isTRUE(discover_species)) args <- c(args, "--discover-species")
  if (isTRUE(discover_genotype)) args <- c(args, "--discover-genotype")
  if (isTRUE(allow_species)) args <- c(args, "--allow-species")
  if (!is.null(max_parameters)) args <- c(args, "--max-parameters", as.character(as.integer(max_parameters)))
  if (!is.null(max_per_group)) args <- c(args, "--max-per-group", as.character(as.integer(max_per_group)))
  args <- .impact_add_values(args, "--outputs", output_files)
  if (!is.null(num_cores)) args <- c(args, "--cores", as.character(as.integer(num_cores)))
  if (!is.null(dssat_exe)) args <- c(args, "--dssat-exe", dssat_exe)
  if (!is.null(dssat_dir)) args <- c(args, "--dssat-dir", dssat_dir)
  if (!is.null(hemp_dir)) args <- c(args, "--hemp-dir", hemp_dir)
  if (isTRUE(keep_existing)) args <- c(args, "--keep-existing")
  if (!isTRUE(write_long)) args <- c(args, "--no-long")
  if (isTRUE(compress_long)) args <- c(args, "--compress-long")
  if (!is.null(effect_tolerance)) args <- c(args, "--effect-tolerance", as.character(as.numeric(effect_tolerance)))
  if (!isTRUE(progress)) args <- c(args, "--no-progress")
  args
}

.impact_read_csv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Run a real-DSSAT one-at-a-time parameter impact atlas.
#'
#' This is the R front end for `run_impact_atlas.py`. It runs real DSSAT spawns,
#' writes the same artifact set as the Python CLI, and returns the generated
#' tables as data frames.
#'
#' @param config Path to a YAML config, or a config list that can be written as
#'   YAML.
#' @param output_dir Output directory. Defaults to
#'   `<results_dir>/<calibrator_name>_impact_atlas`.
#' @param experiments Optional experiment IDs.
#' @param groups Optional parameter groups to sweep.
#' @param levels Parameter levels to run; usually `c("low", "high")`.
#' @param active_only Sweep only config parameters with `active: true`.
#' @param discover_cultivar Discover cultivar coefficients from the crop `.CUL`
#'   file.
#' @param discover_ecotype Discover ecotype coefficients from the crop `.ECO`
#'   file.
#' @param discover_species Discover numeric tokens from the crop `.SPE` file.
#' @param discover_genotype Discover cultivar, ecotype, and species genotype
#'   parameters.
#' @param allow_species Explicitly allow `.SPE` species edits for species sweeps.
#' @param max_parameters Optional total parameter limit for smoke tests.
#' @param max_per_group Optional per-group parameter limit for smoke tests.
#' @param output_files Optional DSSAT `*.OUT` files to collect.
#' @param num_cores Override `calibrator.num_cores`.
#' @param dssat_exe,dssat_dir,hemp_dir Optional path overrides.
#' @param keep_existing Keep an existing output directory.
#' @param write_long Write the large `outputs_long.csv` table.
#' @param compress_long Write `outputs_long.csv.gz`.
#' @param effect_tolerance Absolute delta threshold for counting changed output
#'   summaries.
#' @param python Python executable.
#' @param script Path to `run_impact_atlas.py`.
#' @param progress Show Python progress output.
#' @return An `impact_atlas_result` list with manifests, raw effects, compact
#'   score/output/parameter summaries, `capability_map`, `parameter_catalog`,
#'   and optionally `output_long`.
#' @export
run_impact_atlas <- function(config, output_dir = NULL,
                             experiments = NULL, groups = NULL,
                             levels = c("low", "high"),
                             active_only = FALSE,
                             discover_cultivar = FALSE,
                             discover_ecotype = FALSE,
                             discover_species = FALSE,
                             discover_genotype = FALSE,
                             allow_species = FALSE,
                             max_parameters = NULL, max_per_group = NULL,
                             output_files = NULL, num_cores = NULL,
                             dssat_exe = NULL, dssat_dir = NULL,
                             hemp_dir = NULL, keep_existing = FALSE,
                             write_long = TRUE, compress_long = FALSE,
                             effect_tolerance = NULL,
                             python = Sys.getenv("DSSATCAL_PYTHON", "python"),
                             script = NULL, progress = TRUE) {
  cfg_obj <- list()
  if (is.list(config)) {
    cfg_obj <- config
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("run_impact_atlas() requires the 'yaml' package when config is a list.")
    }
    config_path <- tempfile(fileext = ".yaml")
    yaml::write_yaml(config[setdiff(names(config), "_config_path")], config_path)
  } else {
    config_path <- as.character(config)
    if (requireNamespace("yaml", quietly = TRUE) && file.exists(config_path)) {
      cfg_obj <- yaml::read_yaml(config_path)
      if (is.null(cfg_obj)) cfg_obj <- list()
    }
  }
  if (is.null(output_dir)) output_dir <- .impact_default_outdir(cfg_obj)
  script <- .impact_script_path(script, config_path)
  args <- .impact_atlas_args(
    script, config_path, output_dir = output_dir, experiments = experiments,
    groups = groups, levels = levels, active_only = active_only,
    discover_cultivar = discover_cultivar,
    discover_ecotype = discover_ecotype,
    discover_species = discover_species,
    discover_genotype = discover_genotype,
    allow_species = allow_species,
    max_parameters = max_parameters,
    max_per_group = max_per_group, output_files = output_files,
    num_cores = num_cores, dssat_exe = dssat_exe, dssat_dir = dssat_dir,
    hemp_dir = hemp_dir, keep_existing = keep_existing, write_long = write_long,
    compress_long = compress_long, effect_tolerance = effect_tolerance,
    progress = progress
  )

  out <- system2(python, args = args,
                 stdout = if (isTRUE(progress)) "" else TRUE,
                 stderr = if (isTRUE(progress)) "" else TRUE)
  status <- if (is.character(out)) attr(out, "status") %||% 0L else as.integer(out)
  if (!is.null(status) && status != 0L) {
    msg <- if (is.character(out) && length(out)) paste(utils::tail(out, 20), collapse = "\n") else ""
    stop(sprintf("Impact atlas command failed with status %s.\n%s", status, msg))
  }

  long_name <- if (isTRUE(compress_long)) "outputs_long.csv.gz" else "outputs_long.csv"
  result <- list(
    output_dir = normalizePath(output_dir, mustWork = FALSE),
    run_manifest = .impact_read_csv(file.path(output_dir, "run_manifest.csv")),
    file_manifest = .impact_read_csv(file.path(output_dir, "file_manifest.csv")),
    output_effects = .impact_read_csv(file.path(output_dir, "parameter_output_effects.csv")),
    score_effects = .impact_read_csv(file.path(output_dir, "score_effects.csv")),
    output_impact_summary = .impact_read_csv(file.path(output_dir, "output_impact_summary.csv")),
    parameter_impact_summary = .impact_read_csv(file.path(output_dir, "parameter_impact_summary.csv")),
    capability_map = .impact_read_csv(file.path(output_dir, "capability_map.csv")),
    parameter_catalog = .impact_read_csv(file.path(output_dir, "parameter_catalog.csv")),
    output_long = if (isTRUE(write_long)) .impact_read_csv(file.path(output_dir, long_name)) else data.frame(),
    command = list(python = python, args = args)
  )
  class(result) <- "impact_atlas_result"
  result
}

#' @export
print.impact_atlas_result <- function(x, ...) {
  n <- nrow(x$run_manifest)
  ok <- if (n) sum(x$run_manifest$status %in% c("success", "cached")) else 0L
  cat(sprintf("Impact atlas result: %s/%s successful spawns\n", ok, n))
  cat(sprintf("Output directory: %s\n", x$output_dir))
  invisible(x)
}

# Configuration loading: environment > YAML > built-in default.
#
# R twin of python/dssatcalibrator/config.py. A run is described by one YAML
# file (see config_hemp.yaml). This module loads it, merges defaults, applies a
# handful of environment overrides, and offers small helpers to enumerate the
# *active* parameters and resolve paths.
#
# Design choice (shared with Python): the config is kept as plain nested lists
# (not bespoke S4/R6 classes) so it round-trips to YAML/JSON unchanged and the
# manifest of a run is just the list.
#
# Seeded from DSSAT_Calibration/config.R (dc_load_config / dc_validate_config).

# Minimal built-in defaults; the YAML supplies the rest. Mirrors DEFAULTS in
# config.py key-for-key so that deleting a YAML key reproduces built-in
# behaviour identically in both languages.
.dssatcal_defaults <- function() {
  list(
    calibrator = list(
      name = "run",
      seed = 42L,
      workdir = "results/_workdir",
      results_dir = "results",
      figures_dir = "figures",
      dssat_exe = "",
      dssat_dir = "C:/DSSAT48",
      num_cores = 0L,
      cache_spawns = TRUE,
      cache_evaluations = TRUE,
      evaluation_cache_dir = "",
      evaluation_cache_salt = "",
      keep_run_dirs = FALSE
    ),
    method = list(
      preset = "C",
      sample = list(engine = "lhs", n = 200L),
      validation = list(scheme = "none")
    ),
    objective = list(weighting = "unified", weights = list(), error_model = list(),
                     likelihood = list(type = "gaussian"),
                     model_discrepancy = list()),
    parameters = list(),
    crops = list(),
    experiments = list(),
    execution = list(backend = "native"),
    templates = list(template_dir = ""),
    gating = list(cultivar = "free", ecotype = "gated", species = "blocked"),
    management_options = list(use_source_planting_date = FALSE),
    weather = list(provider = "file", gap_fill = "none", horizon = 0L,
                   cache_dir = "weather_cache"),
    soil = list(provider = "file", source = "ssurgo", cache_dir = "soil_cache"),
    forecast = list(active = FALSE, variables = list("LAID"), n_ensemble = 0L,
                    anchor_continuity = TRUE, decay_days = 21L),
    diagnostics = list(active = FALSE),
    observation_sources = list(),
    fusion = list(conflict_resolution = "keep_all", source_priority = list()),
    assimilation = list(
      active = FALSE,
      mode = "recalibration",
      allow_uncoupled = FALSE,
      recalibration = list(
        engine = "glue",
        recal_sample_size = 100L,
        warm_start = TRUE,
        update_frequency = "on_observation"
      ),
      enkf = list(n_ensemble = 50L, inflation = 1.05,
                  state_variables = list("LAID", "CWAD")),
      forcing = list(min_confidence = 0.8, smoothing = TRUE)
    ),
    sparse = list(
      delta_from_analog = list(active = FALSE),
      hierarchical_priors = list(active = FALSE),
      trait_priors = list(active = FALSE),
      identifiability_gate = list(active = FALSE),
      observation_design = list(active = FALSE)
    )
  )
}

# Environment overrides: ENV name -> c(section, key)
.dssatcal_env_overrides <- list(
  DSSATCAL_DSSAT_EXE   = c("calibrator", "dssat_exe"),
  DSSATCAL_DSSAT_DIR   = c("calibrator", "dssat_dir"),
  DSSATCAL_NUM_CORES   = c("calibrator", "num_cores"),
  DSSATCAL_WORKDIR     = c("calibrator", "workdir"),
  DSSATCAL_RESULTS_DIR = c("calibrator", "results_dir"),
  DSSAT_TEMPLATE_DIR   = c("templates", "template_dir")
)

# Recursive deep-merge of named lists; `over` wins on scalars, dicts merge.
# Mirrors config.py:_deep_merge.
.deep_merge <- function(base, over) {
  out <- base
  for (k in names(over)) {
    v <- over[[k]]
    if (is.list(v) && !is.null(out[[k]]) && is.list(out[[k]]) &&
        !is.null(names(v)) && !is.null(names(out[[k]]))) {
      out[[k]] <- .deep_merge(out[[k]], v)
    } else {
      out[[k]] <- v
    }
  }
  out
}

#' Load a config YAML, merge defaults, apply env overrides, and validate.
#'
#' @param path Path to a YAML configuration file.
#' @param validate Logical; run [validate_config()] before returning (default
#'   TRUE). Set FALSE to load a config that is still being assembled.
#' @return A nested list (the merged config), with `_config_path` recorded.
#' @export
load_config <- function(path, validate = TRUE) {
  user <- yaml::read_yaml(path)
  if (is.null(user)) user <- list()
  cfg <- .deep_merge(.dssatcal_defaults(), user)

  for (env in names(.dssatcal_env_overrides)) {
    sk <- .dssatcal_env_overrides[[env]]
    val <- Sys.getenv(env, unset = "")
    if (nzchar(val)) {
      sec <- sk[1]; key <- sk[2]
      if (is.null(cfg[[sec]])) cfg[[sec]] <- list()
      cfg[[sec]][[key]] <- if (key == "num_cores") as.integer(val) else val
    }
  }

  cfg[["_config_path"]] <- normalizePath(path, mustWork = FALSE)
  if (validate) validate_config(cfg)
  cfg
}

#' Validate a merged config, stopping with ALL problems listed.
#'
#' R twin of config.py:validate_config. Catches unknown engine/preset/weighting
#' vocabulary, inverted or non-numeric parameter bounds, start values outside
#' their bounds, unknown prior distributions, and an empty active-parameter set.
#' Returns `cfg` invisibly on success.
#' @export
validate_config <- function(cfg) {
  # Allowed vocabularies — mirror config.py. Kept local so sourcing all R/ files
  # into one environment can't clash with engine constants (e.g. orchestrator's
  # `.PRESETS` stage map).
  PRESETS            <- c("A", "B", "C", "D")
  WEIGHTING_MODES    <- c("unified", "sigma", "user", "count_scale", "agmip_wls")
  CV_SCHEMES         <- c("none", "loeo", "year", "site", "random")
  PRIOR_DISTS        <- c("uniform", "normal", "lognormal", "triangular")
  GATING_LEVELS      <- c("free", "gated", "blocked")
  EXECUTION_BACKENDS <- c("native", "dssatengine")
  ASSIMILATION_MODES <- c("recalibration", "enkf", "forcing")
  BAYES_ENGINES      <- c("glue", "smc_pf", "mcmc", "dream", "es_mda", "bayesopt",
                          "abc_smc", "history", "none", "")
  OPTIMIZER_ENGINES  <- c("nelder_mead", "neldermead", "nm", "diffevo", "de",
                          "cmaes", "cma_es", "cma", "none", "")
  PARAMETER_SCOPES   <- c("global", "shared", "pooled", "pool",
                          "experiment", "experiments", "per_experiment", "per-experiment",
                          "experiment_specific", "experiment-specific", "local")

  errors <- character(0)
  is_num <- function(x) is.numeric(x) && length(x) == 1L

  preset <- toupper(as.character(.cfg_get(cfg$method, "preset", "C")))
  if (!preset %in% PRESETS)
    errors <- c(errors, sprintf("method.preset '%s' is not one of %s.",
                                preset, paste(PRESETS, collapse = ", ")))

  scheme <- tolower(as.character(.cfg_get(.cfg_get(cfg$method, "validation", list()), "scheme", "none")))
  if (!scheme %in% CV_SCHEMES)
    errors <- c(errors, sprintf("method.validation.scheme '%s' is not one of %s.",
                                scheme, paste(CV_SCHEMES, collapse = ", ")))

  weighting <- tolower(as.character(.cfg_get(cfg$objective, "weighting", "unified")))
  if (!weighting %in% WEIGHTING_MODES)
    errors <- c(errors, sprintf("objective.weighting '%s' is not one of %s.",
                                weighting, paste(WEIGHTING_MODES, collapse = ", ")))

  backend <- tolower(as.character(.cfg_get(cfg$execution, "backend", "native")))
  if (!backend %in% EXECUTION_BACKENDS)
    errors <- c(errors, sprintf("execution.backend '%s' is not one of %s.",
                                backend, paste(EXECUTION_BACKENDS, collapse = ", ")))

  mode <- tolower(as.character(.cfg_get(cfg$assimilation, "mode", "recalibration")))
  if (!mode %in% ASSIMILATION_MODES)
    errors <- c(errors, sprintf("assimilation.mode '%s' is not one of %s.",
                                mode, paste(ASSIMILATION_MODES, collapse = ", ")))

  be <- tolower(as.character(.cfg_get(.cfg_get(cfg$method, "bayesian", list()), "engine", "glue")))
  if (!be %in% BAYES_ENGINES)
    errors <- c(errors, sprintf("method.bayesian.engine '%s' is not one of %s.",
                                be, paste(setdiff(BAYES_ENGINES, ""), collapse = ", ")))
  oe <- tolower(as.character(.cfg_get(.cfg_get(cfg$method, "optimizer", list()), "engine", "none")))
  if (!oe %in% OPTIMIZER_ENGINES)
    errors <- c(errors, sprintf("method.optimizer.engine '%s' is not one of %s.",
                                oe, paste(setdiff(OPTIMIZER_ENGINES, ""), collapse = ", ")))

  for (lvl in c("cultivar", "ecotype", "species")) {
    g <- tolower(as.character(.cfg_get(cfg$gating, lvl, "free")))
    if (!g %in% GATING_LEVELS)
      errors <- c(errors, sprintf("gating.%s '%s' is not one of %s.",
                                  lvl, g, paste(GATING_LEVELS, collapse = ", ")))
  }

  cores <- .cfg_get(cfg$calibrator, "num_cores", 0L)
  if (!(is_num(cores) && cores >= 0))
    errors <- c(errors, sprintf("calibrator.num_cores must be an integer >= 0 (got %s).",
                                format(cores)))

  n_active <- 0L
  params_block <- cfg$parameters
  if (!is.null(params_block)) {
    for (group in names(params_block)) {
      params <- params_block[[group]]
      if (!is.list(params) || is.null(names(params))) next
      for (name in names(params)) {
        spec <- params[[name]]
        tag <- sprintf("parameters.%s.%s", group, name)
        if (!is.list(spec)) {
          errors <- c(errors, sprintf("%s must be a mapping.", tag)); next
        }
        lo <- spec$min; hi <- spec$max
        if (!is_num(lo) || !is_num(hi)) {
          errors <- c(errors, sprintf("%s: min/max must both be numeric.", tag))
        } else if (lo >= hi) {
          errors <- c(errors, sprintf("%s: min (%s) must be < max (%s).",
                                      tag, format(lo), format(hi)))
        }
        start <- spec$start
        if (!is.null(start) && is_num(start) && is_num(lo) && is_num(hi) &&
            !(lo <= start && start <= hi))
          errors <- c(errors, sprintf("%s: start (%s) is outside [min=%s, max=%s].",
                                      tag, format(start), format(lo), format(hi)))
        prior <- spec$prior
        if (is.list(prior)) {
          dist <- tolower(as.character(.cfg_get(prior, "dist", "uniform")))
          if (!dist %in% PRIOR_DISTS)
            errors <- c(errors, sprintf("%s: prior.dist '%s' is not one of %s.",
                                        tag, dist, paste(PRIOR_DISTS, collapse = ", ")))
        }
        scope <- tolower(as.character(.cfg_get(spec, "scope", .cfg_get(spec, "pooling", "global"))))
        if (!scope %in% PARAMETER_SCOPES) {
          errors <- c(errors, sprintf("%s: scope '%s' is not one of %s.",
                                      tag, scope, paste(PARAMETER_SCOPES, collapse = ", ")))
        }
        if (scope %in% c("experiment", "experiments", "per_experiment", "per-experiment",
                         "experiment_specific", "experiment-specific", "local") &&
            length(.cfg_get(cfg, "experiments", list())) == 0) {
          errors <- c(errors, sprintf("%s: scope '%s' requires at least one configured experiment.",
                                      tag, scope))
        }
        if (isTRUE(spec$active)) n_active <- n_active + 1L
      }
    }
  }
  if (n_active == 0L)
    errors <- c(errors, paste0("No active parameters: at least one parameter must ",
                               "have `active: true` to calibrate."))

  if (length(errors)) {
    stop(sprintf("Invalid configuration (%d problem%s):\n  - %s",
                 length(errors), if (length(errors) != 1L) "s" else "",
                 paste(errors, collapse = "\n  - ")), call. = FALSE)
  }
  invisible(cfg)
}

#' Flatten the `parameters` block to the list of ACTIVE parameter specs.
#'
#' Each returned spec carries `group`, `name`, `min`, `max`, `start`, plus any
#' extra keys (`prior`, `role`, `type`, `dssat`, ...). Mirrors
#' config.py:active_parameters.
#' @export
active_parameters <- function(cfg) {
  out <- list()
  params_block <- cfg$parameters
  if (is.null(params_block)) return(out)
  for (group in names(params_block)) {
    params <- params_block[[group]]
    if (!is.list(params) || is.null(names(params))) next
    for (name in names(params)) {
      spec <- params[[name]]
      if (!is.list(spec) || !isTRUE(spec$active)) next
      rec <- c(list(group = group, name = name), spec)
      out[[length(out) + 1L]] <- rec
    }
  }
  out
}

#' Every declared parameter (active or not). Mirrors config.py:all_parameters.
#' @export
all_parameters <- function(cfg) {
  out <- list()
  params_block <- cfg$parameters
  if (is.null(params_block)) return(out)
  for (group in names(params_block)) {
    params <- params_block[[group]]
    if (!is.list(params) || is.null(names(params))) next
    for (name in names(params)) {
      spec <- params[[name]]
      if (is.list(spec)) {
        out[[length(out) + 1L]] <- c(list(group = group, name = name), spec)
      }
    }
  }
  out
}

#' Return the crop block matching a 2-letter DSSAT code (first crop if one).
#' Mirrors config.py:crop_for.
#' @export
crop_for <- function(cfg, code) {
  crops <- cfg$crops
  if (is.null(crops) || length(crops) == 0) return(list())
  for (c in crops) {
    if (!is.null(c$code) && c$code == code) return(c)
  }
  crops[[1]]
}

#' Resolve the DSSAT executable path. Mirrors config.py:resolve_exe.
#' @export
resolve_exe <- function(cfg) {
  exe <- cfg$calibrator$dssat_exe
  if (!is.null(exe) && nzchar(exe)) return(exe)
  file.path(cfg$calibrator$dssat_dir, "DSCSM048.EXE")
}

#' Resolve the DSSAT48 install layout used by every spawn.
#' Mirrors config.py:resolve_dssat_paths.
#' @export
resolve_dssat_paths <- function(cfg) {
  root <- cfg$calibrator$dssat_dir
  list(
    root = root,
    exe = resolve_exe(cfg),
    genotype = file.path(root, "Genotype"),
    weather = file.path(root, "Weather"),
    soil = file.path(root, "Soil")
  )
}

# Sibling DSSAT_Gridded_Run_Tutorial/dssat_templates relative to this package.
.workspace_template_dir <- function() {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  workspaces <- unique(c(dirname(cwd), cwd, dirname(dirname(cwd))))
  for (workspace in workspaces) {
    candidate <- file.path(workspace, "DSSAT_Gridded_Run_Tutorial", "dssat_templates")
    if (dir.exists(candidate)) return(candidate)
  }
  file.path(workspaces[1], "DSSAT_Gridded_Run_Tutorial", "dssat_templates")
}

#' Resolve the shared `dssat_templates` directory.
#'
#' Precedence: DSSAT_TEMPLATE_DIR / templates.template_dir / top-level
#' template_dir / the sibling gridded tutorial template folder. Returns NULL
#' when none configured/discoverable unless `required = TRUE`. Mirrors
#' config.py:resolve_template_dir.
#' @export
resolve_template_dir <- function(cfg = NULL, required = FALSE) {
  if (is.null(cfg)) cfg <- list()
  configured <- Sys.getenv("DSSAT_TEMPLATE_DIR", unset = "")
  if (!nzchar(configured)) {
    configured <- if (!is.null(cfg$templates$template_dir)) cfg$templates$template_dir else ""
  }
  if (!nzchar(configured) && !is.null(cfg$template_dir)) configured <- cfg$template_dir
  if (nzchar(configured)) return(configured)

  default <- .workspace_template_dir()
  if (dir.exists(default)) return(default)
  if (required) {
    stop("Shared DSSAT template directory not found. Set DSSAT_TEMPLATE_DIR ",
         "or templates.template_dir.")
  }
  NULL
}

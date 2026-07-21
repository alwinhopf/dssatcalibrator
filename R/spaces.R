# Parameter search space built from the config's active parameters.
# R twin of python/dssatcalibrator/spaces.py (ParameterSpace dataclass).
#
# Represented as an S3 list of class "parameter_space" with fields
# names/low/high/start/specs, mirroring the Python dataclass.

#' Build a ParameterSpace from a config.
#'
#' @param cfg A merged config (see [load_config]).
#' @return An object of class `parameter_space` with `names`, `low`, `high`,
#'   `start` (numeric vectors) and `specs` (list of active specs).
#' @export
parameter_space_from_config <- function(cfg) {
  specs <- expand_parameter_specs(cfg, active_parameters(cfg))
  if (length(specs) == 0) {
    stop("No active parameters in config (set active: true on some).")
  }
  nm <- vapply(specs, function(s) s$name, character(1))
  low <- vapply(specs, function(s) as.numeric(s$min), numeric(1))
  high <- vapply(specs, function(s) as.numeric(s$max), numeric(1))
  start <- vapply(seq_along(specs), function(i) {
    s <- specs[[i]]
    if (!is.null(s$start)) as.numeric(s$start) else 0.5 * (low[i] + high[i])
  }, numeric(1))
  start <- pmin(pmax(start, low), high)
  structure(
    list(names = nm, low = low, high = high, start = start, specs = specs),
    class = "parameter_space"
  )
}

.GLOBAL_SCOPES <- c("", "global", "shared", "pooled", "pool")
.EXPERIMENT_SCOPES <- c("experiment", "experiments", "per_experiment", "per-experiment",
                        "experiment_specific", "experiment-specific", "local")
.CULTIVAR_SCOPES <- c("cultivar", "cultivars", "per_cultivar", "per-cultivar",
                      "cultivar_specific", "cultivar-specific")

.default_scope_for_group <- function(cfg, group) {
  defaults <- .cfg_get(cfg, "parameter_defaults", list())
  by_group <- .cfg_get(defaults, "scope_by_group", list())
  if (!is.null(by_group[[group]])) return(as.character(by_group[[group]]))
  as.character(.cfg_get(defaults, "scope", "global"))
}

.scope_of <- function(cfg, spec) {
  default <- .default_scope_for_group(cfg, .cfg_get(spec, "group", NULL))
  raw <- tolower(as.character(.cfg_get(spec, "scope", .cfg_get(spec, "pooling", default))))
  if (raw %in% .EXPERIMENT_SCOPES) return("experiment")
  if (raw %in% .CULTIVAR_SCOPES) return("cultivar")
  if (raw %in% .GLOBAL_SCOPES) return("global")
  stop(sprintf("Unknown parameter scope '%s' for %s.%s; use 'global', 'experiment', or 'cultivar'.",
               raw, .cfg_get(spec, "group", "?"), .cfg_get(spec, "name", "?")))
}

.scoped_name <- function(base, exp_id = NULL) {
  if (is.null(exp_id)) base else sprintf("%s__%s", base, exp_id)
}

.cultivars_for_scope <- function(cfg, spec) {
  explicit <- unlist(.cfg_get(spec, "cultivars",
                              .cfg_get(spec, "cultivar_anchors",
                                       .cfg_get(spec, "cultivar_codes", list()))))
  if (length(explicit)) return(unique(as.character(explicit[nzchar(explicit)])))
  cultivars <- character(0)
  for (crop in .cfg_get(cfg, "crops", list())) {
    calibrated <- unlist(.cfg_get(crop, "calibration_cultivars", list()))
    if (length(calibrated)) {
      cultivars <- c(cultivars, as.character(calibrated))
    } else {
      cultivars <- c(cultivars, as.character(unlist(.cfg_get(crop, "cultivar_anchors", list()))))
      anchor <- .cfg_get(crop, "cultivar_anchor", NULL)
      if (!is.null(anchor)) cultivars <- c(cultivars, as.character(anchor))
    }
  }
  unique(cultivars[nzchar(cultivars)])
}

.apply_cultivar_overrides <- function(rec, cultivar) {
  aliases <- list(start = c("start_by_cultivar", "starts_by_cultivar"),
                  min = c("min_by_cultivar", "mins_by_cultivar"),
                  max = c("max_by_cultivar", "maxs_by_cultivar"))
  for (field in names(aliases)) {
    for (key in aliases[[field]]) {
      values <- rec[[key]]
      if (!is.null(values[[cultivar]])) {
        rec[[field]] <- values[[cultivar]]
        break
      }
    }
  }
  rec
}

#' Expand active parameter specs into optimizer dimensions.
#'
#' A parameter is global by default. With `scope: experiment`, one optimizer
#' dimension is created per configured experiment, while `base_name` preserves
#' the original DSSAT coefficient name for the spawn writer.
#' @export
expand_parameter_specs <- function(cfg, specs) {
  experiments <- unlist(.cfg_get(cfg, "experiments", list()))
  out <- list()
  for (spec in specs) {
    base <- spec$name
    scope <- .scope_of(cfg, spec)
    if (scope == "experiment") {
      if (length(experiments) == 0) {
        stop(sprintf("Parameter %s.%s has scope=experiment but no experiments are configured.",
                     spec$group, base))
      }
      for (exp_id in experiments) {
        rec <- spec
        rec$base_name <- base
        rec$name <- .scoped_name(base, exp_id)
        rec$scope <- "experiment"
        rec$exp_id <- as.character(exp_id)
        out[[length(out) + 1L]] <- rec
      }
    } else if (scope == "cultivar") {
      cultivars <- .cultivars_for_scope(cfg, spec)
      if (!length(cultivars)) {
        stop(sprintf("Parameter %s.%s has scope=cultivar but no cultivars are configured.",
                     spec$group, base))
      }
      for (cultivar in cultivars) {
        rec <- spec
        rec$base_name <- base
        rec$name <- .scoped_name(base, cultivar)
        rec$scope <- "cultivar"
        rec$cultivar <- cultivar
        rec <- .apply_cultivar_overrides(rec, cultivar)
        out[[length(out) + 1L]] <- rec
      }
    } else {
      rec <- spec
      rec$base_name <- base
      rec$scope <- "global"
      out[[length(out) + 1L]] <- rec
    }
  }
  out
}

#' Number of dimensions (active parameters). Mirrors ParameterSpace.ndim.
#' @export
ps_ndim <- function(space) length(space$names)

#' Map a parameter vector (native units) to a named theta list.
#' Mirrors ParameterSpace.to_theta.
#' @export
ps_to_theta <- function(space, vector) {
  vals <- as.numeric(vector)
  setNames(as.list(vals), space$names)
}

#' Clip a parameter vector to the [low, high] box. Mirrors ParameterSpace.clip.
#' @export
ps_clip <- function(space, vector) {
  pmin(pmax(as.numeric(vector), space$low), space$high)
}

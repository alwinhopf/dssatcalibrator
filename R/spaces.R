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
  specs <- active_parameters(cfg)
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

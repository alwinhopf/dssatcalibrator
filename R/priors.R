# Prior distributions for parameters — sampling and log-density.
# R twin of python/dssatcalibrator/priors.py. Seeded from
# DSSAT_Calibration/priors.R (regularization / shrinkage idiom).
#
# Supported distributions (all truncated to [min, max]):
#   uniform (default) | normal | lognormal | triangular
#
# Parity note: log-density functions are deterministic and are checked exactly
# against the Python implementation. Sampling uses base-R RNG streams, which
# cannot match NumPy bit-for-bit, so sampling parity is verified statistically
# (matching moments) per the agreed policy.

.prior_bounds <- function(spec) c(as.numeric(spec$min), as.numeric(spec$max))

.prior_center <- function(spec) {
  b <- .prior_bounds(spec); lo <- b[1]; hi <- b[2]
  prior <- spec$prior
  if (is.null(prior)) prior <- list()
  if (!is.null(prior$mean)) return(as.numeric(prior$mean))
  if (!is.null(prior$mode)) return(as.numeric(prior$mode))
  if (!is.null(spec$start)) return(as.numeric(spec$start))
  0.5 * (lo + hi)
}

.prior_dist_name <- function(spec) {
  prior <- spec$prior
  if (is.null(prior) || is.null(prior$dist)) return("uniform")
  tolower(as.character(prior$dist))
}

# Inverse-CDF draw from a triangular(lo, mode, hi).
.rtriangular <- function(n, lo, mode, hi) {
  u <- runif(n)
  fc <- if (hi > lo) (mode - lo) / (hi - lo) else 0.5
  out <- numeric(n)
  left <- u < fc
  out[left]  <- lo + sqrt(u[left] * (hi - lo) * (mode - lo))
  out[!left] <- hi - sqrt((1 - u[!left]) * (hi - lo) * (hi - mode))
  out
}

#' Draw `n` samples for one parameter from its (truncated) prior.
#' Mirrors priors.py:sample_one. `rng` is accepted for signature parity and
#' ignored (base-R RNG is used; seed via set.seed()).
#' @export
sample_one <- function(spec, n, rng = NULL) {
  b <- .prior_bounds(spec); lo <- b[1]; hi <- b[2]
  dist <- .prior_dist_name(spec)
  prior <- spec$prior; if (is.null(prior)) prior <- list()

  if (dist == "uniform") {
    return(runif(n, lo, hi))
  }
  if (dist == "normal") {
    mu <- .prior_center(spec)
    sd <- if (!is.null(prior$sd)) as.numeric(prior$sd) else 0.25 * (hi - lo)
    # Truncated normal via inverse CDF on [lo, hi].
    plo <- pnorm(lo, mu, sd); phi <- pnorm(hi, mu, sd)
    u <- runif(n, plo, phi)
    return(qnorm(u, mu, sd))
  }
  if (dist == "lognormal") {
    center <- max(.prior_center(spec), 1e-9)
    sigma <- if (!is.null(prior$sd)) as.numeric(prior$sd) else 0.5
    draws <- rlnorm(n, meanlog = log(center), sdlog = sigma)
    return(pmin(pmax(draws, lo), hi))
  }
  if (dist == "triangular") {
    mode <- min(max(.prior_center(spec), lo), hi)
    return(.rtriangular(n, lo, mode, hi))
  }
  stop(sprintf("unknown prior dist '%s' for parameter '%s'", dist, spec$name))
}

#' Log prior density of one parameter value (-Inf outside the bounds).
#' Mirrors priors.py:log_prior_one (and matches scipy's truncnorm/lognorm/triang
#' parameterisations).
#' @export
log_prior_one <- function(spec, value) {
  b <- .prior_bounds(spec); lo <- b[1]; hi <- b[2]
  if (!(lo <= value && value <= hi)) return(-Inf)
  dist <- .prior_dist_name(spec)
  prior <- spec$prior; if (is.null(prior)) prior <- list()

  if (dist == "uniform") return(0.0)

  if (dist == "normal") {
    mu <- .prior_center(spec)
    sd <- if (!is.null(prior$sd)) as.numeric(prior$sd) else 0.25 * (hi - lo)
    xi <- (value - mu) / sd
    logZ <- log(pnorm(hi, mu, sd) - pnorm(lo, mu, sd))
    # scipy truncnorm.logpdf = norm.logpdf(xi) - log(scale) - log(Z)
    return(dnorm(xi, log = TRUE) - log(sd) - logZ)
  }
  if (dist == "lognormal") {
    center <- max(.prior_center(spec), 1e-9)
    sigma <- if (!is.null(prior$sd)) as.numeric(prior$sd) else 0.5
    return(dlnorm(value, meanlog = log(center), sdlog = sigma, log = TRUE))
  }
  if (dist == "triangular") {
    mode <- min(max(.prior_center(spec), lo), hi)
    width <- hi - lo
    if (width <= 0) return(-Inf)
    if (value <= mode) {
      dens <- if (mode > lo) 2 * (value - lo) / (width * (mode - lo)) else 0
    } else {
      dens <- if (hi > mode) 2 * (hi - value) / (width * (hi - mode)) else 0
    }
    return(if (dens > 0) log(dens) else -Inf)
  }
  stop(sprintf("unknown prior dist '%s' for parameter '%s'", dist, spec$name))
}

#' Draw `n` parameter sets, one column per active parameter.
#' Returns a data.frame in native units. Mirrors priors.py:sample_prior_design.
#' @export
sample_prior_design <- function(space, n, rng = NULL) {
  cols <- lapply(space$specs, function(s) sample_one(s, n, rng))
  names(cols) <- space$names
  as.data.frame(cols, check.names = FALSE)
}

#' Total log prior of a parameter vector = sum over independent parameters.
#' Mirrors priors.py:log_prior_vec.
#' @export
log_prior_vec <- function(space, theta) {
  total <- 0.0
  for (s in space$specs) {
    total <- total + log_prior_one(s, as.numeric(theta[[s$name]]))
    if (is.infinite(total) && total < 0) break
  }
  total
}

#' TRUE if any active parameter declares a non-uniform prior.
#' Mirrors priors.py:has_informative_prior.
#' @export
has_informative_prior <- function(space) {
  any(vapply(space$specs, function(s) .prior_dist_name(s) != "uniform", logical(1)))
}

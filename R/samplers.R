# Design-of-experiment samplers over a ParameterSpace.
# R twin of python/dssatcalibrator/samplers.py.
#
# LHS via the 'lhs' package, Sobol via 'randtoolbox', plus plain Monte-Carlo and
# a coarse grid. Returns a data.frame in native parameter units, one column per
# parameter. (RNG streams differ from NumPy/SciPy, so sampling parity is
# statistical per the agreed policy; the scaling/þstructure is identical.)

#' Draw `n` parameter sets from `space` using the chosen design engine.
#' engine: "lhs" | "sobol" | "montecarlo" | "grid". When `include_start`, the
#' configured start point is prepended as the first row. Mirrors samplers.py:sample.
#' @export
sample_design <- function(space, n, engine = "lhs", seed = 42, include_start = TRUE) {
  d <- ps_ndim(space)
  engine <- tolower(engine)
  set.seed(seed)

  unit <- switch(engine,
    "lhs" = {
      if (!requireNamespace("lhs", quietly = TRUE)) stop("engine 'lhs' needs the 'lhs' package.")
      lhs::randomLHS(n, d)
    },
    "sobol" = {
      if (!requireNamespace("randtoolbox", quietly = TRUE)) stop("engine 'sobol' needs 'randtoolbox'.")
      m <- randtoolbox::sobol(n, dim = d, scrambling = 1, seed = seed)
      if (is.null(dim(m))) matrix(m, ncol = d) else m
    },
    "montecarlo" = matrix(runif(n * d), nrow = n, ncol = d),
    "grid" = {
      per <- max(2L, as.integer(round(n^(1.0 / d))))
      axes <- replicate(d, seq(0, 1, length.out = per), simplify = FALSE)
      g <- as.matrix(expand.grid(axes))
      g[seq_len(min(n, nrow(g))), , drop = FALSE]
    },
    stop(sprintf("unknown sampler engine: %s", engine))
  )

  # scale unit hypercube to native units: low + unit*(high-low)
  scaled <- sweep(sweep(unit, 2, (space$high - space$low), `*`), 2, space$low, `+`)
  df <- as.data.frame(scaled)
  names(df) <- space$names

  if (include_start) {
    start_row <- as.data.frame(as.list(space$start)); names(start_row) <- space$names
    df <- rbind(start_row, df)
  }
  df
}

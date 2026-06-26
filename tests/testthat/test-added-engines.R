# Offline correctness tests for the R twins of the added engines (CMA-ES, DREAM,
# ES-MDA). A synthetic linear-Gaussian problem with a known optimum is solved
# through each engine's real code path (no DSSAT). Mirrors
# tests/test_added_engines.py. Bayesian optimisation needs DiceKriging and is
# covered in Python; here it is skipped unless the package is installed.

NAMES <- c("p0", "p1", "p2", "p3")
TARGET <- c(1.5, -2.0, 0.7, 3.1)

make_space <- function() {
  specs <- lapply(seq_along(NAMES), function(i)
    list(name = NAMES[i], min = -5, max = 5, start = 0, prior = list(dist = "uniform")))
  list(names = NAMES, low = rep(-5, 4), high = rep(5, 4), start = rep(0, 4), specs = specs)
}

make_scorer <- function(seed = 0, n_obs = 10, sigma = 0.2) {
  set.seed(seed)
  A <- matrix(rnorm(n_obs * length(NAMES)), n_obs, length(NAMES))
  obs_vec <- as.numeric(A %*% TARGET)
  function(thetas) lapply(thetas, function(t) {
    x <- as.numeric(unlist(t[NAMES]))
    sim <- as.numeric(A %*% x); resid <- sim - obs_vec
    df <- data.frame(exp_id = "E", treatment = 1L, dssat = paste0("V", seq_len(n_obs)),
                     user_var = "v", date = NA, obs = obs_vec, sim = sim,
                     sigma = sigma, weight = 1, resid = resid, kind = "scalar",
                     stringsAsFactors = FALSE)
    chi2 <- sum((resid / sigma)^2)
    list(score = chi2 / n_obs, loglik = -0.5 * chi2, residuals = df, per_var = list())
  })
}

err <- function(theta) sqrt(sum((as.numeric(unlist(theta[NAMES])) - TARGET)^2))
cfg_for <- function(engine, ...) list(calibrator = list(seed = 7, num_cores = 1),
                                      method = list(bayesian = c(list(engine = engine), list(...))))

test_that("CMA-ES (R) recovers the optimum", {
  sp <- make_space(); sr <- make_scorer()
  score_batch <- function(ths) vapply(sr(ths), function(r) r$score, numeric(1))
  res <- run_optimizer(sp, score_batch, method = "cmaes", seed = 1, maxiter = 60, popsize = 12)
  expect_lt(err(res$best_theta), 0.15)
})

test_that("DREAM (R) recovers the optimum and yields a posterior cloud", {
  r <- run_dream(cfg_for("dream", n_generations = 120, burn_in = 60),
                 make_scorer(), make_space(), progress = FALSE)
  expect_lt(err(r$best_theta), 0.5)
  expect_gt(nrow(r$design), 10)
})

test_that("ES-MDA (R) recovers the optimum", {
  r <- run_es_mda(cfg_for("es_mda", ensemble_size = 40, iterations = 6),
                  make_scorer(), make_space(), progress = FALSE)
  expect_lt(err(r$best_theta), 0.4)
})

test_that("the added engines resolve through the registry", {
  expect_equal(.resolve_estimator(list(bayesian = list(engine = "dream"))), "dream")
  expect_equal(.resolve_estimator(list(bayesian = list(engine = "es_mda"))), "es_mda")
  expect_equal(.resolve_estimator(list(bayesian = list(engine = "bayesopt"))), "bayesopt")
  expect_equal(.resolve_estimator(list(bayesian = list(engine = "none"),
                                       optimizer = list(engine = "cmaes"))), "optimizer")
})

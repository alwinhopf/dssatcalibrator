# Offline correctness tests for the R twins of the added engines (CMA-ES, DREAM,
# ES-MDA, history matching, ABC-SMC). A synthetic linear-Gaussian problem with a known optimum is solved
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

all_invalid_scorer <- function(thetas) lapply(thetas, function(t) list(
  score = Inf, loglik = -Inf, residuals = data.frame(), per_var = list()
))

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

test_that("history matching (R) keeps a finite NROY region", {
  r <- run_history_matching(cfg_for("history", n = 96, waves = 2, implausibility_cutoff = 4.0),
                            make_scorer(), make_space(), progress = FALSE)
  expect_gt(nrow(r$design), 20)
  expect_gt(r$ess, 0)
  expect_true(is.finite(r$best$score))
  expect_lt(err(r$best_theta), 4.0)
})

test_that("ABC-SMC (R) accepts a final population", {
  r <- run_abc_smc(cfg_for("abc_smc", n_particles = 40, waves = 2, oversample = 2,
                           threshold_quantile = 0.4),
                   make_scorer(), make_space(), progress = FALSE)
  expect_gt(nrow(r$design), 40)
  expect_gt(r$ess, 0)
  expect_equal(length(r$thresholds), 2)
  expect_lte(tail(r$thresholds, 1), r$thresholds[1])
  expect_lt(err(r$best_theta), 3.0)
})

test_that("sparse Bayesian engines reject all-invalid candidates", {
  expect_error(
    run_history_matching(cfg_for("history", n = 8, waves = 2),
                         all_invalid_scorer, make_space(), progress = FALSE),
    "History matching found no valid candidates"
  )
  expect_error(
    run_abc_smc(cfg_for("abc_smc", n_particles = 8, waves = 2, oversample = 1),
                all_invalid_scorer, make_space(), progress = FALSE),
    "ABC-SMC found no valid candidates"
  )
})

test_that("the added engines resolve through the registry", {
  expect_equal(.resolve_estimator(list(bayesian = list(engine = "dream"))), "dream")
  expect_equal(.resolve_estimator(list(bayesian = list(engine = "es_mda"))), "es_mda")
  expect_equal(.resolve_estimator(list(bayesian = list(engine = "bayesopt"))), "bayesopt")
  expect_equal(.resolve_estimator(list(bayesian = list(engine = "history"))), "history")
  expect_equal(.resolve_estimator(list(bayesian = list(engine = "abc_smc"))), "abc_smc")
  expect_equal(.resolve_estimator(list(bayesian = list(engine = "none"),
                                       optimizer = list(engine = "cmaes"))), "optimizer")
})

test_that("Sobol sensitivity (R) ranks the influential parameter", {
  sp <- list(
    names = c("a", "b"), low = c(0, 0), high = c(10, 10), start = c(5, 5),
    specs = list(list(name = "a", min = 0, max = 10, start = 5),
                 list(name = "b", min = 0, max = 10, start = 5))
  )
  scorer <- function(thetas) lapply(thetas, function(theta)
    list(score = 5 * theta$a + 0.001 * theta$b))
  first <- run_sensitivity(sp, scorer, method = "sobol", n_base = 128, seed = 7)
  second <- run_sensitivity(sp, scorer, method = "sobol", n_base = 128, seed = 7)
  expect_equal(first$ranking, second$ranking)
  expect_equal(first$ranking$parameter[1], "a")
  expect_equal(first$n_eval, 128 * (2 * length(sp$names)))
})

test_that("DEoptim differential evolution executes the real optional backend", {
  skip_if_not_installed("DEoptim")
  sp <- make_space()
  score_batch <- function(ths) vapply(ths, function(t) sum((unlist(t[NAMES]) - TARGET)^2), numeric(1))
  res <- run_optimizer(sp, score_batch, method = "diffevo", seed = 4,
                       maxiter = 25, popsize = 10)
  expect_true(is.finite(res$best_score))
  expect_lt(err(res$best_theta), 1.0)
})

test_that("mco NSGA-II executes and returns a finite Pareto population", {
  skip_if_not_installed("mco")
  sp <- make_space()
  evaluate <- function(ths) lapply(ths, function(t) {
    x <- unlist(t[NAMES])
    list(first = sum((x - TARGET)^2), second = sum((x + TARGET)^2))
  })
  res <- run_nsga2(evaluate, sp, c("first", "second"), pop_size = 8,
                   n_gen = 3, seed = 5)
  expect_equal(ncol(res$F), 2)
  expect_true(nrow(res$F) > 0 && all(is.finite(res$F)))
})

test_that("GP and RF surrogate backends execute with finite designs", {
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("ranger")
  sp <- make_space(); scorer <- make_scorer()
  for (backend in c("gp", "rf")) {
    cfg <- list(calibrator = list(seed = 8),
                method = list(surrogate = list(engine = backend, n_train = 16,
                                                n_candidates = 40, top_k = 3)))
    res <- run_surrogate(cfg, sp, scorer, progress = FALSE)
    expect_equal(res$info$model, backend)
    expect_true(nrow(res$design) >= 7 && any(is.finite(res$design$score)))
  }
})

test_that("DiceKriging Bayesian optimisation executes", {
  skip_if_not_installed("DiceKriging")
  sp <- make_space(); scorer <- make_scorer()
  cfg <- cfg_for("bayesopt", n_init = 10, n_iter = 1, batch_size = 2)
  res <- run_bayesopt(cfg, scorer, sp, progress = FALSE)
  expect_true(is.finite(res$best$score))
  expect_equal(nrow(res$design), 12)
})

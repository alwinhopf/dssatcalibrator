# Pluggable calibration engines. R twin of python/dssatcalibrator/engines/*.py.
#
# Deterministic engines (GLUE post-processing, AgMIP selection, optimizer wrap,
# Morris design math) port exactly; the stochastic engines (MCMC, SMC-PF,
# NSGA-II, surrogate, EnKF) are faithful ports validated statistically, since
# RNG streams and third-party solver internals cannot match Python bit-for-bit.

# ---- GLUE / Monte-Carlo (preset C default) — DETERMINISTIC ----------------

#' Turn an evaluated design into GLUE weights + a behavioural set.
#' Mirrors engines/glue.py:run_glue. `design` is a data.frame with the parameter
#' columns plus `score` and `loglik`.
#' @export
run_glue <- function(design, param_names, cfg, space = NULL) {
  d <- design
  rownames(d) <- NULL
  ll <- as.numeric(d$loglik)

  if (!is.null(space) && has_informative_prior(space)) {
    lp <- vapply(seq_len(nrow(d)), function(i) {
      theta <- setNames(as.list(as.numeric(d[i, param_names])), param_names)
      log_prior_vec(space, theta)
    }, numeric(1))
    ll <- ll + lp
  }

  finite <- is.finite(ll)
  w <- numeric(nrow(d))
  if (any(finite)) {
    m <- max(ll[finite])
    w[finite] <- exp(ll[finite] - m)
    total <- sum(w)
    if (total > 0) w <- w / total
  }
  d$weight <- w
  ess <- if (sum(w) > 0) 1.0 / sum(w^2) else 0.0

  q <- as.numeric(.cfg_get(.cfg_get(.cfg_get(cfg, "method", list()), "bayesian", list()),
                           "behavioural_quantile", 0.1))
  score_for_min <- ifelse(is.finite(d$score), d$score, Inf)
  valid_scores <- d$score[is.finite(d$score)]
  threshold <- if (length(valid_scores)) as.numeric(quantile(valid_scores, q, names = FALSE, type = 7)) else Inf
  behavioural <- d[is.finite(d$score) & d$score <= threshold, , drop = FALSE]

  best_id0 <- if (length(valid_scores)) (which.min(score_for_min) - 1L) else 0L
  best_theta <- setNames(as.list(as.numeric(d[best_id0 + 1L, param_names])), param_names)

  structure(list(design = d, behavioural = behavioural, best_theta = best_theta,
                 best_sample_id = best_id0, threshold = threshold, ess = ess),
            class = "glue_result")
}

#' Weighted posterior mean / sd / quantiles per parameter (for reporting).
#' Mirrors engines/glue.py:posterior_summary.
#' @export
posterior_summary <- function(glue, param_names) {
  d <- glue$design; w <- d$weight
  wquantile <- function(x, probs) {
    if (sum(w) <= 0) return(as.numeric(quantile(x, probs, names = FALSE, type = 7)))
    o <- order(x); xs <- x[o]; ws <- w[o] / sum(w[o])
    approx(cumsum(ws), xs, xout = probs, rule = 2, ties = "ordered")$y
  }
  rows <- lapply(param_names, function(n) {
    x <- as.numeric(d[[n]])
    if (sum(w) > 0) {
      mean_ <- sum(w * x); var_ <- sum(w * (x - mean_)^2); sd_ <- sqrt(max(var_, 0))
    } else { mean_ <- mean(x); sd_ <- sd(x) }
    qs <- wquantile(x, c(0.05, 0.95))
    data.frame(parameter = n, best = glue$best_theta[[n]], post_mean = mean_, post_sd = sd_,
               p05 = qs[1], p95 = qs[2],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---- classical optimisers (preset B) --------------------------------------

.OPT_FAIL <- 1e12

#' Minimise `score_batch` over `space`; return the best parameter set.
#' Mirrors engines/optimizers.py:run_optimizer. nelder_mead -> stats::optim;
#' diffevo -> DEoptim.
#' @export
run_optimizer <- function(space, score_batch, method = "diffevo", seed = 42,
                          maxiter = NULL, popsize = 15, restarts = 1, tol = 1e-4,
                          progress = FALSE) {
  st <- new.env()
  st$n <- 0L; st$best_theta <- NULL; st$best_score <- Inf; st$hist <- list()
  note <- function(theta, score) {
    if (score < st$best_score) { st$best_score <- score; st$best_theta <- theta }
    st$hist[[length(st$hist) + 1L]] <- list(iter = st$n, score = st$best_score)
  }
  cost_one <- function(x) {
    theta <- ps_to_theta(space, ps_clip(space, x))
    s <- score_batch(list(theta))[[1]]
    s <- if (is.finite(s)) as.numeric(s) else .OPT_FAIL
    st$n <- st$n + 1L; note(theta, s); s
  }
  method <- tolower(method)
  set.seed(seed)
  if (method %in% c("nelder_mead", "neldermead", "nm")) {
    starts <- list(as.numeric(space$start))
    for (i in seq_len(max(0L, restarts - 1L))) starts[[length(starts) + 1L]] <- runif(ps_ndim(space), space$low, space$high)
    for (x0 in starts) {
      stats::optim(x0, cost_one, method = "Nelder-Mead",
                   control = list(maxit = maxiter %||% (200L * ps_ndim(space)),
                                  reltol = tol))
    }
  } else if (method %in% c("diffevo", "differential_evolution", "de")) {
    if (!requireNamespace("DEoptim", quietly = TRUE)) stop("diffevo needs the 'DEoptim' package.")
    cost <- function(x) cost_one(x)
    DEoptim::DEoptim(cost, lower = space$low, upper = space$high,
                     control = DEoptim::DEoptim.control(
                       itermax = maxiter %||% 30L, NP = popsize * ps_ndim(space),
                       trace = FALSE))
  } else if (method %in% c("cmaes", "cma_es", "cma")) {
    score_pop <- function(thetas) {
      sc <- vapply(score_batch(thetas), function(s) if (is.finite(s)) as.numeric(s) else .OPT_FAIL,
                   numeric(1))
      st$n <- st$n + length(thetas)
      jb <- which.min(sc); note(thetas[[jb]], sc[jb]); sc
    }
    .cma_es_r(space, score_pop, maxiter = maxiter,
              popsize = if (!is.null(popsize) && popsize > 4) popsize else NULL)
  } else stop(sprintf("unknown optimizer method '%s' (use nelder_mead | diffevo | cmaes)", method))

  if (is.null(st$best_theta)) { st$best_theta <- ps_to_theta(space, space$start); st$best_score <- Inf }
  structure(list(best_theta = st$best_theta, best_score = st$best_score, method = method,
                 n_eval = st$n, history = st$hist), class = "optimizer_result")
}

# Minimal CMA-ES on the unit cube, mapped to bounds. Mirrors optimizers.py:_cma_es.
.cma_es_r <- function(space, score_pop, seed = NULL, maxiter = NULL, popsize = NULL) {
  n <- ps_ndim(space); low <- space$low; high <- space$high
  span <- ifelse(high > low, high - low, 1)
  mean_ <- pmin(pmax((space$start - low) / span, 0), 1); sigma <- 0.3
  lam <- if (!is.null(popsize)) as.integer(popsize) else 4L + as.integer(3 * log(n))
  mu <- lam %/% 2L
  w <- log(mu + 0.5) - log(seq_len(mu)); w <- w / sum(w)
  mueff <- 1 / sum(w^2)
  cc <- (4 + mueff / n) / (n + 4 + 2 * mueff / n)
  cs <- (mueff + 2) / (n + mueff + 5)
  c1 <- 2 / ((n + 1.3)^2 + mueff)
  cmu <- min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((n + 2)^2 + mueff))
  damps <- 1 + 2 * max(0, sqrt((mueff - 1) / (n + 1)) - 1) + cs
  chiN <- sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n^2))
  pc <- numeric(n); ps_ <- numeric(n); B <- diag(n); D <- rep(1, n); C <- diag(n)
  n_gen <- if (!is.null(maxiter)) as.integer(maxiter) else 100L + 50L * n
  to_theta <- function(xu) ps_to_theta(space, low + pmin(pmax(xu, 0), 1) * span)

  for (gen in seq_len(n_gen)) {
    z <- matrix(rnorm(lam * n), lam, n)
    y <- z %*% t(B %*% diag(D, n))
    x <- pmin(pmax(matrix(mean_, lam, n, byrow = TRUE) + sigma * y, 0), 1)
    scores <- score_pop(lapply(seq_len(lam), function(k) to_theta(x[k, ])))
    ord <- order(scores); xsel <- x[ord[seq_len(mu)], , drop = FALSE]
    old_mean <- mean_; mean_ <- as.numeric(w %*% xsel)
    y_w <- (mean_ - old_mean) / sigma
    Cinv <- B %*% diag(1 / D, n) %*% t(B)
    ps_ <- (1 - cs) * ps_ + sqrt(cs * (2 - cs) * mueff) * as.numeric(Cinv %*% y_w)
    hsig <- (sqrt(sum(ps_^2)) / sqrt(1 - (1 - cs)^(2 * gen)) / chiN) < (1.4 + 2 / (n + 1))
    pc <- (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mueff) * y_w
    yk <- sweep(xsel, 2, old_mean) / sigma
    C <- (1 - c1 - cmu) * C + c1 * (outer(pc, pc) + (1 - hsig) * cc * (2 - cc) * C) +
      cmu * (t(yk * w) %*% yk)
    sigma <- sigma * exp((cs / damps) * (sqrt(sum(ps_^2)) / chiN - 1))
    sigma <- min(max(sigma, 1e-4), 1.0)
    C <- (C + t(C)) / 2
    eg <- tryCatch(eigen(C, symmetric = TRUE), error = function(e) NULL)
    if (is.null(eg)) { B <- diag(n); D <- rep(1, n); C <- diag(n) }
    else { D <- sqrt(pmax(eg$values, 1e-14)); B <- eg$vectors }
    if (sigma < 1e-3 && max(D) * sigma < 1e-3) break
  }
  invisible(NULL)
}

# ---- AgMIP stepwise selection ---------------------------------------------

.ic_criterion <- function(name, loglik, k, n) {
  if (!is.finite(loglik)) return(Inf)
  if (name == "aicc") {
    aic <- 2 * k - 2 * loglik; denom <- max(n - k - 1, 1)
    return(aic + (2 * k * (k + 1)) / denom)
  }
  k * log(max(n, 1)) - 2 * loglik
}

.subspace <- function(space_full, names) {
  keep <- Filter(function(s) s$name %in% names, space_full$specs)
  structure(list(
    names = vapply(keep, function(s) s$name, character(1)),
    low = vapply(keep, function(s) as.numeric(s$min), numeric(1)),
    high = vapply(keep, function(s) as.numeric(s$max), numeric(1)),
    start = vapply(keep, function(s) if (!is.null(s$start)) as.numeric(s$start) else 0.5 * (as.numeric(s$min) + as.numeric(s$max)), numeric(1)),
    specs = keep), class = "parameter_space")
}

#' AgMIP stepwise BIC/AICc parameter selection. Mirrors engines/selection.py.
#' @export
stepwise_select <- function(space_full, score_results, criterion = "bic",
                            optimizer = "nelder_mead", optimizer_restarts = 2,
                            maxiter = NULL, min_delta = 0.0, seed = 42, progress = FALSE) {
  criterion <- tolower(criterion)
  start_full <- setNames(
    lapply(space_full$specs, function(s) if (!is.null(s$start)) as.numeric(s$start) else 0.5 * (as.numeric(s$min) + as.numeric(s$max))),
    vapply(space_full$specs, function(s) s$name, character(1)))

  obligatory <- vapply(space_full$specs, function(s) if (!is.null(s$role) && s$role == "obligatory") s$name else NA_character_, character(1))
  obligatory <- obligatory[!is.na(obligatory)]
  candidates <- setdiff(space_full$names, obligatory)
  if (length(obligatory) == 0) { obligatory <- space_full$names[1]; candidates <- setdiff(space_full$names, obligatory) }

  fit <- function(names) {
    held <- start_full[setdiff(space_full$names, names)]
    sub <- .subspace(space_full, names)
    score_batch <- function(subset_thetas) {
      full <- lapply(subset_thetas, function(st) modifyList(held, st))
      vapply(score_results(full), function(r) r$score, numeric(1))
    }
    opt <- run_optimizer(sub, score_batch, method = optimizer, seed = seed,
                         restarts = optimizer_restarts, maxiter = maxiter)
    best_full <- modifyList(held, opt$best_theta)
    res <- score_results(list(best_full))[[1]]
    n <- nrow(res$residuals)
    list(best_full = best_full, crit = .ic_criterion(criterion, res$loglik, length(names), n), n = n)
  }

  selected <- as.list(obligatory)
  f0 <- fit(unlist(selected)); best_full <- f0$best_full; best_crit <- f0$crit; n_obs <- f0$n
  history <- list(list(step = 0, added = paste(obligatory, collapse = "+"),
                       value = round(best_crit, 3), k = length(selected), n = n_obs))
  remaining <- as.list(candidates); step <- 0L
  while (length(remaining) > 0) {
    step <- step + 1L
    trials <- lapply(remaining, function(cand) { f <- fit(c(unlist(selected), cand)); list(cand = cand, crit = f$crit, theta = f$best_full) })
    best_trial <- trials[[which.min(vapply(trials, function(t) t$crit, numeric(1)))]]
    if (best_trial$crit < best_crit - min_delta) {
      selected[[length(selected) + 1L]] <- best_trial$cand
      remaining <- remaining[vapply(remaining, function(r) r != best_trial$cand, logical(1))]
      best_crit <- best_trial$crit; best_full <- best_trial$theta
      history[[length(history) + 1L]] <- list(step = step, added = best_trial$cand,
                                              value = round(best_crit, 3), k = length(selected), n = n_obs)
    } else break
  }
  structure(list(selected = unlist(selected), obligatory = obligatory, criterion = criterion,
                 best_theta = best_full, history = history), class = "selection_result")
}

# ---- sensitivity screening (Morris pure; Sobol via 'sensitivity') ---------

#' Morris elementary-effects screening. Mirrors engines/sensitivity.py:run_morris.
#' @export
run_morris <- function(space, score_results, trajectories = 10, levels = 4, seed = 42, progress = FALSE) {
  k <- ps_ndim(space); set.seed(seed)
  delta <- levels / (2.0 * (levels - 1))
  grid <- seq(0, 1, length.out = levels)
  bases <- grid[grid <= 1.0 - delta + 1e-9]
  traj_pts <- list(); traj_order <- list()
  for (t in seq_len(trajectories)) {
    x <- sample(bases, k, replace = TRUE); order <- sample.int(k)
    pts <- list(x); cur <- x
    for (j in order) { cur <- cur; cur[j] <- cur[j] + delta; pts[[length(pts) + 1L]] <- cur }
    traj_pts[[t]] <- do.call(rbind, pts); traj_order[[t]] <- order
  }
  all_unit <- do.call(rbind, traj_pts)
  native <- sweep(sweep(all_unit, 2, (space$high - space$low), `*`), 2, space$low, `+`)
  thetas <- lapply(seq_len(nrow(native)), function(i) ps_to_theta(space, native[i, ]))
  results <- score_results(thetas)
  scores <- vapply(results, function(r) if (is.finite(r$score)) r$score else NA_real_, numeric(1))
  ee <- setNames(vector("list", k), space$names)
  m <- k + 1L
  for (t in seq_len(trajectories)) {
    y <- scores[((t - 1L) * m + 1L):(t * m)]
    ord <- traj_order[[t]]
    for (step in seq_along(ord)) {
      j <- ord[step]; d <- (y[step + 1L] - y[step]) / delta
      if (is.finite(d)) ee[[space$names[j]]] <- c(ee[[space$names[j]]], d)
    }
  }
  rows <- lapply(space$names, function(name) {
    arr <- ee[[name]]
    data.frame(parameter = name,
               mu_star = if (length(arr)) mean(abs(arr)) else NA_real_,
               mu = if (length(arr)) mean(arr) else NA_real_,
               sigma = if (length(arr)) sd_pop(arr) else NA_real_, stringsAsFactors = FALSE)
  })
  ranking <- do.call(rbind, rows)
  ranking <- ranking[order(-ranking$mu_star), ]; rownames(ranking) <- NULL
  structure(list(method = "morris", ranking = ranking, n_eval = length(thetas)), class = "sensitivity_result")
}

# population sd (np.std) for Morris sigma parity
sd_pop <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) return(NA_real_); sqrt(mean((x - mean(x))^2)) }

#' Pick influential parameters from a ranking table. Mirrors influential_params.
#' @export
influential_params <- function(ranking, keep = NULL, rel_threshold = 0.1) {
  metric <- if ("mu_star" %in% names(ranking)) "mu_star" else "ST"
  r <- ranking[order(-ranking[[metric]]), ]
  if (!is.null(keep)) return(head(r$parameter, keep))
  top <- if (nrow(r)) as.numeric(r[[metric]][1]) else 0.0
  if (top <= 0) return(r$parameter)
  r$parameter[r[[metric]] >= rel_threshold * top]
}

#' Dispatch to Morris (Sobol requires the 'sensitivity' package). Mirrors run_sensitivity.
#' @export
run_sensitivity <- function(space, score_results, method = "morris", ...) {
  method <- tolower(method); args <- list(...)
  if (method == "morris") {
    return(run_morris(space, score_results,
                      trajectories = as.integer(args$trajectories %||% 10),
                      levels = as.integer(args$levels %||% 4),
                      seed = as.integer(args$seed %||% 42)))
  }
  stop(sprintf("sensitivity method '%s' not available in R port yet (use morris)", method))
}

#' Share of output variance per discrete factor (one-way ANOVA). Mirrors anova_variance_share.
#' @export
anova_variance_share <- function(design, factor_cols, response_col) {
  y <- as.numeric(design[[response_col]]); total_ss <- sum((y - mean(y))^2)
  rows <- lapply(factor_cols, function(f) {
    between <- 0.0
    for (lev in unique(design[[f]])) {
      gy <- y[design[[f]] == lev]; between <- between + length(gy) * (mean(gy) - mean(y))^2
    }
    data.frame(factor = f, var_share = if (total_ss > 0) between / total_ss else NA_real_, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); out[order(-out$var_share), ]
}

# ---- NSGA-II (via 'mco'), surrogate (via DiceKriging/ranger) --------------

#' NSGA-II multi-objective Pareto front. Mirrors engines/nsga2.py (uses 'mco').
#' @export
run_nsga2 <- function(evaluate_batch, space, objective_vars, pop_size = 16, n_gen = 5, seed = 42) {
  if (!requireNamespace("mco", quietly = TRUE)) stop("run_nsga2 needs the 'mco' package.")
  big <- 1e6
  fn <- function(x) {
    theta <- ps_to_theta(space, x)
    pv <- evaluate_batch(list(theta))[[1]]
    vapply(objective_vars, function(v) as.numeric(pv[[v]] %||% big), numeric(1))
  }
  set.seed(seed)
  res <- mco::nsga2(fn, idim = ps_ndim(space), odim = length(objective_vars),
                    lower.bounds = space$low, upper.bounds = space$high,
                    popsize = pop_size, generations = n_gen)
  X <- res$par; F <- res$value
  if (is.null(dim(X))) X <- matrix(X, nrow = 1); if (is.null(dim(F))) F <- matrix(F, nrow = 1)
  structure(list(objective_vars = objective_vars, X = X, F = F, param_names = space$names),
            class = "nsga2_result")
}

#' Surrogate/emulator-accelerated calibration. Mirrors engines/surrogate.py
#' (GP via DiceKriging, RF via ranger).
#' @export
run_surrogate <- function(cfg, space, score_results, progress = TRUE) {
  scfg <- .cfg_get(.cfg_get(cfg, "method", list()), "surrogate", list())
  model <- tolower(as.character(.cfg_get(scfg, "engine", .cfg_get(scfg, "model", "gp"))))
  n_train <- as.integer(.cfg_get(scfg, "n_train", 64))
  n_candidates <- as.integer(.cfg_get(scfg, "n_candidates", 5000))
  top_k <- as.integer(.cfg_get(scfg, "top_k", 10))
  seed <- as.integer(.cfg_get(cfg$calibrator, "seed", 42))
  lo <- space$low; hi <- space$high; span <- ifelse((hi - lo) == 0, 1.0, hi - lo)

  train <- sample_design(space, n_train, engine = "lhs", seed = seed, include_start = TRUE)
  train_thetas <- lapply(seq_len(nrow(train)), function(i) ps_to_theta(space, as.numeric(train[i, ])))
  train_res <- score_results(train_thetas)
  Xn <- list(); y <- numeric(0); kept <- list()
  for (i in seq_along(train_thetas)) {
    r <- train_res[[i]]; th <- train_thetas[[i]]
    if (is.finite(r$score)) {
      Xn[[length(Xn) + 1L]] <- (as.numeric(unlist(th[space$names])) - lo) / span
      y <- c(y, r$score); kept[[length(kept) + 1L]] <- list(theta = th, res = r)
    }
  }
  if (length(kept) < 4) stop("surrogate: too few successful training runs to fit an emulator")
  Xn <- do.call(rbind, Xn)
  cand <- sample_design(space, n_candidates, engine = "lhs", seed = seed + 1L, include_start = FALSE)
  cand_native <- as.matrix(cand)
  cand_norm <- sweep(sweep(cand_native, 2, lo, `-`), 2, span, `/`)
  if (model == "rf") {
    if (!requireNamespace("ranger", quietly = TRUE)) stop("surrogate rf needs 'ranger'.")
    est <- ranger::ranger(x = Xn, y = y, num.trees = 300, seed = seed)
    pred <- predict(est, data = cand_norm)$predictions
  } else {
    if (!requireNamespace("DiceKriging", quietly = TRUE)) stop("surrogate gp needs 'DiceKriging'.")
    est <- DiceKriging::km(design = as.data.frame(Xn), response = y, covtype = "gauss",
                           control = list(trace = FALSE))
    pred <- DiceKriging::predict(est, newdata = as.data.frame(cand_norm), type = "UK")$mean
  }
  order_idx <- order(pred)[seq_len(top_k)]
  val_thetas <- lapply(order_idx, function(i) ps_to_theta(space, cand_native[i, ]))
  val_res <- score_results(val_thetas)
  rows <- list(); obj_results <- list(); sid <- 0L
  all_pairs <- c(kept, lapply(seq_along(val_thetas), function(i) list(theta = val_thetas[[i]], res = val_res[[i]])))
  for (pair in all_pairs) {
    obj_results[[as.character(sid)]] <- pair$res
    rows[[length(rows) + 1L]] <- c(list(sample_id = sid), pair$theta,
                                   list(score = pair$res$score, loglik = pair$res$loglik,
                                        n_obs = nrow(pair$res$residuals)))
    sid <- sid + 1L
  }
  design <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  info <- list(model = model, n_train = length(y), n_validated = top_k,
               best_predicted_score = as.numeric(pred[order_idx[1]]))
  structure(list(design = design, obj_results = obj_results, info = info), class = "surrogate_result")
}

# ---- MCMC (adaptive random-walk Metropolis) -------------------------------

#' Adaptive random-walk Metropolis posterior. Mirrors engines/mcmc.py:run_mcmc.
#' @export
run_mcmc <- function(cfg, score_results, space, progress = TRUE) {
  bcfg <- .cfg_get(.cfg_get(cfg, "method", list()), "bayesian", list())
  n_walkers <- as.integer(.cfg_get(bcfg, "n_walkers", max(2L * ps_ndim(space), 8L)))
  n_steps <- as.integer(.cfg_get(bcfg, "n_steps", 400L))
  burn <- as.integer(.cfg_get(bcfg, "burn_in", n_steps %/% 2L))
  thin <- max(1L, as.integer(.cfg_get(bcfg, "thin", 1L)))
  target <- as.numeric(.cfg_get(bcfg, "target_accept", 0.234))
  scale <- as.numeric(.cfg_get(bcfg, "proposal_scale", 0.1))
  adapt_every <- as.integer(.cfg_get(bcfg, "adapt_interval", 20L))
  seed <- as.integer(.cfg_get(cfg$calibrator, "seed", 42)); set.seed(seed)
  ranges <- space$high - space$low; names_ <- space$names
  vec <- function(theta) as.numeric(unlist(theta[names_]))

  init <- sample_prior_design(space, n_walkers)
  cur_theta <- lapply(seq_len(n_walkers), function(w) ps_to_theta(space, as.numeric(init[w, ])))
  init_res <- score_results(cur_theta)
  cur_lp <- vapply(cur_theta, function(t) log_prior_vec(space, t), numeric(1))
  cur_logpost <- vapply(seq_len(n_walkers), function(w) cur_lp[w] + (if (is.finite(init_res[[w]]$loglik)) init_res[[w]]$loglik else -1e300), numeric(1))
  cur_res <- init_res
  initial_design <- do.call(rbind, lapply(seq_len(n_walkers), function(w) as.data.frame(c(list(sample_id = w - 1L), cur_theta[[w]]), stringsAsFactors = FALSE)))

  chain_rows <- list(); samples <- list(); accepts <- 0L; proposals <- 0L
  for (step in seq_len(n_steps) - 1L) {
    sd <- scale * ranges
    prop <- lapply(seq_len(n_walkers), function(w) ps_to_theta(space, vec(cur_theta[[w]]) + rnorm(ps_ndim(space), 0, sd)))
    lp_prop <- vapply(prop, function(t) log_prior_vec(space, t), numeric(1))
    idx_in <- which(is.finite(lp_prop))
    res_in <- if (length(idx_in)) score_results(prop[idx_in]) else list()
    logpost_prop <- rep(-Inf, n_walkers); res_prop <- vector("list", n_walkers)
    for (kk in seq_along(idx_in)) {
      w <- idx_in[kk]; r <- res_in[[kk]]
      ll <- if (is.finite(r$loglik)) r$loglik else -1e300
      logpost_prop[w] <- lp_prop[w] + ll; res_prop[[w]] <- r
    }
    for (w in seq_len(n_walkers)) {
      proposals <- proposals + 1L
      if (log(runif(1)) < (logpost_prop[w] - cur_logpost[w])) {
        cur_theta[[w]] <- prop[[w]]; cur_logpost[w] <- logpost_prop[w]; cur_res[[w]] <- res_prop[[w]]; accepts <- accepts + 1L
      }
      chain_rows[[length(chain_rows) + 1L]] <- c(list(step = step, walker = w - 1L, logpost = cur_logpost[w]), cur_theta[[w]])
      if (step >= burn && ((step - burn) %% thin == 0)) samples[[length(samples) + 1L]] <- list(theta = cur_theta[[w]], res = cur_res[[w]])
    }
    if (step < burn && (step + 1L) %% adapt_every == 0 && proposals > 0) {
      ar <- accepts / proposals; scale <- min(max(scale * exp(ar - target), 1e-3), 1.0)
    }
  }
  if (length(samples) == 0) samples <- lapply(seq_len(n_walkers), function(w) list(theta = cur_theta[[w]], res = cur_res[[w]]))
  rows <- list(); obj_results <- list()
  for (sid in seq_along(samples)) {
    s <- samples[[sid]]; obj_results[[as.character(sid - 1L)]] <- s$res
    rows[[length(rows) + 1L]] <- c(list(sample_id = sid - 1L), s$theta, list(score = s$res$score, loglik = s$res$loglik, n_obs = nrow(s$res$residuals)))
  }
  design <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  design$weight <- 1.0 / nrow(design)
  best_sample_id <- which.min(ifelse(is.finite(design$score), design$score, Inf)) - 1L
  best_theta <- setNames(as.list(as.numeric(design[best_sample_id + 1L, names_])), names_)
  best <- obj_results[[as.character(best_sample_id)]]
  q <- as.numeric(.cfg_get(bcfg, "behavioural_quantile", 0.1))
  valid <- design$score[is.finite(design$score)]
  threshold <- if (length(valid)) as.numeric(quantile(valid, q, names = FALSE, type = 7)) else Inf
  behavioural <- design[is.finite(design$score) & design$score <= threshold, , drop = FALSE]
  structure(list(design = design, behavioural = behavioural, best_theta = best_theta,
                 best_sample_id = best_sample_id, threshold = threshold, ess = nrow(design),
                 obj_results = obj_results, best = best, acceptance = accepts / max(proposals, 1L),
                 chain = do.call(rbind, lapply(chain_rows, function(r) as.data.frame(r, stringsAsFactors = FALSE))),
                 initial_design = initial_design), class = "mcmc_result")
}

# ---- assimilation prototypes (EnKF / forcing / recalibration) -------------

#' Systematic resampling indices (1-based). Mirrors smc_pf.py:_systematic_resample.
#' @export
systematic_resample <- function(weights) {
  n <- length(weights)
  positions <- (runif(1) + (seq_len(n) - 1L)) / n
  cumsum_ <- cumsum(weights); cumsum_[n] <- 1.0
  idx <- findInterval(positions, cumsum_, left.open = FALSE) + 1L
  pmin(pmax(idx, 1L), n)
}

#' Ensemble Kalman Filter update step (UNCOUPLED PROTOTYPE). Mirrors engines/enkf.py.
#' @export
enkf_assimilate <- function(cfg, ensemble_states, obs_var, obs_value, obs_sigma) {
  ecfg <- .cfg_get(.cfg_get(cfg, "assimilation", list()), "enkf", list())
  n_ensemble <- as.integer(.cfg_get(ecfg, "n_ensemble", 50L))
  inflation <- as.numeric(.cfg_get(ecfg, "inflation", 1.05))
  state_vars <- unlist(.cfg_get(ecfg, "state_variables", list("LAID", "CWAD")))
  if (!(obs_var %in% state_vars)) return(ensemble_states)
  j <- match(obs_var, state_vars)
  x_mean <- colMeans(ensemble_states)
  X <- sweep(ensemble_states, 2, x_mean, `-`) * inflation
  H_x <- ensemble_states[, j]; HP <- X[, j]
  HPHt <- sd_pop(H_x)^2 + obs_sigma^2
  if (HPHt == 0) return(ensemble_states)
  K <- as.numeric(t(X) %*% HP) / ((n_ensemble - 1) * HPHt)
  obs_perturbed <- obs_value + obs_sigma * rnorm(nrow(ensemble_states))
  innovation <- obs_perturbed - H_x
  ensemble_states + outer(innovation, K)
}

#' Direct state replacement (UNCOUPLED PROTOTYPE). Mirrors engines/forcing.py.
#' @export
forcing_apply <- function(cfg, model_state, observation) {
  fcfg <- .cfg_get(.cfg_get(cfg, "assimilation", list()), "forcing", list())
  min_conf <- as.numeric(.cfg_get(fcfg, "min_confidence", 0.8))
  smoothing <- isTRUE(.cfg_get(fcfg, "smoothing", TRUE))
  updated <- model_state
  var <- observation$variable; val <- observation$value
  conf <- observation$confidence %||% 1.0
  if (conf >= min_conf) {
    if (smoothing && !is.null(updated[[var]])) updated[[var]] <- conf * val + (1.0 - conf) * updated[[var]]
    else updated[[var]] <- val
  }
  updated
}

# ===========================================================================
# Added engines (R twins): DREAM (DE-MC), ES-MDA, Bayesian optimisation.
# CMA-ES lives inside run_optimizer (method = "cmaes"). These mirror the Python
# engines in python/dssatcalibrator/engines/{dream,es_mda,bayesopt}.py.
# ===========================================================================

# Build the shared posterior-style result (same shape as run_mcmc's mcmc_result)
# from a list of list(theta=, res=) samples.
.posterior_result <- function(samples, names_, bcfg, acceptance = NA_real_,
                              chain_rows = list(), initial_design = NULL) {
  rows <- list()
  obj_results <- list()
  for (sid in seq_along(samples)) {
    s <- samples[[sid]]
    obj_results[[as.character(sid - 1L)]] <- s$res
    rows[[length(rows) + 1L]] <- c(list(sample_id = sid - 1L), s$theta,
                                   list(score = s$res$score, loglik = s$res$loglik,
                                        n_obs = nrow(s$res$residuals)))
  }
  design <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  design$weight <- 1.0 / nrow(design)
  best_sample_id <- which.min(ifelse(is.finite(design$score), design$score, Inf)) - 1L
  best_theta <- setNames(as.list(as.numeric(design[best_sample_id + 1L, names_])), names_)
  best <- obj_results[[as.character(best_sample_id)]]
  q <- as.numeric(.cfg_get(bcfg, "behavioural_quantile", 0.1))
  valid <- design$score[is.finite(design$score)]
  threshold <- if (length(valid)) as.numeric(quantile(valid, q, names = FALSE, type = 7)) else Inf
  behavioural <- design[is.finite(design$score) & design$score <= threshold, , drop = FALSE]
  chain <- if (length(chain_rows))
    do.call(rbind, lapply(chain_rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  else data.frame()
  structure(list(design = design, behavioural = behavioural, best_theta = best_theta,
                 best_sample_id = best_sample_id, threshold = threshold, ess = nrow(design),
                 obj_results = obj_results, best = best, acceptance = acceptance,
                 chain = chain, initial_design = initial_design), class = "mcmc_result")
}

#' DREAM (DE-MC) — robust Bayesian posterior. Mirrors engines/dream.py.
#' @export
run_dream <- function(cfg, score_results, space, progress = TRUE) {
  bcfg <- .cfg_get(.cfg_get(cfg, "method", list()), "bayesian", list())
  d <- ps_ndim(space); names_ <- space$names
  n_chains <- as.integer(.cfg_get(bcfg, "n_chains", max(2L * d, 6L)))
  if (n_chains < 4L) n_chains <- 4L
  n_gen <- as.integer(.cfg_get(bcfg, "n_generations", .cfg_get(bcfg, "n_steps", 400L)))
  burn <- as.integer(.cfg_get(bcfg, "burn_in", n_gen %/% 2L))
  thin <- max(1L, as.integer(.cfg_get(bcfg, "thin", 1L)))
  snooker <- as.numeric(.cfg_get(bcfg, "snooker", 0.1))
  eps <- as.numeric(.cfg_get(bcfg, "eps", 1e-4))
  set.seed(as.integer(.cfg_get(cfg$calibrator, "seed", 42L)))
  ranges <- space$high - space$low
  gamma_default <- 2.38 / sqrt(2 * d)
  vec <- function(theta) as.numeric(unlist(theta[names_]))

  init <- sample_prior_design(space, n_chains)
  cur_theta <- lapply(seq_len(n_chains), function(c) ps_to_theta(space, as.numeric(init[c, ])))
  cur_res <- score_results(cur_theta)
  cur_lp <- vapply(cur_theta, function(t) log_prior_vec(space, t), numeric(1))
  cur_lpost <- vapply(seq_len(n_chains), function(c)
    cur_lp[c] + (if (is.finite(cur_res[[c]]$loglik)) cur_res[[c]]$loglik else -1e300), numeric(1))
  initial_design <- do.call(rbind, lapply(seq_len(n_chains), function(c)
    as.data.frame(c(list(sample_id = c - 1L), cur_theta[[c]]), stringsAsFactors = FALSE)))

  samples <- list(); accepts <- 0L; proposals <- 0L
  for (gen in seq_len(n_gen)) {
    cur_vecs <- lapply(cur_theta, vec)
    prop <- vector("list", n_chains)
    for (c in seq_len(n_chains)) {
      others <- setdiff(seq_len(n_chains), c)
      ab <- sample(others, 2L)
      gamma <- if (runif(1) < snooker) 1.0 else gamma_default
      jump <- gamma * (cur_vecs[[ab[1]]] - cur_vecs[[ab[2]]]) + eps * ranges * rnorm(d)
      prop[[c]] <- ps_to_theta(space, ps_clip(space, cur_vecs[[c]] + jump))
    }
    lp_prop <- vapply(prop, function(t) log_prior_vec(space, t), numeric(1))
    idx_in <- which(is.finite(lp_prop))
    res_in <- if (length(idx_in)) score_results(prop[idx_in]) else list()
    lpost_prop <- rep(-Inf, n_chains); res_prop <- vector("list", n_chains)
    for (j in seq_along(idx_in)) {
      c <- idx_in[j]; r <- res_in[[j]]
      ll <- if (is.finite(r$loglik)) r$loglik else -1e300
      lpost_prop[c] <- lp_prop[c] + ll; res_prop[[c]] <- r
    }
    for (c in seq_len(n_chains)) {
      proposals <- proposals + 1L
      if (log(runif(1)) < (lpost_prop[c] - cur_lpost[c])) {
        cur_theta[[c]] <- prop[[c]]; cur_lpost[c] <- lpost_prop[c]; cur_res[[c]] <- res_prop[[c]]
        accepts <- accepts + 1L
      }
      if (gen > burn && ((gen - burn - 1L) %% thin == 0L))
        samples[[length(samples) + 1L]] <- list(theta = cur_theta[[c]], res = cur_res[[c]])
    }
  }
  if (!length(samples))
    samples <- lapply(seq_len(n_chains), function(c) list(theta = cur_theta[[c]], res = cur_res[[c]]))
  .posterior_result(samples, names_, bcfg, acceptance = accepts / max(proposals, 1L),
                    initial_design = initial_design)
}

#' ES-MDA — Ensemble Smoother with Multiple Data Assimilation. Mirrors engines/es_mda.py.
#' @export
run_es_mda <- function(cfg, score_results, space, progress = TRUE) {
  bcfg <- .cfg_get(.cfg_get(cfg, "method", list()), "bayesian", list())
  d <- ps_ndim(space); names_ <- space$names
  ne <- as.integer(.cfg_get(bcfg, "ensemble_size", max(4L * d, 24L)))
  na <- as.integer(.cfg_get(bcfg, "iterations", 4L))
  set.seed(as.integer(.cfg_get(cfg$calibrator, "seed", 42L)))
  alpha <- as.numeric(na)

  init <- sample_prior_design(space, ne)
  ens <- as.matrix(init); colnames(ens) <- names_
  results <- score_results(lapply(seq_len(ne), function(i) ps_to_theta(space, ens[i, ])))

  obs_vectors <- function(results) {
    per <- lapply(results, function(r) {
      rd <- r$residuals
      if (is.null(rd) || nrow(rd) == 0) return(list())
      keys <- paste(rd$exp_id, rd$treatment, rd$dssat,
                    ifelse(is.na(rd$date), "NA", as.character(rd$date)), sep = "|")
      setNames(lapply(seq_len(nrow(rd)), function(i)
        c(sim = rd$sim[i], obs = rd$obs[i], sigma = rd$sigma[i])), keys)
    })
    common <- NULL
    for (m in per) if (length(m)) common <- if (is.null(common)) names(m) else intersect(common, names(m))
    common <- if (is.null(common)) character(0) else sort(common)
    if (!length(common)) return(NULL)
    ref <- per[[which(vapply(per, function(m) length(m) > 0, logical(1)))[1]]]
    d_obs <- vapply(common, function(k) ref[[k]]["obs"], numeric(1))
    sigma <- vapply(common, function(k) ref[[k]]["sigma"], numeric(1))
    d_sim <- matrix(NA_real_, nrow = length(per), ncol = length(common))
    for (i in seq_along(per)) for (j in seq_along(common))
      if (!is.null(per[[i]][[common[j]]])) d_sim[i, j] <- per[[i]][[common[j]]]["sim"]
    list(d_obs = d_obs, sigma = sigma, d_sim = d_sim)
  }

  for (it in seq_len(na)) {
    ov <- obs_vectors(results)
    if (is.null(ov)) break
    d_obs <- ov$d_obs; sigma <- ov$sigma; d_sim <- ov$d_sim; nd <- length(d_obs)
    bad <- !is.finite(d_sim)
    if (any(bad)) d_sim[bad] <- d_obs[col(d_sim)[bad]]
    theta_mean <- colMeans(ens); d_mean <- colMeans(d_sim)
    Ta <- sweep(ens, 2, theta_mean); Da <- sweep(d_sim, 2, d_mean)
    C_td <- (t(Ta) %*% Da) / (ne - 1)
    C_dd <- (t(Da) %*% Da) / (ne - 1)
    R <- alpha * diag(sigma^2, nd)
    K <- C_td %*% MASS_ginv(C_dd + R)
    pert <- matrix(rep(d_obs, each = ne), nrow = ne) +
      sqrt(alpha) * matrix(rep(sigma, each = ne), nrow = ne) * matrix(rnorm(ne * nd), ne, nd)
    ens <- ens + (pert - d_sim) %*% t(K)
    ens <- t(apply(ens, 1, function(x) pmin(pmax(x, space$low), space$high)))
    results <- score_results(lapply(seq_len(ne), function(i) ps_to_theta(space, ens[i, ])))
  }
  samples <- lapply(seq_len(ne), function(i) list(theta = ps_to_theta(space, ens[i, ]), res = results[[i]]))
  .posterior_result(samples, names_, bcfg, initial_design =
    do.call(rbind, lapply(seq_len(ne), function(i)
      as.data.frame(c(list(sample_id = i - 1L), ps_to_theta(space, as.numeric(init[i, ]))),
                    stringsAsFactors = FALSE))))
}

# Moore-Penrose pseudo-inverse without importing MASS (small, robust via SVD).
MASS_ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  s <- svd(X); pos <- s$d > max(tol * s$d[1], 0)
  s$v[, pos, drop = FALSE] %*% ((1 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
}

#' Bayesian optimisation (GP + Expected Improvement). Mirrors engines/bayesopt.py.
#' Requires the 'DiceKriging' package for the GP surrogate.
#' @export
run_bayesopt <- function(cfg, score_results, space, progress = TRUE) {
  if (!requireNamespace("DiceKriging", quietly = TRUE))
    stop("Bayesian optimisation (R) needs the 'DiceKriging' package.")
  bcfg <- .cfg_get(.cfg_get(cfg, "method", list()), "bayesian", list())
  d <- ps_ndim(space); names_ <- space$names
  n_init <- as.integer(.cfg_get(bcfg, "n_init", max(2L * d, 10L)))
  n_iter <- as.integer(.cfg_get(bcfg, "n_iter", 20L))
  batch <- max(1L, as.integer(.cfg_get(bcfg, "batch_size", 4L)))
  xi_frac <- as.numeric(.cfg_get(bcfg, "xi", 0.01))
  set.seed(as.integer(.cfg_get(cfg$calibrator, "seed", 42L)))
  low <- space$low; high <- space$high; span <- ifelse(high > low, high - low, 1)

  init <- sample_prior_design(space, n_init); X <- as.matrix(init)
  results <- score_results(lapply(seq_len(n_init), function(i) ps_to_theta(space, X[i, ])))
  y <- vapply(results, function(r) if (is.finite(r$score)) r$score else 1e12, numeric(1))

  for (it in seq_len(n_iter)) {
    Xu <- sweep(sweep(X, 2, low), 2, span, "/")
    ym <- mean(y); ys <- stats::sd(y); if (!is.finite(ys) || ys == 0) ys <- 1
    m <- suppressWarnings(DiceKriging::km(design = as.data.frame(Xu), response = (y - ym) / ys,
                  covtype = "matern5_2", control = list(trace = FALSE), nugget = 1e-6))
    cand <- matrix(runif(2000 * d), 2000, d)
    pr <- DiceKriging::predict(m, newdata = as.data.frame(cand), type = "UK",
                               checkNames = FALSE)
    mu <- pr$mean * ys + ym; sd_ <- pr$sd * ys
    best_y <- min(y); imp <- best_y - mu - xi_frac * (max(y) - best_y + 1e-9)
    z <- imp / pmax(sd_, 1e-12)
    ei <- imp * pnorm(z) + pmax(sd_, 1e-12) * dnorm(z)
    picks <- order(ei, decreasing = TRUE)[seq_len(batch)]
    newX <- sweep(sweep(cand[picks, , drop = FALSE], 2, span, "*"), 2, low, "+")
    new_res <- score_results(lapply(seq_len(nrow(newX)), function(i) ps_to_theta(space, newX[i, ])))
    X <- rbind(X, newX); y <- c(y, vapply(new_res, function(r) if (is.finite(r$score)) r$score else 1e12, numeric(1)))
    results <- c(results, new_res)
  }
  rows <- list(); obj_results <- list()
  for (sid in seq_along(results)) {
    obj_results[[as.character(sid - 1L)]] <- results[[sid]]
    rows[[length(rows) + 1L]] <- c(list(sample_id = sid - 1L), ps_to_theta(space, X[sid, ]),
                                   list(score = y[sid], loglik = results[[sid]]$loglik,
                                        n_obs = nrow(results[[sid]]$residuals)))
  }
  design <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  best_id0 <- which.min(y) - 1L; design$weight <- 0; design$weight[best_id0 + 1L] <- 1
  best_theta <- ps_to_theta(space, X[best_id0 + 1L, ])
  structure(list(design = design, obj_results = obj_results, best_theta = best_theta,
                 best_sample_id = best_id0, best = obj_results[[as.character(best_id0)]],
                 info = list(n_eval = length(results))), class = "bayesopt_result")
}

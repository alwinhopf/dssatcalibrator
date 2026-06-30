# Multi-variable, multi-experiment objective: align sim vs observed and score.
# R twin of python/dssatcalibrator/objective.py.
#
# Two alignment paths (CONCEPT.md §8):
#   * scalars / phenology come from Evaluate.OUT (DSSAT pairs sim & measured).
#   * time-series are matched from PlantGro.OUT to FileT/CSV obs by (treatment,
#     date), averaging replicates.
#
# Scoring exposes a minimisation `score` and a maximisation `loglik`
# (sigma-weighted Gaussian) so optimisers/GLUE and the Bayesian engines share
# one residual table. Per-variable metrics: RMSE, nRMSE%, MBE, Willmott d, EF, R2.

#' @export
variable_maps <- function(cfg) {
  eng <- if (!is.null(cfg$engine)) cfg$engine else list()
  ts <- if (!is.null(eng$timeseries_outputs)) eng$timeseries_outputs else list()
  sc <- if (!is.null(eng$scalar_outputs)) eng$scalar_outputs else list()
  inv <- function(m) {
    if (length(m) == 0) return(list())
    setNames(as.list(names(m)), unlist(m, use.names = FALSE))
  }
  list(ts = ts, sc = sc, ts_inv = inv(ts), sc_inv = inv(sc))
}

#' DSSAT variables present in observations but NOT scorable.
#' Mirrors objective.py:unmatched_variables.
#' @export
unmatched_variables <- function(obs_table, cfg) {
  vm <- variable_maps(cfg)
  known <- unique(c(unlist(vm$ts, use.names = FALSE), unlist(vm$sc, use.names = FALSE)))
  if (is.null(obs_table) || nrow(obs_table) == 0 || !"variable" %in% names(obs_table)) {
    return(character(0))
  }
  present <- unique(obs_table$variable[!is.na(obs_table$variable)])
  sort(setdiff(present, known))
}

#' RMSE / nRMSE% / MBE / Willmott d / modelling efficiency EF / R2 + n.
#' Mirrors objective.py:metrics exactly.
#' @export
metrics <- function(obs, sim) {
  obs <- as.numeric(obs); sim <- as.numeric(sim)
  ok <- is.finite(obs) & is.finite(sim)
  o <- obs[ok]; s <- sim[ok]
  n <- length(o)
  base <- list(n = n, RMSE = NA_real_, nRMSE_pct = NA_real_, MBE = NA_real_,
               d = NA_real_, EF = NA_real_, R2 = NA_real_)
  if (n == 0) return(base)
  rmse <- sqrt(mean((s - o)^2))
  ob <- mean(o)
  base$RMSE <- rmse
  base$nRMSE_pct <- if (ob != 0) 100 * rmse / ob else NA_real_
  base$MBE <- mean(s - o)
  d_den <- sum((abs(s - ob) + abs(o - ob))^2)
  base$d <- if (d_den > 0) 1 - sum((s - o)^2) / d_den else NA_real_
  o_var <- sum((o - ob)^2)
  base$EF <- if (o_var > 0) 1 - sum((o - s)^2) / o_var else NA_real_
  if (n > 1 && sd(o) > 0 && sd(s) > 0) base$R2 <- cor(o, s)^2
  base
}

.obj_sigma <- function(user_var, obs_value, cfg) {
  em <- if (!is.null(cfg$objective$error_model)) cfg$objective$error_model else list()
  spec <- em[[user_var]]
  if (is.null(spec)) {
    sigma <- max(abs(0.10 * obs_value), 1e-6)   # default: 10% relative
  } else {
    type <- if (!is.null(spec$type)) as.character(spec$type) else "relative"
    sigma <- if (type == "relative") max(abs(as.numeric(spec$value) * obs_value), 1e-6) else as.numeric(spec$value)
  }
  disc_cfg <- if (!is.null(cfg$objective$model_discrepancy)) cfg$objective$model_discrepancy else list()
  disc <- as.numeric(.cfg_get(disc_cfg, "default", .cfg_get(disc_cfg, "value", 0.0)))
  variables <- .cfg_get(disc_cfg, "variables", list())
  relative <- .cfg_get(disc_cfg, "relative", list())
  if (!is.null(variables[[user_var]])) disc <- as.numeric(variables[[user_var]])
  if (!is.null(relative[[user_var]])) disc <- max(disc, abs(as.numeric(relative[[user_var]]) * obs_value))
  if (is.finite(disc) && disc > 0) sigma <- sqrt(sigma^2 + disc^2)
  max(as.numeric(sigma), 1e-12)
}

.standardized_loss <- function(z, cfg) {
  lcfg <- .cfg_get(.cfg_get(cfg, "objective", list()), "likelihood", list(type = "gaussian"))
  if (is.character(lcfg)) {
    kind <- tolower(lcfg)
    lcfg <- list()
  } else {
    kind <- tolower(as.character(.cfg_get(lcfg, "type", "gaussian")))
  }
  z <- as.numeric(z)
  if (kind %in% c("student_t", "student-t", "t")) {
    nu <- max(as.numeric(.cfg_get(lcfg, "df", .cfg_get(lcfg, "nu", 4.0))), 1.01)
    return((nu + 1.0) * log1p((z^2) / nu))
  }
  if (kind == "huber") {
    delta <- max(as.numeric(.cfg_get(lcfg, "delta", 1.5)), 1e-9)
    az <- abs(z)
    return(ifelse(az <= delta, z^2, 2.0 * delta * az - delta^2))
  }
  z^2
}

#' Down-weight dense time-series for serial correlation (obs_autocorr: true).
#' Mirrors objective.py:_downweight_autocorr (AR(1) effective-sample-size factor).
.downweight_autocorr <- function(df) {
  ts <- df[df$kind == "timeseries", , drop = FALSE]
  if (nrow(ts) == 0) return(df)
  key <- paste(ts$exp_id, ts$user_var, ts$treatment, sep = "\r")
  for (k in unique(key)) {
    idx <- which(key == k)
    g <- ts[idx, , drop = FALSE]
    if (nrow(g) < 3) next
    g <- g[order(g$date), ]
    x <- as.numeric(g$obs); x <- x - mean(x)
    denom <- sum(x * x)
    if (denom <= 0) next
    rho <- sum(x[-length(x)] * x[-1]) / denom
    rho <- min(max(rho, 0.0), 0.99)
    factor <- (1.0 - rho) / (1.0 + rho)
    rows <- as.integer(rownames(g))
    df[rows, "weight"] <- df[rows, "weight"] * factor
  }
  df
}

#' Assemble the residual table (one row per matched observation).
#' Mirrors objective.py:build_residuals. `results` is a named list of
#' per-experiment result objects each carrying `$evaluate` and `$plantgro`
#' data.frames.
#' @export
build_residuals <- function(results, obs_table, cfg) {
  vm <- variable_maps(cfg)
  ts_inv <- vm$ts_inv; sc_inv <- vm$sc_inv; sc_map <- vm$sc
  rows <- list()
  add <- function(r) rows[[length(rows) + 1L]] <<- r
  pheno <- c("ADAP", "EDAP", "MDAP")
  seen_scalar <- character(0)

  # scalars / phenology from Evaluate.OUT
  for (exp in names(results)) {
    ev <- results[[exp]]$evaluate
    if (is.null(ev) || nrow(ev) == 0) next
    for (i in seq_len(nrow(ev))) {
      r <- ev[i, ]
      base_var <- r$variable
      if (!base_var %in% names(sc_inv)) next
      sim <- r$sim; meas <- r$meas
      if (is.na(sim) || is.na(meas)) {
        if (!is.na(meas) && is.na(sim)) {
          penalty_sim <- meas + 1000.0
          kind <- if (base_var %in% pheno) "phenology" else "scalar"
          add(data.frame(exp_id = exp, treatment = as.integer(r$treatment),
                         user_var = sc_inv[[base_var]], dssat = base_var, kind = kind,
                         date = as.Date(NA), obs = as.numeric(meas),
                         sim = as.numeric(penalty_sim), stringsAsFactors = FALSE))
        }
        next
      }
      kind <- if (base_var %in% pheno) "phenology" else "scalar"
      add(data.frame(exp_id = exp, treatment = as.integer(r$treatment),
                     user_var = sc_inv[[base_var]], dssat = base_var, kind = kind,
                     date = as.Date(NA), obs = as.numeric(meas), sim = as.numeric(sim),
                     stringsAsFactors = FALSE))
      seen_scalar <- c(seen_scalar, paste(exp, as.integer(r$treatment), sc_inv[[base_var]], kind, sep = "\r"))
    }
  }

  # CSV / fused scalar observations matched to Evaluate.OUT simulated values.
  if (!is.null(obs_table) && nrow(obs_table) > 0) {
    scalar_obs <- obs_table[obs_table$kind %in% c("scalar", "phenology"), , drop = FALSE]
    for (exp in names(results)) {
      ev <- results[[exp]]$evaluate
      if (is.null(ev) || nrow(ev) == 0) next
      o <- scalar_obs[scalar_obs$exp_id == exp, , drop = FALSE]
      if (nrow(o) == 0) next
      agg <- aggregate(value ~ treatment + variable + kind, data = o, FUN = mean)
      for (i in seq_len(nrow(agg))) {
        r <- agg[i, ]
        user_var <- r$variable
        dssat_var <- if (!is.null(sc_map[[user_var]])) sc_map[[user_var]] else user_var
        if (!(dssat_var %in% ev$variable)) next
        key <- paste(exp, as.integer(r$treatment), user_var, r$kind, sep = "\r")
        if (key %in% seen_scalar) next
        sub <- ev[ev$treatment == r$treatment & ev$variable == dssat_var, , drop = FALSE]
        if (nrow(sub) == 0 || is.na(sub$sim[1]) || is.na(r$value)) next
        add(data.frame(exp_id = exp, treatment = as.integer(r$treatment),
                       user_var = user_var, dssat = dssat_var, kind = r$kind,
                       date = as.Date(NA), obs = as.numeric(r$value),
                       sim = as.numeric(sub$sim[1]), stringsAsFactors = FALSE))
        seen_scalar <- c(seen_scalar, key)
      }
    }
  }

  # time-series from PlantGro matched to FileT/CSV obs
  if (!is.null(obs_table) && nrow(obs_table) > 0) {
    ts_obs <- obs_table[obs_table$kind == "timeseries", , drop = FALSE]
    for (exp in names(results)) {
      pg <- results[[exp]]$plantgro
      if (is.null(pg) || nrow(pg) == 0) next
      o <- ts_obs[ts_obs$exp_id == exp, , drop = FALSE]
      if (nrow(o) == 0) next
      agg <- aggregate(value ~ treatment + date + variable, data = o, FUN = mean)
      for (i in seq_len(nrow(agg))) {
        r <- agg[i, ]
        col <- r$variable
        if (!col %in% names(pg)) next
        sub <- pg[pg$treatment == r$treatment & pg$date == r$date, , drop = FALSE]
        if (nrow(sub) == 0) next
        simv <- sub[[col]][1]
        if (is.na(simv)) next
        uv <- if (!is.null(ts_inv[[col]])) ts_inv[[col]] else col
        add(data.frame(exp_id = exp, treatment = as.integer(r$treatment),
                       user_var = uv, dssat = col, kind = "timeseries",
                       date = r$date, obs = as.numeric(r$value), sim = as.numeric(simv),
                       stringsAsFactors = FALSE))
      }
    }
  }

  if (length(rows) == 0) return(data.frame())
  df <- do.call(rbind, rows)
  rownames(df) <- seq_len(nrow(df))

  wts <- if (!is.null(cfg$objective$weights)) cfg$objective$weights else list()
  df$sigma <- mapply(function(uv, ov) .obj_sigma(uv, ov, cfg), df$user_var, df$obs)
  df$weight <- vapply(df$user_var, function(v) {
    if (!is.null(wts[[v]])) as.numeric(wts[[v]]) else 1.0
  }, numeric(1))

  df$resid <- df$sim - df$obs
  if (isTRUE(cfg$objective$obs_autocorr)) df <- .downweight_autocorr(df)
  df
}

#' Score a set of per-experiment spawn results against observations.
#' Returns a list with `score` (minimise), `loglik` (maximise), `residuals`,
#' `per_var`, `per_exp_var`. Mirrors objective.py:score.
#' @export
score <- function(results, obs_table, cfg) {
  resid <- build_residuals(results, obs_table, cfg)
  if (nrow(resid) == 0) {
    return(structure(list(score = Inf, loglik = -Inf, residuals = resid,
                          per_var = list(), per_exp_var = data.frame()),
                     class = "objective_result"))
  }
  weighting <- if (!is.null(cfg$objective$weighting)) cfg$objective$weighting else "unified"
  wts <- if (!is.null(cfg$objective$weights)) cfg$objective$weights else list()
  w <- function(uv) if (!is.null(wts[[uv]])) as.numeric(wts[[uv]]) else 1.0

  resid$`_loss` <- .standardized_loss(resid$resid / resid$sigma, cfg)
  loglik <- -0.5 * sum(resid$`_loss` * resid$weight)

  by_var <- split(resid, resid$user_var)
  if (weighting == "sigma") {
    sc <- sum(resid$`_loss` * resid$weight)
  } else if (weighting == "user") {
    sc <- 0.0
    for (uv in names(by_var)) {
      g <- by_var[[uv]]
      sc <- sc + w(uv) * mean(g$`_loss`)
    }
  } else if (weighting == "count_scale") {
    per <- vapply(names(by_var), function(uv) {
      g <- by_var[[uv]]; w(uv) * mean(g$`_loss`)
    }, numeric(1))
    sc <- if (length(per)) mean(per) else Inf
  } else {  # "unified" (default) and "agmip_wls"
    sc <- 0.0
    for (uv in names(by_var)) {
      g <- by_var[[uv]]
      mse <- mean(g$`_loss`)
      sc <- sc + w(uv) * mse
    }
  }

  per_var <- lapply(by_var, function(g) metrics(g$obs, g$sim))
  pev_rows <- list()
  by_ev <- split(resid, list(resid$exp_id, resid$user_var), drop = TRUE)
  for (nm in names(by_ev)) {
    g <- by_ev[[nm]]
    m <- metrics(g$obs, g$sim)
    pev_rows[[length(pev_rows) + 1L]] <- cbind(
      data.frame(exp_id = g$exp_id[1], user_var = g$user_var[1], stringsAsFactors = FALSE),
      as.data.frame(m)
    )
  }
  pev <- if (length(pev_rows)) do.call(rbind, pev_rows) else data.frame()

  structure(list(score = as.numeric(sc), loglik = as.numeric(loglik),
                 residuals = resid, per_var = per_var, per_exp_var = pev),
            class = "objective_result")
}

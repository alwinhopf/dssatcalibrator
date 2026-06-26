# In-season forecast, identifiability/adequacy diagnostics, and new-crop
# scaffolding. R twins of python/dssatcalibrator/{forecast,diagnostics,scaffold}.py.
# The pure numeric helpers (percentiles, anchor, diagnostics) port exactly.

# ---- forecast (LAI nowcast) -----------------------------------------------

#' Daily percentiles of `variable` across an ensemble of PlantGro curves.
#' Mirrors forecast.py:ensemble_percentiles.
#' @export
ensemble_percentiles <- function(curves, variable = "LAID", quantiles = c(0.1, 0.5, 0.9)) {
  frames <- list()
  for (c in curves) {
    if (is.null(c) || nrow(c) == 0 || !(variable %in% names(c)) || !("date" %in% names(c))) next
    sub <- c[, c("date", variable)]; sub <- sub[stats::complete.cases(sub), ]
    frames[[length(frames) + 1L]] <- sub
  }
  if (length(frames) == 0) return(data.frame(date = as.Date(character(0)), p10 = numeric(0),
                                             p50 = numeric(0), p90 = numeric(0), mean = numeric(0), n = integer(0)))
  allc <- do.call(rbind, frames)
  dates <- sort(unique(allc$date))
  do.call(rbind, lapply(dates, function(d) {
    v <- allc[[variable]][allc$date == d]
    data.frame(date = d,
               p10 = as.numeric(quantile(v, quantiles[1], names = FALSE, type = 7)),
               p50 = as.numeric(quantile(v, quantiles[2], names = FALSE, type = 7)),
               p90 = as.numeric(quantile(v, quantiles[3], names = FALSE, type = 7)),
               mean = mean(v), n = length(v), stringsAsFactors = FALSE)
  }))
}

#' Shift a forecast to start from the last observation, decaying over `decay_days`.
#' Mirrors forecast.py:anchor_correction.
#' @export
anchor_correction <- function(forecast, last_obs_value, last_obs_date, decay_days = 21,
                              mode = "additive", cols = c("p10", "p50", "p90", "mean")) {
  out <- forecast
  if (nrow(out) == 0 || is.na(last_obs_value)) return(out)
  last_obs_date <- as.Date(last_obs_date)
  anchor_row <- out[out$date == last_obs_date, , drop = FALSE]
  if (nrow(anchor_row) == 0) {
    prior <- out[out$date <= last_obs_date, , drop = FALSE]
    if (nrow(prior) == 0) return(out)
    anchor_row <- prior[nrow(prior), , drop = FALSE]
  }
  sim_at_anchor <- as.numeric(anchor_row$p50[1])
  if (mode == "multiplicative") { if (sim_at_anchor == 0) return(out); full <- last_obs_value / sim_at_anchor }
  else full <- last_obs_value - sim_at_anchor
  days_past <- pmax(as.numeric(out$date - last_obs_date), 0)
  weight <- ifelse(days_past <= decay_days, 1.0 - days_past / max(decay_days, 1), 0.0)
  weight <- ifelse(out$date < last_obs_date, 0.0, weight)
  for (c in cols) {
    if (!(c %in% names(out))) next
    if (mode == "multiplicative") out[[paste0(c, "_adj")]] <- out[[c]] * (1.0 + (full - 1.0) * weight)
    else out[[paste0(c, "_adj")]] <- out[[c]] + full * weight
  }
  out$anchor_weight <- weight
  out
}

#' Forecast spread (P90-P10) vs lead time. Mirrors forecast.py:lead_time_table.
#' @export
lead_time_table <- function(forecast, last_obs_date) {
  if (nrow(forecast) == 0) return(data.frame(lead_days = integer(0), date = as.Date(character(0)),
                                            p50 = numeric(0), spread = numeric(0), rel_spread = numeric(0)))
  last_obs_date <- as.Date(last_obs_date)
  out <- forecast[forecast$date >= last_obs_date, , drop = FALSE]
  out$lead_days <- as.integer(out$date - last_obs_date)
  out$spread <- out$p90 - out$p10
  out$rel_spread <- out$spread / ifelse(out$p50 == 0, NA, out$p50)
  out[, c("lead_days", "date", "p50", "spread", "rel_spread")]
}

.behavioural_thetas <- function(result, n) {
  space <- result$space; design <- result$design
  if (n <= 0 || is.null(design) || nrow(design) == 0) return(list(result$best_theta))
  cols <- intersect(space$names, names(design))
  if (length(cols) == 0) return(list(result$best_theta))
  ranked <- if ("score" %in% names(design)) design[order(design$score), ] else design
  k <- min(n, nrow(ranked))
  thetas <- lapply(seq_len(k), function(i) setNames(as.list(as.numeric(ranked[i, cols])), cols))
  c(list(result$best_theta), thetas)[seq_len(max(n, 1))]
}

#' Propagate the calibrated ensemble forward into per-experiment LAI forecasts.
#' Mirrors forecast.py:forecast_lai. (Used by the LAI assimilation/nowcast mode.)
#' @export
forecast_lai <- function(cfg, result, last_obs = NULL, variable = "LAID") {
  fcfg <- .cfg_get(cfg, "forecast", list())
  n_ens <- as.integer(.cfg_get(fcfg, "n_ensemble", 0))
  anchor <- isTRUE(.cfg_get(fcfg, "anchor_continuity", TRUE))
  decay <- as.integer(.cfg_get(fcfg, "decay_days", 21))
  last_obs <- last_obs %||% list()
  thetas <- .behavioural_thetas(result, n_ens)
  out <- list()
  for (exp in result$experiments) {
    curves <- list()
    for (th in thetas) {
      spawns <- spawn_results_for(cfg, th, exp)
      pg <- spawns[[exp]]$plantgro
      if (nrow(pg) > 0) curves[[length(curves) + 1L]] <- pg
    }
    fc <- ensemble_percentiles(curves, variable = variable)
    if (anchor && !is.null(last_obs[[exp]])) {
      lo <- last_obs[[exp]]; fc <- anchor_correction(fc, lo[[2]], lo[[1]], decay_days = decay)
    }
    out[[exp]] <- fc
  }
  out
}

# ---- diagnostics ----------------------------------------------------------

.prior_std <- function(spec) {
  prior <- spec$prior; if (is.null(prior)) prior <- list()
  dist <- tolower(as.character(.cfg_get(prior, "dist", "uniform")))
  lo <- as.numeric(spec$min); hi <- as.numeric(spec$max)
  if (dist %in% c("normal", "lognormal") && !is.null(prior$sd)) return(as.numeric(prior$sd))
  (hi - lo) / sqrt(12.0)
}

#' Per-parameter identifiability (posterior vs prior width + collinearity).
#' Mirrors diagnostics.py:identifiability.
#' @export
identifiability <- function(result, behavioural_quantile = 0.1) {
  cfg <- result$cfg; space <- result$space; design <- result$design
  specs <- setNames(active_parameters(cfg), vapply(active_parameters(cfg), function(s) s$name, character(1)))
  names_ <- intersect(space$names, if (!is.null(design)) names(design) else character(0))
  if (is.null(design) || nrow(design) == 0 || length(names_) == 0) {
    return(data.frame(parameter = character(0), posterior_sd = numeric(0), prior_sd = numeric(0),
                      sd_ratio = numeric(0), max_abs_corr = numeric(0), identifiable = logical(0)))
  }
  d <- design
  if ("score" %in% names(d) && nrow(d) > 5) {
    k <- max(5L, as.integer(ceiling(behavioural_quantile * nrow(d))))
    d <- d[order(d$score), ][seq_len(k), ]
  }
  sub <- as.data.frame(lapply(d[names_], as.numeric))
  corr <- abs(suppressWarnings(cor(sub))); diag(corr) <- NA
  rows <- lapply(names_, function(n) {
    post_sd <- if (nrow(sub) > 1) sd(sub[[n]]) else NA_real_
    pri_sd <- if (!is.null(specs[[n]])) .prior_std(specs[[n]]) else NA_real_
    ratio <- if (!is.na(pri_sd) && pri_sd > 0) post_sd / pri_sd else NA_real_
    max_corr <- if (n %in% colnames(corr)) suppressWarnings(max(corr[, n], na.rm = TRUE)) else NA_real_
    if (is.infinite(max_corr)) max_corr <- NA_real_
    data.frame(parameter = n, posterior_sd = post_sd, prior_sd = pri_sd, sd_ratio = ratio,
               max_abs_corr = max_corr, identifiable = isTRUE(is.finite(ratio) && ratio < 0.6),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); out[order(out$sd_ratio), ]
}

#' Per-variable structural-adequacy check on the best fit.
#' Mirrors diagnostics.py:structural_adequacy.
#' @export
structural_adequacy <- function(result, ef_floor = 0.0, nrmse_ceiling = 50.0) {
  per_var <- result$best$per_var %||% list()
  rows <- lapply(names(per_var), function(var) {
    m <- per_var[[var]]; ef <- m$EF %||% NA_real_; nrmse <- m$nRMSE_pct %||% NA_real_
    reasons <- character(0)
    if (is.finite(ef) && ef < ef_floor) reasons <- c(reasons, sprintf("EF=%.2f<%s", ef, ef_floor))
    if (is.finite(nrmse) && nrmse > nrmse_ceiling) reasons <- c(reasons, sprintf("nRMSE=%.0f%%>%s%%", nrmse, nrmse_ceiling))
    data.frame(variable = var, EF = ef, nRMSE_pct = nrmse, n = m$n %||% NA_real_,
               flag = length(reasons) > 0, reason = paste(reasons, collapse = "; "), stringsAsFactors = FALSE)
  })
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# ---- scaffold a new crop --------------------------------------------------

.PHENOLOGY_HINT <- c("CSDL", "PPSEN", "EM-FL", "FL-SD", "SD-PM", "FL-SH", "FL-LF",
                     "PL-EM", "PLEM", "P1", "P2", "P3", "P4", "P5", "P1V", "P1D",
                     "PHINT", "EM-V1", "PHTHRS")

.scaffold_role <- function(name) if (toupper(name) %in% .PHENOLOGY_HINT) "obligatory" else "candidate"

#' Clone analog genotype files and emit a starter parameter block.
#' Mirrors scaffold.py:scaffold_crop.
#' @export
scaffold_crop <- function(dssat_dir, analog_stem, new_stem, new_code, source_anchor,
                          new_anchor = NULL, out_dir = NULL, spread = 0.3, copy_spe = TRUE) {
  geno <- file.path(dssat_dir, "Genotype")
  if (is.null(out_dir)) out_dir <- resolve_template_dir(required = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  exts <- c("CUL", "ECO", if (copy_spe) "SPE")
  files <- list()
  for (e in exts) {
    src <- file.path(geno, sprintf("%s.%s", analog_stem, e))
    if (file.exists(src)) { dst <- file.path(out_dir, sprintf("%s.%s", new_stem, e)); file.copy(src, dst, overwrite = TRUE); files[[e]] <- dst }
  }
  cul <- file.path(out_dir, sprintf("%s.CUL", new_stem))
  fmap <- cultivar_field_map(cul); starts <- read_cultivar_values(cul, source_anchor); bounds <- read_cul_calibration_bounds(cul)
  coeffs <- list()
  for (name in names(fmap)) {
    start <- starts[[name]]; if (is.null(start) || is.na(start)) next
    if (!is.null(bounds[[name]])) { lo <- bounds[[name]]$min; hi <- bounds[[name]]$max }
    else { lo <- start * (1 - spread); hi <- start * (1 + spread); if (lo > hi) { tmp <- lo; lo <- hi; hi <- tmp } }
    coeffs[[name]] <- list(min = round(lo, 4), max = round(hi, 4), start = round(start, 4), role = .scaffold_role(name))
  }
  list(files = files, out_dir = out_dir, parameters_yaml = .emit_parameters_yaml(coeffs),
       coefficients = coeffs, new_anchor = new_anchor %||% source_anchor)
}

.emit_parameters_yaml <- function(coeffs) {
  lines <- c("parameters:", "  genetic_cultivar:")
  for (name in names(coeffs)) {
    c <- coeffs[[name]]
    active <- if (c$role == "obligatory") "true" else "false"
    sd <- round((c$max - c$min) / 6.0, 4); if (sd == 0) sd <- 1.0
    lines <- c(lines, sprintf('    "%s": { active: %s, role: %s, min: %s, max: %s, start: %s, prior: {dist: normal, sd: %s} }',
                              name, active, c$role, c$min, c$max, c$start, sd))
  }
  paste0(paste(lines, collapse = "\n"), "\n")
}

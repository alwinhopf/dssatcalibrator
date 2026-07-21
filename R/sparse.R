# Sparse-data cultivar/species calibration helpers.
# R twin of python/dssatcalibrator/sparse.py plus the lightweight history and
# ABC-SMC engines. Everything here is opt-in and dependency-light.

.PHENOLOGY_HINT_SPARSE <- c("CSDL", "PPSEN", "EM-FL", "FL-SD", "SD-PM", "FL-SH",
                            "FL-LF", "PL-EM", "PLEM", "P1", "P2", "P3", "P4",
                            "P5", "P1V", "P1D", "PHINT", "EM-V1", "PHTHRS",
                            "ADAP", "MDAP")

.param_names_active <- function(cfg) vapply(active_parameters(cfg), function(s) s$name, character(1))

.set_normal_prior <- function(spec, center, sd, overwrite) {
  if (isTRUE(overwrite) || is.null(spec$prior)) {
    spec$prior <- list(dist = "normal", mean = as.numeric(center), sd = max(as.numeric(sd), 1e-9))
  }
  spec
}

.center_bounds <- function(spec, center, half_width) {
  lo0 <- as.numeric(spec$min); hi0 <- as.numeric(spec$max)
  lo <- max(lo0, center - half_width); hi <- min(hi0, center + half_width)
  if (lo >= hi) c(lo0, hi0) else c(lo, hi)
}

#' Apply analog-centered delta calibration priors/bounds.
#' @export
apply_delta_from_analog <- function(cfg) {
  out <- cfg
  dc <- .cfg_get(.cfg_get(out, "sparse", list()), "delta_from_analog", list())
  if (!isTRUE(.cfg_get(dc, "active", FALSE))) return(out)
  analog <- .cfg_get(dc, "theta", .cfg_get(dc, "analog_theta", list()))
  scale <- as.numeric(.cfg_get(dc, "relative_width", .cfg_get(dc, "scale", 0.20)))
  min_width <- as.numeric(.cfg_get(dc, "min_width", 1e-6))
  widths <- .cfg_get(dc, "widths", list())
  overwrite <- isTRUE(.cfg_get(dc, "overwrite_prior", TRUE))
  for (group in names(out$parameters)) {
    params <- out$parameters[[group]]; if (!is.list(params)) next
    for (name in names(params)) {
      if (is.null(analog[[name]])) next
      spec <- params[[name]]
      center <- as.numeric(analog[[name]])
      hw <- as.numeric(.cfg_get(widths, name, max(abs(center) * scale, min_width)))
      b <- .center_bounds(spec, center, hw)
      spec$min <- b[1]; spec$max <- b[2]
      spec$start <- min(max(center, b[1]), b[2])
      spec$analog <- center
      spec <- .set_normal_prior(spec, center, hw / 2.0, overwrite)
      out$parameters[[group]][[name]] <- spec
    }
  }
  out
}

#' Apply empirical-Bayes hierarchical priors from analog cultivars/species.
#' @export
apply_hierarchical_priors <- function(cfg) {
  out <- cfg
  hc <- .cfg_get(.cfg_get(out, "sparse", list()), "hierarchical_priors", list())
  if (!isTRUE(.cfg_get(hc, "active", FALSE))) return(out)
  params <- .cfg_get(hc, "parameters", list())
  pop <- .cfg_get(hc, "population", list())
  if (length(pop)) {
    nms <- unique(unlist(lapply(pop, names)))
    for (nm in nms) {
      vals <- as.numeric(unlist(lapply(pop, function(row) row[[nm]])))
      vals <- vals[is.finite(vals)]
      if (length(vals)) {
        if (is.null(params[[nm]])) params[[nm]] <- list()
        if (is.null(params[[nm]]$mean)) params[[nm]]$mean <- mean(vals)
        if (is.null(params[[nm]]$sd)) params[[nm]]$sd <- if (length(vals) > 1) sd(vals) else 0.1 * max(abs(mean(vals)), 1.0)
      }
    }
  }
  bounds_sd <- as.numeric(.cfg_get(hc, "bounds_sd", 3.0))
  overwrite <- isTRUE(.cfg_get(hc, "overwrite_prior", FALSE))
  for (group in names(out$parameters)) {
    params_group <- out$parameters[[group]]; if (!is.list(params_group)) next
    for (name in names(params_group)) {
      if (is.null(params[[name]])) next
      spec <- params_group[[name]]
      pc <- params[[name]]
      center <- as.numeric(.cfg_get(pc, "mean", .cfg_get(pc, "mu", .cfg_get(spec, "start", 0.5 * (spec$min + spec$max)))))
      sdev <- as.numeric(.cfg_get(pc, "sd", .cfg_get(pc, "sigma", max(0.1 * abs(center), 1.0))))
      if (isTRUE(.cfg_get(hc, "shrink_bounds", FALSE))) {
        b <- .center_bounds(spec, center, bounds_sd * sdev)
        spec$min <- b[1]; spec$max <- b[2]
      }
      if (isTRUE(.cfg_get(hc, "set_start", TRUE))) spec$start <- min(max(center, spec$min), spec$max)
      spec$hierarchical_mean <- center; spec$hierarchical_sd <- sdev
      spec <- .set_normal_prior(spec, center, sdev, overwrite)
      out$parameters[[group]][[name]] <- spec
    }
  }
  out
}

.trait_center <- function(name, spec, traits) {
  start <- as.numeric(.cfg_get(spec, "start", 0.5 * (spec$min + spec$max)))
  up <- toupper(name)
  if (!is.null(traits$maturity_days) && up %in% .PHENOLOGY_HINT_SPARSE) {
    delta <- as.numeric(traits$maturity_days) - 120.0
    strength <- if (up %in% c("SD-PM", "P5", "PHTHRS")) 0.35 else 0.15
    return(start + strength * delta)
  }
  if (!is.null(traits$photoperiod_sensitivity) && up %in% c("PPSEN", "CSDL", "P1D")) {
    return(start * (1.0 + 0.25 * as.numeric(traits$photoperiod_sensitivity)))
  }
  NULL
}

#' Apply trait-informed priors.
#' @export
apply_trait_priors <- function(cfg) {
  out <- cfg
  tc <- .cfg_get(.cfg_get(out, "sparse", list()), "trait_priors", list())
  if (!isTRUE(.cfg_get(tc, "active", FALSE))) return(out)
  traits <- .cfg_get(tc, "traits", .cfg_get(.cfg_get(out, "sparse", list()), "traits", list()))
  rules <- .cfg_get(tc, "rules", list())
  overwrite <- isTRUE(.cfg_get(tc, "overwrite_prior", FALSE))
  for (group in names(out$parameters)) {
    params <- out$parameters[[group]]; if (!is.list(params)) next
    for (name in names(params)) {
      spec <- params[[name]]
      center <- NULL; sdev <- NULL
      rule <- rules[[name]]
      if (!is.null(rule)) {
        if (!is.null(rule$value)) center <- as.numeric(rule$value)
        else if (!is.null(rule$trait) && !is.null(traits[[rule$trait]]))
          center <- as.numeric(.cfg_get(rule, "intercept", 0.0)) + as.numeric(.cfg_get(rule, "slope", 1.0)) * as.numeric(traits[[rule$trait]])
        sdev <- rule$sd
      }
      if (is.null(center)) center <- .trait_center(name, spec, traits)
      if (is.null(center)) next
      center <- min(max(center, spec$min), spec$max)
      span <- as.numeric(spec$max - spec$min)
      if (is.null(sdev)) sdev <- max(0.20 * span, 1e-6)
      if (isTRUE(.cfg_get(tc, "set_start", TRUE))) spec$start <- center
      spec$trait_prior_center <- center
      spec <- .set_normal_prior(spec, center, sdev, overwrite)
      out$parameters[[group]][[name]] <- spec
    }
  }
  out
}

#' Apply all opt-in sparse-data config transforms once.
#' @export
apply_sparse_config <- function(cfg) {
  if (isTRUE(.cfg_get(cfg, "_sparse_applied", FALSE))) return(cfg)
  out <- apply_delta_from_analog(cfg)
  out <- apply_hierarchical_priors(out)
  out <- apply_trait_priors(out)
  out[["_sparse_applied"]] <- TRUE
  out
}

#' Return a small robust first-pass calibration config.
#' @export
make_quick_dirty_config <- function(cfg, n = NULL, engine = "lhs", estimator = "glue") {
  out <- apply_sparse_config(cfg)
  if (is.null(n)) n <- 40L
  out$method$preset <- "C"
  out$method$sample <- list(engine = engine, n = as.integer(n))
  out$method$bayesian <- list(engine = estimator, behavioural_quantile = 0.25)
  out$method$staged <- NULL
  if (is.null(out$objective$likelihood)) out$objective$likelihood <- list(type = "huber", delta = 2.0)
  out
}

#' Run the quick first-pass calibration.
#' @export
calibrate_quick <- function(cfg, n = NULL, progress = TRUE) calibrate(make_quick_dirty_config(cfg, n = n), progress = progress)

.apply_active_subset_sparse <- function(cfg, keep) {
  out <- cfg; keep <- unlist(keep)
  for (group in names(out$parameters)) {
    params <- out$parameters[[group]]; if (!is.list(params)) next
    for (name in names(params)) if (is.list(params[[name]]) && isTRUE(params[[name]]$active))
      out$parameters[[group]][[name]]$active <- name %in% keep
  }
  out
}

.stage_keep <- function(cfg, stage) {
  active <- .param_names_active(cfg)
  keep <- unlist(.cfg_get(stage, "params", list()))
  groups <- unlist(.cfg_get(stage, "groups", list()))
  roles <- unlist(.cfg_get(stage, "roles", list()))
  for (s in active_parameters(cfg)) {
    if (s$group %in% groups || .cfg_get(s, "role", "") %in% roles) keep <- unique(c(keep, s$name))
  }
  lname <- tolower(as.character(.cfg_get(stage, "name", "")))
  if (length(keep) == 0 && grepl("phen", lname)) {
    keep <- vapply(active_parameters(cfg), function(s) {
      if (toupper(s$name) %in% .PHENOLOGY_HINT_SPARSE || identical(.cfg_get(s, "role", ""), "obligatory")) s$name else NA_character_
    }, character(1))
    keep <- keep[!is.na(keep)]
  }
  if (length(keep) == 0) keep <- active
  active[active %in% keep]
}

.focus_objective <- function(cfg, variables) {
  if (length(variables) == 0) return(cfg)
  out <- cfg
  known <- unique(c(names(.cfg_get(.cfg_get(out, "engine", list()), "timeseries_outputs", list())),
                    names(.cfg_get(.cfg_get(out, "engine", list()), "scalar_outputs", list()))))
  weights <- .cfg_get(.cfg_get(out, "objective", list()), "weights", list())
  for (v in known) weights[[v]] <- if (v %in% variables) 1.0 else 0.0
  out$objective$weights <- weights
  out
}

.sparse_sampler_engine <- function(engine) {
  engine <- tolower(as.character(engine))
  if (engine == "lhs" && !requireNamespace("lhs", quietly = TRUE)) "montecarlo" else engine
}

#' Default sparse-data stages.
#' @export
default_stages <- function(cfg) list(
  list(name = "phenology", roles = list("obligatory"), variables = list("phenology", "anthesis", "maturity")),
  list(name = "canopy_growth", roles = list("candidate"), variables = list("lai", "biomass")),
  list(name = "yield", groups = list("genetic_cultivar", "management"), variables = list("grain_yield", "yield"))
)

#' Build per-stage configs for staged sparse calibration.
#' @export
build_staged_configs <- function(cfg) {
  base <- apply_sparse_config(cfg)
  scfg <- .cfg_get(.cfg_get(base, "method", list()), "staged", list())
  stages <- .cfg_get(scfg, "stages", default_stages(base))
  out <- list()
  for (stage in stages) {
    keep <- .stage_keep(base, stage)
    st <- .apply_active_subset_sparse(base, keep)
    st <- .focus_objective(st, unlist(.cfg_get(stage, "variables", list())))
    st$method$staged <- list(active = FALSE)
    if (!is.null(stage$n)) st$method$sample$n <- as.integer(stage$n)
    if (!is.null(stage$engine)) st$method$bayesian$engine <- stage$engine
    out[[length(out) + 1L]] <- list(name = as.character(.cfg_get(stage, "name", paste0("stage_", length(out) + 1L))), cfg = st)
  }
  out
}

#' Run sequential staged calibration.
#' @export
calibrate_staged <- function(cfg, progress = TRUE) {
  stages <- build_staged_configs(cfg)
  results <- list(); prev_theta <- list()
  for (st in stages) {
    scfg <- st$cfg
    for (group in names(scfg$parameters)) {
      params <- scfg$parameters[[group]]; if (!is.list(params)) next
      for (name in names(params)) if (!is.null(prev_theta[[name]])) scfg$parameters[[group]][[name]]$start <- as.numeric(prev_theta[[name]])
    }
    if (isTRUE(progress)) message(sprintf("[staged] %s: %d active parameter(s)", st$name, length(.param_names_active(scfg))))
    res <- calibrate(scfg, progress = progress)
    results[[length(results) + 1L]] <- list(stage = st$name, result = res, active = .param_names_active(scfg))
    prev_theta <- modifyList(prev_theta, res$best_theta)
  }
  final <- results[[length(results)]]$result
  final$extras$stages <- results
  final$best_theta <- prev_theta
  final
}

.implausibility <- function(result, mode = "max_z") {
  resid <- result$residuals
  if (is.null(resid) || nrow(resid) == 0) return(Inf)
  if (mode == "score") return(if (is.finite(result$score)) sqrt(result$score) else Inf)
  max(abs(resid$resid / resid$sigma), na.rm = TRUE)
}

#' Run Bayesian history matching.
#' @export
run_history_matching <- function(cfg, score_results, space, progress = TRUE) {
  hcfg <- .cfg_get(.cfg_get(cfg, "method", list()), "bayesian", list())
  waves <- as.integer(.cfg_get(hcfg, "waves", 3)); n <- as.integer(.cfg_get(hcfg, "n", .cfg_get(hcfg, "n_per_wave", 128)))
  cutoff <- as.numeric(.cfg_get(hcfg, "implausibility_cutoff", 3.0))
  mode <- tolower(as.character(.cfg_get(hcfg, "implausibility", "max_z")))
  sampler <- .sparse_sampler_engine(.cfg_get(hcfg, "sampler", "lhs"))
  seed <- as.integer(.cfg_get(cfg$calibrator, "seed", 42))
  low <- space$low; high <- space$high
  rows <- list(); obj_results <- list(); sid <- 0L; wave_info <- list()
  for (wave in seq_len(waves) - 1L) {
    unit <- sample_design(space, n = n, engine = sampler, seed = seed + wave, include_start = FALSE)
    u <- sweep(sweep(as.matrix(unit), 2, space$low, `-`), 2, (space$high - space$low), `/`)
    vals <- sweep(sweep(u, 2, (high - low), `*`), 2, low, `+`)
    design <- as.data.frame(vals); names(design) <- space$names
    if (wave == 0L) {
      start_row <- as.data.frame(as.list(space$start), check.names = FALSE)
      names(start_row) <- space$names
      design <- rbind(start_row, design)
    }
    names(design) <- space$names
    thetas <- lapply(seq_len(nrow(design)), function(i) ps_to_theta(space, as.numeric(design[i, ])))
    results <- score_results(thetas)
    if (!any(vapply(results, function(r) is.finite(r$score), logical(1)))) {
      stop(paste(
        "History matching found no valid candidates: every evaluated score is non-finite.",
        "Inspect the spawn manifest and per-run DSSAT errors."
      ), call. = FALSE)
    }
    impl <- vapply(results, .implausibility, numeric(1), mode = mode)
    nroy <- impl <= cutoff
    if (!any(nroy) && any(is.finite(impl))) nroy <- impl <= as.numeric(quantile(impl[is.finite(impl)], 0.9, names = FALSE))
    for (i in seq_along(thetas)) {
      obj_results[[as.character(sid)]] <- results[[i]]
      rows[[length(rows) + 1L]] <- c(list(sample_id = sid, wave = wave), thetas[[i]],
                                     list(score = results[[i]]$score, loglik = results[[i]]$loglik,
                                          n_obs = nrow(results[[i]]$residuals),
                                          implausibility = impl[i], nroy = nroy[i]))
      sid <- sid + 1L
    }
    wave_info[[length(wave_info) + 1L]] <- list(wave = wave, n = nrow(design), nroy = sum(nroy), low = low, high = high)
    if (isTRUE(progress)) message(sprintf("  history wave %d/%d: NROY %d/%d", wave + 1L, waves, sum(nroy), nrow(design)))
    if (!any(nroy)) break
    nr <- design[nroy, , drop = FALSE]
    span <- pmax(apply(nr, 2, max) - apply(nr, 2, min), 1e-12)
    low <- pmax(space$low, apply(nr, 2, min) - 0.10 * span)
    high <- pmin(space$high, apply(nr, 2, max) + 0.10 * span)
  }
  design <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  for (cc in c("sample_id", "wave", "score", "loglik", "n_obs", "implausibility")) {
    if (cc %in% names(design)) design[[cc]] <- as.numeric(design[[cc]])
  }
  if ("nroy" %in% names(design)) design$nroy <- as.logical(design$nroy)
  design$weight <- 0.0
  behavioural <- design[as.logical(design$nroy), , drop = FALSE]
  if (nrow(behavioural) > 0) design$weight[as.logical(design$nroy)] <- 1.0 / nrow(behavioural)
  valid <- which(is.finite(design$score))
  if (!length(valid)) {
    stop(paste(
      "History matching found no valid candidates: every evaluated score is non-finite.",
      "Inspect the spawn manifest and per-run DSSAT errors."
    ), call. = FALSE)
  }
  best_sample_id <- design$sample_id[valid[which.min(design$score[valid])]]
  best_theta <- as.list(design[design$sample_id == best_sample_id, space$names, drop = FALSE][1, ])
  best <- obj_results[[as.character(best_sample_id)]]
  structure(list(design = design, behavioural = behavioural, best_theta = best_theta,
                 best_sample_id = best_sample_id, threshold = cutoff,
                 ess = nrow(behavioural), obj_results = obj_results, best = best,
                 waves = wave_info), class = "history_result")
}

#' Run ABC-SMC.
#' @export
run_abc_smc <- function(cfg, score_results, space, progress = TRUE) {
  bcfg <- .cfg_get(.cfg_get(cfg, "method", list()), "bayesian", list())
  n_particles <- as.integer(.cfg_get(bcfg, "n_particles", 128))
  waves <- as.integer(.cfg_get(bcfg, "waves", 4))
  oversample <- as.integer(.cfg_get(bcfg, "oversample", 3))
  q <- as.numeric(.cfg_get(bcfg, "threshold_quantile", 0.5))
  sampler <- .sparse_sampler_engine(.cfg_get(bcfg, "sampler", "lhs"))
  seed <- as.integer(.cfg_get(cfg$calibrator, "seed", 42))
  set.seed(seed)
  init <- sample_design(space, n = n_particles * oversample, engine = sampler, seed = seed, include_start = TRUE)
  accepted <- lapply(seq_len(nrow(init)), function(i) ps_to_theta(space, as.numeric(init[i, ])))
  accepted_scores <- rep(Inf, length(accepted))
  rows <- list(); obj_results <- list(); sid <- 0L; thresholds <- numeric(0)
  for (wave in seq_len(waves) - 1L) {
    if (wave == 0L) {
      candidates <- accepted
    } else {
      sdv <- pmax((space$high - space$low) * 0.05, 1e-9)
      candidates <- lapply(seq_len(n_particles * oversample), function(i) {
        parent <- accepted[[sample(seq_along(accepted), 1)]]
        vals <- pmin(pmax(as.numeric(unlist(parent[space$names])) + rnorm(length(space$names), 0, sdv), space$low), space$high)
        ps_to_theta(space, vals)
      })
    }
    results <- score_results(candidates)
    scores <- vapply(results, function(r) if (is.finite(r$score)) r$score else Inf, numeric(1))
    finite <- scores[is.finite(scores)]
    if (!length(finite)) {
      stop(paste(
        "ABC-SMC found no valid candidates: every evaluated score is non-finite.",
        "Inspect the spawn manifest and per-run DSSAT errors."
      ), call. = FALSE)
    }
    eps <- as.numeric(quantile(finite, q, names = FALSE))
    if (length(thresholds)) eps <- min(eps, tail(thresholds, 1))
    keep <- which(scores <= eps)
    if (length(keep) < max(8L, n_particles %/% 4L)) keep <- order(scores)[seq_len(min(length(scores), n_particles))]
    keep <- keep[order(scores[keep])][seq_len(min(length(keep), n_particles))]
    accepted <- candidates[keep]; accepted_scores <- scores[keep]; thresholds <- c(thresholds, eps)
    keep_set <- keep
    for (i in seq_along(candidates)) {
      obj_results[[as.character(sid)]] <- results[[i]]
      rows[[length(rows) + 1L]] <- c(list(sample_id = sid, wave = wave), candidates[[i]],
                                     list(score = results[[i]]$score, loglik = results[[i]]$loglik,
                                          n_obs = nrow(results[[i]]$residuals),
                                          threshold = eps, accepted = i %in% keep_set))
      sid <- sid + 1L
    }
    if (isTRUE(progress)) message(sprintf("  ABC wave %d/%d: accepted %d/%d", wave + 1L, waves, length(keep), length(candidates)))
  }
  design <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  for (cc in c("sample_id", "wave", "score", "loglik", "n_obs", "threshold")) {
    if (cc %in% names(design)) design[[cc]] <- as.numeric(design[[cc]])
  }
  if ("accepted" %in% names(design)) design$accepted <- as.logical(design$accepted)
  design$weight <- 0.0
  final <- design[design$wave == max(design$wave) & as.logical(design$accepted), , drop = FALSE]
  if (nrow(final) > 0) design$weight[design$wave == max(design$wave) & as.logical(design$accepted)] <- 1.0 / nrow(final)
  valid <- which(is.finite(design$score))
  if (!length(valid)) {
    stop(paste(
      "ABC-SMC found no valid candidates: every evaluated score is non-finite.",
      "Inspect the spawn manifest and per-run DSSAT errors."
    ), call. = FALSE)
  }
  best_sample_id <- design$sample_id[valid[which.min(design$score[valid])]]
  best_theta <- as.list(design[design$sample_id == best_sample_id, space$names, drop = FALSE][1, ])
  best <- obj_results[[as.character(best_sample_id)]]
  structure(list(design = design, behavioural = final, best_theta = best_theta,
                 best_sample_id = best_sample_id, threshold = tail(thresholds, 1),
                 ess = nrow(final), obj_results = obj_results, best = best,
                 thresholds = thresholds, initial_design = init), class = "abc_smc_result")
}

#' Freeze weakly identified parameters for a follow-up run.
#' @export
apply_identifiability_gate <- function(cfg, result, max_sd_ratio = NULL, protect_roles = c("obligatory")) {
  gate <- .cfg_get(.cfg_get(cfg, "sparse", list()), "identifiability_gate", list())
  if (is.null(max_sd_ratio)) max_sd_ratio <- as.numeric(.cfg_get(gate, "max_sd_ratio", 0.8))
  diag <- identifiability(result)
  weak <- if (nrow(diag)) diag$parameter[diag$sd_ratio >= max_sd_ratio] else character(0)
  out <- cfg; frozen <- character(0)
  for (group in names(out$parameters)) {
    params <- out$parameters[[group]]; if (!is.list(params)) next
    for (name in names(params)) if (name %in% weak && !(params[[name]]$role %in% protect_roles)) {
      out$parameters[[group]][[name]]$active <- FALSE; frozen <- c(frozen, name)
    }
  }
  list(cfg = out, diagnostics = diag, frozen = frozen)
}

#' Rank candidate next observations by a posterior-spread proxy.
#' @export
recommend_observations <- function(result, candidates = NULL) {
  cfg <- result$cfg
  if (is.null(candidates)) {
    vars <- sort(unique(c(names(.cfg_get(.cfg_get(cfg, "engine", list()), "timeseries_outputs", list())),
                          names(.cfg_get(.cfg_get(cfg, "engine", list()), "scalar_outputs", list())))))
    candidates <- lapply(vars, function(v) list(variable = v))
  }
  rows <- list()
  for (cand in candidates) {
    uv <- cand$variable; sims <- numeric(0); sigmas <- numeric(0)
    for (ores in result$obj_results) {
      resid <- ores$residuals
      if (is.null(resid) || nrow(resid) == 0 || !(uv %in% resid$user_var)) next
      g <- resid[resid$user_var == uv, , drop = FALSE]
      sims <- c(sims, mean(g$sim)); sigmas <- c(sigmas, mean(g$sigma))
    }
    utility <- if (length(sims) >= 2) var(sims) / max(mean(sigmas)^2, 1e-12) else 1.0
    rows[[length(rows) + 1L]] <- as.data.frame(c(cand, list(utility = utility, available_samples = length(sims))), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out[order(-out$utility), , drop = FALSE]
}

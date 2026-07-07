# Orchestration: evaluate a sampled design and drive a calibration engine.
# R twin of python/dssatcalibrator/orchestrator.py.
#
# Flow (preset C / GLUE):
#   config -> ParameterSpace -> sample design -> spawn (sample x experiment) in
#   parallel -> parse -> score -> GLUE post-process (weights, behavioural, best).

# Build the experiment/observation/treatment setup shared by every engine.
.setup <- function(cfg) {
  space <- parameter_space_from_config(cfg)
  crops <- .cfg_get(cfg, "crops", list())
  first_code <- if (length(crops)) (.cfg_get(crops[[1]], "code", "HM")) else "HM"
  crop <- crop_for(cfg, first_code)
  exe <- resolve_exe(cfg)
  specs <- active_parameters(cfg)
  hemp_dir <- cfg$source$hemp_dir
  run_root <- file.path(cfg$calibrator$workdir, cfg$calibrator$name)
  dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

  src_block <- .cfg_get(cfg, "source", list())
  if (!is.null(src_block$table)) {
    obs <- observations(src_block$table)
  } else if (length(.cfg_get(cfg, "observation_sources", list()))) {
    obs <- observations_from_sources(cfg, unlist(.cfg_get(cfg, "experiments", list())))
  } else {
    src <- .cfg_get(src_block, "observations", "dssat")
    obs <- if (identical(src, "dssat"))
      observations_from_dssat(hemp_dir, unlist(.cfg_get(cfg, "experiments", list())), crop_ext = crop$code)
    else observations_from_csv(src)
  }
  experiments <- intersect(unlist(.cfg_get(cfg, "experiments", list())), obs_experiments(obs))
  treatments <- setNames(lapply(experiments, function(e)
    parse_treatments(file.path(hemp_dir, sprintf("%s.%s", e, crop$filex_ext)))), experiments)

  if (isTRUE(.cfg_get(.cfg_get(cfg, "management_options", list()), "use_source_planting_date", FALSE)) &&
      is.null(cfg[["_planting_dates"]])) {
    cfg[["_planting_dates"]] <- planting_dates(obs)
  }
  list(space = space, crop = crop, exe = exe, specs = specs, run_root = run_root,
       obs = obs, experiments = experiments, treatments = treatments, cfg = cfg)
}

#' Load, fuse, and return observations from all configured active sources.
#' Mirrors Observations.from_sources.
#' @export
observations_from_sources <- function(cfg, experiments) {
  sources <- build_sources(cfg)
  if (length(sources) == 0) {
    crops <- .cfg_get(cfg, "crops", list())
    code <- if (length(crops)) .cfg_get(crops[[1]], "code", "HM") else "HM"
    crop <- crop_for(cfg, code)
    df <- observations_from_dssat(cfg$source$hemp_dir, experiments, crop_ext = crop$code)$table
    for (col in SCHEMA_EXTENDED) if (!(col %in% names(df))) {
      df[[col]] <- if (col == "source") "field_measurements" else if (col == "quality_flag") 0L else NA
    }
    return(observations(df[SCHEMA_EXTENDED]))
  }
  fuser <- observation_fuser(sources, cfg)
  frames <- list()
  for (exp in experiments) {
    df <- fuser_collect(fuser, exp, list(as.Date("1970-01-01"), as.Date("2099-12-31")))
    if (nrow(df) > 0) frames[[length(frames) + 1L]] <- df
  }
  observations(if (length(frames)) do.call(rbind, frames) else .empty_extended())
}

# Score many parameter sets in ONE parallel batch — the shared engine primitive.
#' @export
evaluate_thetas <- function(cfg, thetas, setup = NULL, n_workers = NULL, progress = FALSE) {
  if (is.null(setup)) setup <- .setup(cfg)
  if (is.null(n_workers)) n_workers <- resolve_cores(.cfg_get(cfg$calibrator, "num_cores", 0))
  cache <- .evaluation_cache_from_setup(cfg, setup$crop, setup$specs, setup$experiments,
                                        setup$treatments, setup$obs$table, setup$exe)
  out <- vector("list", length(thetas))
  miss_by_key <- list()
  miss_keys <- character(0)
  for (ti in seq_along(thetas)) {
    key <- if (isTRUE(cache$enabled)) .eval_cache_key(cache, thetas[[ti]], setup$experiments) else paste0("no-cache-", ti)
    cached <- if (isTRUE(cache$enabled)) .eval_cache_get(cache, key) else NULL
    if (!is.null(cached)) {
      out[[ti]] <- cached
    } else {
      if (is.null(miss_by_key[[key]])) {
        miss_by_key[[key]] <- list(theta = thetas[[ti]], indices = integer(0))
        miss_keys <- c(miss_keys, key)
      }
      miss_by_key[[key]]$indices <- c(miss_by_key[[key]]$indices, ti)
    }
  }

  jobs <- list(); idx <- list()
  for (mi in seq_along(miss_keys)) {
    theta <- miss_by_key[[miss_keys[[mi]]]]$theta
    for (exp in setup$experiments) {
      jobs[[length(jobs) + 1L]] <- list(theta = theta, exp_id = exp, cfg = cfg, crop = setup$crop,
                                        param_specs = setup$specs, run_root = setup$run_root,
                                        treatments = setup$treatments[[exp]], exe = setup$exe)
      idx[[length(idx) + 1L]] <- list(mi = mi, exp = exp)
    }
  }
  results <- if (length(jobs)) run_many(jobs, n_workers = n_workers) else list()
  per_miss <- lapply(seq_along(miss_keys), function(i) list())
  for (k in seq_along(idx)) {
    mi <- idx[[k]]$mi; exp <- idx[[k]]$exp
    per_miss[[mi]][[exp]] <- results[[k]]
  }
  for (mi in seq_along(miss_keys)) {
    key <- miss_keys[[mi]]
    scored <- score(per_miss[[mi]], setup$obs$table, cfg)
    if (isTRUE(cache$enabled)) .eval_cache_put(cache, key, scored)
    for (ti in miss_by_key[[key]]$indices) out[[ti]] <- scored
  }
  list(results = out, setup = setup)
}

# Run every (sample x experiment) spawn, score per sample.
#' @export
evaluate_design <- function(cfg, samples, progress = TRUE) {
  setup <- .setup(cfg)
  n_workers <- resolve_cores(.cfg_get(cfg$calibrator, "num_cores", 0))
  rows <- list(); obj_results <- list()
  best_score <- Inf; best_obj_res <- NULL; best_sid <- NULL
  cache <- .evaluation_cache_from_setup(cfg, setup$crop, setup$specs, setup$experiments,
                                        setup$treatments, setup$obs$table, setup$exe)
  miss_by_key <- list(); miss_keys <- character(0)
  record_objective <- function(sid, theta, o) {
    if (o$score < best_score) {
      best_score <<- o$score; best_obj_res <<- o; best_sid <<- as.integer(sid)
    }
    rows[[length(rows) + 1L]] <<- c(list(sample_id = as.integer(sid)), theta,
                                    list(score = o$score, loglik = o$loglik, n_obs = nrow(o$residuals)))
  }
  for (sid in seq_len(nrow(samples)) - 1L) {
    theta <- ps_to_theta(setup$space, as.numeric(samples[sid + 1L, ]))
    key <- if (isTRUE(cache$enabled)) .eval_cache_key(cache, theta, setup$experiments) else paste0("no-cache-", sid)
    cached <- if (isTRUE(cache$enabled)) .eval_cache_get(cache, key) else NULL
    if (!is.null(cached)) {
      record_objective(sid, theta, cached)
    } else {
      if (is.null(miss_by_key[[key]])) {
        miss_by_key[[key]] <- list(theta = theta, sample_ids = integer(0))
        miss_keys <- c(miss_keys, key)
      }
      miss_by_key[[key]]$sample_ids <- c(miss_by_key[[key]]$sample_ids, sid)
    }
  }

  jobs <- list(); idx <- list()
  for (mi in seq_along(miss_keys)) {
    theta <- miss_by_key[[miss_keys[[mi]]]]$theta
    for (exp in setup$experiments) {
      jobs[[length(jobs) + 1L]] <- list(theta = theta, exp_id = exp, cfg = cfg, crop = setup$crop,
                                        param_specs = setup$specs, run_root = setup$run_root,
                                        treatments = setup$treatments[[exp]], exe = setup$exe)
      idx[[length(idx) + 1L]] <- list(mi = mi, exp = exp)
    }
  }
  results <- if (length(jobs)) run_many(jobs, n_workers = n_workers) else list()
  per_miss <- list()
  for (k in seq_along(idx)) {
    mi <- as.character(idx[[k]]$mi); exp <- idx[[k]]$exp
    if (is.null(per_miss[[mi]])) per_miss[[mi]] <- list()
    per_miss[[mi]][[exp]] <- results[[k]]
  }
  for (mi_name in names(per_miss)) {
    mi <- as.integer(mi_name)
    key <- miss_keys[[mi]]
    o <- score(per_miss[[mi_name]], setup$obs$table, cfg)
    if (isTRUE(cache$enabled)) .eval_cache_put(cache, key, o)
    for (sid in miss_by_key[[key]]$sample_ids) {
      theta <- ps_to_theta(setup$space, as.numeric(samples[as.integer(sid) + 1L, ]))
      record_objective(sid, theta, o)
    }
  }
  if (!is.null(best_sid)) obj_results[[as.character(best_sid)]] <- best_obj_res
  design <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE, check.names = FALSE)))
  design <- design[order(design$sample_id), ]; rownames(design) <- NULL
  list(design = design, obj_results = obj_results, space = setup$space, obs = setup$obs, experiments = setup$experiments)
}

# Preset -> per-stage engine defaults (CONCEPT.md §14a).
.PRESETS <- list(
  A = list(sample = list(engine = "lhs"), bayesian = list(engine = "smc_pf"), sensitivity = list(engine = "morris")),
  B = list(optimizer = list(engine = "diffevo"), bayesian = list(engine = "none"), sensitivity = list(engine = "morris")),
  C = list(sample = list(engine = "lhs"), bayesian = list(engine = "glue")),
  D = list(sample = list(engine = "sobol"), bayesian = list(engine = "mcmc"), sensitivity = list(engine = "morris"))
)

.resolve_method <- function(cfg) {
  method <- .cfg_get(cfg, "method", list())
  preset <- toupper(as.character(.cfg_get(method, "preset", "C")))
  if (!is.null(.PRESETS[[preset]])) {
    for (stage in names(.PRESETS[[preset]])) {
      method[[stage]] <- modifyList(.PRESETS[[preset]][[stage]], .cfg_get(method, stage, list()))
    }
  }
  if (is.null(method$sample)) method$sample <- list(engine = "lhs", n = 200L)
  if (is.null(method$bayesian)) method$bayesian <- list(engine = "glue")
  method
}

.stage_on <- function(block) isTRUE(.cfg_get(block, "active", FALSE))

.apply_staging <- function(cfg) {
  st <- .cfg_get(.cfg_get(cfg, "method", list()), "staging", list())
  fg <- unlist(.cfg_get(st, "freeze_groups", list())); fp <- unlist(.cfg_get(st, "freeze_params", list()))
  if (length(fg) == 0 && length(fp) == 0) return(cfg)
  for (group in names(cfg$parameters)) {
    params <- cfg$parameters[[group]]
    if (!is.list(params)) next
    for (name in names(params)) {
      spec <- params[[name]]
      if (is.list(spec) && isTRUE(spec$active) && (group %in% fg || name %in% fp)) {
        cfg$parameters[[group]][[name]]$active <- FALSE
      }
    }
  }
  cfg
}

.apply_active_subset <- function(cfg, keep) {
  keep <- unlist(keep)
  for (group in names(cfg$parameters)) {
    params <- cfg$parameters[[group]]
    if (!is.list(params)) next
    for (name in names(params)) {
      spec <- params[[name]]
      if (is.list(spec) && isTRUE(spec$active)) cfg$parameters[[group]][[name]]$active <- name %in% keep
    }
  }
  cfg
}

.results_scorer <- function(cfg, setup, n_workers) {
  function(thetas) evaluate_thetas(cfg, thetas, setup = setup, n_workers = n_workers)$results
}

# --- Estimator registry (mirrors orchestrator.py) ---------------------------
# Each estimator has one signature and returns a calibration_result, so adding a
# main estimator is a function + one registry entry (no edit to `calibrate`).
.OPTIMIZER_ALIASES <- c("nelder_mead", "diffevo", "neldermead", "nm", "de",
                        "cmaes", "cma_es", "cma")

# Pick the estimator key from the resolved `method` block. Precedence:
# surrogate accelerator > explicit bayesian engine > optimiser > GLUE default.
.resolve_estimator <- function(method) {
  if (.stage_on(method$surrogate)) return("surrogate")
  bayes <- tolower(as.character(.cfg_get(method$bayesian, "engine", "glue")))
  if (bayes %in% names(.ESTIMATOR_REGISTRY) && !(bayes %in% c("none", "", "surrogate", "optimizer")))
    return(bayes)
  opt <- tolower(as.character(.cfg_get(method$optimizer, "engine", "none")))
  if (bayes %in% c("none", "") && opt %in% .OPTIMIZER_ALIASES) return("optimizer")
  "glue"
}

.estimate_surrogate <- function(work_cfg, space, setup, method, seed, n_workers, extras, progress) {
  obs <- setup$obs; experiments <- setup$experiments
  scorer <- .results_scorer(work_cfg, setup, n_workers)
  sur <- run_surrogate(work_cfg, space, scorer, progress = progress)
  glue <- run_glue(sur$design, space$names, work_cfg, space = space)
  best <- sur$obj_results[[as.character(glue$best_sample_id)]]
  extras$surrogate_info <- sur$info; extras$engine <- "surrogate"
  .calib_result(work_cfg, space, obs, experiments, glue$design, sur$obj_results, glue$best_theta, best, glue, extras)
}

.estimate_smc_pf <- function(work_cfg, space, setup, method, seed, n_workers, extras, progress) {
  obs <- setup$obs; experiments <- setup$experiments
  smc <- run_smc_pf(work_cfg, progress = progress)
  extras$initial_design <- smc$initial_design; extras$engine <- "smc_pf"
  .calib_result(work_cfg, space, obs, experiments, smc$design, smc$obj_results, smc$best_theta, smc$best, smc, extras)
}

.estimate_mcmc <- function(work_cfg, space, setup, method, seed, n_workers, extras, progress) {
  obs <- setup$obs; experiments <- setup$experiments
  scorer <- .results_scorer(work_cfg, setup, n_workers)
  mc <- run_mcmc(work_cfg, scorer, space, progress = progress)
  extras$initial_design <- mc$initial_design; extras$engine <- "mcmc"; extras$acceptance <- mc$acceptance
  .calib_result(work_cfg, space, obs, experiments, mc$design, mc$obj_results, mc$best_theta, mc$best, mc, extras)
}

.estimate_optimizer <- function(work_cfg, space, setup, method, seed, n_workers, extras, progress) {
  obs <- setup$obs; experiments <- setup$experiments
  ocfg <- method$optimizer
  opt <- tolower(as.character(.cfg_get(ocfg, "engine", "none")))
  scorer <- .results_scorer(work_cfg, setup, n_workers)
  score_batch <- function(ths) vapply(scorer(ths), function(r) r$score, numeric(1))
  ores <- run_optimizer(space, score_batch, method = opt, seed = seed, maxiter = ocfg$maxiter,
                        popsize = as.integer(.cfg_get(ocfg, "popsize", 15)),
                        restarts = as.integer(.cfg_get(ocfg, "restarts", 1)), progress = progress)
  best <- scorer(list(ores$best_theta))[[1]]
  design <- as.data.frame(c(list(sample_id = 0L), ores$best_theta,
                            list(score = best$score, loglik = best$loglik, n_obs = nrow(best$residuals), weight = 1.0)),
                          stringsAsFactors = FALSE, check.names = FALSE)
  extras$engine <- "optimizer"; extras$optimizer_history <- ores$history; extras$n_eval <- ores$n_eval
  .calib_result(work_cfg, space, obs, experiments, design, list(`0` = best), ores$best_theta, best, NULL, extras)
}

.estimate_glue <- function(work_cfg, space, setup, method, seed, n_workers, extras, progress) {
  n <- as.integer(.cfg_get(method$sample, "n", 200))
  engine <- .cfg_get(method$sample, "engine", "lhs")
  samples <- sample_design(space, n = n, engine = engine, seed = seed, include_start = TRUE)
  ed <- evaluate_design(work_cfg, samples, progress = progress)
  glue <- run_glue(ed$design, ed$space$names, work_cfg, space = ed$space)
  best <- ed$obj_results[[as.character(glue$best_sample_id)]]
  extras$engine <- "glue"
  .calib_result(work_cfg, ed$space, ed$obs, ed$experiments, glue$design, ed$obj_results, glue$best_theta, best, glue, extras)
}

.estimate_dream <- function(work_cfg, space, setup, method, seed, n_workers, extras, progress) {
  obs <- setup$obs; experiments <- setup$experiments
  scorer <- .results_scorer(work_cfg, setup, n_workers)
  mc <- run_dream(work_cfg, scorer, space, progress = progress)
  extras$initial_design <- mc$initial_design; extras$engine <- "dream"; extras$acceptance <- mc$acceptance
  .calib_result(work_cfg, space, obs, experiments, mc$design, mc$obj_results, mc$best_theta, mc$best, mc, extras)
}

.estimate_es_mda <- function(work_cfg, space, setup, method, seed, n_workers, extras, progress) {
  obs <- setup$obs; experiments <- setup$experiments
  scorer <- .results_scorer(work_cfg, setup, n_workers)
  es <- run_es_mda(work_cfg, scorer, space, progress = progress)
  extras$initial_design <- es$initial_design; extras$engine <- "es_mda"
  .calib_result(work_cfg, space, obs, experiments, es$design, es$obj_results, es$best_theta, es$best, es, extras)
}

.estimate_bayesopt <- function(work_cfg, space, setup, method, seed, n_workers, extras, progress) {
  obs <- setup$obs; experiments <- setup$experiments
  scorer <- .results_scorer(work_cfg, setup, n_workers)
  bo <- run_bayesopt(work_cfg, scorer, space, progress = progress)
  extras$engine <- "bayesopt"; extras$bayesopt_info <- bo$info
  .calib_result(work_cfg, space, obs, experiments, bo$design, bo$obj_results, bo$best_theta, bo$best, NULL, extras)
}

.estimate_history <- function(work_cfg, space, setup, method, seed, n_workers, extras, progress) {
  obs <- setup$obs; experiments <- setup$experiments
  scorer <- .results_scorer(work_cfg, setup, n_workers)
  hm <- run_history_matching(work_cfg, scorer, space, progress = progress)
  extras$engine <- "history"; extras$history_waves <- hm$waves
  .calib_result(work_cfg, space, obs, experiments, hm$design, hm$obj_results, hm$best_theta, hm$best, hm, extras)
}

.estimate_abc_smc <- function(work_cfg, space, setup, method, seed, n_workers, extras, progress) {
  obs <- setup$obs; experiments <- setup$experiments
  scorer <- .results_scorer(work_cfg, setup, n_workers)
  abc <- run_abc_smc(work_cfg, scorer, space, progress = progress)
  extras$engine <- "abc_smc"; extras$thresholds <- abc$thresholds; extras$initial_design <- abc$initial_design
  .calib_result(work_cfg, space, obs, experiments, abc$design, abc$obj_results, abc$best_theta, abc$best, abc, extras)
}

.ESTIMATOR_REGISTRY <- list(
  surrogate = .estimate_surrogate,
  smc_pf    = .estimate_smc_pf,
  mcmc      = .estimate_mcmc,
  dream     = .estimate_dream,
  es_mda    = .estimate_es_mda,
  bayesopt  = .estimate_bayesopt,
  history   = .estimate_history,
  abc_smc   = .estimate_abc_smc,
  optimizer = .estimate_optimizer,
  glue      = .estimate_glue
)

#' Run the configured calibration pipeline. Mirrors orchestrator.py:calibrate.
#' Steps: resolve preset -> [sensitivity] -> [select] -> estimate (glue|smc_pf|
#' mcmc|optimizer|surrogate) -> [NSGA-II add-on].
#' @export
calibrate <- function(cfg, progress = TRUE) {
  if (!isTRUE(.cfg_get(cfg, "_sparse_applied", FALSE))) cfg <- apply_sparse_config(cfg)
  if (isTRUE(.cfg_get(.cfg_get(.cfg_get(cfg, "method", list()), "staged", list()), "active", FALSE))) {
    return(calibrate_staged(cfg, progress = progress))
  }
  cfg <- .apply_staging(cfg)
  method <- .resolve_method(cfg)
  seed <- as.integer(.cfg_get(cfg$calibrator, "seed", 42))
  n_workers <- resolve_cores(.cfg_get(cfg$calibrator, "num_cores", 0))
  extras <- list()
  work_cfg <- cfg
  setup <- .setup(work_cfg); space <- setup$space

  sens_block <- method$sensitivity
  if (.stage_on(sens_block)) {
    scorer <- .results_scorer(work_cfg, setup, n_workers)
    sens <- run_sensitivity(space, scorer, method = .cfg_get(sens_block, "engine", "morris"),
                            trajectories = .cfg_get(sens_block, "trajectories", 10),
                            seed = seed)
    extras$sensitivity <- sens$ranking
    if (isTRUE(.cfg_get(sens_block, "auto_activate", FALSE))) {
      keep <- influential_params(sens$ranking, keep = sens_block$keep,
                                 rel_threshold = as.numeric(.cfg_get(sens_block, "rel_threshold", 0.1)))
      work_cfg <- .apply_active_subset(work_cfg, keep); setup <- .setup(work_cfg); space <- setup$space
    }
  }

  sel_block <- method$select
  if (.stage_on(sel_block)) {
    crit <- if (grepl("aicc", tolower(as.character(.cfg_get(sel_block, "engine", ""))))) "aicc" else "bic"
    scorer <- .results_scorer(work_cfg, setup, n_workers)
    sel <- stepwise_select(space, scorer, criterion = crit,
                           optimizer = .cfg_get(sel_block, "optimizer", "nelder_mead"),
                           optimizer_restarts = as.integer(.cfg_get(sel_block, "restarts", 2)),
                           maxiter = sel_block$maxiter, seed = seed)
    extras$selection <- sel
    work_cfg <- .apply_active_subset(work_cfg, sel$selected); setup <- .setup(work_cfg); space <- setup$space
  }

  if (tolower(as.character(.cfg_get(.cfg_get(cfg, "objective", list()), "weighting", ""))) == "agmip_wls") {
    work_cfg <- .agmip_reweight(work_cfg, setup, n_workers, seed, progress)
  }

  estimator <- .resolve_estimator(method)
  result <- .ESTIMATOR_REGISTRY[[estimator]](work_cfg, space, setup, method,
                                             seed, n_workers, extras, progress)

  mo <- method$multiobjective
  if (!is.null(mo) && identical(.cfg_get(mo, "engine", NULL), "nsga2")) {
    obj_vars <- .cfg_get(mo, "variables", sort(names(result$best$per_var)))
    eval_batch <- function(thetas) {
      res <- evaluate_thetas(work_cfg, thetas, setup = setup, n_workers = n_workers)$results
      lapply(res, function(r) setNames(lapply(obj_vars, function(v) (r$per_var[[v]]$nRMSE_pct) %||% 1e6), obj_vars))
    }
    result$nsga2 <- run_nsga2(eval_batch, space, obj_vars,
                              pop_size = as.integer(.cfg_get(mo, "pop_size", 16)),
                              n_gen = as.integer(.cfg_get(mo, "n_gen", 5)), seed = seed)
  }
  result
}

.calib_result <- function(cfg, space, obs, experiments, design, obj_results, best_theta, best, glue, extras) {
  structure(list(cfg = cfg, space = space, obs = obs, experiments = experiments, design = design,
                 obj_results = obj_results, best_theta = best_theta, best = best, glue = glue,
                 nsga2 = NULL, extras = extras), class = "calibration_result")
}

.agmip_reweight <- function(cfg, setup, n_workers, seed, progress) {
  space <- setup$space
  samples <- sample_design(space, n = as.integer(.cfg_get(cfg$calibrator, "wls_probe_n", 40)),
                           engine = "lhs", seed = seed, include_start = TRUE)
  thetas <- lapply(seq_len(nrow(samples)), function(i) ps_to_theta(space, as.numeric(samples[i, ])))
  results <- evaluate_thetas(cfg, thetas, setup = setup, n_workers = n_workers)$results
  scores <- vapply(results, function(r) if (is.finite(r$score)) r$score else Inf, numeric(1))
  best <- results[[which.min(scores)]]
  weights <- .cfg_get(.cfg_get(cfg, "objective", list()), "weights", list())
  if (nrow(best$residuals) > 0) {
    for (uv in unique(best$residuals$user_var)) {
      g <- best$residuals[best$residuals$user_var == uv, ]
      v <- if (nrow(g) > 1) var(g$resid) else g$resid[1]^2
      weights[[uv]] <- if (v > 0) 1.0 / v else 1.0
    }
  }
  if (is.null(cfg$objective)) cfg$objective <- list()
  cfg$objective$weights <- weights; cfg$objective$weighting <- "unified"
  cfg
}

# ---- SMC particle filter (needs _setup/run_many; lives here) --------------

#' SMC particle filter + Metropolis-Hastings move (preset A). Mirrors smc_pf.py.
#' @export
run_smc_pf <- function(cfg, progress = TRUE) {
  setup <- .setup(cfg)
  space <- setup$space; crop <- setup$crop; exe <- setup$exe; specs <- setup$specs
  run_root <- setup$run_root; obs <- setup$obs; experiments <- setup$experiments; treatments <- setup$treatments
  n_workers <- resolve_cores(.cfg_get(cfg$calibrator, "num_cores", 0))
  method <- .cfg_get(cfg, "method", list()); bcfg <- .cfg_get(method, "bayesian", list())
  n_particles <- as.integer(.cfg_get(bcfg, "n_particles", 200))
  ess_frac <- as.numeric(.cfg_get(bcfg, "ess_frac", 0.5))
  mutation_scale <- as.numeric(.cfg_get(bcfg, "mutation_scale", 0.02))
  move_kernel <- as.character(.cfg_get(bcfg, "move_kernel", "adaptive"))
  kernel_floor <- as.numeric(.cfg_get(bcfg, "kernel_floor", 0.01))
  kernel_scale <- .cfg_get(bcfg, "kernel_scale", NULL)
  seed <- as.integer(.cfg_get(cfg$calibrator, "seed", 42)); set.seed(seed)

  sample_engine <- .cfg_get(.cfg_get(method, "sample", list()), "engine", "lhs")
  if (has_informative_prior(space)) {
    prior_draws <- sample_prior_design(space, n_particles)
    start_row <- as.data.frame(as.list(space$start), check.names = FALSE); names(start_row) <- space$names
    samples <- rbind(start_row, prior_draws)
  } else {
    samples <- sample_design(space, n = n_particles, engine = sample_engine, seed = seed, include_start = TRUE)
  }
  n_particles <- nrow(samples)
  initial_design <- do.call(rbind, lapply(seq_len(n_particles) - 1L, function(sid)
    as.data.frame(c(list(sample_id = sid), ps_to_theta(space, as.numeric(samples[sid + 1L, ]))), stringsAsFactors = FALSE, check.names = FALSE)))

  jobs <- list(); idx <- list()
  for (sid in seq_len(n_particles) - 1L) {
    theta <- ps_to_theta(space, as.numeric(samples[sid + 1L, ]))
    for (exp in experiments) {
      jobs[[length(jobs) + 1L]] <- list(theta = theta, exp_id = exp, cfg = cfg, crop = crop,
                                        param_specs = specs, run_root = run_root, treatments = treatments[[exp]], exe = exe)
      idx[[length(idx) + 1L]] <- list(sid = sid, exp = exp)
    }
  }
  results <- run_many(jobs, n_workers = n_workers)
  particles <- vector("list", n_particles)
  for (sid in seq_len(n_particles) - 1L) {
    p_results <- list()
    for (k in seq_along(idx)) if (idx[[k]]$sid == sid) p_results[[idx[[k]]$exp]] <- results[[k]]
    theta <- ps_to_theta(space, as.numeric(samples[sid + 1L, ]))
    resid_df <- build_residuals(p_results, obs$table, cfg)
    particles[[sid + 1L]] <- list(theta = theta, results = p_results, residuals = resid_df,
                                  loglik = if (nrow(resid_df) == 0) -1e10 else 0.0)
  }

  unique_dates <- character(0); scalar_vars <- character(0)
  for (p in particles) {
    if (nrow(p$residuals) == 0) next
    for (i in seq_len(nrow(p$residuals))) {
      r <- p$residuals[i, ]
      if (is.na(r$date)) scalar_vars <- union(scalar_vars, r$user_var)
      else unique_dates <- union(unique_dates, as.character(r$date))
    }
  }
  steps <- c(
    lapply(sort(unique_dates), function(d) list(date = as.Date(d), var = NULL, label = d)),
    lapply(sort(scalar_vars), function(v) list(date = as.Date(NA), var = v, label = paste0("scalar:", v)))
  )
  loglik_contrib <- function(resid_df, step) {
    if (nrow(resid_df) == 0) return(-1e10)
    mask <- if (is.null(step$var)) (!is.na(resid_df$date) & resid_df$date == step$date)
            else (is.na(resid_df$date) & resid_df$user_var == step$var)
    sub <- resid_df[mask, , drop = FALSE]
    if (nrow(sub) == 0) return(0.0)
    -0.5 * sum(((sub$resid / sub$sigma)^2) * sub$weight)
  }
  loglik_accum <- vapply(particles, function(p) p$loglik, numeric(1))
  n_dim <- length(space$names); default_c <- 2.38 / sqrt(max(n_dim, 1))
  weights <- rep(1 / n_particles, n_particles); ess_trace <- list()

  move_sd <- function() {
    sd <- setNames(numeric(length(space$names)), space$names)
    if (move_kernel == "adaptive") {
      mat <- do.call(rbind, lapply(particles, function(p) as.numeric(unlist(p$theta[space$names]))))
      wmean <- apply(mat, 2, function(col) sum(col * weights) / sum(weights))
      wstd <- sqrt(vapply(seq_len(ncol(mat)), function(j) sum(weights * (mat[, j] - wmean[j])^2) / sum(weights), numeric(1)))
      c_ <- if (!is.null(kernel_scale)) as.numeric(kernel_scale) else default_c
      for (j in seq_along(space$names)) sd[space$names[j]] <- max(c_ * wstd[j], kernel_floor * (space$high[j] - space$low[j]))
    } else {
      for (j in seq_along(space$names)) sd[space$names[j]] <- mutation_scale * (space$high[j] - space$low[j])
    }
    sd
  }
  perturb <- function(parent, sd) {
    mutated <- parent
    for (name in names(parent)) {
      spec <- Find(function(s) s$name == name, space$specs)
      lo <- as.numeric(spec$min); hi <- as.numeric(spec$max)
      mutated[[name]] <- min(max(parent[[name]] + rnorm(1, 0, sd[[name]]), lo), hi)
    }
    mutated
  }

  for (step_idx in seq_along(steps)) {
    step <- steps[[step_idx]]
    for (i in seq_len(n_particles)) if (loglik_accum[i] > -1e9) loglik_accum[i] <- loglik_accum[i] + loglik_contrib(particles[[i]]$residuals, step)
    max_ll <- max(loglik_accum)
    if (max_ll == -1e10 || !is.finite(max_ll)) weights <- rep(1 / n_particles, n_particles)
    else { w <- exp(loglik_accum - max_ll); w[!is.finite(w)] <- 0; total <- sum(w); weights <- if (total > 0) w / total else rep(1 / n_particles, n_particles) }
    ess <- 1.0 / sum(weights^2)
    ess_trace[[length(ess_trace) + 1L]] <- list(step = step_idx, label = step$label, ess = ess, n = n_particles, resampled = FALSE)
    if (ess < n_particles * ess_frac && step_idx < length(steps)) {
      ess_trace[[length(ess_trace)]]$resampled <- TRUE
      sd <- move_sd(); resampled_idx <- systematic_resample(weights)
      mutated_thetas <- lapply(seq_len(n_particles), function(i) perturb(particles[[resampled_idx[i]]]$theta, sd))
      mjobs <- list(); midx <- list()
      for (i in seq_len(n_particles)) for (exp in experiments) {
        mjobs[[length(mjobs) + 1L]] <- list(theta = mutated_thetas[[i]], exp_id = exp, cfg = cfg, crop = crop,
                                            param_specs = specs, run_root = run_root, treatments = treatments[[exp]], exe = exe)
        midx[[length(midx) + 1L]] <- list(i = i, exp = exp)
      }
      mres <- run_many(mjobs, n_workers = n_workers)
      new_particles <- vector("list", n_particles)
      for (i in seq_len(n_particles)) {
        p_mut <- list()
        for (k in seq_along(midx)) if (midx[[k]]$i == i) p_mut[[midx[[k]]$exp]] <- mres[[k]]
        succeeded <- all(vapply(p_mut, function(r) r$status %in% c("success", "cached"), logical(1)))
        accepted <- FALSE
        if (succeeded) {
          mut_resids <- build_residuals(p_mut, obs$table, cfg)
          ll_mut <- 0; ll_par <- 0
          for (s_idx in seq_len(step_idx)) { ll_mut <- ll_mut + loglik_contrib(mut_resids, steps[[s_idx]]); ll_par <- ll_par + loglik_contrib(particles[[resampled_idx[i]]]$residuals, steps[[s_idx]]) }
          lp_mut <- log_prior_vec(space, mutated_thetas[[i]]); lp_par <- log_prior_vec(space, particles[[resampled_idx[i]]]$theta)
          alpha <- exp((ll_mut + lp_mut) - (ll_par + lp_par))
          if (!is.nan(alpha) && runif(1) < alpha) {
            accepted <- TRUE
            new_particles[[i]] <- list(theta = mutated_thetas[[i]], results = p_mut, residuals = mut_resids, loglik = ll_mut)
          }
        }
        if (!accepted) {
          par <- particles[[resampled_idx[i]]]
          new_particles[[i]] <- list(theta = par$theta, results = par$results, residuals = par$residuals, loglik = loglik_accum[resampled_idx[i]])
        }
      }
      particles <- new_particles; loglik_accum <- rep(0, n_particles)
    }
  }

  rows <- list(); obj_results <- list()
  for (i in seq_len(n_particles)) {
    o <- score(particles[[i]]$results, obs$table, cfg); obj_results[[as.character(i - 1L)]] <- o
    rows[[length(rows) + 1L]] <- c(list(sample_id = i - 1L), particles[[i]]$theta,
                                   list(score = o$score, loglik = o$loglik, n_obs = nrow(o$residuals)))
  }
  design <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE, check.names = FALSE)))
  design$weight <- weights
  best_sample_id <- which.min(ifelse(is.finite(design$score), design$score, Inf)) - 1L
  best_theta <- particles[[best_sample_id + 1L]]$theta; best <- obj_results[[as.character(best_sample_id)]]
  q <- as.numeric(.cfg_get(bcfg, "behavioural_quantile", 0.1))
  valid <- design$score[is.finite(design$score)]
  threshold <- if (length(valid)) as.numeric(quantile(valid, q, names = FALSE, type = 7)) else Inf
  behavioural <- design[is.finite(design$score) & design$score <= threshold, , drop = FALSE]
  structure(list(design = design, behavioural = behavioural, best_theta = best_theta,
                 best_sample_id = best_sample_id, threshold = threshold, ess = 1.0 / sum(weights^2),
                 obj_results = obj_results, best = best, ess_trace = ess_trace, initial_design = initial_design),
            class = "smc_result")
}

# ---- in-season recalibration + assimilation drivers -----------------------

#' Mid-season parameter re-estimation (coupled). Mirrors recalibration.py.
#' @export
recalibrate <- function(cfg, obs_df, current_date, warm_start_theta = NULL) {
  ts_limit <- as.Date(current_date)
  filtered <- obs_df[is.na(obs_df$date) | obs_df$date <= ts_limit, , drop = FALSE]
  cfg2 <- cfg
  cfg2$source <- modifyList(.cfg_get(cfg2, "source", list()), list(table = filtered))
  recal_cfg <- .cfg_get(.cfg_get(cfg, "assimilation", list()), "recalibration", list())
  recal_n <- .cfg_get(recal_cfg, "recal_sample_size", 100)
  if (!is.null(warm_start_theta) && isTRUE(.cfg_get(recal_cfg, "warm_start", TRUE))) {
    for (group in names(cfg2$parameters)) {
      params <- cfg2$parameters[[group]]; if (!is.list(params)) next
      for (name in names(params)) if (is.list(params[[name]]) && !is.null(warm_start_theta[[name]]))
        cfg2$parameters[[group]][[name]]$start <- as.numeric(warm_start_theta[[name]])
    }
  }
  if (is.null(cfg2$method)) cfg2$method <- list()
  if (is.null(cfg2$method$sample)) cfg2$method$sample <- list()
  cfg2$method$sample$n <- recal_n
  engine <- .cfg_get(recal_cfg, "engine", NULL)
  if (!is.null(engine) && tolower(as.character(engine)) != "none") {
    if (is.null(cfg2$method$bayesian)) cfg2$method$bayesian <- list()
    cfg2$method$bayesian$engine <- engine
  }
  calibrate(cfg2, progress = FALSE)$best_theta
}

#' Run in-season data assimilation. modes: recalibration | enkf | forcing.
#' Mirrors orchestrator.py:assimilate (enkf/forcing gated behind allow_uncoupled).
#' @export
assimilate <- function(cfg, progress = TRUE) {
  assim_cfg <- .cfg_get(cfg, "assimilation", list())
  mode <- .cfg_get(assim_cfg, "mode", "recalibration")
  if (!(mode %in% c("recalibration", "enkf", "forcing"))) {
    stop(sprintf("Unknown assimilation mode '%s'. Expected: recalibration | enkf | forcing.", mode))
  }
  if (mode %in% c("enkf", "forcing") && !isTRUE(.cfg_get(assim_cfg, "allow_uncoupled", FALSE))) {
    stop(sprintf("Assimilation mode '%s' is an UNCOUPLED prototype; set assimilation.allow_uncoupled: true to run it, or use mode: recalibration.", mode))
  }
  setup <- .setup(cfg); obs <- setup$obs
  if (mode == "recalibration") {
    obs_df <- obs$table
    valid_dates <- sort(unique(obs_df$date[!is.na(obs_df$date)]))
    if (length(valid_dates) == 0) return(list())
    warm_start <- isTRUE(.cfg_get(.cfg_get(assim_cfg, "recalibration", list()), "warm_start", TRUE))
    trace <- list(); prev_theta <- NULL
    for (d in valid_dates) {
      best_theta <- recalibrate(cfg, obs_df, d, warm_start_theta = if (warm_start) prev_theta else NULL)
      trace[[length(trace) + 1L]] <- list(date = d, theta = best_theta); prev_theta <- best_theta
    }
    return(list(mode = mode, trace = trace, final_theta = if (length(trace)) trace[[length(trace)]]$theta else NULL))
  } else if (mode == "forcing") {
    obs_df <- obs$table; dated <- obs_df[!is.na(obs_df$date), ]; dated <- dated[order(dated$date), ]
    state_history <- list(); current_state <- list()
    for (i in seq_len(nrow(dated))) {
      r <- dated[i, ]
      current_state <- forcing_apply(cfg, current_state, list(variable = r$variable, value = r$value,
                                                              confidence = if (is.na(r$weight)) 1.0 else r$weight))
      state_history[[length(state_history) + 1L]] <- list(date = r$date, variable = r$variable, state = current_state)
    }
    return(list(mode = mode, state_history = state_history, final_state = current_state))
  } else {
    obs_df <- obs$table
    ecfg <- .cfg_get(assim_cfg, "enkf", list())
    state_vars <- unlist(.cfg_get(ecfg, "state_variables", list("LAID", "CWAD")))
    n_ens <- as.integer(.cfg_get(ecfg, "n_ensemble", 50)); n_vars <- length(state_vars)
    set.seed(42); ensemble <- matrix(rnorm(n_ens * n_vars, 1.0, 0.2), nrow = n_ens)
    dated <- obs_df[!is.na(obs_df$date), ]; dated <- dated[order(dated$date), ]
    filter_history <- list()
    for (i in seq_len(nrow(dated))) {
      r <- dated[i, ]
      if (r$variable %in% state_vars) {
        obs_sig <- if (is.na(r$sigma)) 0.1 else r$sigma
        ensemble <- enkf_assimilate(cfg, ensemble, r$variable, r$value, obs_sig)
        filter_history[[length(filter_history) + 1L]] <- list(date = r$date, variable = r$variable,
                                                              mean_state = colMeans(ensemble), std_state = apply(ensemble, 2, sd_pop))
      }
    }
    return(list(mode = mode, filter_history = filter_history, final_ensemble_mean = colMeans(ensemble)))
  }
}

#' Calibration then in-season assimilation. Mirrors orchestrator.py:combined_mode.
#' @export
combined_mode <- function(cfg, progress = TRUE) {
  cal_result <- calibrate(cfg, progress = progress)
  cfg_assim <- cfg
  for (group in names(cfg_assim$parameters)) {
    params <- cfg_assim$parameters[[group]]; if (!is.list(params)) next
    for (name in names(params)) if (is.list(params[[name]]) && !is.null(cal_result$best_theta[[name]]))
      cfg_assim$parameters[[group]][[name]]$start <- as.numeric(cal_result$best_theta[[name]])
  }
  list(calibration = cal_result, assimilation = assimilate(cfg_assim, progress = progress))
}

# Score one theta across a set of experiments (used by combine_runs / validate_cv).
.score_theta <- function(theta, experiments, cfg, crop, specs, run_root, treatments, exe, obs, n_workers) {
  cache <- .evaluation_cache_from_setup(cfg, crop, specs, experiments, treatments[experiments],
                                        obs$table, exe)
  if (isTRUE(cache$enabled)) {
    key <- .eval_cache_key(cache, theta, experiments)
    cached <- .eval_cache_get(cache, key)
    if (!is.null(cached)) return(cached)
  }
  jobs <- lapply(experiments, function(exp) list(theta = theta, exp_id = exp, cfg = cfg, crop = crop,
                                                 param_specs = specs, run_root = run_root,
                                                 treatments = treatments[[exp]], exe = exe))
  results <- run_many(jobs, n_workers = n_workers)
  rmap <- setNames(results, experiments)
  scored <- score(rmap, obs$table, cfg)
  if (isTRUE(cache$enabled)) .eval_cache_put(cache, key, scored)
  scored
}

#' Return {exp_id: SpawnResult} for a theta (cached spawns are instant).
#' Mirrors orchestrator.py:spawn_results_for.
#' @export
spawn_results_for <- function(cfg, theta, experiments = NULL) {
  setup <- .setup(cfg)
  experiments <- experiments %||% setup$experiments
  out <- list()
  for (e in experiments) {
    out[[e]] <- spawn_and_run(theta, exp_id = e, cfg = cfg, crop = setup$crop, param_specs = setup$specs,
                              run_root = setup$run_root, treatments = setup$treatments[[e]], exe = setup$exe)
  }
  out
}

#' Combine results from multiple completed calibration directories.
#' Mirrors orchestrator.py:combine_runs.
#' @export
combine_runs <- function(cfg, run_dirs) {
  setup <- .setup(cfg); space <- setup$space
  n_workers <- resolve_cores(.cfg_get(cfg$calibrator, "num_cores", 0))
  dfs <- lapply(run_dirs, function(rdir) {
    p <- file.path(rdir, "design.csv")
    if (!file.exists(p)) stop(sprintf("No design.csv found in %s", rdir))
    utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  })
  combined <- do.call(rbind, dfs)
  param_cols <- intersect(space$names, names(combined))
  if (length(param_cols) == 0) stop("No active parameters found in design.csv files matching config")
  key <- apply(combined[param_cols], 1, function(r) paste(round(as.numeric(r), 6), collapse = "_"))
  combined <- combined[!duplicated(key), ]; rownames(combined) <- NULL
  combined$sample_id <- seq_len(nrow(combined)) - 1L
  glue <- run_glue(combined, space$names, cfg, space = space)
  best <- .score_theta(glue$best_theta, setup$experiments, cfg = cfg, crop = setup$crop, specs = setup$specs,
                       run_root = setup$run_root, treatments = setup$treatments, exe = setup$exe,
                       obs = setup$obs, n_workers = n_workers)
  .calib_result(cfg, space, setup$obs, setup$experiments, glue$design,
                setNames(list(best), as.character(glue$best_sample_id)), glue$best_theta, best, glue, list())
}

.year_key <- function(exp) { d <- gsub("[^0-9]", "", exp); if (nchar(d) >= 2) substr(d, 1, 2) else exp }
.site_key <- function(exp) { s <- gsub("[0-9]", "", exp); if (nzchar(s)) s else exp }

.make_folds <- function(experiments, scheme, seed) {
  exps <- experiments
  if (scheme %in% c("loeo", "none", "")) return(lapply(exps, function(e) list(label = e, held = e)))
  if (scheme == "year") {
    groups <- split(exps, vapply(exps, .year_key, character(1)))
    return(lapply(names(groups), function(k) list(label = paste0("year_", k), held = groups[[k]])))
  }
  if (scheme == "site") {
    groups <- split(exps, vapply(exps, .site_key, character(1)))
    return(lapply(names(groups), function(k) list(label = paste0("site_", k), held = groups[[k]])))
  }
  if (scheme == "random") {
    set.seed(seed); shuffled <- sample(exps); k <- min(length(shuffled), 5L)
    folds <- lapply(seq_len(k) - 1L, function(i) shuffled[seq(i + 1L, length(shuffled), by = k)])
    return(lapply(seq_along(folds), function(i) list(label = paste0("fold_", i - 1L), held = folds[[i]])))
  }
  stop(sprintf("Unknown validation scheme '%s'. Use loeo | year | site | random.", scheme))
}

#' Generalised cross-validation. Mirrors orchestrator.py:validate_cv.
#' @export
validate_cv <- function(cfg, scheme = NULL, progress = FALSE) {
  setup <- .setup(cfg); space <- setup$space
  n_workers <- resolve_cores(.cfg_get(cfg$calibrator, "num_cores", 0))
  method <- .cfg_get(cfg, "method", list())
  n <- as.integer(.cfg_get(.cfg_get(method, "sample", list()), "n", 100))
  seed <- as.integer(.cfg_get(cfg$calibrator, "seed", 42))
  scheme <- scheme %||% .cfg_get(.cfg_get(method, "validation", list()), "scheme", "loeo")
  rows <- list()
  for (fold in .make_folds(setup$experiments, scheme, seed)) {
    held <- fold$held; train <- setdiff(setup$experiments, held)
    if (length(train) == 0 || length(held) == 0) next
    cfg_train <- cfg; cfg_train$experiments <- as.list(train)
    samples <- sample_design(space, n = n, engine = .cfg_get(.cfg_get(method, "sample", list()), "engine", "lhs"),
                             seed = seed, include_start = TRUE)
    ed <- evaluate_design(cfg_train, samples, progress = progress)
    glue <- run_glue(ed$design, space$names, cfg_train, space = space)
    cal <- ed$obj_results[[as.character(glue$best_sample_id)]]
    cfg_held <- cfg; cfg_held$experiments <- as.list(held)
    ev <- .score_theta(glue$best_theta, held, cfg = cfg_held, crop = setup$crop, specs = setup$specs,
                       run_root = setup$run_root, treatments = setup$treatments, exe = setup$exe, obs = setup$obs, n_workers = n_workers)
    for (uv in names(cal$per_var)) rows[[length(rows) + 1L]] <- as.data.frame(c(list(fold = fold$label, held_out = paste(held, collapse = ","), split = "calibration", variable = uv), cal$per_var[[uv]]), stringsAsFactors = FALSE)
    for (uv in names(ev$per_var)) rows[[length(rows) + 1L]] <- as.data.frame(c(list(fold = fold$label, held_out = paste(held, collapse = ","), split = "evaluation", variable = uv), ev$per_var[[uv]]), stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

#' Leave-one-environment-out CV (wrapper over validate_cv).
#' @export
validate_loeo <- function(cfg, progress = FALSE) validate_cv(cfg, scheme = "loeo", progress = progress)

#' Operational in-season nowcast: (re)calibrate up to a date, persist, forecast.
#' Mirrors orchestrator.py:nowcast. This is the LAI-assimilation entry point.
#' @export
nowcast <- function(cfg, as_of_date, progress = TRUE) {
  name <- cfg$calibrator$name
  out_dir <- file.path(.cfg_get(cfg$calibrator, "results_dir", "results"), name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  state_path <- file.path(out_dir, "nowcast_state.json")
  prev_theta <- NULL
  if (file.exists(state_path) && requireNamespace("jsonlite", quietly = TRUE)) {
    prev_theta <- tryCatch(jsonlite::fromJSON(state_path)$theta, error = function(e) NULL)
  }
  obs_all <- .setup(cfg)$obs$table
  ts <- as.Date(as_of_date)
  filtered <- obs_all[is.na(obs_all$date) | obs_all$date <= ts, , drop = FALSE]
  work <- cfg; work$source <- modifyList(.cfg_get(work, "source", list()), list(table = filtered))
  if (!is.null(prev_theta)) {
    for (group in names(work$parameters)) {
      params <- work$parameters[[group]]; if (!is.list(params)) next
      for (nm in names(params)) if (is.list(params[[nm]]) && !is.null(prev_theta[[nm]])) work$parameters[[group]][[nm]]$start <- as.numeric(prev_theta[[nm]])
    }
  }
  recal_n <- .cfg_get(.cfg_get(.cfg_get(cfg, "assimilation", list()), "recalibration", list()), "recal_sample_size", NULL)
  if (!is.null(recal_n)) { if (is.null(work$method)) work$method <- list(); if (is.null(work$method$sample)) work$method$sample <- list(); work$method$sample$n <- as.integer(recal_n) }
  result <- calibrate(work, progress = progress)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    writeLines(jsonlite::toJSON(list(as_of = as.character(ts), theta = result$best_theta), auto_unbox = TRUE), state_path)
  }
  last_obs <- list()
  lai <- filtered[filtered$variable == "LAID" & !is.na(filtered$date), , drop = FALSE]
  for (exp in unique(lai$exp_id)) {
    g <- lai[lai$exp_id == exp, ]; g <- g[order(g$date), ]; r <- g[nrow(g), ]
    last_obs[[exp]] <- list(as.Date(r$date), as.numeric(r$value))
  }
  forecasts <- list()
  if (isTRUE(.cfg_get(.cfg_get(cfg, "forecast", list()), "active", FALSE))) {
    for (var in unlist(.cfg_get(.cfg_get(cfg, "forecast", list()), "variables", list("LAID")))) {
      forecasts[[var]] <- forecast_lai(work, result, last_obs = last_obs, variable = var)
    }
  }
  list(as_of = as.character(ts), best_theta = result$best_theta, result = result, forecast = forecasts, last_obs = last_obs)
}

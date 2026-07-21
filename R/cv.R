# Cross-validation framework for dssatcalibrator.
# R twin of python/dssatcalibrator/cv.py. This is the report/final-parameter
# layer over the established .make_folds(), calibrate(), and evaluate_thetas()
# primitives in orchestrator.R.

#' Parse the cross-validation configuration block
#'
#' @param cfg A complete calibrator configuration list.
#' @return A validated list of CV settings.
#' @export
parse_cv_config <- function(cfg) {
  cv <- .cfg_get(cfg, "cross_validation", list())
  strategy <- as.character(.cfg_get(cv, "strategy", "leave_one_out"))
  valid <- c("leave_one_out", "k_fold", "leave_site_out", "temporal_forward")
  if (!(strategy %in% valid)) {
    stop(sprintf(
      "Unknown cross-validation strategy '%s'. Expected one of: %s.",
      strategy, paste(sort(valid), collapse = ", ")
    ), call. = FALSE)
  }
  k <- as.integer(.cfg_get(cv, "k", 5L))
  if (strategy == "k_fold" && (is.na(k) || k < 2L)) {
    stop("cross_validation.k must be at least 2 for k_fold.", call. = FALSE)
  }
  final_theta <- as.character(.cfg_get(
    cv, "final_theta", .cfg_get(cv, "final_model", "full_refit")
  ))
  if (!(final_theta %in% c("full_refit", "best_fold", "ensemble_mean"))) {
    stop(paste(
      "cross_validation.final_theta must be full_refit, best_fold,",
      "or ensemble_mean."
    ), call. = FALSE)
  }
  scheme <- switch(
    strategy,
    leave_one_out = "loeo",
    k_fold = "random",
    leave_site_out = "site",
    temporal_forward = "year"
  )
  overfit_threshold <- as.numeric(.cfg_get(cv, "overfit_threshold", 1.5))
  if (!is.finite(overfit_threshold) || overfit_threshold < 1.2) {
    stop("cross_validation.overfit_threshold must be finite and at least 1.2.",
         call. = FALSE)
  }
  list(
    enabled = isTRUE(.cfg_get(cv, "enabled", FALSE)),
    strategy = strategy,
    scheme = scheme,
    k = k,
    seed = as.integer(.cfg_get(cv, "seed", .cfg_get(cfg$calibrator, "seed", 42L))),
    report = isTRUE(.cfg_get(cv, "report", TRUE)),
    overfit_threshold = overfit_threshold,
    final_theta = final_theta
  )
}

.cv_overfit_ratio <- function(cal_rmse, val_rmse) {
  cal_rmse <- as.numeric(cal_rmse); val_rmse <- as.numeric(val_rmse)
  if (!length(cal_rmse) || !length(val_rmse) || !is.finite(cal_rmse) ||
      !is.finite(val_rmse) || cal_rmse == 0) return(NA_real_)
  val_rmse / cal_rmse
}

.cv_metric <- function(objective, variable, metric) {
  if (is.null(objective) || is.null(objective$per_var) ||
      is.null(objective$per_var[[variable]]) ||
      is.null(objective$per_var[[variable]][[metric]])) return(NA_real_)
  as.numeric(objective$per_var[[variable]][[metric]])
}

.cv_mean <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

.aggregate_cv_results <- function(folds, cv_cfg) {
  all_vars <- sort(unique(unlist(lapply(folds, function(f) {
    c(names(.cfg_get(f$cal_obj, "per_var", list())),
      names(.cfg_get(f$val_obj, "per_var", list())))
  }))))
  rows <- lapply(folds, function(f) {
    row <- list(
      fold = as.integer(f$fold),
      strategy = cv_cfg$strategy,
      train_experiments = paste(f$train_exps, collapse = ","),
      test_experiments = paste(f$test_exps, collapse = ","),
      cal_score = as.numeric(f$cal_obj$score),
      val_score = as.numeric(f$val_obj$score)
    )
    ratios <- numeric(0)
    for (v in all_vars) {
      cal_rmse <- .cv_metric(f$cal_obj, v, "RMSE")
      val_rmse <- .cv_metric(f$val_obj, v, "RMSE")
      row[[paste0("cal_RMSE_", v)]] <- cal_rmse
      row[[paste0("val_RMSE_", v)]] <- val_rmse
      row[[paste0("cal_d_", v)]] <- .cv_metric(f$cal_obj, v, "d")
      row[[paste0("val_d_", v)]] <- .cv_metric(f$val_obj, v, "d")
      ratio <- .cv_overfit_ratio(cal_rmse, val_rmse)
      if (is.finite(ratio)) ratios <- c(ratios, ratio)
    }
    row$overfit_ratio <- if (length(ratios)) max(ratios) else NA_real_
    as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  })
  report_df <- do.call(rbind, rows)

  per_var <- setNames(lapply(all_vars, function(v) {
    cal_mean <- .cv_mean(report_df[[paste0("cal_RMSE_", v)]])
    val_mean <- .cv_mean(report_df[[paste0("val_RMSE_", v)]])
    list(cal_RMSE_mean = cal_mean, val_RMSE_mean = val_mean,
         overfit_ratio = .cv_overfit_ratio(cal_mean, val_mean))
  }), all_vars)
  max_ratio <- if (all(is.na(report_df$overfit_ratio))) NA_real_ else
    max(report_df$overfit_ratio, na.rm = TRUE)
  recommendation <- if (is.na(max_ratio) || max_ratio <= 1.2) {
    "good"
  } else if (max_ratio <= cv_cfg$overfit_threshold) {
    "mild_overfit"
  } else {
    "overfit"
  }
  finite_val <- which(is.finite(report_df$val_score))
  worst_fold <- if (length(finite_val)) {
    as.integer(report_df$fold[finite_val[which.max(report_df$val_score[finite_val])]])
  } else NULL
  summary <- list(
    cal_score_mean = .cv_mean(report_df$cal_score),
    val_score_mean = .cv_mean(report_df$val_score),
    overfit_ratio_mean = .cv_mean(report_df$overfit_ratio),
    overfit_ratio_max = max_ratio,
    per_variable = per_var,
    recommendation = recommendation,
    worst_fold = worst_fold
  )
  list(summary = summary, report_df = report_df)
}

.select_cv_final_theta <- function(cfg, folds, progress) {
  method <- parse_cv_config(cfg)$final_theta
  if (method == "full_refit") {
    if (isTRUE(progress)) message("Running full refit for final theta...")
    return(list(theta = calibrate(cfg, progress = progress)$best_theta, method = method))
  }
  if (!length(folds)) return(list(theta = list(), method = method))
  if (method == "best_fold") {
    scores <- vapply(folds, function(f) as.numeric(f$val_obj$score), numeric(1))
    scores[!is.finite(scores)] <- Inf
    if (!any(is.finite(scores))) {
      stop("Cross-validation produced no finite validation scores.", call. = FALSE)
    }
    return(list(theta = folds[[which.min(scores)]]$best_theta, method = method))
  }
  params <- names(folds[[1]]$best_theta)
  if (length(folds) > 1L && any(vapply(folds[-1], function(f) {
    !setequal(names(f$best_theta), params)
  }, logical(1)))) {
    stop("Cross-validation folds have inconsistent parameter sets.", call. = FALSE)
  }
  theta <- setNames(lapply(params, function(p) {
    mean(vapply(folds, function(f) as.numeric(f$best_theta[[p]]), numeric(1)))
  }), params)
  list(theta = theta, method = method)
}

#' Run cross-validation and produce overfit diagnostics
#'
#' @param cfg Configuration list with `cross_validation.enabled: true`.
#' @param progress Whether to show fold progress.
#' @return A `cv_result` list matching Python's `CVResult` fields.
#' @export
run_cross_validation <- function(cfg, progress = TRUE) {
  cv_cfg <- parse_cv_config(cfg)
  if (!cv_cfg$enabled) stop("Cross-validation is not enabled in config.", call. = FALSE)
  setup <- .setup(cfg)
  experiments <- as.character(setup$experiments)
  if (length(experiments) < 2L) {
    stop("Cross-validation requires at least two active experiments.", call. = FALSE)
  }
  folds_def <- .make_folds(
    experiments, cv_cfg$scheme, cv_cfg$seed, k = cv_cfg$k
  )
  if (cv_cfg$strategy == "temporal_forward") {
    keys <- vapply(folds_def, function(f) {
      if (length(f$held)) .year_key(as.character(f$held[[1]])) else f$label
    }, character(1))
    folds_def <- folds_def[order(keys)]
    accumulated <- character(0)
    fold_jobs <- list()
    for (fold_def in folds_def) {
      held <- as.character(fold_def$held)
      if (length(accumulated)) {
        fold_jobs[[length(fold_jobs) + 1L]] <- list(
          label = fold_def$label, held = held, train = accumulated
        )
      }
      accumulated <- c(accumulated, held)
    }
  } else {
    fold_jobs <- lapply(folds_def, function(fold_def) {
      held <- as.character(fold_def$held)
      list(label = fold_def$label, held = held,
           train = experiments[!(experiments %in% held)])
    })
  }
  n_workers <- resolve_cores(.cfg_get(cfg$calibrator, "num_cores", 0))

  folds <- list()
  for (fold_def in fold_jobs) {
    held <- fold_def$held
    train <- fold_def$train
    if (!length(train) || !length(held)) next
    i <- length(folds) + 1L
    if (isTRUE(progress)) {
      message(sprintf("Running Fold %d/%d: %s", i, length(fold_jobs), fold_def$label))
    }
    train_cfg <- cfg
    train_cfg$experiments <- train
    train_res <- calibrate(train_cfg, progress = progress)

    test_cfg <- cfg
    test_cfg$experiments <- held
    test_eval <- evaluate_thetas(
      test_cfg, list(train_res$best_theta), setup = .setup(test_cfg),
      n_workers = n_workers, progress = progress
    )$results
    if (!length(test_eval)) {
      stop(sprintf("Cross-validation fold '%s' returned no validation result.",
                   fold_def$label), call. = FALSE)
    }
    folds[[i]] <- list(
      fold = i,
      label = fold_def$label,
      train_exps = train,
      test_exps = held,
      best_theta = train_res$best_theta,
      cal_obj = train_res$best,
      val_obj = test_eval[[1]]
    )
  }
  if (!length(folds)) {
    stop("Cross-validation produced no non-empty train/test folds.", call. = FALSE)
  }
  final <- .select_cv_final_theta(cfg, folds, progress)
  aggregated <- .aggregate_cv_results(folds, cv_cfg)
  structure(list(
    strategy = cv_cfg$strategy,
    folds = folds,
    summary = aggregated$summary,
    report_df = aggregated$report_df,
    final_theta = final$theta,
    final_theta_method = final$method
  ), class = "cv_result")
}

#' Write a cross-validation report to disk
#'
#' @param result The `cv_result` returned by `run_cross_validation()`.
#' @param outdir Directory for `cv_report.csv` and `cv_summary.json`.
#' @return Invisibly, named paths to the two report files.
#' @export
write_cv_report <- function(result, outdir) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to write a CV report.", call. = FALSE)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  csv_file <- file.path(outdir, "cv_report.csv")
  json_file <- file.path(outdir, "cv_summary.json")
  utils::write.csv(result$report_df, csv_file, row.names = FALSE)
  writeLines(jsonlite::toJSON(
    result$summary, auto_unbox = TRUE, pretty = TRUE, na = "null"
  ), json_file, useBytes = TRUE)
  invisible(list(csv = csv_file, json = json_file))
}

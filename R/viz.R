# Reporting & figures. R twin of python/dssatcalibrator/viz.py.
# summary_fit_table is pure; make_report writes the same CSV tables and (when
# ggplot2 is available) the diagnostic figures the CLI emits.

#' Per-variable best-fit metrics table. Mirrors viz.py:summary_fit_table.
#' @export
summary_fit_table <- function(result) {
  per_var <- result$best$per_var %||% list()
  if (length(per_var) == 0) return(data.frame())
  rows <- lapply(names(per_var), function(v) {
    m <- per_var[[v]]
    data.frame(variable = v, n = m$n %||% NA, RMSE = m$RMSE %||% NA, nRMSE_pct = m$nRMSE_pct %||% NA,
               MBE = m$MBE %||% NA, d = m$d %||% NA, EF = m$EF %||% NA, R2 = m$R2 %||% NA,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Write the calibration report (CSV tables + optional figures), returning the
#' written paths keyed by name. Mirrors viz.py:make_report.
#' @export
make_report <- function(result, outdir, best_spawns = NULL, figdir = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  paths <- list()

  design_path <- file.path(outdir, "design.csv")
  utils::write.csv(result$design, design_path, row.names = FALSE)
  paths$design <- design_path

  fit_path <- file.path(outdir, "fit_summary.csv")
  utils::write.csv(summary_fit_table(result), fit_path, row.names = FALSE)
  paths$fit_summary <- fit_path

  best_path <- file.path(outdir, "best_theta.csv")
  utils::write.csv(data.frame(parameter = names(result$best_theta),
                              value = unlist(result$best_theta), stringsAsFactors = FALSE),
                   best_path, row.names = FALSE)
  paths$best_theta <- best_path

  if (!is.null(result$glue)) {
    post <- tryCatch(posterior_summary(result$glue, result$space$names), error = function(e) NULL)
    if (!is.null(post)) {
      post_path <- file.path(outdir, "posterior_summary.csv")
      utils::write.csv(post, post_path, row.names = FALSE)
      paths$posterior_summary <- post_path
    }
  }

  if (!is.null(figdir) && requireNamespace("ggplot2", quietly = TRUE)) {
    dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
    resid <- result$best$residuals
    if (!is.null(resid) && nrow(resid) > 0) {
      p <- ggplot2::ggplot(resid, ggplot2::aes(x = obs, y = sim, colour = user_var)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::labs(title = "Observed vs simulated (best fit)", x = "Observed", y = "Simulated") +
        ggplot2::theme_minimal()
      ov <- file.path(figdir, "obs_vs_sim.png")
      suppressMessages(ggplot2::ggsave(ov, p, width = 6, height = 5, dpi = 120))
      paths$obs_vs_sim <- ov
    }
  }
  paths
}

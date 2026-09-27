# Optional adapters reuse the shared execution, calibration and scoring engine.
model_platform <- function(cfg) {
  value <- gsub("-", "_", tolower(as.character(.cfg_get(
    .cfg_get(cfg, "model", list()), "platform", .cfg_get(cfg, "platform", "dssat")))))
  aliases <- c(dssat_csm = "dssat", apsim = "apsim_ng", apsimx = "apsim_ng",
               apsim_next_generation = "apsim_ng", apsim_nextgen = "apsim_ng")
  if (value %in% names(aliases)) unname(aliases[[value]]) else value
}

.model_adapter <- function(cfg) {
  platform <- model_platform(cfg)
  if (identical(platform, "dssat")) return(NULL)
  package <- .cfg_get(.cfg_get(cfg, "model", list()), "adapter_package", "cropmodelcalibrator")
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Model platform '", platform, "' needs the installed adapter package '", package, "'.")
  }
  getExportedValue(package, "get_model_adapter")(platform)
}

# Stable helper for model adapters; keep parameter scope handling in one place.
effective_parameters <- function(theta, param_specs, exp_id, cultivars = NULL) {
  .effective_theta(theta, param_specs, exp_id, cultivars)
}

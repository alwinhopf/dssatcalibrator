# Optional weather/soil acquisition for synthesized / new-site experiments.
# R twin of python/dssatcalibrator/{acquisition,weather}.py.
#
# Acquisition is delegated to the shared `dssatutils` package (the same provider
# the rest of the workspace uses), so the calibrator carries no download code of
# its own. These wrappers map config -> a dssatutils process_* call.

.dssatutils_required <- function() {
  if (!requireNamespace("dssatutils", quietly = TRUE)) {
    stop("weather/soil acquisition requires dssatutils. Install dssatcalibrator's ",
         "'acquire' extra or the dssatutils R package.")
  }
}

#' Acquire a single-site soil profile via dssatutils and write a .SOL.
#' Mirrors acquisition.py:acquire_soil_profile (provider from soil.source).
#' @export
acquire_soil_profile <- function(cfg, site_id, lat, lon, out_path) {
  .dssatutils_required()
  source <- tolower(as.character(.cfg_get(.cfg_get(cfg, "soil", list()), "source", "ssurgo")))
  fn_name <- paste0("process_soils_", source)
  if (!exists(fn_name, where = asNamespace("dssatutils"))) {
    stop(sprintf("dssatutils has no soil provider '%s'", source))
  }
  fn <- get(fn_name, envir = asNamespace("dssatutils"))
  fn(site_id = site_id, lat = lat, lon = lon, out_path = out_path)
  out_path
}

#' Acquire daily weather via dssatutils and write a .WTH for [start, end].
#' Mirrors weather.py:acquire_wth (provider from weather.provider).
#' @export
acquire_wth <- function(cfg, station, lat, lon, start, end, out_path) {
  .dssatutils_required()
  provider <- tolower(as.character(.cfg_get(.cfg_get(cfg, "weather", list()), "provider", "nasapower")))
  fn_name <- paste0("process_weather_", provider)
  if (!exists(fn_name, where = asNamespace("dssatutils"))) {
    stop(sprintf("dssatutils has no weather provider '%s'", provider))
  }
  fn <- get(fn_name, envir = asNamespace("dssatutils"))
  fn(station = station, lat = lat, lon = lon, start = start, end = end, out_path = out_path)
  out_path
}

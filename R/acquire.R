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
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("weather/soil acquisition requires the 'sf' package.")
  }
}

.make_single_point_sf <- function(id, lat, lon) {
  lat_f <- as.numeric(lat); lon_f <- as.numeric(lon)
  if (is.na(lat_f) || is.na(lon_f) || abs(lat_f) > 90 || abs(lon_f) > 180) {
    stop(sprintf("Invalid site coordinates: lat=%s, lon=%s", lat, lon))
  }
  df <- data.frame(ID = as.character(id), LAT = lat_f, LONG = lon_f, stringsAsFactors = FALSE)
  sf::st_as_sf(df, coords = c("LONG", "LAT"), crs = 4326, remove = FALSE)
}

#' Acquire a single-site soil profile via dssatutils and write a .SOL.
#' Mirrors acquisition.py:acquire_soil_profile (provider from soil.source).
#' @export
acquire_soil_profile <- function(cfg, site_id, lat, lon, out_path) {
  scfg <- .cfg_get(cfg, "soil", list())
  provider <- tolower(as.character(.cfg_get(scfg, "provider", "file")))
  if (provider %in% c("", "file", "none")) return(NULL)
  if (provider != "dssatutils") {
    stop("soil.provider must be 'file' or 'dssatutils'.")
  }
  .dssatutils_required()

  source <- tolower(as.character(.cfg_get(scfg, "source", "ssurgo")))
  fn_name <- paste0("process_soils_", source)
  if (!exists(fn_name, where = asNamespace("dssatutils"))) {
    stop(sprintf("dssatutils has no soil provider '%s'", source))
  }

  pts <- .make_single_point_sf(site_id, lat, lon)
  cache_dir <- .cfg_get(scfg, "cache_dir", tempfile("soil_cache"))
  sol_dir <- file.path(cache_dir, paste0(source, "_individual_SOL"))
  dir.create(sol_dir, recursive = TRUE, showWarnings = FALSE)
  map_csv <- file.path(cache_dir, paste0(source, "_soil_map.csv"))
  n_cores <- as.integer(.cfg_get(scfg, "n_cores", 1L))

  fn <- get(fn_name, envir = asNamespace("dssatutils"))
  if (source == "ssurgo") {
    fn(grid_points = pts, output_dir_csv = map_csv, output_dir_individual = sol_dir,
       n_cores = n_cores, id_col = "ID", lat_col = "LAT", long_col = "LONG")
  } else {
    fn(pts, map_csv, sol_dir, n_cores, "ID", "LAT", "LONG")
  }

  direct <- file.path(sol_dir, paste0(site_id, ".SOL"))
  src <- if (file.exists(direct)) direct else NULL
  if (is.null(src) && file.exists(map_csv)) {
    mapping <- utils::read.csv(map_csv, stringsAsFactors = FALSE)
    if ("ID" %in% names(mapping)) {
      sub <- mapping[as.character(mapping$ID) == as.character(site_id), , drop = FALSE]
      if (nrow(sub) > 0) {
        for (col in c("SOIL_ID", "soil_id", "SOURCE_SOIL_ID")) {
          if (col %in% names(sub) && !is.na(sub[[col]][1]) && nzchar(as.character(sub[[col]][1]))) {
            cand <- file.path(sol_dir, paste0(sub[[col]][1], ".SOL"))
            if (file.exists(cand)) { src <- cand; break }
          }
        }
      }
    }
  }

  if (is.null(src)) {
    stop(sprintf("dssatutils did not produce a .SOL for site '%s' in %s", site_id, sol_dir))
  }

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  file.copy(src, out_path, overwrite = TRUE)
  out_path
}

#' Acquire daily weather via dssatutils and write a .WTH for [start, end].
#' Mirrors weather.py:acquire_wth (provider from weather.provider).
#' @export
acquire_wth <- function(cfg, station, lat, lon, start, end, out_path) {
  wcfg <- .cfg_get(cfg, "weather", list())
  provider <- tolower(as.character(.cfg_get(wcfg, "provider", "file")))
  if (provider %in% c("", "file", "none")) return(NULL)

  if (provider == "dssatutils") {
    provider <- tolower(as.character(.cfg_get(wcfg, "dssatutils_provider", "nasapower")))
  }
  if (provider == "nasa_power") provider <- "nasapower"

  fn_name <- paste0("process_weather_", provider)
  if (!exists(fn_name, where = asNamespace("dssatutils"))) {
    stop(sprintf("dssatutils has no weather provider '%s'", provider))
  }
  .dssatutils_required()

  pts <- .make_single_point_sf(station, lat, lon)
  cache_dir <- .cfg_get(wcfg, "cache_dir", tempfile("weather_cache"))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(cache_dir, "dssatutils_weather.log")
  start_year <- as.integer(substr(as.character(start), 1, 4))
  end_year <- as.integer(substr(as.character(end), 1, 4))
  n_cores <- as.integer(.cfg_get(wcfg, "n_cores", 1L))

  fn <- get(fn_name, envir = asNamespace("dssatutils"))
  fn(shapefile = pts, start_year = start_year, end_year = end_year,
     output_dir = cache_dir, id_col = "ID", lat_col = "LAT", lon_col = "LONG",
     n_cores = n_cores, log_file = log_file)

  cand <- file.path(cache_dir, paste0(station, ".WTH"))
  if (!file.exists(cand)) {
    # Try finding any .WTH in cache_dir
    found <- list.files(cache_dir, pattern = "\\.WTH$", full.names = TRUE, ignore.case = TRUE)
    if (length(found)) cand <- found[1]
  }

  if (!file.exists(cand)) {
    stop(sprintf("dssatutils did not produce %s.WTH in %s", station, cache_dir))
  }

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  file.copy(cand, out_path, overwrite = TRUE)
  out_path
}

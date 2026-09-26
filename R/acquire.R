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

.get_utils_provider <- function(name) {
  ns <- asNamespace("dssatutils")
  if (!exists(name, envir = ns, inherits = FALSE)) stop("Unknown dssatutils provider: ", name)
  get(name, envir = ns, inherits = FALSE)
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

  pts <- .make_single_point_sf(site_id, lat, lon)
  cache_dir <- .cfg_get(scfg, "cache_dir", "soil_cache")
  sol_dir <- file.path(cache_dir, paste0(source, "_individual_SOL"))
  dir.create(sol_dir, recursive = TRUE, showWarnings = FALSE)
  map_csv <- file.path(cache_dir, paste0(source, "_soil_map.csv"))
  n_cores <- as.integer(.cfg_get(scfg, "n_cores", 1L))

  fn <- .get_utils_provider(fn_name)
  if (source == "ssurgo") {
    fn(grid_points = pts, output_dir_csv = map_csv, output_dir_individual = sol_dir,
       n_cores = n_cores, id_col = "ID", lat_col = "LAT", long_col = "LONG")
  } else if (source == "soilgrids") {
    source_sol <- scfg$source_sol_file %||% scfg$external_soil_file
    if (is.null(source_sol) || !nzchar(source_sol)) stop("soil.source: soilgrids requires soil.source_sol_file or soil.external_soil_file")
    fn(grid_points = pts, source_sol_file = source_sol, output_csv_path = map_csv,
       output_sol_dir = sol_dir, id_col = "ID")
  } else if (source == "soilgrids_online") {
    mode <- toupper(.cfg_get(scfg, "soilgrids_mode", "REST"))
    if (!mode %in% c("REST", "VRT")) stop("soilgrids_mode must be REST or VRT")
    fn(gridfile = pts, soilfile_csv_path = map_csv, output_sol_dir = sol_dir,
       id_col = "ID", use_rest_api = mode == "REST")
  } else {
    stop("soil.source must be one of: ssurgo, soilgrids, soilgrids_online")
  }

  direct <- file.path(sol_dir, paste0(site_id, ".SOL"))
  src <- if (file.exists(direct)) direct else NULL
  if (is.null(src) && file.exists(map_csv)) {
    mapping <- utils::read.csv(map_csv, stringsAsFactors = FALSE, colClasses = "character")
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
  if (!file.copy(src, out_path, overwrite = TRUE)) stop("Could not copy acquired soil to ", out_path)
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
  if (provider %in% c("nasa_power", "nasa-power")) provider <- "nasapower"

  .dssatutils_required()
  fn_name <- paste0("process_weather_", provider)
  .dssatutils_required()

  pts <- .make_single_point_sf(station, lat, lon)
  cache_dir <- .cfg_get(wcfg, "cache_dir", "weather_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(cache_dir, "dssatutils_weather.log")
  start_year <- as.integer(substr(as.character(start), 1, 4))
  end_year <- as.integer(substr(as.character(end), 1, 4))
  n_cores <- as.integer(.cfg_get(wcfg, "n_cores", 1L))

  fn <- .get_utils_provider(fn_name)
  fn(shapefile = pts, start_year = start_year, end_year = end_year,
     output_dir = cache_dir, id_col = "ID", lat_col = "LAT", lon_col = "LONG",
     n_cores = n_cores, log_file = log_file)

  cand <- file.path(cache_dir, paste0(station, ".WTH"))

  if (!file.exists(cand) || file.info(cand)$size == 0) {
    stop(sprintf("dssatutils did not produce %s.WTH in %s", station, cache_dir))
  }

  if (!dssatutils::is_wth_valid(cand, required_columns = c("SRAD", "TMAX", "TMIN", "RAIN"),
                               start_date = as.character(start), end_date = as.character(end))) {
    stop("Acquired weather is invalid or incomplete for station ", station)
  }
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(cand, out_path, overwrite = TRUE)) stop("Could not copy acquired weather to ", out_path)
  out_path
}

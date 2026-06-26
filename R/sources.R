# Pluggable observation sources. R twin of python/dssatcalibrator/sources/*.py.
#
# Modelled as S3 objects (one class per adapter) dispatching the generics
# src_fetch / src_error_model / src_variable_mapping / src_quality_filter.
# Function names and error-model formulas mirror the Python adapters exactly.

# small helper: nested config lookup with default
.cfg_get <- function(block, key, default = NULL) {
  if (is.null(block) || is.null(block[[key]])) default else block[[key]]
}

# ---- generics --------------------------------------------------------------

#' @export
src_fetch <- function(src, experiment, date_range, ...) UseMethod("src_fetch")
#' @export
src_error_model <- function(src, variable, value, metadata) UseMethod("src_error_model")
#' @export
src_variable_mapping <- function(src) UseMethod("src_variable_mapping")
#' @export
src_quality_filter <- function(src, df) UseMethod("src_quality_filter")
#' @export
src_quality_filter.default <- function(src, df) df

.new_source <- function(adapter, name, source_type, config) {
  structure(list(config = config, name = name, source_type = source_type),
            class = c(adapter, "obs_source"))
}

# Generic CSV fetch shared by the satellite/uav/iot/farm adapters: read CSV,
# parse dates, filter to (experiment, date_range). Returns the filtered frame.
.read_source_csv <- function(src, experiment, date_range) {
  path <- .cfg_get(src$config, "data_path")
  if (is.null(path) || !nzchar(path)) return(NULL)
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  df$date <- as.Date(as.character(df$date))
  df[df$exp_id == experiment & df$date >= date_range[[1]] & df$date <= date_range[[2]], , drop = FALSE]
}

# ---- field measurements ----------------------------------------------------

#' @export
source_field_measurements <- function(config)
  .new_source("field_measurements", "field_measurements", "field", config)

#' @export
src_fetch.field_measurements <- function(src, experiment, date_range, ...) {
  hemp_dir <- .cfg_get(src$config, "hemp_dir", .cfg_get(src$config, "data_path"))
  if (is.null(hemp_dir)) return(.empty_schema())
  crop_ext <- .cfg_get(src$config, "crop_ext", "HM")
  fa <- file.path(hemp_dir, sprintf("%s.%sA", experiment, crop_ext))
  ft <- file.path(hemp_dir, sprintf("%s.%sT", experiment, crop_ext))
  frames <- list()
  if (file.exists(fa)) frames[[length(frames) + 1L]] <- read_filea(fa, experiment)
  if (file.exists(ft)) frames[[length(frames) + 1L]] <- read_filet(ft, experiment)
  if (length(frames) == 0) return(.empty_schema())
  df <- do.call(rbind, frames)
  in_range <- is.na(df$date) | (df$date >= date_range[[1]] & df$date <= date_range[[2]])
  df <- df[in_range, , drop = FALSE]
  df$source <- src$name; df$quality_flag <- 0L; df$spatial_res_m <- NA_real_
  df$sigma <- mapply(function(v, x) src_error_model(src, v, x, list()),
                     df$variable, df$value)
  df
}

#' @export
src_error_model.field_measurements <- function(src, variable, value, metadata) {
  defaults <- list(LAID = c("relative", 0.15), CWAD = c("relative", 0.12),
                   HWAM = c("relative", 0.08), GSTD = c("absolute", 1.0),
                   ADAT = c("absolute", 3.0), MDAT = c("absolute", 3.0))
  cfg_models <- .cfg_get(src$config, "error_model", list())
  if (!is.null(cfg_models[[variable]])) {
    kind <- .cfg_get(cfg_models[[variable]], "type", "relative")
    val <- .cfg_get(cfg_models[[variable]], "value", 0.15)
  } else if (!is.null(defaults[[variable]])) {
    kind <- defaults[[variable]][1]; val <- as.numeric(defaults[[variable]][2])
  } else {
    kind <- "relative"; val <- 0.15
  }
  if (kind == "relative") return(max(abs(val * value), 1e-6))
  as.numeric(val)
}

#' @export
src_variable_mapping.field_measurements <- function(src) list()

# ---- UAV multispectral -----------------------------------------------------

#' @export
source_uav_multispectral <- function(config)
  .new_source("uav_multispectral", "uav_multispectral", "uav", config)

#' @export
src_fetch.uav_multispectral <- function(src, experiment, date_range, ...) {
  df <- .read_source_csv(src, experiment, date_range)
  if (is.null(df) || nrow(df) == 0) return(.empty_extended())
  rows <- lapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    meta <- list(flight_quality = if (!is.null(r$flight_quality)) r$flight_quality else "good")
    .ext_row(r$exp_id, r$treatment, r$variable, "timeseries", r$date, as.numeric(r$value),
             src_error_model(src, r$variable, as.numeric(r$value), meta),
             1.0, src$name, 0L, 0.05)
  })
  do.call(rbind, rows)
}

#' @export
src_error_model.uav_multispectral <- function(src, variable, value, metadata) {
  fq <- .cfg_get(metadata, "flight_quality", "good")
  cfg_models <- .cfg_get(src$config, "error_model", list())
  if (!is.null(cfg_models[[variable]])) {
    val <- .cfg_get(cfg_models[[variable]], "value", 0.15)
  } else {
    base <- list(LAID = 0.4, canopy_cover = 0.05, canopy_height = 0.03)
    val <- if (!is.null(base[[variable]])) base[[variable]] else 0.15 * abs(value)
  }
  if (identical(fq, "poor")) val <- val * 1.5
  as.numeric(val)
}

#' @export
src_variable_mapping.uav_multispectral <- function(src)
  list(uav_lai = "LAID", uav_canopy_cover = "canopy_cover", uav_canopy_height = "canopy_height")

# ---- IoT: soil moisture ----------------------------------------------------

#' @export
source_soil_moisture_iot <- function(config)
  .new_source("soil_moisture_iot", "soil_moisture_iot", "iot", config)

#' @export
src_fetch.soil_moisture_iot <- function(src, experiment, date_range, ...) {
  df <- .read_source_csv(src, experiment, date_range)
  if (is.null(df) || nrow(df) == 0) return(.empty_extended())
  meta <- list(sensor_type = .cfg_get(src$config, "sensor_type", "capacitance"),
               calibration_status = .cfg_get(src$config, "calibration_status", "factory"))
  rows <- lapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    .ext_row(r$exp_id, r$treatment, "SW", "timeseries", r$date, as.numeric(r$value),
             src_error_model(src, "SW", as.numeric(r$value), meta), 1.0, src$name, 0L, NA_real_)
  })
  do.call(rbind, rows)
}

#' @export
src_error_model.soil_moisture_iot <- function(src, variable, value, metadata) {
  sensor <- .cfg_get(metadata, "sensor_type", "capacitance")
  cal <- .cfg_get(metadata, "calibration_status", "factory")
  base <- c(capacitance = 0.04, tdr = 0.02, tensiometer = 0.03)
  b <- if (!is.na(base[sensor])) base[[sensor]] else 0.04
  if (identical(cal, "field_calibrated")) b <- b * 0.6
  as.numeric(b)
}

#' @export
src_variable_mapping.soil_moisture_iot <- function(src) list(soil_moisture = "SW")

# ---- IoT: canopy temperature -----------------------------------------------

#' @export
source_canopy_temperature <- function(config)
  .new_source("canopy_temperature", "canopy_temperature", "iot", config)

#' @export
src_fetch.canopy_temperature <- function(src, experiment, date_range, ...) {
  df <- .read_source_csv(src, experiment, date_range)
  if (is.null(df) || nrow(df) == 0) return(.empty_extended())
  rows <- lapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    .ext_row(r$exp_id, r$treatment, "TMEAN", "timeseries", r$date, as.numeric(r$value),
             src_error_model(src, "TMEAN", as.numeric(r$value), list()), 1.0, src$name, 0L, NA_real_)
  })
  do.call(rbind, rows)
}

#' @export
src_error_model.canopy_temperature <- function(src, variable, value, metadata) {
  as.numeric(.cfg_get(.cfg_get(src$config, "error_model", list()), "value", 1.0))
}

#' @export
src_variable_mapping.canopy_temperature <- function(src) list(canopy_temp = "TMEAN")

# ---- satellite: Sentinel-2 LAI ---------------------------------------------

# Optional LAI observation operator: value' = scale*value + offset (identity by default).
.apply_obs_operator <- function(config, value) {
  op <- .cfg_get(config, "obs_operator", list())
  as.numeric(.cfg_get(op, "scale", 1.0)) * value + as.numeric(.cfg_get(op, "offset", 0.0))
}

#' @export
source_sentinel2_lai <- function(config)
  .new_source("sentinel2_lai", "sentinel2_lai", "satellite", config)

#' @export
src_fetch.sentinel2_lai <- function(src, experiment, date_range, ...) {
  df <- .read_source_csv(src, experiment, date_range)
  if (is.null(df) || nrow(df) == 0) return(.empty_extended())
  max_cloud <- as.numeric(.cfg_get(src$config, "max_cloud_fraction", 1.0))
  rows <- list()
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    cloud <- as.numeric(if (!is.null(r$cloud_fraction)) r$cloud_fraction else 0.0)
    if (cloud > max_cloud) next
    val <- .apply_obs_operator(src$config, as.numeric(r$value))
    qf <- as.integer(if (!is.null(r$quality_flag)) r$quality_flag else 0L)
    rows[[length(rows) + 1L]] <- .ext_row(r$exp_id, r$treatment, "LAID", "timeseries",
      r$date, val, src_error_model(src, "LAID", val, list(cloud_fraction = cloud)),
      1.0, src$name, qf, 10.0)
  }
  if (length(rows) == 0) return(.empty_extended())
  do.call(rbind, rows)
}

#' @export
src_quality_filter.sentinel2_lai <- function(src, df) {
  if (nrow(df) == 0 || !("quality_flag" %in% names(df))) return(df)
  if (!isTRUE(.cfg_get(src$config, "drop_bad_quality", TRUE))) return(df)
  qf <- df$quality_flag; qf[is.na(qf)] <- 0
  df[qf == 0, , drop = FALSE]
}

#' @export
src_error_model.sentinel2_lai <- function(src, variable, value, metadata) {
  em <- .cfg_get(src$config, "error_model", list())
  base_rmse <- as.numeric(.cfg_get(em, "base_rmse", 0.7))
  sat_lai <- as.numeric(.cfg_get(em, "saturation_lai", 4.0))
  if (value > sat_lai) base_rmse <- base_rmse * (1.0 + 0.3 * (value - sat_lai))
  cloud_factor <- 1.0 + as.numeric(.cfg_get(metadata, "cloud_fraction", 0.0)) * 0.5
  base_rmse * cloud_factor
}

#' @export
src_variable_mapping.sentinel2_lai <- function(src) list(sentinel_lai = "LAID")

# ---- satellite: MODIS LAI --------------------------------------------------

#' @export
source_modis_lai <- function(config)
  .new_source("modis_lai", "modis_lai", "satellite", config)

#' @export
src_fetch.modis_lai <- function(src, experiment, date_range, ...) {
  df <- .read_source_csv(src, experiment, date_range)
  if (is.null(df) || nrow(df) == 0) return(.empty_extended())
  max_qc <- as.integer(.cfg_get(src$config, "max_qc_flag", 99L))
  rows <- list()
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    qc <- as.integer(if (!is.null(r$qc_flag)) r$qc_flag else 0L)
    if (qc > max_qc) next
    val <- .apply_obs_operator(src$config, as.numeric(r$value))
    rows[[length(rows) + 1L]] <- .ext_row(r$exp_id, r$treatment, "LAID", "timeseries",
      r$date, val, src_error_model(src, "LAID", val, list(qc_flag = qc)),
      0.5, src$name, qc, 250.0)
  }
  if (length(rows) == 0) return(.empty_extended())
  do.call(rbind, rows)
}

#' @export
src_error_model.modis_lai <- function(src, variable, value, metadata) {
  base_rmse <- as.numeric(.cfg_get(.cfg_get(src$config, "error_model", list()), "base_rmse", 0.66))
  qc <- as.numeric(.cfg_get(metadata, "qc_flag", 0))
  base_rmse * (1.0 + 0.5 * qc)
}

#' @export
src_variable_mapping.modis_lai <- function(src) list(modis_lai = "LAID")

# ---- farm software: phenology ----------------------------------------------

#' @export
source_farm_phenology <- function(config)
  .new_source("farm_phenology", "farm_phenology", "farm_software", config)

#' @export
src_fetch.farm_phenology <- function(src, experiment, date_range, ...) {
  df <- .read_source_csv(src, experiment, date_range)
  if (is.null(df) || nrow(df) == 0) return(.empty_extended())
  prec <- .cfg_get(src$config, "date_precision", "exact")
  rows <- lapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]; var <- r$variable
    is_gstd <- var == "GSTD"
    val_for_sigma <- if (is_gstd) as.numeric(r$value) else 1.0
    .ext_row(r$exp_id, r$treatment, var,
             if (var %in% c("ADAT", "MDAT")) "phenology" else "timeseries",
             r$date, as.numeric(r$value),
             src_error_model(src, var, val_for_sigma, list(date_precision = prec)),
             1.0, src$name, 0L, NA_real_)
  })
  do.call(rbind, rows)
}

#' @export
src_error_model.farm_phenology <- function(src, variable, value, metadata) {
  precision <- .cfg_get(metadata, "date_precision", "exact")
  if (variable %in% c("GSTD", "growth_stage")) {
    m <- c(exact = 1.0, weekly = 2.0, biweekly = 3.0)
  } else {
    m <- c(exact = 2.0, weekly = 5.0, biweekly = 7.0)
  }
  as.numeric(if (!is.na(m[precision])) m[[precision]] else m[["exact"]])
}

#' @export
src_variable_mapping.farm_phenology <- function(src)
  list(growth_stage = "GSTD", anthesis_date = "ADAT", maturity_date = "MDAT")

# ---- farm software: management events --------------------------------------

#' @export
source_farm_management <- function(config)
  .new_source("farm_management", "farm_management", "farm_software", config)

#' @export
src_fetch.farm_management <- function(src, experiment, date_range, ...) {
  df <- .read_source_csv(src, experiment, date_range)
  if (is.null(df) || nrow(df) == 0) return(.empty_extended())
  rows <- lapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    .ext_row(r$exp_id, r$treatment, r$variable, "management_constraint", r$date,
             as.numeric(r$value), src_error_model(src, r$variable, as.numeric(r$value), list()),
             1.0, src$name, 0L, NA_real_)
  })
  do.call(rbind, rows)
}

#' @export
src_error_model.farm_management <- function(src, variable, value, metadata) {
  if (grepl("date", tolower(variable))) return(1.0)
  abs(value * 0.05)
}

#' @export
src_variable_mapping.farm_management <- function(src) list()

# ---- registry --------------------------------------------------------------

#' Registry mapping adapter names to source constructors. Mirrors
#' sources/registry.py:ADAPTER_REGISTRY.
#' @export
ADAPTER_REGISTRY <- list(
  sentinel2_lai      = source_sentinel2_lai,
  modis_lai          = source_modis_lai,
  farm_phenology     = source_farm_phenology,
  farm_management    = source_farm_management,
  field_measurements = source_field_measurements,
  soil_moisture_iot  = source_soil_moisture_iot,
  canopy_temperature = source_canopy_temperature,
  uav_multispectral  = source_uav_multispectral
)

#' Instantiate all active observation sources from config.
#' Mirrors sources/registry.py:build_sources.
#' @export
build_sources <- function(cfg) {
  sources <- list()
  blocks <- .cfg_get(cfg, "observation_sources", list())
  for (name in names(blocks)) {
    block <- blocks[[name]]
    if (!isTRUE(.cfg_get(block, "active", FALSE))) next
    adapter_name <- .cfg_get(block, "adapter", name)
    ctor <- ADAPTER_REGISTRY[[adapter_name]]
    if (is.null(ctor)) {
      stop(sprintf("Unknown adapter: %s. Available: %s",
                   adapter_name, paste(names(ADAPTER_REGISTRY), collapse = ", ")))
    }
    sources[[length(sources) + 1L]] <- ctor(block)
  }
  sources
}

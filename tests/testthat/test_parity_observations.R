# Cross-language parity tests for Phase 2: observations, sources, fusion.
# Loads golden fixtures generated from the Python implementation.

fixture_dir <- file.path("..", "fixtures")
if (!dir.exists(fixture_dir)) fixture_dir <- file.path("tests", "fixtures")
read_fix <- function(name) jsonlite::fromJSON(file.path(fixture_dir, name),
                                              simplifyVector = FALSE)

test_that("dssat_io::yyddd_to_date matches Python for all codes", {
  gold <- read_fix("yyddd_dates.json")
  for (code in names(gold)) {
    got <- yyddd_to_date(as.numeric(code))
    exp <- gold[[code]]
    if (is.null(exp)) {
      expect_true(is.na(got), info = sprintf("code %s expected NA", code))
    } else {
      expect_equal(format(got, "%Y-%m-%d"), exp, info = sprintf("code %s", code))
    }
  }
})

test_that("source error models match Python", {
  gold <- read_fix("source_error_models.json")
  ctor_for <- function(label) {
    if (startsWith(label, "field")) return(source_field_measurements(list()))
    if (startsWith(label, "uav")) return(source_uav_multispectral(list()))
    if (startsWith(label, "iot_sw")) return(source_soil_moisture_iot(list()))
    if (startsWith(label, "iot_tmean")) return(source_canopy_temperature(list()))
    if (startsWith(label, "sentinel")) return(source_sentinel2_lai(list()))
    if (startsWith(label, "modis")) return(source_modis_lai(list()))
    if (startsWith(label, "farm")) return(source_farm_phenology(list()))
    if (startsWith(label, "mgmt")) return(source_farm_management(list()))
    stop("no ctor for ", label)
  }
  for (label in names(gold)) {
    if (startsWith(label, "_")) next
    c <- gold[[label]]
    src <- ctor_for(label)
    got <- src_error_model(src, c$variable, c$value, c$metadata)
    expect_equal(got, as.numeric(c$sigma), tolerance = 1e-12, info = label)
  }
  # observation operator (satellite)
  op <- gold[["_obs_operator"]]
  expect_equal(.apply_obs_operator(list(), 3.0), as.numeric(op$identity), tolerance = 1e-12)
  expect_equal(.apply_obs_operator(list(obs_operator = list(scale = 0.9, offset = 0.2)), 3.0),
               as.numeric(op$scaled), tolerance = 1e-12)
})

test_that("inverse-variance fusion merge matches Python", {
  gold <- read_fix("fusion_inverse_variance.json")
  rows <- gold$input
  df <- do.call(rbind, lapply(rows, function(r) {
    data.frame(exp_id = r$exp_id, treatment = as.integer(r$treatment),
               variable = r$variable, kind = r$kind, date = as.Date(r$date),
               value = as.numeric(r$value), sigma = as.numeric(r$sigma),
               weight = as.numeric(r$weight), source = r$source,
               quality_flag = as.integer(r$quality_flag),
               spatial_res_m = as.numeric(r$spatial_res_m), stringsAsFactors = FALSE)
  }))
  fuser <- observation_fuser(list(), list(fusion = list(conflict_resolution = "inverse_variance")))
  merged <- fuser_resolve_conflicts(fuser, df)
  merged <- merged[order(merged$variable, merged$date), ]

  for (em in gold$merged) {
    row <- merged[merged$variable == em$variable &
                    format(merged$date, "%Y-%m-%d") == em$date, ]
    expect_equal(nrow(row), 1L, info = em$variable)
    expect_equal(row$value[1], as.numeric(em$value), tolerance = 1e-9, info = em$variable)
    expect_equal(row$sigma[1], as.numeric(em$sigma), tolerance = 1e-9, info = em$variable)
    expect_equal(row$source[1], em$source, info = em$variable)
  }
})

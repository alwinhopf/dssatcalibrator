# Cross-language parity tests for Phase 3: dssat_io parsers + spawn helpers.

fixture_dir <- file.path("..", "fixtures")
if (!dir.exists(fixture_dir)) fixture_dir <- file.path("tests", "fixtures")
read_fix <- function(name) jsonlite::fromJSON(file.path(fixture_dir, name),
                                              simplifyVector = FALSE)

test_that("parse_plantgro matches Python", {
  gold <- read_fix("parsers.json")$plantgro
  pg <- parse_plantgro(file.path(fixture_dir, "PlantGro.OUT"))
  expect_equal(nrow(pg), gold$nrow)
  expect_equal(as.numeric(pg$LAID), as.numeric(unlist(gold$LAID)), tolerance = 1e-9)
  expect_equal(as.numeric(pg$CWAD), as.numeric(unlist(gold$CWAD)), tolerance = 1e-9)
  expect_equal(as.integer(pg$treatment), as.integer(unlist(gold$treatment)))
  expect_equal(format(pg$date, "%Y-%m-%d"), unlist(gold$dates))
})

test_that("parse_evaluate matches Python", {
  gold <- read_fix("parsers.json")$evaluate
  ev <- parse_evaluate(file.path(fixture_dir, "Evaluate.OUT"))
  ev <- ev[order(ev$treatment, ev$variable), ]
  expect_equal(nrow(ev), length(gold))
  for (i in seq_along(gold)) {
    g <- gold[[i]]
    expect_equal(ev$treatment[i], as.integer(g$treatment), info = g$variable)
    expect_equal(ev$variable[i], g$variable)
    expect_equal(ev$sim[i], as.numeric(g$sim), tolerance = 1e-9, info = g$variable)
    expect_equal(ev$meas[i], as.numeric(g$meas), tolerance = 1e-9, info = g$variable)
  }
})

test_that("parse_summary matches Python", {
  gold <- read_fix("parsers.json")$summary
  sm <- parse_summary(file.path(fixture_dir, "Summary.OUT"))
  expect_equal(nrow(sm), gold$nrow)
  expect_equal(as.integer(sm$RUNNO), as.integer(unlist(gold$RUNNO)))
  expect_equal(as.numeric(sm$CWAM), as.numeric(unlist(gold$CWAM)), tolerance = 1e-9)
  expect_equal(as.numeric(sm$HWAM), as.numeric(unlist(gold$HWAM)), tolerance = 1e-9)
})

test_that("theta_hash matches Python (identical SHA-256 of identical blob)", {
  skip_if_not_installed("digest")
  gold <- read_fix("spawn_helpers.json")
  for (k in names(gold$thetas)) {
    expect_equal(theta_hash(gold$thetas[[k]]), gold$theta_hash[[k]], info = k)
  }
})

test_that("theta_hash supports FileX code values", {
  skip_if_not_installed("digest")
  expect_equal(theta_hash(list(irrig_code = "IR004", x = 1.0)), "06a431c288780c62")
  expect_equal(theta_hash(list(x = 1L, irrig_code = "IR004")), "06a431c288780c62")
  expect_false(theta_hash(list(irrig_code = "IR004", x = 1.0000001)) == "06a431c288780c62")
  expect_false(theta_hash(list(irrig_code = "IR005", x = 1.0)) == "06a431c288780c62")
})

test_that("write_dssbatch and normalize_treatments match Python", {
  gold <- read_fix("spawn_helpers.json")
  d <- tempfile(); dir.create(d)
  p <- write_dssbatch(d, "EXP0001.HMX", c(1, 2, 10))
  got <- paste0(paste(readLines(p, warn = FALSE), collapse = "\n"), "\n")
  expect_identical(got, gold$write_dssbatch)
  expect_equal(normalize_treatments(c(3, 1, 1, 10, 3)), as.integer(unlist(gold$normalize$dedup_order)))
  expect_equal(normalize_treatments(c(5)), as.integer(unlist(gold$normalize$single)))
})

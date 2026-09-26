test_that("spawn provenance changes for every simulation input", {
  root <- tempfile(); dir.create(root)
  paths <- list(root = root, soil = file.path(root, "Soil"),
                weather = file.path(root, "Weather"), genotype = file.path(root, "Genotype"))
  for (p in paths[-1]) dir.create(p)
  crop <- list(code = "MZ", genotype_stem = "MZCER048", model = "MZCER048")
  files <- c(file.path(root, "E1.MZX"), file.path(root, "fake.exe"),
             file.path(root, "DSSATPRO.V48"), file.path(root, "E1.MZA"),
             file.path(paths$genotype, paste0(crop$genotype_stem, c(".CUL", ".ECO", ".SPE"))),
             file.path(paths$soil, "CUSTOM.SOL"), file.path(paths$weather, "SECOND.WTH"))
  for (p in files) writeLines("original", p)
  cfg <- list()
  provenance <- function() .spawn_provenance(cfg, crop, list(), files[1],
                                            paths$genotype, paths, files[2], 1L, list())
  baseline <- provenance()
  for (p in files) {
    writeLines("modified", p)
    expect_false(identical(baseline, provenance()), info = p)
    writeLines("original", p)
  }
  cfg[["_planting_dates"]] <- list(E1 = "2001-04-01")
  expect_false(identical(baseline, provenance()))
  cfg <- list()
  crop$model <- "OTHER"
  expect_false(identical(baseline, provenance()))
})

test_that("spawn cache requires successful complete output and matching manifest", {
  root <- tempfile(); dir.create(root)
  paths <- list(root = root, soil = file.path(root, "Soil"),
                weather = file.path(root, "Weather"), genotype = file.path(root, "Genotype"))
  for (p in paths[-1]) dir.create(p)
  filex <- file.path(root, "E1.MZX")
  writeLines(c("*TREATMENTS", "@N R O C TNAME", " 1 1 1 0 first", "*CULTIVARS"), filex)
  exe <- file.path(root, "fake.exe"); writeLines("fake", exe)
  cfg <- list(source = list(hemp_dir = root), calibrator = list(cache_spawns = TRUE))
  crop <- list(code = "MZ", genotype_stem = "MZCER048", filex_ext = "MZX", model = "MZCER048")
  env <- new.env(parent = environment(spawn_and_run))
  env$resolve_dssat_paths <- function(cfg) paths
  env$.execution_backend <- function(cfg) "native"
  calls <- 0L; fail <- FALSE; complete <- TRUE
  env$.run_native_dssat <- function(run_dir, ...) {
    calls <<- calls + 1L
    if (fail) return("injected failure")
    writeLines("new output", file.path(run_dir, "PlantGro.OUT"))
    ""
  }
  env$.collect_core_outputs <- function(...) list(
    plantgro = if (complete) data.frame(treatment = 1L, LAID = 2) else data.frame(),
    evaluate = data.frame())
  spawn <- spawn_and_run; environment(spawn) <- env
  run <- function() spawn(list(x = 1), "E1", cfg, crop,
                           list(list(name = "x", group = "unused")),
                           file.path(root, "runs"), exe = exe)
  first <- run()
  expect_equal(first$status, "success")
  expect_equal(run()$status, "cached")
  expect_equal(calls, 1L)
  writeLines("weather change", file.path(paths$weather, "NEW.WTH"))
  second <- run()
  expect_equal(second$status, "success")
  expect_false(identical(first$run_dir, second$run_dir))
  expect_equal(calls, 2L)
  complete <- FALSE
  expect_equal(run()$status, "error")
  expect_false(file.exists(file.path(second$run_dir, "spawn_manifest.rds")))
  complete <- TRUE
  expect_equal(run()$status, "success")
  cfg$calibrator$cache_spawns <- FALSE; fail <- TRUE
  expect_equal(run()$status, "error")
  expect_false(file.exists(file.path(second$run_dir, "spawn_manifest.rds")))
  expect_false(file.exists(file.path(second$run_dir, "PlantGro.OUT")))
  cfg$calibrator$cache_spawns <- TRUE; fail <- FALSE
  expect_equal(run()$status, "success")
})

test_that("cache completeness fails closed", {
  expect_false(.spawn_outputs_complete(data.frame(), 1L))
  expect_false(.spawn_outputs_complete(data.frame(LAID = 1), 1L))
  expect_false(.spawn_outputs_complete(data.frame(treatment = 1L), c(1L, 2L)))
  expect_false(.spawn_outputs_complete(data.frame(treatment = c(1L, 2L)), 1L))
  expect_true(.spawn_outputs_complete(data.frame(treatment = c(2L, 1L)), c(1L, 2L)))
})

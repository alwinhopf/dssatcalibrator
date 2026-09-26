test_that("rank diagnostics identify separated chains and preserve undefined results", {
  set.seed(8)
  d <- expand.grid(step = 0:499, walker = 0:3)
  d$moving <- rnorm(nrow(d)) + d$walker * 10
  d$constant <- 1
  result <- chain_diagnostics(d, "moving")
  expect_gt(result$rhat, 2)
  expect_lt(result$ess, 20)
  expect_true(is.nan(chain_diagnostics(d, c("constant", "moving"))$rhat))
  expect_true(is.nan(chain_diagnostics(d, c("moving", "constant"))$rhat))
})

test_that("R SoilGrids acquisition passes the provider's named arguments", {
  root <- tempfile(); dir.create(root)
  env <- new.env(parent = environment(acquire_soil_profile))
  env$.dssatutils_required <- function() NULL
  env$.make_single_point_sf <- function(...) data.frame(ID = "0001")
  env$.get_utils_provider <- function(name) {
    expect_equal(name, "process_soils_soilgrids")
    function(grid_points, source_sol_file, output_csv_path, output_sol_dir, id_col) {
      expect_equal(source_sol_file, "source.SOL")
      expect_equal(id_col, "ID")
      writeLines("mock soil", file.path(output_sol_dir, "0001.SOL"))
    }
  }
  acquire <- acquire_soil_profile; environment(acquire) <- env
  cfg <- list(soil = list(provider = "dssatutils", source = "soilgrids",
                          source_sol_file = "source.SOL", cache_dir = root))
  out <- file.path(root, "output.SOL")
  expect_equal(acquire(cfg, "0001", 0, 0, out), out)
  expect_equal(readLines(out), "mock soil")
})

test_that("R weather acquisition never substitutes another cached station", {
  root <- tempfile(); dir.create(root); writeLines("other station", file.path(root, "OTHER.WTH"))
  env <- new.env(parent = environment(acquire_wth))
  env$.dssatutils_required <- function() NULL
  env$.make_single_point_sf <- function(...) data.frame(ID = "MISSING")
  env$.get_utils_provider <- function(name) function(...) NULL
  acquire <- acquire_wth; environment(acquire) <- env
  cfg <- list(weather = list(provider = "dssatutils", cache_dir = root))
  expect_error(acquire(cfg, "MISSING", 0, 0, "2020-01-01", "2020-01-02", file.path(root, "out.WTH")),
               "did not produce MISSING.WTH")
})

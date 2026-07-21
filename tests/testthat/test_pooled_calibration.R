pooled_cfg <- function() {
  list(
    experiments = list("E1", "E2"),
    parameters = list(
      genetic_species = list(
        SPEP = list(active = TRUE, min = 1, max = 3, start = 2, spe_key = "SPEP")
      ),
      genetic_cultivar = list(
        CULP = list(active = TRUE, min = 10, max = 20, start = 15, scope = "experiment")
      ),
      genetic_ecotype = list(
        ECOP = list(active = TRUE, min = 100, max = 200, start = 150, pooling = "per_experiment")
      )
    )
  )
}

test_that("experiment-scoped parameters expand per experiment", {
  cfg <- pooled_cfg()
  specs <- expand_parameter_specs(cfg, active_parameters(cfg))

  expect_identical(vapply(specs, function(s) s$name, character(1)),
                   c("SPEP", "CULP__E1", "CULP__E2", "ECOP__E1", "ECOP__E2"))
  expect_identical(vapply(specs, function(s) s$base_name, character(1)),
                   c("SPEP", "CULP", "CULP", "ECOP", "ECOP"))
  expect_identical(vapply(specs, function(s) s$scope, character(1)),
                   c("global", "experiment", "experiment", "experiment", "experiment"))
  expect_identical(vapply(specs, function(s) .cfg_get(s, "exp_id", NA_character_), character(1)),
                   c(NA_character_, "E1", "E2", "E1", "E2"))

  space <- parameter_space_from_config(cfg)
  expect_identical(space$names, c("SPEP", "CULP__E1", "CULP__E2", "ECOP__E1", "ECOP__E2"))
  expect_equal(as.numeric(space$start), c(2, 15, 15, 150, 150))
})

test_that("spawn partition uses only values for the current experiment", {
  cfg <- pooled_cfg()
  specs <- expand_parameter_specs(cfg, active_parameters(cfg))
  theta <- list(SPEP = 2.5, CULP__E1 = 11, CULP__E2 = 19, ECOP__E1 = 111, ECOP__E2 = 199)
  spawn_env <- environment(spawn_and_run)
  partition_theta <- get(".partition_theta", envir = spawn_env)
  effective_theta <- get(".effective_theta", envir = spawn_env)

  groups <- partition_theta(theta, specs, exp_id = "E1")
  expect_equal(groups$genetic_species, list(SPEP = 2.5))
  expect_equal(groups$genetic_cultivar, list(CULP = 11))
  expect_equal(groups$genetic_ecotype, list(ECOP = 111))

  effective <- effective_theta(theta, specs, exp_id = "E1")
  expect_equal(effective, list(SPEP = 2.5, CULP = 11, ECOP = 111))

  e2_groups <- partition_theta(theta, specs, exp_id = "E2")
  expect_equal(e2_groups$genetic_cultivar, list(CULP = 19))
  expect_equal(e2_groups$genetic_ecotype, list(ECOP = 199))
})

test_that("FileX overrides merge global and experiment records", {
  spawn_env <- environment(spawn_and_run)
  filex_overrides_for <- get(".filex_overrides_for", envir = spawn_env)
  cfg <- list(filex_overrides = list(
    all = list(list(section = "FIELDS", field = "ID_SOIL", value = "BASE")),
    CNKU2101 = list(list(section = "FIELDS", field = "WSTA", value = "CNKU2101"))
  ))

  cnku <- filex_overrides_for(cfg, "CNKU2101")
  expect_equal(length(cnku), 2L)
  expect_equal(cnku[[1]]$field, "ID_SOIL")
  expect_equal(cnku[[2]]$field, "WSTA")
  expect_equal(length(filex_overrides_for(cfg, "YUKU2101")), 1L)
})

test_that("validate_config accepts parameter scope and rejects unknown scope", {
  cfg <- pooled_cfg()
  expect_identical(validate_config(cfg), invisible(cfg))

  bad <- pooled_cfg()
  bad$parameters$genetic_cultivar$CULP$scope <- "plot"
  expect_error(validate_config(bad), "scope 'plot'", fixed = TRUE)
})

test_that("cultivar-scoped parameters expand, filter, and preserve fixed values", {
  cfg <- pooled_cfg()
  cfg$crops <- list(list(cultivar_anchor = "IB1", cultivar_anchors = list("IB1", "IB2")))
  cfg$parameters$genetic_cultivar$CULP$scope <- "cultivar"
  cfg$parameters$genetic_cultivar$CULP$start_by_cultivar <- list(IB1 = 12, IB2 = 18)
  cfg$parameters$genetic_ecotype$FIXED <- list(active = FALSE, fixed = TRUE,
                                               min = 1, max = 5, start = 3)
  specs <- expand_parameter_specs(cfg, active_parameters(cfg))
  expect_true(all(c("CULP__IB1", "CULP__IB2") %in%
                    vapply(specs, function(s) s$name, character(1))))
  starts <- setNames(vapply(specs, function(s) as.numeric(s$start), numeric(1)),
                     vapply(specs, function(s) s$name, character(1)))
  expect_equal(starts[c("CULP__IB1", "CULP__IB2")], c(CULP__IB1 = 12, CULP__IB2 = 18))

  fixed <- expand_parameter_specs(cfg, fixed_parameters(cfg))
  theta <- list(CULP__IB1 = 11, CULP__IB2 = 19)
  groups <- .partition_theta(theta, c(specs, fixed), exp_id = "E1", cultivars = "IB1")
  expect_equal(groups$genetic_cultivar_by_cultivar$IB1, list(CULP = 11))
  expect_null(groups$genetic_cultivar_by_cultivar$IB2)
  expect_equal(groups$genetic_ecotype$FIXED, 3)
})

test_that("treatment selections use distinct cache directories", {
  expect_equal(.treatment_run_key(c(2, 1, 2)), "T1-2")
  expect_equal(.treatment_run_key(2), "T2")
  expect_null(.treatment_run_key(NULL))
})

# Cross-language parity for the deterministic engine post-processing (GLUE).
# (Stochastic engines — MCMC, SMC-PF, NSGA-II, surrogate — are validated
# statistically rather than by fixture, per the agreed parity policy.)

fixture_dir <- file.path("..", "fixtures")
if (!dir.exists(fixture_dir)) fixture_dir <- file.path("tests", "fixtures")
read_fix <- function(name) jsonlite::fromJSON(file.path(fixture_dir, name), simplifyVector = TRUE)

test_that("run_glue matches Python (weights, ESS, threshold, best, behavioural)", {
  g <- read_fix("glue.json")
  inp <- g$input
  design <- data.frame(
    sample_id = inp$sample_id,
    P1 = as.numeric(inp$P1), P5 = as.numeric(inp$P5),
    score = vapply(inp$score, function(v) if (is.null(v)) Inf else as.numeric(v), numeric(1)),
    loglik = vapply(inp$loglik, function(v) if (is.null(v)) -Inf else as.numeric(v), numeric(1)),
    stringsAsFactors = FALSE
  )
  res <- run_glue(design, c("P1", "P5"),
                  list(method = list(bayesian = list(behavioural_quantile = g$behavioural_quantile))),
                  space = NULL)
  expect_equal(as.numeric(res$design$weight), as.numeric(g$weight), tolerance = 1e-12)
  expect_equal(res$ess, g$ess, tolerance = 1e-9)
  expect_equal(res$threshold, g$threshold, tolerance = 1e-9)
  expect_equal(res$best_sample_id, g$best_sample_id)
  expect_equal(nrow(res$behavioural), g$n_behavioural)
  expect_equal(as.numeric(res$best_theta$P1), as.numeric(g$best_theta$P1), tolerance = 1e-9)
})

#!/usr/bin/env Rscript
# Run the R parity suite WITHOUT devtools or installing the package.
#
#   Rscript run_r_tests.R
#
# It installs the few packages the parity tests need (if missing), sources the
# package's R/ files into the session, and runs testthat over tests/testthat.
# (The heavier engine packages — lhs, DEoptim, mco, DiceKriging, ranger, etc. —
# are only needed to RUN those engines, not for the parity tests.)

needed <- c("testthat", "jsonlite", "yaml", "digest")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  message("Installing: ", paste(miss, collapse = ", "))
  install.packages(miss, repos = "https://cloud.r-project.org")
}

suppressMessages(library(testthat))

# Source every package R file into the global environment so the test files can
# see the functions (no build/install step required).
r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
invisible(lapply(r_files, function(f) sys.source(f, envir = globalenv())))

cat(sprintf("Sourced %d R modules; running parity tests...\n\n", length(r_files)))
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)

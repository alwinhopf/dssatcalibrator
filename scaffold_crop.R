#!/usr/bin/env Rscript
# CLI: scaffold a new crop/cultivar from an analog DSSAT module (R twin of
# scaffold_crop.py).
#
#   Rscript scaffold_crop.R --dssat-dir C:/DSSAT48 \
#       --analog-stem SBGRO048 --new-stem QUGRO048 --new-code QU \
#       --source-anchor IB0001 --out-dir templates/quinoa

suppressMessages(library(dssatcalibrator))

.parse_args <- function(argv) {
  out <- list(); bool_flags <- c("no-spe"); i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (startsWith(a, "--")) {
      key <- substring(a, 3)
      if (key %in% bool_flags) { out[[key]] <- TRUE; i <- i + 1L }
      else { out[[key]] <- argv[i + 1L]; i <- i + 2L }
    } else i <- i + 1L
  }
  out
}

a <- .parse_args(commandArgs(trailingOnly = TRUE))
for (req in c("dssat-dir", "analog-stem", "new-stem", "new-code", "source-anchor")) {
  if (is.null(a[[req]])) stop(sprintf("missing required --%s", req))
}

res <- scaffold_crop(
  dssat_dir = a[["dssat-dir"]], analog_stem = a[["analog-stem"]], new_stem = a[["new-stem"]],
  new_code = a[["new-code"]], source_anchor = a[["source-anchor"]], new_anchor = a[["new-anchor"]],
  out_dir = a[["out-dir"]], spread = as.numeric(a$spread %||% 0.3), copy_spe = !isTRUE(a[["no-spe"]])
)

cat("Cloned genotype files:\n")
for (ext in names(res$files)) cat(sprintf("  .%s: %s\n", ext, res$files[[ext]]))
block_path <- file.path(res$out_dir, "parameters_block.yaml")
writeLines(res$parameters_yaml, block_path)
cat(sprintf("\nStarter parameter block -> %s\n", block_path))
cat("Paste it into your config and review every bound against literature.\n\n")
cat(res$parameters_yaml)

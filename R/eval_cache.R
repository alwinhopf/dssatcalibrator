# Persistent objective-result cache for expensive DSSAT evaluations.
#
# This is the R twin of python/dssatcalibrator/eval_cache.py. The cache stores a
# scored objective for a whole theta across the requested experiments, above the
# existing per-spawn DSSAT-output cache.

.EVAL_CACHE_SCHEMA_VERSION <- 1L

.eval_norm <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) {
    out <- x
    for (nm in names(out)) {
      if (inherits(out[[nm]], c("Date", "POSIXct", "POSIXt"))) out[[nm]] <- as.character(out[[nm]])
    }
    return(out)
  }
  if (is.list(x)) {
    if (!is.null(names(x))) {
      nms <- sort(names(x), method = "radix")
      return(setNames(lapply(nms, function(nm) .eval_norm(x[[nm]])), nms))
    }
    return(lapply(x, .eval_norm))
  }
  if (is.numeric(x)) {
    out <- round(as.numeric(x), 12)
    out[is.nan(out)] <- NA_real_
    return(out)
  }
  if (inherits(x, c("Date", "POSIXct", "POSIXt"))) return(as.character(x))
  x
}

.eval_digest <- function(x) {
  digest::digest(.eval_norm(x), algo = "sha1", serialize = TRUE)
}

.eval_file_fingerprint <- function(path) {
  if (is.null(path) || is.na(path) || !file.exists(path)) return(NULL)
  info <- file.info(path)
  if (is.na(info$size) || isTRUE(info$isdir)) return(NULL)
  list(
    path = normalizePath(path, mustWork = FALSE),
    size = as.numeric(info$size),
    sha1 = digest::digest(file = path, algo = "sha1", serialize = FALSE)
  )
}

.eval_cfg_fingerprint <- function(cfg) {
  keep <- list(
    parameters = .cfg_get(cfg, "parameters", list()),
    crops = .cfg_get(cfg, "crops", list()),
    source = .cfg_get(cfg, "source", list()),
    engine = .cfg_get(cfg, "engine", list()),
    objective = .cfg_get(cfg, "objective", list()),
    execution = .cfg_get(cfg, "execution", list()),
    templates = .cfg_get(cfg, "templates", list()),
    gating = .cfg_get(cfg, "gating", list()),
    management_options = .cfg_get(cfg, "management_options", list()),
    weather = .cfg_get(cfg, "weather", list()),
    soil = .cfg_get(cfg, "soil", list()),
    observation_sources = .cfg_get(cfg, "observation_sources", list()),
    fusion = .cfg_get(cfg, "fusion", list()),
    experiments = .cfg_get(cfg, "experiments", list()),
    cache_salt = .cfg_get(.cfg_get(cfg, "calibrator", list()), "evaluation_cache_salt", "")
  )
  .eval_digest(keep)
}

.eval_cache_enabled <- function(cfg) {
  isTRUE(.cfg_get(.cfg_get(cfg, "calibrator", list()), "cache_evaluations", TRUE)) &&
    requireNamespace("digest", quietly = TRUE)
}

.eval_cache_dir <- function(cfg) {
  ccfg <- .cfg_get(cfg, "calibrator", list())
  run_root <- file.path(.cfg_get(ccfg, "workdir", "results/_workdir"), .cfg_get(ccfg, "name", "run"))
  raw <- .cfg_get(ccfg, "evaluation_cache_dir", .cfg_get(ccfg, "eval_cache_dir", ""))
  if (!is.null(raw) && nzchar(as.character(raw))) {
    root <- as.character(raw)
    if (!grepl("^([A-Za-z]:)?[\\\\/]", root)) root <- file.path(run_root, root)
    return(root)
  }
  file.path(run_root, "evaluation_cache")
}

.evaluation_cache_from_setup <- function(cfg, crop, specs, experiments, treatments, obs_table, exe) {
  if (!.eval_cache_enabled(cfg)) {
    return(list(enabled = FALSE, root = NULL, context = NULL))
  }

  hemp_dir <- .cfg_get(.cfg_get(cfg, "source", list()), "hemp_dir", "")
  input_files <- list(exe = .eval_file_fingerprint(exe))
  input_files$genotype <- tryCatch({
    dssat_paths <- resolve_dssat_paths(cfg)
    geno_dir <- dssat_paths$genotype
    stem <- .cfg_get(crop, "genotype_stem", "")
    setNames(lapply(c("CUL", "ECO", "SPE"), function(ext) {
      .eval_file_fingerprint(file.path(geno_dir, sprintf("%s.%s", stem, ext)))
    }), c("CUL", "ECO", "SPE"))
  }, error = function(e) list(error = conditionMessage(e)))

  filex_ext <- .cfg_get(crop, "filex_ext", "")
  code <- .cfg_get(crop, "code", "")
  input_files$experiments <- setNames(lapply(experiments, function(exp) {
    list(
      filex = .eval_file_fingerprint(file.path(hemp_dir, sprintf("%s.%s", exp, filex_ext))),
      filea = .eval_file_fingerprint(file.path(hemp_dir, sprintf("%s.%sA", exp, code))),
      filet = .eval_file_fingerprint(file.path(hemp_dir, sprintf("%s.%sT", exp, code)))
    )
  }), experiments)

  context <- list(
    schema = .EVAL_CACHE_SCHEMA_VERSION,
    cfg = .eval_cfg_fingerprint(cfg),
    crop = crop,
    specs = specs,
    treatments = treatments,
    obs = .eval_digest(obs_table),
    inputs = input_files
  )
  list(enabled = TRUE, root = .eval_cache_dir(cfg), context = context)
}

.eval_cache_key <- function(cache, theta, experiments) {
  digest::digest(list(
    context = cache$context,
    experiments = as.character(experiments),
    theta = .eval_norm(theta)
  ), algo = "sha1", serialize = TRUE)
}

.eval_cache_path <- function(cache, key) {
  file.path(cache$root, substr(key, 1L, 2L), paste0(key, ".rds"))
}

.eval_cache_get <- function(cache, key) {
  if (!isTRUE(cache$enabled)) return(NULL)
  path <- .eval_cache_path(cache, key)
  if (!file.exists(path)) return(NULL)
  tryCatch({
    payload <- readRDS(path)
    if (!identical(payload$schema, .EVAL_CACHE_SCHEMA_VERSION)) return(NULL)
    payload$result
  }, error = function(e) NULL)
}

.eval_cache_put <- function(cache, key, result) {
  if (!isTRUE(cache$enabled)) return(invisible(FALSE))
  path <- .eval_cache_path(cache, key)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(list(schema = .EVAL_CACHE_SCHEMA_VERSION, result = result), tmp)
  ok <- file.rename(tmp, path)
  if (!isTRUE(ok)) {
    if (file.exists(path)) unlink(path)
    ok <- file.rename(tmp, path)
  }
  if (!isTRUE(ok)) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(isTRUE(ok))
}

# Spawn = deterministic function of (theta, experiment) -> a DSSAT run + outputs.
# R twin of python/dssatcalibrator/spawn.py.
#
# A spawn copies the base FileX and genotype files into an isolated run dir,
# writes perturbed coefficients (and optionally management / initial conditions),
# runs dscsm048, and parses PlantGro.OUT + Evaluate.OUT.

GENETIC_GROUPS <- c("genetic_cultivar")
BATCH_FILE <- "DSSBatch.V48"

# Python-repr-compatible number formatting for theta_hash blob parity.
.py_num <- function(x) {
  v <- round(as.numeric(x), 6)
  if (abs(v - round(v)) < 1e-9 && abs(v) < 1e15) return(sprintf("%.1f", v))
  s <- formatC(v, format = "f", digits = 6)
  s <- sub("0+$", "", s); s <- sub("\\.$", ".0", s)
  s
}

#' Stable 10-char hash of a (rounded) theta. Mirrors spawn.py:theta_hash.
#' Builds the identical JSON blob (sorted keys, ", "/": " separators) and SHA-1s it.
#' @export
theta_hash <- function(theta) {
  keys <- sort(names(theta), method = "radix")
  parts <- vapply(keys, function(k) sprintf('"%s": %s', k, .py_num(theta[[k]])), character(1))
  blob <- paste0("{", paste(parts, collapse = ", "), "}")
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("theta_hash requires the 'digest' package.")
  }
  substr(digest::digest(blob, algo = "sha1", serialize = FALSE), 1, 10)
}

#' Return the treatment numbers from a FileX `*TREATMENTS` section.
#' Mirrors spawn.py:parse_treatments.
#' @export
parse_treatments <- function(filex_path) {
  lines <- readLines(filex_path, warn = FALSE)
  start <- which(startsWith(lines, "*TREATMENTS"))
  if (length(start) == 0) return(1L)
  start <- start[1]
  trts <- integer(0)
  if (start < length(lines)) {
    for (ln in lines[(start + 1L):length(lines)]) {
      if (startsWith(ln, "*")) break
      if (grepl("^\\s*\\d", ln) && !startsWith(trimws(ln, which = "left"), "@")) {
        trts <- c(trts, as.integer(strsplit(trimws(ln), "\\s+")[[1]][1]))
      }
    }
  }
  if (length(trts) == 0) 1L else trts
}

#' Write a DSSAT B-mode batch file listing the requested treatments.
#' Mirrors spawn.py:write_dssbatch (byte-identical layout).
#' @export
write_dssbatch <- function(run_dir, filex_name, treatments) {
  header <- paste0("$BATCH(CALIB)\n!\n@FILEX", strrep(" ", 86),
                   "TRTNO     RP     SQ     OP     CO\n")
  rows <- paste0(vapply(treatments, function(t)
    sprintf("%-93s%6d      1      0      1      0\n", filex_name, as.integer(t)), character(1)),
    collapse = "")
  p <- file.path(run_dir, BATCH_FILE)
  cat(header, rows, file = p, sep = "")
  p
}

.execution_backend <- function(cfg) {
  backend <- tolower(as.character(.cfg_get(.cfg_get(cfg, "execution", list()), "backend", "native")))
  if (!(backend %in% c("native", "dssatengine"))) {
    stop("execution.backend must be 'native' or 'dssatengine'.")
  }
  backend
}

#' Normalize a treatment selection (dedupe, positive, order-preserving).
#' Mirrors spawn.py:_normalize_treatments (native path).
#' @export
normalize_treatments <- function(treatments, backend = "native") {
  if (backend == "dssatengine" && requireNamespace("dssatengine", quietly = TRUE)) {
    return(dssatengine::normalize_treatment_list(1, 1, treatment_list = treatments))
  }
  seen <- integer(0); out <- integer(0)
  for (value in treatments) {
    trt <- as.integer(value)
    if (trt < 1) stop("Treatment IDs must be positive integers.")
    if (!(trt %in% seen)) { seen <- c(seen, trt); out <- c(out, trt) }
  }
  if (length(out) == 0) stop("No valid treatments selected.")
  out
}

# Split a flat theta into per-group update lists using the param specs.
.partition_theta <- function(theta, param_specs) {
  group_of <- list()
  for (p in param_specs) group_of[[p$name]] <- p$group
  groups <- list()
  for (name in names(theta)) {
    g <- if (!is.null(group_of[[name]])) group_of[[name]] else "genetic_cultivar"
    if (is.null(groups[[g]])) groups[[g]] <- list()
    groups[[g]][[name]] <- theta[[name]]
  }
  groups
}

# Native DSSAT subprocess. Returns "" on success, else an error/timeout message.
.run_native_dssat <- function(run_dir, exe, model, timeout) {
  old <- setwd(run_dir); on.exit(setwd(old), add = TRUE)
  out <- tryCatch(
    system2(exe, args = c(model, "B", BATCH_FILE), stdout = TRUE, stderr = TRUE,
            timeout = timeout),
    error = function(e) structure(character(0), status = 1L, msg = conditionMessage(e)))
  status <- attr(out, "status")
  if (length(out) > 0) {
    writeLines(out, file.path(run_dir, "dssat_B_stdout_stderr.log"))
  }
  if (!is.null(status) && status != 0) {
    tail <- if (length(out)) paste(utils::tail(out, 12), collapse = " | ") else "<no stdout/stderr captured>"
    return(sprintf("DSSAT exited with status %s. Tail: %s", status, tail))
  }
  ""
}

#' SpawnResult constructor (S3). Mirrors spawn.py:SpawnResult.
#' @export
spawn_result <- function(status, run_dir, theta, plantgro = data.frame(),
                         evaluate = data.frame(), message = "") {
  structure(list(status = status, run_dir = run_dir, theta = theta,
                 plantgro = plantgro, evaluate = evaluate, message = message),
            class = "spawn_result")
}

#' Materialize and run one spawn; return parsed PlantGro + Evaluate tables.
#' Mirrors spawn.py:spawn_and_run. (DSSAT execution requires the binary; the
#' file-staging and parsing logic mirror the Python path.)
#' @export
spawn_and_run <- function(theta, exp_id, cfg, crop, param_specs, run_root,
                          treatments = NULL, exe, timeout = 600) {
  backend <- .execution_backend(cfg)
  dssat_paths <- resolve_dssat_paths(cfg)
  hemp_dir <- cfg$source$hemp_dir
  geno_dir <- dssat_paths$genotype
  stem <- crop$genotype_stem
  ext <- crop$filex_ext
  code <- crop$code

  run_dir <- file.path(run_root, exp_id, paste0("s_", theta_hash(theta)))
  pg_path <- file.path(run_dir, "PlantGro.OUT")

  if (isTRUE(.cfg_get(cfg$calibrator, "cache_spawns", TRUE)) &&
      file.exists(pg_path) && file.info(pg_path)$size > 0) {
    return(spawn_result("cached", run_dir, theta,
                        plantgro = parse_plantgro(pg_path),
                        evaluate = parse_evaluate(file.path(run_dir, "Evaluate.OUT"))))
  }

  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  for (e in c("CUL", "ECO", "SPE")) {
    src <- file.path(geno_dir, sprintf("%s.%s", stem, e))
    if (file.exists(src)) file.copy(src, file.path(run_dir, sprintf("%s.%s", stem, e)), overwrite = TRUE)
  }

  filex_name <- sprintf("%s.%s", exp_id, ext)
  file.copy(file.path(hemp_dir, filex_name), file.path(run_dir, filex_name), overwrite = TRUE)
  for (obs_ext in c(sprintf("%sA", code), sprintf("%sT", code))) {
    src <- file.path(hemp_dir, sprintf("%s.%s", exp_id, obs_ext))
    if (file.exists(src)) file.copy(src, file.path(run_dir, sprintf("%s.%s", exp_id, obs_ext)), overwrite = TRUE)
  }

  groups <- .partition_theta(theta, param_specs)
  cul_updates <- list()
  for (g in GENETIC_GROUPS) cul_updates <- modifyList(cul_updates, .cfg_get(groups, g, list()))
  if (length(cul_updates) > 0) {
    edit_cultivar(file.path(run_dir, sprintf("%s.CUL", stem)), crop$cultivar_anchor, cul_updates)
  }
  eco_updates <- .cfg_get(groups, "genetic_ecotype", list())
  if (length(eco_updates) > 0) {
    edit_ecotype(file.path(run_dir, sprintf("%s.ECO", stem)), crop$ecotype, eco_updates)
  }

  spe_updates <- .cfg_get(groups, "genetic_species", list())
  if (length(spe_updates) > 0 &&
      tolower(as.character(.cfg_get(.cfg_get(cfg, "gating", list()), "species", "blocked"))) == "free") {
    spe_file <- file.path(run_dir, sprintf("%s.SPE", stem))
    if (file.exists(spe_file)) {
      updates <- list()
      for (name in names(spe_updates)) {
        spec <- Find(function(s) s$name == name, param_specs)
        key <- if (!is.null(spec) && !is.null(spec$spe_key)) spec$spe_key else name
        updates[[key]] <- spe_updates[[name]]
      }
      edit_species(spe_file, updates)
    }
  }

  mgt_updates <- .cfg_get(groups, "management", list())
  init_updates <- .cfg_get(groups, "initial_conditions", list())
  mgt_fields <- list()
  for (name in names(mgt_updates)) {
    spec <- Find(function(s) s$name == name, param_specs)
    if (!is.null(spec) && !is.null(spec$dssat)) mgt_fields[[spec$dssat]] <- mgt_updates[[name]]
  }
  pdate <- .cfg_get(.cfg_get(cfg, "_planting_dates", list()), exp_id, NULL)
  if (!is.null(pdate)) {
    ts <- as.Date(pdate)
    yy <- as.integer(format(ts, "%y")); doy <- as.integer(format(ts, "%j"))
    mgt_fields[["PDATE"]] <- as.integer(sprintf("%02d%03d", yy, doy))
  }
  if (length(mgt_fields) > 0 || length(init_updates) > 0) {
    edit_filex(file.path(run_dir, filex_name), mgt_fields, init_updates)
  }

  # (soil/weather acquisition + editing mirror the Python path; acquisition is
  # delegated to dssatutils when configured. Omitted providers default to DSSAT's
  # central files.)
  soil_updates <- .cfg_get(groups, "soil", list())
  weather_updates <- .cfg_get(groups, "weather", list())
  if (length(soil_updates) > 0 || length(weather_updates) > 0) {
    fields <- parse_fields(file.path(run_dir, filex_name))
    if (length(soil_updates) > 0 && !is.null(fields$id_soil)) {
      local_sol <- file.path(run_dir, "SOIL.SOL")
      src_sol <- if (file.exists(local_sol)) local_sol else file.path(dssat_paths$soil, "SOIL.SOL")
      pid <- fields$id_soil
      if (file.exists(src_sol)) {
        layer_mults <- list(); profile_sets <- list()
        for (name in names(soil_updates)) {
          spec <- Find(function(s) s$name == name, param_specs)
          if (is.null(spec) || is.null(spec$dssat)) next
          if (!is.null(spec$op) && spec$op == "set") profile_sets[[spec$dssat]] <- soil_updates[[name]]
          else layer_mults[[spec$dssat]] <- soil_updates[[name]]
        }
        if (src_sol != local_sol) cat(extract_soil_profile(src_sol, pid), file = local_sol)
        edit_soil(local_sol, pid, layer_mults = layer_mults, profile_sets = profile_sets)
      }
    }
    if (length(weather_updates) > 0 && !is.null(fields$wsta)) {
      wsta <- fields$wsta
      local_wth <- file.path(run_dir, sprintf("%s.WTH", wsta))
      src_wth <- if (file.exists(local_wth)) local_wth else file.path(dssat_paths$weather, sprintf("%s.WTH", wsta))
      if (file.exists(src_wth)) {
        ops <- list()
        for (name in names(weather_updates)) {
          spec <- Find(function(s) s$name == name, param_specs)
          if (is.null(spec) || is.null(spec$dssat)) next
          ops[[spec$dssat]] <- list(if (!is.null(spec$op)) spec$op else "mult", weather_updates[[name]])
        }
        if (src_wth != local_wth) file.copy(src_wth, local_wth, overwrite = TRUE)
        edit_weather(local_wth, ops)
      }
    }
  }

  if (is.null(treatments)) treatments <- parse_treatments(file.path(run_dir, filex_name))
  treatments <- normalize_treatments(treatments, backend)
  write_dssbatch(run_dir, filex_name, treatments)

  run_error <- .run_native_dssat(run_dir, exe, crop$model, timeout)
  if (nzchar(run_error)) return(spawn_result("error", run_dir, theta, message = run_error))

  if (!file.exists(pg_path) || file.info(pg_path)$size == 0) {
    return(spawn_result("error", run_dir, theta, message = "no PlantGro.OUT produced"))
  }
  pg <- parse_plantgro(pg_path)
  ev <- parse_evaluate(file.path(run_dir, "Evaluate.OUT"))
  spawn_result("success", run_dir, theta, plantgro = pg, evaluate = ev)
}

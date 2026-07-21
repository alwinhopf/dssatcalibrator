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

.py_json_value <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 1L && !is.na(v) && !isTRUE(is.logical(x))) return(.py_num(x))
  jsonlite::toJSON(as.character(x), auto_unbox = TRUE)
}

#' Stable 10-char hash of a (rounded) theta. Mirrors spawn.py:theta_hash.
#' Builds the identical JSON blob (sorted keys, ", "/": " separators) and SHA-1s it.
#' @export
theta_hash <- function(theta) {
  keys <- sort(names(theta), method = "radix")
  parts <- vapply(keys, function(k) sprintf('"%s": %s', k, .py_json_value(theta[[k]])), character(1))
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

#' Return cultivar codes listed in a FileX `*CULTIVARS` section.
#' Mirrors spawn.py:parse_cultivars.
#' @export
parse_cultivars <- function(filex_path) {
  lines <- readLines(filex_path, warn = FALSE)
  start <- which(startsWith(lines, "*CULTIVARS"))
  if (!length(start)) return(character(0))
  header <- NULL; out <- character(0)
  if (start[1] >= length(lines)) return(out)
  for (ln in lines[(start[1] + 1L):length(lines)]) {
    if (startsWith(ln, "*")) break
    if (startsWith(trimws(ln, which = "left"), "@")) {
      header <- strsplit(sub("^@", "", trimws(ln)), "\\s+")[[1]]
      next
    }
    if (!is.null(header) && grepl("^\\s*\\d", ln)) {
      values <- strsplit(trimws(ln), "\\s+")[[1]]
      row <- setNames(as.list(values[seq_len(min(length(values), length(header)))]),
                      header[seq_len(min(length(values), length(header)))])
      code <- .cfg_get(row, "INGENO", .cfg_get(row, "CULTIVAR", .cfg_get(row, "VAR#", NULL)))
      if (!is.null(code) && !(code %in% out)) out <- c(out, as.character(code))
    }
  }
  out
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
.spec_applies <- function(spec, exp_id = NULL, cultivars = NULL) {
  scope <- .cfg_get(spec, "scope", "global")
  if (identical(scope, "experiment")) {
    return(is.null(exp_id) || as.character(.cfg_get(spec, "exp_id", "")) == as.character(exp_id))
  }
  if (identical(scope, "cultivar")) {
    return(is.null(cultivars) || as.character(.cfg_get(spec, "cultivar", "")) %in% as.character(cultivars))
  }
  TRUE
}

.effective_theta <- function(theta, param_specs, exp_id = NULL, cultivars = NULL) {
  out <- list()
  matched <- character(0)
  for (spec in param_specs) {
    name <- spec$name
    if (!.spec_applies(spec, exp_id, cultivars)) next
    value <- theta[[name]]
    if (is.null(value) && isTRUE(spec$fixed)) value <- spec$start
    if (is.null(value)) next
    base <- .cfg_get(spec, "base_name", name)
    if (identical(.cfg_get(spec, "scope", "global"), "cultivar")) {
      base <- sprintf("%s__%s", base, .cfg_get(spec, "cultivar", ""))
    }
    out[[base]] <- value
    matched <- c(matched, name)
  }
  if (length(matched) == 0 && length(param_specs) == 0) theta else out
}

.partition_theta <- function(theta, param_specs, exp_id = NULL, cultivars = NULL) {
  groups <- list()
  matched <- character(0)
  for (spec in param_specs) {
    name <- spec$name
    if (!.spec_applies(spec, exp_id, cultivars)) next
    value <- theta[[name]]
    if (is.null(value) && isTRUE(spec$fixed)) value <- spec$start
    if (is.null(value)) next
    g <- .cfg_get(spec, "group", "genetic_cultivar")
    base <- .cfg_get(spec, "base_name", name)
    if (identical(.cfg_get(spec, "scope", "global"), "cultivar")) {
      scoped_group <- paste0(g, "_by_cultivar")
      anchor <- as.character(.cfg_get(spec, "cultivar", ""))
      if (is.null(groups[[scoped_group]])) groups[[scoped_group]] <- list()
      if (is.null(groups[[scoped_group]][[anchor]])) groups[[scoped_group]][[anchor]] <- list()
      groups[[scoped_group]][[anchor]][[base]] <- value
    } else {
      if (is.null(groups[[g]])) groups[[g]] <- list()
      groups[[g]][[base]] <- value
    }
    matched <- c(matched, name)
  }
  if (length(matched) == 0 && length(param_specs) == 0) {
    for (name in names(theta)) {
      if (is.null(groups[["genetic_cultivar"]])) groups[["genetic_cultivar"]] <- list()
      groups[["genetic_cultivar"]][[name]] <- theta[[name]]
    }
  }
  groups
}

.cultivar_ecotype_map <- function(crop) {
  mapping <- .cfg_get(crop, "cultivar_ecotypes", list())
  mapping <- as.list(mapping)
  anchor <- .cfg_get(crop, "cultivar_anchor", NULL)
  ecotype <- .cfg_get(crop, "ecotype", NULL)
  if (!is.null(anchor) && !is.null(ecotype) && is.null(mapping[[as.character(anchor)]])) {
    mapping[[as.character(anchor)]] <- as.character(ecotype)
  }
  mapping
}

.treatment_run_key <- function(treatments) {
  if (is.null(treatments)) return(NULL)
  vals <- sort(unique(as.integer(treatments)))
  if (!length(vals)) NULL else paste0("T", paste(vals, collapse = "-"))
}

.stamp_single_treatment <- function(outputs, treatments) {
  vals <- sort(unique(as.integer(treatments %||% integer(0))))
  if (length(vals) != 1L) return(outputs)
  trt <- vals[1]
  lapply(outputs, function(df) {
    if (!is.data.frame(df) || !nrow(df)) return(df)
    if ("treatment" %in% names(df)) {
      missing <- is.na(df$treatment)
      df$treatment[missing] <- trt
      df$treatment <- as.integer(df$treatment)
    } else {
      df$treatment <- trt
    }
    df
  })
}

.collect_core_outputs <- function(run_dir, treatments = NULL) {
  outputs <- list(
    plantgro = parse_plantgro(file.path(run_dir, "PlantGro.OUT")),
    evaluate = parse_evaluate(file.path(run_dir, "Evaluate.OUT")),
    summary = parse_summary(file.path(run_dir, "Summary.OUT"))
  )
  .stamp_single_treatment(outputs, treatments)
}

.filex_overrides_for <- function(cfg, exp_id) {
  block <- .cfg_get(cfg, "filex_overrides", list())
  updates <- list()
  add_updates <- function(records) {
    if (is.null(records)) return()
    for (rec in records) {
      if (is.list(rec)) updates[[paste0("override_", length(updates) + 1L)]] <<- rec
    }
  }
  add_updates(.cfg_get(block, "all", list()))
  add_updates(.cfg_get(block, exp_id, list()))
  updates
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
    writeLines(out, "dssat_B_stdout_stderr.log")
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
                         evaluate = data.frame(), outputs = list(), message = "") {
  structure(list(status = status, run_dir = run_dir, theta = theta,
                 plantgro = plantgro, evaluate = evaluate, outputs = outputs,
                 message = message),
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
  filex_name <- sprintf("%s.%s", exp_id, ext)
  source_filex <- file.path(hemp_dir, filex_name)
  exp_cultivars <- parse_cultivars(source_filex)

  effective_theta <- .effective_theta(theta, param_specs, exp_id, exp_cultivars)
  treatment_key <- .treatment_run_key(treatments)
  run_dir <- file.path(run_root, exp_id)
  if (!is.null(treatment_key)) run_dir <- file.path(run_dir, treatment_key)
  run_dir <- file.path(run_dir, paste0("s_", theta_hash(effective_theta)))
  pg_path <- file.path(run_dir, "PlantGro.OUT")

  if (isTRUE(.cfg_get(cfg$calibrator, "cache_spawns", TRUE)) &&
      file.exists(pg_path) && file.info(pg_path)$size > 0) {
    outputs <- .collect_core_outputs(run_dir, treatments)
    return(spawn_result("cached", run_dir, theta,
                        plantgro = outputs$plantgro,
                        evaluate = outputs$evaluate,
                        outputs = outputs))
  }

  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  for (profile_name in c("DSSATPRO.L48", "DSSATPRO.V48", "DSSATPRO.v48", "DSCSM048.CTR")) {
    src <- file.path(dssat_paths$root, profile_name)
    if (file.exists(src)) file.copy(src, file.path(run_dir, profile_name), overwrite = TRUE)
  }

  for (e in c("CUL", "ECO", "SPE")) {
    src <- file.path(geno_dir, sprintf("%s.%s", stem, e))
    if (file.exists(src)) file.copy(src, file.path(run_dir, sprintf("%s.%s", stem, e)), overwrite = TRUE)
  }

  if (!file.exists(source_filex)) stop(sprintf("FileX not found: %s", source_filex))
  if (!file.copy(source_filex, file.path(run_dir, filex_name), overwrite = TRUE)) {
    stop(sprintf("Could not copy FileX into run directory: %s", source_filex))
  }
  filex_overrides <- .filex_overrides_for(cfg, exp_id)
  if (length(filex_overrides) > 0) {
    edit_filex(file.path(run_dir, filex_name), filex_overrides, list())
  }
  for (obs_ext in c(sprintf("%sA", code), sprintf("%sT", code))) {
    src <- file.path(hemp_dir, sprintf("%s.%s", exp_id, obs_ext))
    if (file.exists(src)) file.copy(src, file.path(run_dir, sprintf("%s.%s", exp_id, obs_ext)), overwrite = TRUE)
  }

  groups <- .partition_theta(theta, param_specs, exp_id, exp_cultivars)
  cul_updates <- list()
  for (g in GENETIC_GROUPS) cul_updates <- modifyList(cul_updates, .cfg_get(groups, g, list()))
  if (length(cul_updates) > 0) {
    anchors <- unlist(.cfg_get(crop, "cultivar_anchors", list(crop$cultivar_anchor)))
    for (anchor in anchors) {
      edit_cultivar(file.path(run_dir, sprintf("%s.CUL", stem)), anchor, cul_updates)
    }
  }
  for (anchor in names(.cfg_get(groups, "genetic_cultivar_by_cultivar", list()))) {
    edit_cultivar(file.path(run_dir, sprintf("%s.CUL", stem)), anchor,
                  groups$genetic_cultivar_by_cultivar[[anchor]])
  }
  eco_updates <- .cfg_get(groups, "genetic_ecotype", list())
  if (length(eco_updates) > 0) {
    edit_ecotype(file.path(run_dir, sprintf("%s.ECO", stem)), crop$ecotype, eco_updates)
  }
  cultivar_ecotypes <- .cultivar_ecotype_map(crop)
  for (anchor in names(.cfg_get(groups, "genetic_ecotype_by_cultivar", list()))) {
    eco_anchor <- cultivar_ecotypes[[anchor]]
    if (is.null(eco_anchor)) {
      stop(sprintf("No ecotype mapping for cultivar '%s'. Add crops[].cultivar_ecotypes.", anchor))
    }
    edit_ecotype(file.path(run_dir, sprintf("%s.ECO", stem)), eco_anchor,
                 groups$genetic_ecotype_by_cultivar[[anchor]])
  }

  spe_updates <- .cfg_get(groups, "genetic_species", list())
  if (length(spe_updates) > 0 &&
      tolower(as.character(.cfg_get(.cfg_get(cfg, "gating", list()), "species", "blocked"))) == "free") {
    spe_file <- file.path(run_dir, sprintf("%s.SPE", stem))
    if (file.exists(spe_file)) {
      updates <- list()
      for (name in names(spe_updates)) {
        spec <- Find(function(s) .cfg_get(s, "base_name", s$name) == name, param_specs)
        key <- if (!is.null(spec) && !is.null(spec$spe_key)) spec$spe_key else name
        if (!is.null(spec) && (!is.null(spec$spe_index) || !is.null(spec$token_index))) {
          updates[[key]] <- list(
            value = spe_updates[[name]],
            index = as.integer(.cfg_get(spec, "spe_index", .cfg_get(spec, "token_index", 0L)))
          )
        } else {
          updates[[key]] <- spe_updates[[name]]
        }
      }
      edit_species(spe_file, updates)
    }
  }

  mgt_updates <- .cfg_get(groups, "management", list())
  init_updates <- .cfg_get(groups, "initial_conditions", list())
  filex_update_from_spec <- function(name, val, spec, default_section) {
    if (is.null(spec)) return(val)
    section <- .cfg_get(spec, "section", .cfg_get(spec, "filex_section", NULL))
    field <- .cfg_get(spec, "field", .cfg_get(spec, "filex_field", .cfg_get(spec, "dssat", NULL)))
    is_soil_water_mult <- identical(default_section, "INITIAL CONDITIONS") &&
      identical(name, "initial_soil_water_mult")
    if (is_soil_water_mult && is.null(field)) field <- "SH2O"
    generic_keys <- c("header_prefix", "row", "treatment", "trt", "trtno",
                      "clip_01", "required", "type", "format")
    if (is.null(section) && !is.null(field)) {
      uses_generic <- any(generic_keys %in% names(spec)) ||
        tolower(as.character(.cfg_get(spec, "op", "set"))) != "set"
      if (identical(default_section, "PLANTING DETAILS") && !uses_generic) return(val)
      out <- list(section = default_section, field = field, value = val,
                  op = .cfg_get(spec, "op", if (is_soil_water_mult) "mult" else "set"))
      if (is_soil_water_mult && is.null(spec$clip_01)) out$clip_01 <- TRUE
      for (key in generic_keys) {
        if (!is.null(spec[[key]])) out[[key]] <- spec[[key]]
      }
      return(out)
    }
    out <- list(section = section %||% default_section,
                field = field %||% name,
                value = val,
                op = .cfg_get(spec, "op", if (is_soil_water_mult) "mult" else "set"))
    if (is_soil_water_mult && is.null(spec$clip_01)) out$clip_01 <- TRUE
    for (key in generic_keys) {
      if (!is.null(spec[[key]])) out[[key]] <- spec[[key]]
    }
    out
  }
  mgt_fields <- list()
  for (name in names(mgt_updates)) {
    spec <- Find(function(s) .cfg_get(s, "base_name", s$name) == name, param_specs)
    if (!is.null(spec)) {
      key <- .cfg_get(spec, "dssat", .cfg_get(spec, "field", .cfg_get(spec, "filex_field", name)))
      mgt_fields[[key]] <- filex_update_from_spec(name, mgt_updates[[name]], spec, "PLANTING DETAILS")
    }
  }
  pdate <- .cfg_get(.cfg_get(cfg, "_planting_dates", list()), exp_id, NULL)
  if (!is.null(pdate)) {
    ts <- as.Date(pdate)
    yy <- as.integer(format(ts, "%y")); doy <- as.integer(format(ts, "%j"))
    mgt_fields[["PDATE"]] <- as.integer(sprintf("%02d%03d", yy, doy))
  }
  if (length(mgt_fields) > 0 || length(init_updates) > 0) {
    init_fields <- list()
    for (name in names(init_updates)) {
      spec <- Find(function(s) .cfg_get(s, "base_name", s$name) == name, param_specs)
      if (!is.null(spec)) {
        key <- .cfg_get(spec, "dssat", .cfg_get(spec, "field", .cfg_get(spec, "filex_field", name)))
        init_fields[[key]] <- filex_update_from_spec(name, init_updates[[name]], spec, "INITIAL CONDITIONS")
      } else {
        init_fields[[name]] <- init_updates[[name]]
      }
    }
    edit_filex(file.path(run_dir, filex_name), mgt_fields, init_fields)
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
      pid <- fields$id_soil
      candidates <- c(if (file.exists(local_sol)) local_sol else character(0),
                      file.path(dssat_paths$soil, "SOIL.SOL"),
                      list.files(dssat_paths$soil, pattern = "[.]SOL$", full.names = TRUE))
      profile_text <- NULL
      for (src_sol in unique(candidates)) {
        if (!file.exists(src_sol)) next
        profile_text <- tryCatch(extract_soil_profile(src_sol, pid), error = function(e) NULL)
        if (!is.null(profile_text)) break
      }
      if (is.null(profile_text)) {
        stop(sprintf("Soil profile '%s' not found in DSSAT soil directory %s", pid, dssat_paths$soil))
      }
      layer_mults <- list(); profile_sets <- list()
      for (name in names(soil_updates)) {
        spec <- Find(function(s) .cfg_get(s, "base_name", s$name) == name, param_specs)
        if (is.null(spec) || is.null(spec$dssat)) next
        if (!is.null(spec$op) && spec$op == "set") profile_sets[[spec$dssat]] <- soil_updates[[name]]
        else layer_mults[[spec$dssat]] <- soil_updates[[name]]
      }
      cat(profile_text, file = local_sol)
      edit_soil(local_sol, pid, layer_mults = layer_mults, profile_sets = profile_sets)
    }
    if (length(weather_updates) > 0 && !is.null(fields$wsta)) {
      wsta <- fields$wsta
      local_wth <- file.path(run_dir, sprintf("%s.WTH", wsta))
      src_wth <- if (file.exists(local_wth)) local_wth else file.path(dssat_paths$weather, sprintf("%s.WTH", wsta))
      if (file.exists(src_wth)) {
        ops <- list()
        for (name in names(weather_updates)) {
          spec <- Find(function(s) .cfg_get(s, "base_name", s$name) == name, param_specs)
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
  outputs <- .collect_core_outputs(run_dir, treatments)
  spawn_result("success", run_dir, theta, plantgro = outputs$plantgro,
               evaluate = outputs$evaluate, outputs = outputs)
}

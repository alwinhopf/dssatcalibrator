# Read observed data into one tidy long-format Observations table.
# R twin of python/dssatcalibrator/observations.py.
#
# Schema (one row per measurement):
#   exp_id | treatment | variable | kind | date | value | sigma | weight
# kind in {timeseries, scalar, phenology}. -99 is treated as missing.

SCHEMA <- c("exp_id", "treatment", "variable", "kind", "date",
            "value", "sigma", "weight")

# Columns whose values are DSSAT date codes (YYDDD) rather than magnitudes.
.DATE_COLS <- c("EDAT", "ADAT", "MDAT", "IDAT", "DRAT", "GDAT", "PD1T", "PDFT",
                "R1", "R3", "R5", "R7", "R8", "TSAT", "HDAT")

.empty_schema <- function() {
  data.frame(exp_id = character(0), treatment = integer(0), variable = character(0),
             kind = character(0), date = as.Date(character(0)), value = numeric(0),
             sigma = numeric(0), weight = numeric(0), stringsAsFactors = FALSE)
}

# Read a DSSAT A/T file into a list of wide data.frames (one per @ header).
# Mirrors observations.py:_read_abt_blocks.
.read_abt_blocks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  blocks <- list(); header <- NULL; rows <- list()
  flush <- function() {
    if (!is.null(header) && length(rows) > 0) {
      ncol <- length(header)
      norm <- lapply(rows, function(r) {
        r <- r[seq_len(min(length(r), ncol))]
        if (length(r) < ncol) r <- c(r, rep(as.character(MISSING), ncol - length(r)))
        r
      })
      m <- do.call(rbind, norm)
      df <- as.data.frame(m, stringsAsFactors = FALSE)
      names(df) <- header
      blocks[[length(blocks) + 1L]] <<- df
    }
    header <<- NULL; rows <<- list()
  }
  for (ln in lines) {
    s <- sub("\\s+$", "", ln)
    if (s == "" || startsWith(s, "!")) next
    if (startsWith(s, "*")) { flush(); next }
    st <- sub("^\\s+", "", s)
    if (startsWith(st, "@")) {
      flush()
      h <- strsplit(sub("^@", "", st), "\\s+")[[1]]
      h <- h[h != ""]
      if (length(h) > 0 && h[1] %in% c("TRNO", "TRT", "TR")) h[1] <- "TRNO"
      # loop-body locals use `<-`; only flush() uses `<<-` to mutate them.
      header <- h; rows <- list()
    } else if (!is.null(header) && grepl("^\\s*[0-9]", s)) {
      toks <- strsplit(sub("^\\s+", "", s), "\\s+")[[1]]
      rows[[length(rows) + 1L]] <- toks[toks != ""]
    }
  }
  flush()
  blocks
}

#' Read a DSSAT FileA (end-of-season averages) into the long schema.
#' Mirrors observations.py:read_filea.
#' @export
read_filea <- function(path, exp_id = NULL) {
  if (is.null(exp_id)) exp_id <- tools::file_path_sans_ext(basename(path))
  out <- list()
  for (wide in .read_abt_blocks(path)) {
    if (!("TRNO" %in% names(wide))) next
    num <- as.data.frame(lapply(wide, function(c) suppressWarnings(as.numeric(c))),
                         stringsAsFactors = FALSE)
    names(num) <- names(wide)
    for (var in setdiff(names(num), "TRNO")) {
      for (i in seq_len(nrow(num))) {
        val <- num[[var]][i]
        if (is.na(val) || abs(val - MISSING) <= 1e-3) next
        if (var %in% .DATE_COLS) {
          out[[length(out) + 1L]] <- data.frame(
            exp_id = exp_id, treatment = as.integer(num$TRNO[i]), variable = var,
            kind = "phenology", date = yyddd_to_date(val), value = as.numeric(val),
            sigma = NA_real_, weight = 1.0, stringsAsFactors = FALSE)
        } else {
          out[[length(out) + 1L]] <- data.frame(
            exp_id = exp_id, treatment = as.integer(num$TRNO[i]), variable = var,
            kind = "scalar", date = as.Date(NA), value = as.numeric(val),
            sigma = NA_real_, weight = 1.0, stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (length(out) == 0) return(.empty_schema())
  do.call(rbind, out)
}

#' Read a DSSAT FileT (in-season time-series; replicate rows preserved).
#' Mirrors observations.py:read_filet.
#' @export
read_filet <- function(path, exp_id = NULL) {
  if (is.null(exp_id)) exp_id <- tools::file_path_sans_ext(basename(path))
  out <- list()
  for (wide in .read_abt_blocks(path)) {
    if (!all(c("TRNO", "DATE") %in% names(wide))) next
    num <- as.data.frame(lapply(wide, function(c) suppressWarnings(as.numeric(c))),
                         stringsAsFactors = FALSE)
    names(num) <- names(wide)
    value_cols <- setdiff(names(wide), c("TRNO", "DATE"))
    for (i in seq_len(nrow(num))) {
      d <- yyddd_to_date(num$DATE[i])
      if (is.na(d)) next
      for (var in value_cols) {
        val <- num[[var]][i]
        if (is.na(val) || abs(val - MISSING) <= 1e-3) next
        out[[length(out) + 1L]] <- data.frame(
          exp_id = exp_id, treatment = as.integer(num$TRNO[i]), variable = var,
          kind = "timeseries", date = d, value = as.numeric(val),
          sigma = NA_real_, weight = 1.0, stringsAsFactors = FALSE)
      }
    }
  }
  if (length(out) == 0) return(.empty_schema())
  do.call(rbind, out)
}

#' Read a user long-format observations CSV and normalise to the schema.
#' Mirrors observations.py:read_csv.
#' @export
read_csv <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- trimws(names(df))
  ren <- c(experiment = "exp_id", exp = "exp_id", trt = "treatment",
           var = "variable", obs = "value", val = "value")
  for (k in names(ren)) if (k %in% names(df)) names(df)[names(df) == k] <- ren[[k]]
  if ("date" %in% names(df)) {
    df$date <- suppressWarnings(as.Date(as.character(df$date)))
  }
  if (!("kind" %in% names(df))) {
    if ("date" %in% names(df)) {
      df$kind <- ifelse(!is.na(df$date), "timeseries", "scalar")
    } else {
      df$kind <- "scalar"
    }
  }
  if (!("sigma" %in% names(df))) df$sigma <- NA_real_
  if (!("weight" %in% names(df))) df$weight <- 1.0
  if (!("date" %in% names(df))) df$date <- as.Date(NA)
  for (col in c("value", "sigma")) {
    if (col %in% names(df)) df[[col]][df[[col]] == MISSING] <- NA
  }
  keep <- intersect(SCHEMA, names(df))
  df[keep]
}

# ----- Observations bundle (S3) --------------------------------------------

#' Construct an Observations bundle from a long-format table.
#' @export
observations <- function(table) {
  structure(list(table = table), class = "observations")
}

#' Load FileA/FileT observations for a list of DSSAT experiment codes.
#' Mirrors Observations.from_dssat.
#' @export
observations_from_dssat <- function(hemp_dir, experiments, crop_ext = "HM") {
  frames <- list()
  for (exp in experiments) {
    fa <- file.path(hemp_dir, sprintf("%s.%sA", exp, crop_ext))
    ft <- file.path(hemp_dir, sprintf("%s.%sT", exp, crop_ext))
    if (file.exists(fa)) frames[[length(frames) + 1L]] <- read_filea(fa, exp)
    if (file.exists(ft)) frames[[length(frames) + 1L]] <- read_filet(ft, exp)
  }
  tbl <- if (length(frames)) do.call(rbind, frames) else .empty_schema()
  observations(tbl)
}

#' @export
observations_from_csv <- function(path) observations(read_csv(path))

#' Return {exp_id: date} from any ingested planting-date rows.
#' Mirrors Observations.planting_dates.
#' @export
planting_dates <- function(obs) {
  tbl <- obs$table
  if (nrow(tbl) == 0 || !("variable" %in% names(tbl))) return(list())
  keys <- c("planting_date", "pdate", "planting", "sowing_date", "sowing")
  m <- tbl[tolower(as.character(tbl$variable)) %in% keys, , drop = FALSE]
  out <- list()
  for (i in seq_len(nrow(m))) {
    if (!is.na(m$date[i])) out[[as.character(m$exp_id[i])]] <- m$date[i]
  }
  out
}

#' Per-experiment x variable count of observations.
#' Mirrors Observations.coverage.
#' @export
coverage <- function(obs) {
  tbl <- obs$table
  if (nrow(tbl) == 0) return(data.frame())
  agg <- aggregate(list(n = seq_len(nrow(tbl))),
                   by = list(exp_id = tbl$exp_id, kind = tbl$kind, variable = tbl$variable),
                   FUN = length)
  agg
}

#' @export
obs_experiments <- function(obs) {
  if (nrow(obs$table) == 0) return(character(0))
  sort(unique(as.character(obs$table$exp_id)))
}

#' @export
obs_variables <- function(obs) {
  if (nrow(obs$table) == 0) return(character(0))
  sort(unique(as.character(obs$table$variable)))
}

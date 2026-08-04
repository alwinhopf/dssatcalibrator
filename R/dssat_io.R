# DSSAT I/O helpers. R twin of python/dssatcalibrator/dssat_io.py.
#
# Phase 2 scope: the shared constants and date conversion used by the
# observation readers. The folder builders / output parsers (write FileX,
# parse Summary/Evaluate/PlantGro) are ported in the Phase 3 spawn layer.

# Sentinel value DSSAT uses for "missing"; treated as NA throughout.
MISSING <- -99.0

#' Convert a DSSAT YYDDD (or YYYYDDD) date code to a Date.
#'
#' Two-digit years < 80 are read as 20xx, otherwise 19xx (DSSAT convention).
#' Returns NA for missing / unparseable codes. Mirrors dssat_io.py:yyddd_to_date.
#' @export
yyddd_to_date <- function(code) {
  v <- suppressWarnings(as.numeric(code))
  if (length(v) != 1L || is.na(v)) return(as.Date(NA))
  n <- as.integer(trunc(v))            # int(float(code)) — truncates toward zero
  s <- as.character(n)
  if (s %in% c("", "-99", "0")) return(as.Date(NA))
  if (nchar(s) <= 5L) {                # YYDDD
    s5 <- sprintf("%05d", n)
    yy <- as.integer(substr(s5, 1, 2))
    doy <- as.integer(substr(s5, 3, 5))
    year <- if (yy < 80) 2000L + yy else 1900L + yy
  } else {                             # YYYYDDD
    s7 <- sprintf("%07d", n)
    year <- as.integer(substr(s7, 1, 4))
    doy <- as.integer(substr(s7, 5, 7))
  }
  # Match Python and DSSAT's offset interpretation: YY366 is accepted for any
  # year and naturally rolls to 1 January of the next year when the encoded
  # year is not a leap year. Values outside the DSSAT 1..366 range are invalid.
  if (is.na(doy) || doy < 1L || doy > 366L) return(as.Date(NA))
  as.Date(sprintf("%04d-01-01", year)) + (doy - 1L)
}

# Coerce all columns to numeric where possible and map -99 -> NA.
# Mirrors dssat_io.py:_to_numeric (np.isclose tol ~1e-3 around -99).
.to_numeric <- function(df) {
  out <- as.data.frame(lapply(df, function(c) suppressWarnings(as.numeric(c))),
                       stringsAsFactors = FALSE)
  names(out) <- names(df)
  for (j in seq_along(out)) {
    col <- out[[j]]
    col[!is.na(col) & abs(col - MISSING) <= 1e-3] <- NA
    out[[j]] <- col
  }
  out
}

# union-rbind: stack frames that may differ in columns (fill missing with NA),
# mirroring pandas.concat's column-union behaviour.
.rbind_union <- function(frames) {
  cols <- unique(unlist(lapply(frames, names)))
  filled <- lapply(frames, function(f) {
    for (c in setdiff(cols, names(f))) f[[c]] <- NA
    f[cols]
  })
  do.call(rbind, filled)
}

#' Parse `PlantGro.OUT` into a tidy daily data.frame (one row per treatment-day).
#' Mirrors dssat_io.py:parse_plantgro.
#' @export
parse_plantgro <- function(path) {
  text <- readLines(path, warn = FALSE)
  runs <- list()
  run_no <- NULL; treatment <- NULL; header <- NULL; rows <- list()

  flush <- function() {
    if (!is.null(header) && length(rows) > 0) {
      m <- do.call(rbind, rows)
      df <- as.data.frame(m, stringsAsFactors = FALSE); names(df) <- header
      df <- .to_numeric(df)
      df$run <- if (is.null(run_no)) NA_integer_ else run_no
      df$treatment <- if (!is.null(treatment)) treatment else run_no
      runs[[length(runs) + 1L]] <<- df
    }
    rows <<- list(); header <<- NULL
  }

  # NOTE: loop-body assignments use `<-` (these are locals of parse_plantgro);
  # only the nested flush() uses `<<-` to mutate them. Using `<<-` here would
  # write to the global scope and leave flush() reading empty state.
  for (ln in text) {
    if (startsWith(ln, "*RUN")) {
      flush()
      m <- regmatches(ln, regexec("\\*RUN\\s+(\\d+)", ln))[[1]]
      run_no <- if (length(m) >= 2) as.integer(m[2]) else NA_integer_
      treatment <- NULL
    } else if (startsWith(trimws(ln, which = "left"), "TREATMENT")) {
      m <- regmatches(ln, regexec("\\s*TREATMENT\\s+(\\d+)", ln))[[1]]
      if (length(m) >= 2) treatment <- as.integer(m[2])
    } else if (startsWith(ln, "@YEAR") || startsWith(ln, "@ YEAR")) {
      header <- strsplit(sub("^@+", "", ln), "\\s+")[[1]]
      header <- header[header != ""]
      rows <- list()
    } else if (!is.null(header) && grepl("^\\s*\\d{4}\\s", ln)) {
      parts <- strsplit(trimws(ln), "\\s+")[[1]]
      if (length(parts) >= length(header)) {
        rows[[length(rows) + 1L]] <- parts[seq_len(length(header))]
      } else if (length(parts) > 0) {
        rows[[length(rows) + 1L]] <- c(parts, rep(as.character(MISSING), length(header) - length(parts)))
      }
    }
  }
  flush()

  if (length(runs) == 0) return(data.frame())
  out <- .rbind_union(runs)
  if (all(c("YEAR", "DOY") %in% names(out))) {
    out$date <- as.Date(sprintf("%d-%03d", as.integer(out$YEAR), as.integer(out$DOY)),
                        format = "%Y-%j")
  }
  out$run <- as.integer(out$run)
  out$treatment <- as.integer(out$treatment)
  out
}

#' Parse `Evaluate.OUT` into a long table: one row per (treatment, variable).
#' Columns: treatment, run, variable, sim, meas. Mirrors dssat_io.py:parse_evaluate.
#' @export
parse_evaluate <- function(path) {
  empty <- data.frame(treatment = integer(0), run = integer(0),
                      variable = character(0), sim = numeric(0), meas = numeric(0),
                      stringsAsFactors = FALSE)
  if (!file.exists(path)) return(empty)
  lines <- readLines(path, warn = FALSE)
  hdr_idx <- which(startsWith(lines, "@RUN"))
  if (length(hdr_idx) == 0) return(empty)
  hdr_idx <- hdr_idx[1]
  header <- strsplit(sub("^@+", "", lines[hdr_idx]), "\\s+")[[1]]
  header <- header[header != ""]
  data_lines <- lines[(hdr_idx + 1L):length(lines)]
  data_lines <- data_lines[grepl("^\\s*\\d", data_lines)]
  if (length(data_lines) == 0) return(empty)
  splitrows <- lapply(data_lines, function(l) strsplit(trimws(l), "\\s+")[[1]])
  ncol <- length(splitrows[[1]])
  cols <- header[seq_len(ncol)]
  wide <- as.data.frame(do.call(rbind, lapply(splitrows, function(r) r[seq_len(ncol)])),
                        stringsAsFactors = FALSE)
  names(wide) <- cols

  id_cols <- c("RUN", "EXCODE", "TN", "RN", "CR")
  sim_cols <- cols[endsWith(cols, "S") &
                     (paste0(substr(cols, 1, nchar(cols) - 1L), "M") %in% cols) &
                     !(cols %in% id_cols)]
  recs <- list()
  for (i in seq_len(nrow(wide))) {
    row <- wide[i, ]
    trt <- suppressWarnings(as.numeric(if (!is.null(row[["TN"]])) row[["TN"]] else NA))
    run <- suppressWarnings(as.numeric(if (!is.null(row[["RUN"]])) row[["RUN"]] else NA))
    for (sc in sim_cols) {
      base <- substr(sc, 1, nchar(sc) - 1L)
      sim <- suppressWarnings(as.numeric(row[[sc]]))
      meas <- suppressWarnings(as.numeric(row[[paste0(base, "M")]]))
      recs[[length(recs) + 1L]] <- data.frame(treatment = trt, run = run, variable = base,
                                               sim = sim, meas = meas, stringsAsFactors = FALSE)
    }
  }
  if (length(recs) == 0) return(empty)
  out <- do.call(rbind, recs)
  out$sim[!is.na(out$sim) & abs(out$sim - MISSING) <= 1e-3] <- NA
  out$meas[!is.na(out$meas) & abs(out$meas - MISSING) <= 1e-3] <- NA
  out$treatment <- as.integer(out$treatment)
  out$run <- as.integer(out$run)
  out
}

# Summary.OUT numeric columns of interest (TNAM/FNAM have spaces, so read by
# left header index for the numeric tail). Mirrors dssat_io.py:_SUMMARY_NUM.
.SUMMARY_NUM <- c("RUNNO", "TRNO", "CWAM", "HWAM", "HWAH", "BWAH", "PWAM", "LAIX",
                  "ADAT", "MDAT", "EDAT", "PDAT", "HDAT")

#' Best-effort parse of `Summary.OUT` numeric columns by header name.
#' Mirrors dssat_io.py:parse_summary.
#' @export
parse_summary <- function(path) {
  if (!file.exists(path)) return(data.frame())
  lines <- readLines(path, warn = FALSE)
  hdr_idx <- which(startsWith(trimws(lines, which = "left"), "@") & grepl("RUNNO", lines))
  if (length(hdr_idx) == 0) return(data.frame())
  hdr_idx <- hdr_idx[1]
  header <- strsplit(sub("^@+", "", lines[hdr_idx]), "\\s+")[[1]]
  header <- header[header != ""]
  recs <- list()
  for (ln in lines[(hdr_idx + 1L):length(lines)]) {
    if (!grepl("^\\s*\\d", ln)) next
    parts <- strsplit(trimws(ln), "\\s+")[[1]]
    row <- list()
    for (col in .SUMMARY_NUM) {
      if (col %in% header) {
        idx <- match(col, header)
        if (idx <= length(parts)) row[[col]] <- parts[idx]
      }
    }
    recs[[length(recs) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  if (length(recs) == 0) return(data.frame())
  df <- .rbind_union(recs)
  .to_numeric(df)
}

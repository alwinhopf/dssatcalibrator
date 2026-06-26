# Write perturbed parameter values into DSSAT input files.
# R twin of python/dssatcalibrator/writers.py.
#
# All writers are column-aware, driven by the `@` header so they tolerate
# DSSAT's variable-width spacing. Column bounds are kept Python-style
# (0-indexed, half-open [lo, hi)); a small set of helpers translate to R's
# 1-indexed inclusive substring/char operations so the byte output matches the
# Python implementation exactly.

# Right-justify a number into `width` chars at the highest precision that fits.
# Mirrors writers.py:_fmt.
.fmt <- function(value, width) {
  for (dec in 4:0) {
    s <- formatC(value, format = "f", digits = dec)
    if (nchar(s) <= width) return(sprintf("%*s", width, s))
  }
  s <- formatC(value, format = "f", digits = 0)
  substr(sprintf("%*s", width, s), 1, width)
}

# Parse a fixed-width cell to numeric, or NA if not a number. Mirrors _parse.
.parse_cell <- function(cell) {
  v <- suppressWarnings(as.numeric(trimws(cell)))
  if (length(v) != 1L) return(NA_real_)
  v
}

# Python line[lo:hi] (0-indexed, half-open), tolerant of short lines.
.py_slice <- function(line, lo, hi) {
  n <- nchar(line)
  if (lo >= n) return("")
  substr(line, lo + 1L, min(hi, n))
}

# Replace chars at Python 0-indexed [lo, hi) with `text` (length hi-lo).
.splice <- function(chars, lo, hi, text) {
  tch <- strsplit(text, "", fixed = TRUE)[[1]]
  n <- hi - lo
  if (length(tch) < n) tch <- c(tch, rep(" ", n - length(tch)))
  tch <- tch[seq_len(n)]
  chars[(lo + 1L):hi] <- tch
  chars
}

.chars <- function(line) strsplit(line, "", fixed = TRUE)[[1]]
.ljust <- function(line, width) if (nchar(line) < width) paste0(line, strrep(" ", width - nchar(line))) else line
.write_lines <- function(lines, path) writeLines(lines, con = path, sep = "\n")

# Token start/end (Python 0-indexed: start0, end0 exclusive) for \S+ matches.
.token_spans <- function(s) {
  m <- gregexpr("\\S+", s)[[1]]
  if (length(m) == 1L && m[1] == -1L) return(list(tok = character(0), start0 = integer(0), end0 = integer(0)))
  lens <- attr(m, "match.length")
  starts0 <- as.integer(m) - 1L
  list(tok = vapply(seq_along(m), function(i) substr(s, m[i], m[i] + lens[i] - 1L), character(1)),
       start0 = starts0, end0 = starts0 + lens)
}

#' {coefficient: c(start_col, end_col)} from the `@VAR#` header (Python-style cols).
#' Mirrors writers.py:cultivar_field_map.
#' @export
cultivar_field_map <- function(cul_path) {
  lines <- readLines(cul_path, warn = FALSE)
  header <- lines[startsWith(lines, "@VAR#")]
  if (length(header) == 0) stop(sprintf("No '@VAR#' header found in %s", cul_path))
  header <- header[1]
  sp <- .token_spans(header)
  eco_i <- match("ECO#", sp$tok)
  coeff_idx <- (eco_i + 1L):length(sp$tok)
  boundaries <- c(sp$end0[eco_i], sp$end0[coeff_idx])
  fmap <- list()
  for (k in seq_along(coeff_idx)) {
    fmap[[sp$tok[coeff_idx[k]]]] <- c(boundaries[k], boundaries[k + 1L])
  }
  fmap
}

#' {coefficient: c(start_col, end_col)} from the `@ECO#` header.
#' Mirrors writers.py:ecotype_field_map.
#' @export
ecotype_field_map <- function(eco_path) {
  lines <- readLines(eco_path, warn = FALSE)
  header <- lines[startsWith(lines, "@ECO#")]
  if (length(header) == 0) stop(sprintf("No '@ECO#' header found in %s", eco_path))
  header <- header[1]
  sp <- .token_spans(header)
  name_i <- which(startsWith(sp$tok, "ECONAME"))[1]
  coeff_idx <- (name_i + 1L):length(sp$tok)
  boundaries <- c(sp$end0[name_i], sp$end0[coeff_idx])
  fmap <- list()
  for (k in seq_along(coeff_idx)) {
    fmap[[sp$tok[coeff_idx[k]]]] <- c(boundaries[k], boundaries[k + 1L])
  }
  fmap
}

#' {token: c(start_col, end_col)} for columns under a DSSAT header line.
#' Mirrors writers.py:parse_header_boundaries.
#' @export
parse_header_boundaries <- function(header) {
  sp <- .token_spans(header)
  fmap <- list()
  for (i in seq_along(sp$tok)) {
    name <- sp$tok[i]
    if (startsWith(name, "@")) name <- substring(name, 2)
    start <- if (i > 1L) sp$end0[i - 1L] else 0L
    fmap[[name]] <- c(start, sp$end0[i])
  }
  fmap
}

# shared row-editor for tabular genotype files (.CUL/.ECO)
.edit_tabular_row <- function(path, anchor_code, updates, fmap, what) {
  lines <- readLines(path, warn = FALSE)
  for (name in names(updates)) {
    if (!(name %in% names(fmap))) {
      stop(sprintf("Coefficient '%s' not a column in %s; known: %s",
                   name, basename(path), paste(sort(names(fmap)), collapse = ", ")))
    }
  }
  target <- which(startsWith(lines, anchor_code) &
                    !startsWith(trimws(lines, which = "left"), "!"))
  if (length(target) == 0) {
    stop(sprintf("Active %s row '%s' not found in %s", what, anchor_code, basename(path)))
  }
  idx <- target[1]
  last_end <- max(vapply(fmap, function(b) b[2], numeric(1)))
  line <- .ljust(lines[idx], last_end)
  chars <- .chars(line)
  for (name in names(fmap)) {
    if (name %in% names(updates)) {
      b <- fmap[[name]]
      chars <- .splice(chars, b[1], b[2], .fmt(as.numeric(updates[[name]]), b[2] - b[1]))
    }
  }
  lines[idx] <- paste(chars, collapse = "")
  .write_lines(lines, path)
  invisible(NULL)
}

#' In-place edit of the `anchor_code` cultivar row in a CROPGRO `.CUL`.
#' Mirrors writers.py:edit_cultivar.
#' @export
edit_cultivar <- function(cul_path, anchor_code, updates) {
  .edit_tabular_row(cul_path, anchor_code, updates, cultivar_field_map(cul_path), "cultivar")
}

#' In-place edit of the `anchor_code` ecotype row in a CROPGRO `.ECO`.
#' Mirrors writers.py:edit_ecotype.
#' @export
edit_ecotype <- function(eco_path, anchor_code, updates) {
  .edit_tabular_row(eco_path, anchor_code, updates, ecotype_field_map(eco_path), "ecotype")
}

#' Read current coefficient values for a cultivar row. Mirrors read_cultivar_values.
#' @export
read_cultivar_values <- function(cul_path, anchor_code) {
  lines <- readLines(cul_path, warn = FALSE)
  fmap <- cultivar_field_map(cul_path)
  line <- lines[startsWith(lines, anchor_code) & !startsWith(trimws(lines, which = "left"), "!")]
  if (length(line) == 0) stop(sprintf("Cultivar '%s' not found in %s", anchor_code, basename(cul_path)))
  line <- line[1]
  out <- list()
  for (name in names(fmap)) {
    b <- fmap[[name]]
    out[[name]] <- if (nchar(line) >= b[2]) .parse_cell(.py_slice(line, b[1], b[2])) else NA_real_
  }
  out
}

#' Read current coefficient values for an ecotype row. Mirrors read_ecotype_values.
#' @export
read_ecotype_values <- function(eco_path, anchor_code) {
  lines <- readLines(eco_path, warn = FALSE)
  fmap <- ecotype_field_map(eco_path)
  line <- lines[startsWith(lines, anchor_code) & !startsWith(trimws(lines, which = "left"), "!")]
  if (length(line) == 0) stop(sprintf("Ecotype '%s' not found in %s", anchor_code, basename(eco_path)))
  line <- line[1]
  out <- list()
  for (name in names(fmap)) {
    b <- fmap[[name]]
    val <- .parse_cell(.py_slice(line, b[1], b[2]))
    if (!is.na(val)) out[[name]] <- val
  }
  out
}

#' Edit management (planting details) and initial conditions in a FileX in-place.
#' Mirrors writers.py:edit_filex.
#' @export
edit_filex <- function(filex_path, mgt_fields = list(), init_updates = list()) {
  lines <- readLines(filex_path, warn = FALSE)
  is_data <- function(ln) { t <- trimws(ln); nzchar(t) && grepl("^[0-9]", t) }

  if (length(mgt_fields) > 0) {
    psec <- which(startsWith(lines, "*PLANTING DETAILS"))
    if (length(psec) > 0) {
      psec <- psec[1]; header_idx <- NA_integer_
      for (i in (psec + 1L):length(lines)) {
        if (startsWith(lines[i], "*")) break
        if (startsWith(trimws(lines[i], which = "left"), "@P")) { header_idx <- i; break }
      }
      if (!is.na(header_idx)) {
        fmap <- parse_header_boundaries(lines[header_idx])
        last_end <- max(vapply(fmap, function(b) b[2], numeric(1)))
        for (i in (header_idx + 1L):length(lines)) {
          ln <- lines[i]
          if (startsWith(ln, "*") || startsWith(trimws(ln, which = "left"), "@")) break
          if (startsWith(trimws(ln, which = "left"), "!") || !nzchar(trimws(ln))) next
          if (is_data(ln)) {
            chars <- .chars(.ljust(ln, last_end))
            for (name in names(mgt_fields)) {
              if (name %in% names(fmap)) {
                b <- fmap[[name]]
                chars <- .splice(chars, b[1], b[2], .fmt(as.numeric(mgt_fields[[name]]), b[2] - b[1]))
              }
            }
            lines[i] <- paste(chars, collapse = "")
          }
        }
      }
    }
  }

  if (length(init_updates) > 0 && !is.null(init_updates[["initial_soil_water_mult"]])) {
    mult <- as.numeric(init_updates[["initial_soil_water_mult"]])
    isec <- which(startsWith(lines, "*INITIAL CONDITIONS"))
    if (length(isec) > 0) {
      isec <- isec[1]; header_idx <- NA_integer_
      for (i in (isec + 1L):length(lines)) {
        if (startsWith(lines[i], "*")) break
        if (startsWith(trimws(lines[i], which = "left"), "@C") && grepl("SH2O", lines[i])) {
          header_idx <- i; break
        }
      }
      if (!is.na(header_idx)) {
        fmap <- parse_header_boundaries(lines[header_idx])
        if ("SH2O" %in% names(fmap)) {
          b <- fmap[["SH2O"]]; lo <- b[1]; hi <- b[2]
          for (i in (header_idx + 1L):length(lines)) {
            ln <- lines[i]
            if (startsWith(ln, "*") || startsWith(trimws(ln, which = "left"), "@")) break
            if (startsWith(trimws(ln, which = "left"), "!") || !nzchar(trimws(ln))) next
            if (is_data(ln)) {
              chars <- .chars(if (nchar(ln) < hi) .ljust(ln, hi) else ln)
              val_str <- trimws(paste(chars[(lo + 1L):hi], collapse = ""))
              val <- suppressWarnings(as.numeric(val_str))
              if (!is.na(val)) {
                new_val <- min(max(val * mult, 0.01), 1.0)
                chars <- .splice(chars, lo, hi, .fmt(new_val, hi - lo))
                lines[i] <- paste(chars, collapse = "")
              }
            }
          }
        }
      }
    }
  }
  .write_lines(lines, filex_path)
  invisible(NULL)
}

#' Return field metadata from the FileX `*FIELDS` section. Mirrors parse_fields.
#' @export
parse_fields <- function(filex_path) {
  norm_token <- function(t) gsub("^@", "", gsub("\\.+$", "", t))
  as_float <- function(v) {
    val <- suppressWarnings(as.numeric(v))
    if (length(val) != 1L || is.na(val)) return(NULL)
    if (val %in% c(-99.0, -999.0)) return(NULL)
    val
  }
  lines <- readLines(filex_path, warn = FALSE)
  sec <- which(startsWith(lines, "*FIELDS"))
  if (length(sec) == 0) return(list())
  fields <- list(); i <- sec[1] + 1L
  while (i <= length(lines) && !startsWith(lines[i], "*")) {
    if (!startsWith(trimws(lines[i], which = "left"), "@L")) { i <- i + 1L; next }
    header <- lines[i]; data <- NULL
    j <- i + 1L
    while (j <= length(lines)) {
      lj <- lines[j]
      if (!startsWith(trimws(lj, which = "left"), "@") && !startsWith(lj, "*") &&
          nzchar(trimws(lj)) && grepl("^[0-9]", trimws(lj))) { data <- lj; break }
      if (startsWith(lj, "*")) break
      j <- j + 1L
    }
    if (is.null(data)) { i <- i + 1L; next }
    htoks <- vapply(strsplit(trimws(header), "\\s+")[[1]], norm_token, character(1))
    dtoks <- strsplit(trimws(data), "\\s+")[[1]]
    if (length(dtoks) >= length(htoks)) {
      for (k in seq_along(htoks)) fields[[htoks[k]]] <- dtoks[k]
    }
    i <- i + 1L
  }
  if (length(fields) == 0) return(list())
  list(wsta = fields[["WSTA"]], id_soil = fields[["ID_SOIL"]], id_field = fields[["ID_FIELD"]],
       lat = as_float(fields[["YCRD"]]), lon = as_float(fields[["XCRD"]]),
       elev = as_float(fields[["ELEV"]]))
}

#' Return the single `*<profile_id>` block from a (multi-profile) `.SOL`.
#' Mirrors writers.py:extract_soil_profile.
#' @export
extract_soil_profile <- function(sol_path, profile_id) {
  lines <- readLines(sol_path, warn = FALSE)
  out <- character(0); capturing <- FALSE
  for (ln in lines) {
    if (startsWith(ln, "*")) {
      rest <- substring(ln, 2)
      toks <- strsplit(trimws(rest), "\\s+")[[1]]
      token <- if (nchar(ln) > 1 && length(toks) > 0 && nzchar(toks[1])) toks[1] else ""
      if (capturing) break
      if (token == profile_id) capturing <- TRUE
    }
    if (capturing) out <- c(out, ln)
  }
  while (length(out) > 0 && (!nzchar(trimws(out[length(out)])) ||
                             startsWith(trimws(out[length(out)], which = "left"), "!"))) {
    out <- out[-length(out)]
  }
  if (length(out) == 0) stop(sprintf("Soil profile '%s' not found in %s", profile_id, basename(sol_path)))
  paste0(paste(out, collapse = "\n"), "\n")
}

#' Edit one soil profile in a `.SOL` in place. Mirrors writers.py:edit_soil.
#' @export
edit_soil <- function(sol_path, profile_id, layer_mults = list(), profile_sets = list()) {
  lines <- readLines(sol_path, warn = FALSE)
  start <- NA_integer_
  for (i in seq_along(lines)) {
    ln <- lines[i]
    if (startsWith(ln, "*") && nchar(ln) > 1) {
      toks <- strsplit(trimws(substring(ln, 2)), "\\s+")[[1]]
      if (length(toks) > 0 && toks[1] == profile_id) { start <- i; break }
    }
  }
  if (is.na(start)) stop(sprintf("Soil profile '%s' not found in %s", profile_id, basename(sol_path)))
  end <- length(lines) + 1L
  if (start < length(lines)) {
    nxt <- which(startsWith(lines[(start + 1L):length(lines)], "*"))
    if (length(nxt) > 0) end <- start + nxt[1]
  }

  if (length(profile_sets) > 0) {
    for (j in start:(end - 1L)) {
      ln <- lines[j]
      if (startsWith(trimws(ln, which = "left"), "@") && grepl("SALB", ln) && !grepl("SDUL", ln)) {
        fmap <- parse_header_boundaries(ln)
        k <- j + 1L; chars <- .chars(lines[k])
        for (name in names(profile_sets)) {
          if (name %in% names(fmap)) {
            b <- fmap[[name]]
            if (length(chars) < b[2]) chars <- c(chars, rep(" ", b[2] - length(chars)))
            chars <- .splice(chars, b[1], b[2], .fmt(as.numeric(profile_sets[[name]]), b[2] - b[1]))
          }
        }
        lines[k] <- paste(chars, collapse = "")
        break
      }
    }
  }

  if (length(layer_mults) > 0) {
    for (j in start:(end - 1L)) {
      ln <- lines[j]
      if (startsWith(trimws(ln, which = "left"), "@") && grepl("SDUL", ln) && grepl("SLB", ln)) {
        fmap <- parse_header_boundaries(ln)
        for (k in (j + 1L):(end - 1L)) {
          row <- lines[k]
          lstrip <- trimws(row, which = "left")
          if (startsWith(lstrip, "@") || startsWith(lstrip, "*")) break
          if (!nzchar(trimws(row)) || startsWith(lstrip, "!")) next
          if (!grepl("^[0-9]", lstrip)) next
          chars <- .chars(row)
          for (name in names(layer_mults)) {
            if (!(name %in% names(fmap))) next
            b <- fmap[[name]]
            cur <- if (length(chars) >= b[2]) .parse_cell(paste(chars[(b[1] + 1L):b[2]], collapse = "")) else NA_real_
            if (is.na(cur) || cur == -99) next
            new <- min(max(cur * as.numeric(layer_mults[[name]]), 0.0), 1.0)
            if (length(chars) < b[2]) chars <- c(chars, rep(" ", b[2] - length(chars)))
            chars <- .splice(chars, b[1], b[2], .fmt(new, b[2] - b[1]))
          }
          lines[k] <- paste(chars, collapse = "")
        }
        break
      }
    }
  }
  .write_lines(lines, sol_path)
  invisible(NULL)
}

#' Edit a `.WTH` in place. `ops` maps a column to list(mode, value) where mode is
#' "mult" (factor) or "off" (additive). Missing (-99) cells untouched.
#' Mirrors writers.py:edit_weather.
#' @export
edit_weather <- function(wth_path, ops) {
  lines <- readLines(wth_path, warn = FALSE)
  hdr_idx <- NA_integer_
  for (i in seq_along(lines)) {
    if (startsWith(trimws(lines[i], which = "left"), "@") && grepl("SRAD", lines[i])) { hdr_idx <- i; break }
  }
  if (is.na(hdr_idx)) return(invisible(NULL))
  fmap <- parse_header_boundaries(lines[hdr_idx])
  fmt_w <- function(value, width) {
    s <- formatC(value, format = "f", digits = 1)
    if (nchar(s) <= width) sprintf("%*s", width, s) else .fmt(value, width)
  }
  for (i in (hdr_idx + 1L):length(lines)) {
    ln <- lines[i]
    lstrip <- trimws(ln, which = "left")
    if (!nzchar(trimws(ln)) || grepl("^[@*!$]", lstrip)) next
    if (!grepl("^[0-9]", lstrip)) next
    chars <- .chars(ln)
    for (col in names(ops)) {
      if (!(col %in% names(fmap))) next
      b <- fmap[[col]]
      cur <- if (length(chars) >= b[2]) .parse_cell(paste(chars[(b[1] + 1L):b[2]], collapse = "")) else NA_real_
      if (is.na(cur) || cur == -99) next
      mode <- ops[[col]][[1]]; val <- as.numeric(ops[[col]][[2]])
      new <- if (mode == "mult") cur * val else cur + val
      if (length(chars) < b[2]) chars <- c(chars, rep(" ", b[2] - length(chars)))
      chars <- .splice(chars, b[1], b[2], fmt_w(new, b[2] - b[1]))
    }
    lines[i] <- paste(chars, collapse = "")
  }
  .write_lines(lines, wth_path)
  invisible(NULL)
}

.NUM_RE <- "-?[0-9]+\\.?[0-9]*([eE][-+]?[0-9]+)?"

#' Best-effort in-place edit of named scalar values in a `.SPE`.
#' Mirrors writers.py:edit_species.
#' @export
edit_species <- function(spe_path, updates) {
  lines <- readLines(spe_path, warn = FALSE)
  for (key in names(updates)) {
    val <- updates[[key]]
    matches <- which(vapply(lines, function(ln) {
      grepl(key, ln, fixed = TRUE) && nzchar(trimws(ln)) &&
        !grepl("^[*@!]", trimws(ln, which = "left")) &&
        grepl(.NUM_RE, ln)
    }, logical(1)))
    if (length(matches) != 1L) {
      stop(sprintf("Species key '%s' matched %d lines in %s (need exactly 1); refine the key.",
                   key, length(matches), basename(spe_path)))
    }
    i <- matches[1]; ln <- lines[i]
    m <- regexpr(.NUM_RE, ln)
    s0 <- as.integer(m); len <- attr(m, "match.length")
    width <- len
    new <- if (width > 0) .fmt(as.numeric(val), width) else formatC(as.numeric(val), format = "f", digits = 3)
    new <- sprintf("%*s", width, trimws(new))
    lines[i] <- paste0(substr(ln, 1, s0 - 1L), new, substring(ln, s0 + len))
  }
  .write_lines(lines, spe_path)
  invisible(NULL)
}

#' Read DSSAT `.CUL` MINIMA/MAXIMA (VAR# 999991/999992) calibration rows into
#' per-coefficient bounds. Mirrors writers.py:read_cul_calibration_bounds.
#' @export
read_cul_calibration_bounds <- function(cul_path) {
  lines <- readLines(cul_path, warn = FALSE)
  fmap <- cultivar_field_map(cul_path)
  rows <- list()
  for (ln in lines) {
    code <- trimws(substr(ln, 1, 6))
    if (code %in% c("999991", "999992")) {
      vals <- list()
      for (name in names(fmap)) {
        b <- fmap[[name]]
        v <- if (nchar(ln) >= b[2]) .parse_cell(.py_slice(ln, b[1], b[2])) else NA_real_
        if (!is.na(v)) vals[[name]] <- v
      }
      rows[[if (code == "999991") "min" else "max"]] <- vals
    }
  }
  if (is.null(rows[["min"]]) || is.null(rows[["max"]])) return(list())
  out <- list()
  for (name in names(fmap)) {
    if (!is.null(rows[["min"]][[name]]) && !is.null(rows[["max"]][[name]])) {
      out[[name]] <- list(min = rows[["min"]][[name]], max = rows[["max"]][[name]])
    }
  }
  out
}

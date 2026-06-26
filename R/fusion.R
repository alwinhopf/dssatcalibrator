# Merge observations from multiple sources with conflict resolution.
# R twin of python/dssatcalibrator/fusion.py.

#' Extended observation schema (adds provenance columns to SCHEMA).
#' @export
SCHEMA_EXTENDED <- c("exp_id", "treatment", "variable", "kind", "date", "value",
                     "sigma", "weight", "source", "quality_flag", "spatial_res_m")

# one extended-schema row as a data.frame (used by the source adapters)
.ext_row <- function(exp_id, treatment, variable, kind, date, value, sigma,
                     weight, source, quality_flag, spatial_res_m) {
  data.frame(exp_id = as.character(exp_id), treatment = as.integer(treatment),
             variable = as.character(variable), kind = as.character(kind),
             date = as.Date(date), value = as.numeric(value), sigma = as.numeric(sigma),
             weight = as.numeric(weight), source = as.character(source),
             quality_flag = as.integer(quality_flag),
             spatial_res_m = as.numeric(spatial_res_m), stringsAsFactors = FALSE)
}

.empty_extended <- function() {
  data.frame(exp_id = character(0), treatment = integer(0), variable = character(0),
             kind = character(0), date = as.Date(character(0)), value = numeric(0),
             sigma = numeric(0), weight = numeric(0), source = character(0),
             quality_flag = integer(0), spatial_res_m = numeric(0),
             stringsAsFactors = FALSE)
}

#' Construct an ObservationFuser over a list of sources.
#' Mirrors fusion.py:ObservationFuser.__init__.
#' @export
observation_fuser <- function(sources, cfg) {
  named <- list()
  for (s in sources) named[[s$name]] <- s
  structure(list(sources = named, cfg = cfg), class = "observation_fuser")
}

#' Gather from all active sources, apply QC, merge. Mirrors ObservationFuser.collect.
#' @export
fuser_collect <- function(fuser, experiment, date_range) {
  frames <- list()
  for (name in names(fuser$sources)) {
    src <- fuser$sources[[name]]
    res <- tryCatch({
      df <- src_fetch(src, experiment, date_range)
      if (nrow(df) == 0) NULL else {
        df <- src_quality_filter(src, df)
        df$source <- name
        if (!("sigma" %in% names(df)) || any(is.na(df$sigma))) {
          df$sigma <- mapply(function(v, x) src_error_model(src, v, x, list()),
                             df$variable, df$value)
        }
        for (col in SCHEMA_EXTENDED) if (!(col %in% names(df))) df[[col]] <- NA
        df[SCHEMA_EXTENDED]
      }
    }, error = function(e) {
      warning(sprintf("Source %s failed to fetch for %s: %s", name, experiment,
                      conditionMessage(e)))
      NULL
    })
    if (!is.null(res)) frames[[length(frames) + 1L]] <- res
  }
  if (length(frames) == 0) return(.empty_extended())
  merged <- do.call(rbind, frames)
  fuser_resolve_conflicts(fuser, merged)
}

#' Handle overlapping observations from different sources.
#' Strategies: keep_all | inverse_variance | priority. Mirrors resolve_conflicts.
#' @export
fuser_resolve_conflicts <- function(fuser, df) {
  if (nrow(df) == 0) return(df)
  strategy <- .cfg_get(.cfg_get(fuser$cfg, "fusion", list()),
                       "conflict_resolution", "keep_all")
  if (strategy == "inverse_variance") return(.inverse_variance_merge(df))
  if (strategy == "priority") return(.priority_merge(fuser, df))
  df  # keep_all (and unknown) pass through
}

# Inverse-variance weighted merge of coincident observations.
.inverse_variance_merge <- function(df) {
  df$date_key <- as.Date(df$date)
  key <- paste(df$exp_id, df$treatment, df$variable, as.character(df$date_key), sep = "\r")
  out <- list()
  for (k in unique(key)) {
    g <- df[key == k, , drop = FALSE]
    if (nrow(g) == 1) {
      out[[length(out) + 1L]] <- g[, setdiff(names(g), "date_key"), drop = FALSE]
      next
    }
    gc <- g[!is.na(g$sigma), , drop = FALSE]
    if (nrow(gc) == 0) {
      out[[length(out) + 1L]] <- g[1, setdiff(names(g), "date_key"), drop = FALSE]
      next
    }
    sigmas <- gc$sigma; sigmas[sigmas == 0] <- 1e-6
    w <- 1.0 / (sigmas^2)
    val <- sum(gc$value * w) / sum(w)
    sig <- 1.0 / sqrt(sum(w))
    row <- gc[1, , drop = FALSE]
    row$value <- val; row$sigma <- sig
    row$source <- paste(sort(unique(gc$source)), collapse = "+")
    row$quality_flag <- as.integer(max(gc$quality_flag))
    row$spatial_res_m <- mean(gc$spatial_res_m)
    out[[length(out) + 1L]] <- row[, setdiff(names(row), "date_key"), drop = FALSE]
  }
  do.call(rbind, out)
}

# Keep only the highest-priority source for coincident measurements.
.priority_merge <- function(fuser, df) {
  priority_list <- .cfg_get(.cfg_get(fuser$cfg, "fusion", list()), "source_priority", list())
  if (length(priority_list) == 0) return(df)
  priority_list <- unlist(priority_list, use.names = FALSE)
  p_map <- setNames(seq_along(priority_list) - 1L, priority_list)
  df$date_key <- as.Date(df$date)
  df$priority_rank <- vapply(df$source, function(x) {
    if (!is.na(p_map[x])) as.integer(p_map[[x]]) else 9999L
  }, integer(1))
  df <- df[order(df$priority_rank), , drop = FALSE]
  dedup_key <- paste(df$exp_id, df$treatment, df$variable, as.character(df$date_key), sep = "\r")
  keep <- !duplicated(dedup_key)
  df <- df[keep, , drop = FALSE]
  df[, setdiff(names(df), c("date_key", "priority_rank")), drop = FALSE]
}

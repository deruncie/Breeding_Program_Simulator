# Monitoring and plotting helpers for readable-state simulations.
utils::globalVariables(c("year", "value", "stage"))

# Parse trait input and return index + phenotype_log label.
bp_resolve_trait <- function(trait = 1L) {
  if (is.numeric(trait)) {
    idx <- as.integer(trait)[1]
    if (is.na(idx) || idx < 1L) {
      stop("trait must be >= 1", call. = FALSE)
    }
    return(list(index = idx, label = paste0("trait", idx)))
  }

  tchr <- as.character(trait)[1]
  if (grepl("^trait[0-9]+$", tchr)) {
    idx <- as.integer(sub("^trait", "", tchr))
    return(list(index = idx, label = tchr))
  }

  stop("trait must be numeric index (e.g. 1) or label (e.g. 'trait1')", call. = FALSE)
}

# Extract one trait vector from a Pop slot matrix/vector.
bp_pop_trait_vector <- function(x, trait_index) {
  if (is.null(x)) return(rep(NA_real_, 0))
  if (is.null(dim(x))) return(as.numeric(x))
  if (ncol(x) < trait_index) return(rep(NA_real_, nrow(x)))
  as.numeric(x[, trait_index])
}

# Map source_cohort_id to its first source id (semicolon-delimited).
bp_first_source_id <- function(source_cohort_id) {
  s <- as.character(source_cohort_id %||% "")
  if (length(s) == 0L || is.na(s) || !nzchar(s) || identical(s, "NA")) {
    return(NA_character_)
  }
  strsplit(s, ";", fixed = TRUE)[[1]][1]
}

# Trace a cohort back to origin, optionally stopping at an origin_stage.
bp_trace_origin_cohort <- function(state, cohort_id, origin_stage = NULL) {
  cid <- as.character(cohort_id)
  seen <- character()

  while (!is.na(cid) && nzchar(cid) && !(cid %in% seen)) {
    seen <- c(seen, cid)
    idx <- match(cid, state$cohorts$cohort_id)
    if (is.na(idx)) return(NA_character_)

    stg <- as.character(state$cohorts$stage[idx])
    if (!is.null(origin_stage) && stg %in% origin_stage) {
      return(cid)
    }

    next_id <- bp_first_source_id(state$cohorts$source_cohort_id[idx])
    if (is.na(next_id)) return(cid)
    cid <- next_id
  }

  NA_character_
}

#' Monitor Cohorts Over Time
#'
#' Adds created/available/origin year metadata to cohort table rows.
#'
#' @param state Program state.
#' @param ticks_per_year Ticks per year conversion.
#' @param origin_stage Optional stage(s) used to anchor origin tracing.
#'
#' @return Cohort metadata `data.frame` with added year/origin columns.
#' @export
bp_monitor_cohorts <- function(state, ticks_per_year = as.integer(round(1 / state$time$dt)), origin_stage = NULL) {
  if (nrow(state$cohorts) == 0L) return(state$cohorts)
  tpy <- as.integer(ticks_per_year)
  out <- state$cohorts

  out$created_year <- as.integer(floor(out$created_tick / tpy) + 1L)
  out$available_year <- as.integer(floor(out$available_tick / tpy) + 1L)

  out$origin_cohort_id <- vapply(out$cohort_id, function(cid) {
    bp_trace_origin_cohort(state, cid, origin_stage = origin_stage)
  }, character(1))

  origin_idx <- match(out$origin_cohort_id, state$cohorts$cohort_id)
  out$origin_stage <- ifelse(is.na(origin_idx), NA_character_, as.character(state$cohorts$stage[origin_idx]))
  out$origin_year <- ifelse(
    is.na(origin_idx),
    NA_integer_,
    as.integer(floor(state$cohorts$created_tick[origin_idx] / tpy) + 1L)
  )

  out
}

# Estimate h2/H2 on a line-mean basis from phenotype_log and genetic values.
# h2 = Var(A_line) / Var(P_line), H2 = Var(G_line) / Var(P_line),
# where A is additive breeding value and G is total genetic value.
bp_estimate_h2 <- function(state, cohort_id, trait_label, trait_index) {
  idx <- match(cohort_id, state$cohorts$cohort_id)
  if (is.na(idx)) return(c(h2 = NA_real_, H2 = NA_real_))
  pop <- state$pops[[cohort_id]]
  if (is.null(pop)) return(c(h2 = NA_real_, H2 = NA_real_))

  ph <- state$phenotype_log
  ph <- ph[ph$cohort_id == cohort_id & ph$trait == trait_label, , drop = FALSE]
  if (nrow(ph) == 0L) return(c(h2 = NA_real_, H2 = NA_real_))

  gv <- bp_pop_trait_vector(pop@gv, trait_index)
  av <- tryCatch(
    bp_pop_trait_vector(AlphaSimR::bv(pop, simParam = state$sim$SP), trait_index),
    error = function(e) rep(NA_real_, length(gv))
  )
  ids <- as.integer(pop@id)
  g <- gv[match(ph$individual_id, ids)]
  a <- av[match(ph$individual_id, ids)]
  y <- as.numeric(ph$phenotype_value)

  keep <- !is.na(g) & !is.na(a) & !is.na(y) & !is.na(ph$individual_id)
  if (sum(keep) < 3L) return(c(h2 = NA_real_, H2 = NA_real_))

  ph2 <- data.frame(
    id = as.integer(ph$individual_id[keep]),
    y = as.numeric(y[keep]),
    a = as.numeric(a[keep]),
    g = as.numeric(g[keep]),
    stringsAsFactors = FALSE
  )

  y_mean <- stats::aggregate(y ~ id, data = ph2, FUN = mean)
  a_mean <- stats::aggregate(a ~ id, data = ph2, FUN = mean)
  g_mean <- stats::aggregate(g ~ id, data = ph2, FUN = mean)
  mm <- Reduce(function(x, y) merge(x, y, by = "id", all = FALSE), list(y_mean, a_mean, g_mean))
  if (nrow(mm) < 3L) return(c(h2 = NA_real_, H2 = NA_real_))

  va <- stats::var(mm$a)
  vg <- stats::var(mm$g)
  vp <- stats::var(mm$y)
  if (!is.finite(va) || !is.finite(vg) || !is.finite(vp) || vp <= 0 || va < 0 || vg < 0) {
    return(c(h2 = NA_real_, H2 = NA_real_))
  }

  h2_line <- min(max(as.numeric(va / vp), 0), 1)
  H2_line <- min(max(as.numeric(vg / vp), 0), 1)
  c(h2 = h2_line, H2 = H2_line)
}

#' Extract Cohort Metrics
#'
#' Computes per-cohort summary metrics (e.g. mean/variance/max genetic value,
#' EBV-GV correlation, and line-mean h2/H2).
#'
#' @param state Program state.
#' @param stages Optional stage filter.
#' @param trait Trait index or label.
#' @param origin_stage Optional origin-stage filter.
#' @param include_inactive Whether to include inactive cohorts.
#' @param ticks_per_year Ticks per year conversion.
#'
#' @return `data.frame` with cohort metadata and metrics.
#' @export
bp_extract_cohort_metrics <- function(
  state,
  stages = NULL,
  trait = 1L,
  origin_stage = NULL,
  include_inactive = TRUE,
  ticks_per_year = as.integer(round(1 / state$time$dt))
) {
  tr <- bp_resolve_trait(trait)
  meta <- bp_monitor_cohorts(state, ticks_per_year = ticks_per_year, origin_stage = origin_stage)

  if (!isTRUE(include_inactive)) {
    meta <- meta[meta$active, , drop = FALSE]
  }
  if (!is.null(stages)) {
    meta <- meta[meta$stage %in% stages, , drop = FALSE]
  }
  if (nrow(meta) == 0L) return(meta)

  rows <- vector("list", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    pop <- state$pops[[row$cohort_id]]
    gv <- bp_pop_trait_vector(pop@gv, tr$index)
    ebv <- bp_pop_trait_vector(pop@ebv, tr$index)

    mean_gv <- if (length(gv) > 0L && any(!is.na(gv))) mean(gv, na.rm = TRUE) else NA_real_
    var_gv <- if (length(gv) > 1L && any(!is.na(gv))) stats::var(gv, na.rm = TRUE) else NA_real_
    max_gv <- if (length(gv) > 0L && any(!is.na(gv))) max(gv, na.rm = TRUE) else NA_real_

    cor_ebv_gv <- NA_real_
    if (length(ebv) == length(gv) && length(gv) > 2L) {
      keep <- !is.na(gv) & !is.na(ebv)
      if (sum(keep) > 2L) {
        cor_ebv_gv <- stats::cor(ebv[keep], gv[keep])
      }
    }

    h <- bp_estimate_h2(state, row$cohort_id, tr$label, tr$index)

    rows[[i]] <- cbind(
      row[, c(
        "cohort_id", "stage", "stream", "cycle_id", "active", "n_ind",
        "created_tick", "available_tick", "created_year", "available_year",
        "origin_cohort_id", "origin_stage", "origin_year"
      )],
      data.frame(
        trait = tr$label,
        mean_gv = mean_gv,
        var_gv = var_gv,
        max_gv = max_gv,
        cor_ebv_gv = cor_ebv_gv,
        h2 = h[["h2"]],
        H2 = h[["H2"]],
        stringsAsFactors = FALSE
      )
    )
  }

  do.call(rbind, rows)
}

#' Summarize Metric by Year
#'
#' Aggregates one metric column by year and stage.
#'
#' @param metrics_df Metrics data frame from [bp_extract_cohort_metrics()].
#' @param metric Metric column name.
#' @param year_col Year column (`"available_year"` or `"origin_year"`).
#' @param stage_col Stage column name.
#' @param fun Summary function (default `mean`).
#' @param na.rm Remove missing values before aggregation.
#'
#' @return Aggregated `data.frame` with `year`, `stage`, and `value`.
#' @export
bp_summarize_metric_by_year <- function(
  metrics_df,
  metric = "mean_gv",
  year_col = c("available_year", "origin_year"),
  stage_col = "stage",
  fun = mean,
  na.rm = TRUE
) {
  year_col <- match.arg(year_col)
  if (!metric %in% names(metrics_df)) stop("metric column not found", call. = FALSE)
  if (!year_col %in% names(metrics_df)) stop("year column not found", call. = FALSE)
  if (!stage_col %in% names(metrics_df)) stop("stage column not found", call. = FALSE)

  df <- metrics_df[, c(year_col, stage_col, metric), drop = FALSE]
  names(df) <- c("year", "stage", "value")
  df <- df[!is.na(df$year) & !is.na(df$stage), , drop = FALSE]
  if (nrow(df) == 0L) return(df)

  stats::aggregate(value ~ year + stage, data = df, FUN = function(x) fun(x, na.rm = na.rm))
}

#' Plot Metric by Year
#'
#' Creates a `ggplot2` line plot for a selected metric by year and stage.
#'
#' @param metrics_df Metrics data frame from [bp_extract_cohort_metrics()].
#' @param metric Metric column name.
#' @param year_col Year column (`"available_year"` or `"origin_year"`).
#' @param stage_col Stage column name.
#' @param aggregate Whether to aggregate cohorts within stage-year.
#'
#' @return A `ggplot` object.
#' @export
bp_plot_metric_by_year <- function(
  metrics_df,
  metric = "mean_gv",
  year_col = c("available_year", "origin_year"),
  stage_col = "stage",
  aggregate = TRUE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting", call. = FALSE)
  }
  year_col <- match.arg(year_col)

  plot_df <- if (isTRUE(aggregate)) {
    bp_summarize_metric_by_year(metrics_df, metric = metric, year_col = year_col, stage_col = stage_col)
  } else {
    tmp <- metrics_df[, c(year_col, stage_col, metric), drop = FALSE]
    names(tmp) <- c("year", "stage", "value")
    tmp
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(x = year, y = value, color = stage, group = stage)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(x = ifelse(year_col == "available_year", "Available Year", "Origin Year"), y = metric)
}

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

bp_select_report_cohorts <- function(state, stage, stream = NULL, use = c("latest_available", "previous_available", "latest", "all_available")) {
  use <- match.arg(use)
  df <- state$cohorts
  if (!is.null(stage)) df <- df[df$stage %in% stage, , drop = FALSE]
  if (!is.null(stream)) df <- df[df$stream %in% stream, , drop = FALSE]
  if (nrow(df) == 0L) return(df)

  if (use %in% c("latest_available", "previous_available", "all_available")) {
    df <- df[df$available_tick <= as.integer(state$time$tick), , drop = FALSE]
  }
  if (nrow(df) == 0L) return(df)
  df <- df[order(df$available_tick, df$created_tick, decreasing = TRUE), , drop = FALSE]

  if (use == "all_available") return(df)
  if (use == "previous_available") {
    if (nrow(df) < 2L) return(df[0, , drop = FALSE])
    return(df[2L, , drop = FALSE])
  }
  df[1L, , drop = FALSE]
}

bp_report_matrix <- function(pop, traits = NULL, slot = c("gv", "pheno", "ebv")) {
  slot <- match.arg(slot)
  if (!methods::is(pop, "Pop") && !methods::is(pop, "RawPop")) {
    stop("Stage reporting requires AlphaSimR Pop-like populations for built-in metrics.", call. = FALSE)
  }
  mat <- switch(slot, gv = pop@gv, pheno = pop@pheno, ebv = pop@ebv)
  if (is.null(mat)) return(matrix(numeric(0), nrow = pop_n_ind(pop), ncol = 0L))
  mat <- as.matrix(mat)
  if (ncol(mat) == 0L) return(mat)
  if (is.null(traits)) return(mat)
  if (is.numeric(traits)) return(mat[, as.integer(traits), drop = FALSE])
  idx <- match(as.character(traits), colnames(mat))
  if (anyNA(idx)) {
    stop("Some requested traits were not found in the population matrix.", call. = FALSE)
  }
  mat[, idx, drop = FALSE]
}

#' Compute Weighted Index Values
#'
#' @param values Numeric matrix-like object with individuals in rows and traits
#'   in columns.
#' @param weights Numeric index weights.
#' @param traits Optional trait indices or column names to use.
#'
#' @return Numeric vector of weighted index values.
#' @export
bp_index_values <- function(values, weights, traits = NULL) {
  mat <- as.matrix(values)
  if (!is.null(traits)) {
    if (is.numeric(traits)) {
      idx <- as.integer(traits)
      if (any(is.na(idx)) || any(idx < 1L) || any(idx > ncol(mat))) {
        stop("traits must be valid column indices.", call. = FALSE)
      }
      mat <- mat[, idx, drop = FALSE]
    } else {
      if (is.null(colnames(mat))) stop("Character traits require column names in values.", call. = FALSE)
      idx <- match(as.character(traits), colnames(mat))
      if (anyNA(idx)) stop("Some requested traits were not found in values.", call. = FALSE)
      mat <- mat[, idx, drop = FALSE]
    }
  }
  w <- as.numeric(weights)
  if (length(w) != ncol(mat)) stop("weights length must match number of selected traits.", call. = FALSE)
  drop(mat %*% matrix(w, ncol = 1L))
}

bp_named_metric <- function(values, mat) {
  values <- as.numeric(values)
  nms <- colnames(mat)
  if (is.null(nms) || any(!nzchar(nms))) nms <- paste0("trait", seq_along(values))
  stats::setNames(values, nms[seq_along(values)])
}

bp_subset_report_baseline <- function(baseline, traits = NULL, value_names = NULL) {
  if (is.null(baseline)) return(NULL)
  n <- length(value_names %||% traits %||% baseline$mean)
  if (n == 0L) return(baseline)

  idx <- rep(NA_integer_, n)
  baseline_names <- names(baseline$mean)
  has_baseline_names <- !is.null(baseline_names) && all(nzchar(baseline_names))
  valid_baseline_idx <- function(i) {
    !is.na(i) &&
      i >= 1L &&
      i <= length(baseline$mean) &&
      is.finite(as.numeric(baseline$mean[[i]])) &&
      is.finite(as.numeric(baseline$sd[[i]])) &&
      as.numeric(baseline$sd[[i]]) > 0
  }
  candidate_idx <- function(j) {
    candidates <- integer(0)
    if (!is.null(value_names) && isTRUE(has_baseline_names)) {
      exact <- which(value_names[[j]] == baseline_names)
      folded <- which(tolower(value_names[[j]]) == tolower(baseline_names))
      candidates <- c(candidates, exact, folded)
    }
    if (is.numeric(traits) && length(traits) == n) {
      trait_idx <- as.integer(traits[[j]])
      if (!is.na(trait_idx) && trait_idx <= length(baseline$mean)) {
        candidates <- c(candidates, trait_idx)
      }
    }
    if (!isTRUE(has_baseline_names) && j <= length(baseline$mean)) {
      candidates <- c(candidates, j)
    }
    candidates <- unique(candidates[!is.na(candidates)])
    for (i in candidates) {
      if (valid_baseline_idx(i)) return(i)
    }
    NA_integer_
  }

  idx <- vapply(seq_len(n), candidate_idx, integer(1))
  if (anyNA(idx)) return(NULL)

  out <- baseline
  out$mean <- baseline$mean[idx]
  out$sd <- baseline$sd[idx]
  names(out$mean) <- names(out$sd) <- value_names %||% names(out$mean)
  if (!is.null(baseline$cov)) {
    out$cov <- as.matrix(baseline$cov)[idx, idx, drop = FALSE]
    rownames(out$cov) <- colnames(out$cov) <- names(out$mean)
  }
  out
}

bp_scale_report_values <- function(values, baseline, scale) {
  scale <- match.arg(scale, c("none", "mean", "var"))
  if (scale == "none") return(values)
  if (is.null(baseline)) return(values)
  nms <- names(values)
  if (is.null(nms)) {
    b_mean <- as.numeric(baseline$mean)[seq_along(values)]
    b_sd <- as.numeric(baseline$sd)[seq_along(values)]
  } else {
    b_mean <- as.numeric(baseline$mean[nms])
    b_sd <- as.numeric(baseline$sd[nms])
    missing_baseline <- is.na(b_mean) | is.na(b_sd)
    if (any(missing_baseline)) {
      pos <- seq_along(values)
      b_mean[missing_baseline] <- as.numeric(baseline$mean)[pos[missing_baseline]]
      b_sd[missing_baseline] <- as.numeric(baseline$sd)[pos[missing_baseline]]
    }
  }
  if (scale == "mean") return((values - b_mean) / b_sd)
  values / (b_sd^2)
}

bp_report_prediction_matrix <- function(pop, traits = NULL) {
  ebv <- bp_report_matrix(pop, traits = traits, slot = "ebv")
  if (ncol(ebv) > 0L && any(!is.na(ebv))) {
    return(ebv)
  }
  ph <- bp_report_matrix(pop, traits = traits, slot = "pheno")
  if (ncol(ph) > 0L && any(!is.na(ph))) {
    return(ph)
  }
  matrix(NA_real_, nrow = pop_n_ind(pop), ncol = 0L)
}

bp_within_family_accuracy <- function(pop, gv, traits = NULL) {
  pred <- bp_report_prediction_matrix(pop, traits = traits)
  if (ncol(pred) == 0L) return(rep(NA_real_, ncol(gv)))

  fam <- paste(pop@mother, pop@father)
  out <- rep(NA_real_, ncol(gv))
  for (j in seq_len(ncol(gv))) {
    p <- as.numeric(pred[, j])
    a <- as.numeric(gv[, j])
    keep <- !is.na(p) & !is.na(a)
    if (sum(keep) < 3L) next
    p <- p[keep]
    a <- a[keep]
    f <- fam[keep]
    if (length(unique(f)) > 1L) {
      p <- stats::resid(stats::lm(p ~ f))
      a <- stats::resid(stats::lm(a ~ f))
    }
    out[[j]] <- suppressWarnings(stats::cor(p, a, use = "complete.obs"))
  }
  out
}

bp_report_one_stage_metric <- function(
  state,
  stage,
  stream = NULL,
  metric = c("meanG", "maxG", "varG", "H2", "accEBV", "wf_accEBV"),
  traits = NULL,
  use = c("latest_available", "previous_available", "latest", "all_available"),
  baseline = NULL,
  synthetic_trait = NULL,
  cfg = NULL,
  synthetic_varE = NULL,
  synthetic_gv_n_trials = 100L,
  synthetic_gv_n_plants = 10L,
  synthetic_gv_seed = NULL,
  custom = NULL,
  prefix = stage
) {
  metric <- match.arg(metric)
  use <- match.arg(use)
  synthetic_def <- bp_resolve_synthetic_traits(state, synthetic_trait)
  if (length(synthetic_def) > 1L) stop("Reporting accepts one synthetic_trait at a time.", call. = FALSE)
  synthetic_def <- if (length(synthetic_def) == 1L) synthetic_def[[1L]] else NULL
  rows <- bp_select_report_cohorts(state, stage = stage, stream = stream, use = use)
  if (nrow(rows) == 0L) {
    if (!is.null(synthetic_def)) return(stats::setNames(NA_real_, synthetic_def$name))
    return(NA_real_)
  }
  metric_scale <- switch(metric, meanG = "mean", maxG = "mean", varG = "var", H2 = "none", accEBV = "none", wf_accEBV = "none")

  one_metric <- function(row) {
    pop <- state$pops[[as.character(row$cohort_id)]]
    if (!is.null(custom) && metric %in% names(custom)) {
      fn <- custom[[metric]]
      if (!is.function(fn)) stop("custom metrics must be functions.", call. = FALSE)
      return(fn(state, pop, row, list(metric = metric, traits = traits, baseline = baseline, prefix = prefix)))
    }

    if (!is.null(synthetic_def)) {
      gv_vec <- bp_get_stored_synthetic_values(pop, synthetic_def$name, "gv", missing_ok = TRUE)
      if (is.null(gv_vec)) {
        gv_vec <- bp_get_synthetic_gv(
          pop = pop,
          synthetic_trait = synthetic_def,
          state = state,
          varE = synthetic_varE %||% cfg$varE %||% NULL,
          n_trials = synthetic_gv_n_trials,
          n_plants_per_trial = synthetic_gv_n_plants,
          seed = synthetic_gv_seed
        )
      }
      ph_vec <- bp_get_stored_synthetic_values(pop, synthetic_def$name, "pheno", missing_ok = TRUE)
      if (is.null(ph_vec)) ph_vec <- bp_synthetic_values(pop, synthetic_def, "pheno")
      ebv_vec <- bp_get_stored_synthetic_values(pop, synthetic_def$name, "ebv", missing_ok = TRUE)
      if (is.null(ebv_vec)) ebv_vec <- bp_synthetic_values(pop, synthetic_def, "ebv")
      accuracy <- function(pred, act, within_family = FALSE) {
        keep <- !is.na(pred) & !is.na(act)
        if (sum(keep) < 3L) return(NA_real_)
        pred <- pred[keep]
        act <- act[keep]
        if (isTRUE(within_family)) {
          fam <- paste(pop@mother, pop@father)[keep]
          if (length(unique(fam)) > 1L) {
            pred <- stats::resid(stats::lm(pred ~ fam))
            act <- stats::resid(stats::lm(act ~ fam))
          }
        }
        if (stats::sd(pred) <= 0 || stats::sd(act) <= 0) return(NA_real_)
        suppressWarnings(stats::cor(pred, act, use = "complete.obs"))
      }
      value <- switch(
        metric,
        meanG = mean(gv_vec, na.rm = TRUE),
        maxG = max(gv_vec, na.rm = TRUE),
        varG = stats::var(gv_vec, na.rm = TRUE),
        H2 = {
          vg <- stats::var(gv_vec, na.rm = TRUE)
          vp <- stats::var(ph_vec, na.rm = TRUE)
          if (is.finite(vg) && is.finite(vp) && vp > 0) vg / vp else NA_real_
        },
        accEBV = accuracy(ebv_vec, gv_vec),
        wf_accEBV = accuracy(ebv_vec, gv_vec, within_family = TRUE)
      )
      value <- stats::setNames(as.numeric(value), synthetic_def$name)
      if (is.null(baseline)) baseline <- state$sim$trait_baselines[["default"]]
      if (
        !is.null(baseline) &&
        (!(synthetic_def$name %in% names(baseline$mean)) ||
          !(synthetic_def$name %in% names(baseline$sd)))
      ) {
        baseline <- NULL
      }
      return(bp_scale_report_values(value, baseline, metric_scale))
    }

    gv <- bp_report_matrix(pop, traits = traits, slot = "gv")
    values <- switch(
      metric,
      meanG = colMeans(gv, na.rm = TRUE),
      maxG = apply(gv, 2L, max, na.rm = TRUE),
      varG = apply(gv, 2L, stats::var, na.rm = TRUE),
      H2 = {
        ph <- bp_report_matrix(pop, traits = traits, slot = "pheno")
        out <- rep(NA_real_, ncol(gv))
        for (j in seq_len(ncol(gv))) {
          vg <- stats::var(gv[, j], na.rm = TRUE)
          vp <- stats::var(ph[, j], na.rm = TRUE)
          out[[j]] <- if (is.finite(vg) && is.finite(vp) && vp > 0) vg / vp else NA_real_
        }
        out
      },
      accEBV = {
        ebv <- bp_report_prediction_matrix(pop, traits = traits)
        out <- rep(NA_real_, ncol(gv))
        if (ncol(ebv) > 0L) {
          for (j in seq_len(ncol(gv))) {
            keep <- !is.na(gv[, j]) & !is.na(ebv[, j])
            out[[j]] <- if (sum(keep) > 2L && stats::sd(gv[keep, j]) > 0 && stats::sd(ebv[keep, j]) > 0) {
              stats::cor(gv[keep, j], ebv[keep, j])
            } else {
              NA_real_
            }
          }
        }
        out
      },
      wf_accEBV = bp_within_family_accuracy(pop, gv, traits = traits),
      stop(sprintf("Unknown metric: %s", metric), call. = FALSE)
    )
    values <- bp_named_metric(values, gv)
    if (is.null(baseline)) baseline <- state$sim$trait_baselines[["default"]]
    baseline <- bp_subset_report_baseline(baseline, traits = traits, value_names = names(values))
    bp_scale_report_values(values, baseline, metric_scale)
  }

  vals <- lapply(seq_len(nrow(rows)), function(i) one_metric(rows[i, , drop = FALSE]))
  if (use != "all_available") return(vals[[1L]])
  vals
}

#' Report Stage Metrics
#'
#' Return reporting values as a flat named list of unnamed numeric scalars.
#' Named `metrics` provide the base output names. Biological trait or synthetic
#' trait labels can be appended predictably, making the result safe to combine
#' incrementally with `c()` before creating a one-row `data.frame`.
#'
#' Synthetic genetic values are materialized once in a local copy of `state`
#' and reused across all requested metrics. The input state is not modified.
#'
#' @param state Program state.
#' @param stage Stage name.
#' @param stream Optional stream filter.
#' @param metrics Character vector containing `meanG`, `maxG`, `varG`, `H2`,
#'   `accEBV`, or `wf_accEBV`. Supply names to set output base names, for
#'   example `c(meanCandidates = "meanG", varCandidates = "varG")`.
#' @param traits Optional biological trait indices or names.
#' @param use Cohort selection rule. `all_available` is not supported because
#'   this function returns one flat reporting row.
#' @param baseline Baseline object. If `NULL`, uses
#'   `state$sim$trait_baselines[["default"]]`.
#' @param synthetic_trait Optional synthetic-trait definition or registered
#'   name. Existing numeric `traits` behavior is unchanged when omitted.
#' @param cfg Optional configuration list. `cfg$varE` is used for nonlinear
#'   synthetic GV integration when `synthetic_varE` is omitted.
#' @param synthetic_varE Optional residual covariance override.
#' @param synthetic_gv_n_trials Number of random TPE trials.
#' @param synthetic_gv_n_plants Number of plant residual draws per trial.
#' @param synthetic_gv_seed Optional Monte Carlo seed.
#' @param custom Optional named list of metric functions.
#' @param append_trait Whether to append trait labels to output names:
#'   `"auto"` appends when a metric returns multiple traits, `"always"`
#'   appends even for one trait, and `"never"` does not append.
#' @param name_fn Optional function `(base_name, metric, stage, trait)` returning
#'   one output name. When supplied, it replaces `append_trait` naming.
#'
#' @return Flat named list of unnamed numeric scalars.
#'
#' @examples
#' \dontrun{
#' biological <- bp_report_stage_metrics(
#'   state,
#'   stage = "Candidates",
#'   metrics = c(meanCandidates = "meanG", varCandidates = "varG"),
#'   traits = 1:2,
#'   append_trait = "always"
#' )
#' index <- bp_report_stage_metrics(
#'   state,
#'   stage = "Candidates",
#'   metrics = c(meanCandidates = "meanG"),
#'   synthetic_trait = "Index",
#'   append_trait = "always"
#' )
#' results_year <- c(results_year, biological, index)
#' }
#' @export
bp_report_stage_metrics <- function(
  state,
  stage,
  stream = NULL,
  metrics = c("meanG", "varG", "accEBV"),
  traits = NULL,
  use = c("latest_available", "previous_available", "latest", "all_available"),
  baseline = NULL,
  synthetic_trait = NULL,
  cfg = NULL,
  synthetic_varE = NULL,
  synthetic_gv_n_trials = 100L,
  synthetic_gv_n_plants = 10L,
  synthetic_gv_seed = NULL,
  custom = NULL,
  append_trait = c("auto", "always", "never"),
  name_fn = NULL
) {
  use <- match.arg(use)
  if (use == "all_available") {
    stop("bp_report_stage_metrics does not support use = 'all_available'.", call. = FALSE)
  }
  append_trait <- match.arg(append_trait)
  valid_metrics <- c("meanG", "maxG", "varG", "H2", "accEBV", "wf_accEBV")
  metric_names <- names(metrics)
  metrics <- as.character(metrics)
  if (length(metrics) == 0L || anyNA(metrics) || any(!metrics %in% valid_metrics)) {
    stop(sprintf("metrics must contain only: %s.", paste(valid_metrics, collapse = ", ")), call. = FALSE)
  }
  if (is.null(metric_names)) metric_names <- metrics
  blank_names <- is.na(metric_names) | !nzchar(metric_names)
  metric_names[blank_names] <- metrics[blank_names]
  if (!is.null(name_fn) && !is.function(name_fn)) stop("name_fn must be a function.", call. = FALSE)

  defs <- bp_resolve_synthetic_traits(state, synthetic_trait)
  if (length(defs) > 1L) stop("bp_report_stage_metrics accepts one synthetic_trait at a time.", call. = FALSE)
  state_for_report <- state
  if (length(defs) == 1L) {
    rows <- bp_select_report_cohorts(state, stage = stage, stream = stream, use = use)
    for (cid in as.character(rows$cohort_id)) {
      state_for_report$pops[[cid]] <- bp_materialize_synthetic_gv(
        pop = state_for_report$pops[[cid]],
        synthetic_traits = defs[[1L]],
        state = state_for_report,
        varE = synthetic_varE %||% cfg$varE %||% NULL,
        n_trials = synthetic_gv_n_trials,
        n_plants_per_trial = synthetic_gv_n_plants,
        seed = synthetic_gv_seed
      )
    }
  }

  values_by_metric <- lapply(metrics, function(metric) {
    bp_report_one_stage_metric(
      state = state_for_report,
      stage = stage,
      stream = stream,
      metric = metric,
      traits = traits,
      use = use,
      baseline = baseline,
      synthetic_trait = synthetic_trait,
      cfg = cfg,
      synthetic_varE = synthetic_varE,
      synthetic_gv_n_trials = synthetic_gv_n_trials,
      synthetic_gv_n_plants = synthetic_gv_n_plants,
      synthetic_gv_seed = synthetic_gv_seed,
      custom = custom,
      prefix = stage
    )
  })

  expected_traits <- if (length(defs) == 1L) {
    defs[[1L]]$name
  } else {
    bp_report_expected_names(state, traits = traits, baseline = baseline)
  }
  out <- list()
  for (i in seq_along(values_by_metric)) {
    values <- values_by_metric[[i]]
    if (length(values) == 1L && is.na(values) && is.null(names(values)) && length(expected_traits) > 0L) {
      values <- stats::setNames(rep(NA_real_, length(expected_traits)), expected_traits)
    }
    trait_names <- names(values)
    if (is.null(trait_names) || any(!nzchar(trait_names))) {
      trait_names <- if (length(expected_traits) == length(values)) expected_traits else paste0("value", seq_along(values))
    }
    add_suffix <- switch(
      append_trait,
      always = rep(TRUE, length(values)),
      never = rep(FALSE, length(values)),
      auto = rep(length(values) > 1L, length(values))
    )
    output_names <- vapply(seq_along(values), function(j) {
      if (!is.null(name_fn)) {
        as.character(name_fn(metric_names[[i]], metrics[[i]], stage, trait_names[[j]]))[[1L]]
      } else if (add_suffix[[j]]) {
        paste(metric_names[[i]], trait_names[[j]], sep = "_")
      } else {
        metric_names[[i]]
      }
    }, character(1))
    if (any(!nzchar(output_names))) stop("Reporting output names must be non-empty.", call. = FALSE)
    if (any(output_names %in% names(out)) || anyDuplicated(output_names)) {
      stop("Reporting output names are duplicated; use append_trait or name_fn to make them unique.", call. = FALSE)
    }
    for (j in seq_along(values)) out[[output_names[[j]]]] <- unname(as.numeric(values[[j]]))
  }
  out
}

bp_report_expected_names <- function(state, traits = NULL, baseline = NULL) {
  if (is.null(baseline) && !is.null(state$sim$trait_baselines)) {
    baseline <- state$sim$trait_baselines[["default"]]
  }
  base_names <- if (!is.null(baseline) && !is.null(baseline$mean)) {
    setdiff(names(baseline$mean), "Index")
  } else {
    character(0)
  }
  if (is.null(traits)) {
    out <- base_names
  } else if (is.numeric(traits)) {
    idx <- as.integer(traits)
    out <- if (length(base_names) >= max(idx, na.rm = TRUE)) base_names[idx] else paste0("trait", idx)
  } else {
    out <- as.character(traits)
  }
  out <- out[!is.na(out) & nzchar(out)]
  out
}

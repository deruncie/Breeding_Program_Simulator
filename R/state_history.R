# State-history, reporting setup, safe mutation, and cfg discovery helpers.

bp_validate_state_dt <- function(state) {
  if (is.null(state$time) || is.null(state$time$dt)) {
    stop("state$time$dt is required.", call. = FALSE)
  }
  dt <- state$time$dt
  if (!is.numeric(dt) || length(dt) != 1L || is.na(dt) || dt <= 0) {
    stop("state$time$dt must be a single positive numeric value.", call. = FALSE)
  }
  as.numeric(dt)
}

#' Convert Years to Ticks
#'
#' Convert a calendar duration in years to simulation ticks using
#' `state$time$dt`.
#'
#' @param state Program state.
#' @param years Single non-negative numeric duration in years.
#' @param require_integer Whether the duration must map to an integer number of
#'   ticks.
#' @param tol Tolerance used when checking integer representability.
#'
#' @return Integer number of ticks.
#' @export
bp_years_to_ticks <- function(state, years, require_integer = TRUE, tol = 1e-8) {
  dt <- bp_validate_state_dt(state)
  if (!is.numeric(years) || length(years) != 1L || is.na(years)) {
    stop("`years` must be a single non-missing numeric value.", call. = FALSE)
  }
  if (years < 0) {
    stop("`years` must be non-negative.", call. = FALSE)
  }
  raw_ticks <- as.numeric(years) / dt
  rounded_ticks <- round(raw_ticks)
  if (isTRUE(require_integer) && abs(raw_ticks - rounded_ticks) > tol) {
    stop(
      sprintf(
        "`years = %s` cannot be represented exactly with state$time$dt = %s. It gives %s ticks.",
        years, dt, raw_ticks
      ),
      call. = FALSE
    )
  }
  as.integer(rounded_ticks)
}

#' Advance Time by Years
#'
#' Advance simulation time by a biological/calendar duration rather than by raw
#' ticks.
#'
#' @inheritParams bp_years_to_ticks
#'
#' @return Updated program state.
#' @export
bp_advance_time_years <- function(state, years, require_integer = TRUE, tol = 1e-8) {
  n_ticks <- bp_years_to_ticks(
    state = state,
    years = years,
    require_integer = require_integer,
    tol = tol
  )
  bp_advance_time(state, n_ticks = n_ticks)
}

#' Ticks per Year
#'
#' @param state Program state.
#'
#' @return Integer number of simulation ticks per calendar year.
#' @export
bp_ticks_per_year <- function(state) {
  dt <- bp_validate_state_dt(state)
  raw <- 1 / dt
  rounded <- round(raw)
  if (abs(raw - rounded) > 1e-8) {
    stop(sprintf("state$time$dt = %s does not divide one year into an integer number of ticks.", dt), call. = FALSE)
  }
  as.integer(rounded)
}

#' Current Tick Fraction Within Year
#'
#' @param state Program state.
#'
#' @return Numeric fraction in `[0, 1)`.
#' @export
bp_tick_fraction <- function(state) {
  tpy <- bp_ticks_per_year(state)
  tick_in_year <- as.integer(round(state$time$tick)) %% tpy
  tick_in_year / tpy
}

#' Test Current Tick Fraction
#'
#' @param state Program state.
#' @param fraction Numeric fraction within the year, e.g. `0`, `0.25`, `0.5`.
#' @param tol Numeric tolerance.
#'
#' @return Logical scalar.
#' @export
bp_is_fraction_tick <- function(state, fraction, tol = 1e-8) {
  if (!is.numeric(fraction) || length(fraction) != 1L || is.na(fraction)) {
    stop("`fraction` must be a single non-missing numeric value.", call. = FALSE)
  }
  abs(bp_tick_fraction(state) - as.numeric(fraction)) <= tol
}

#' Reset Cost Log
#'
#' @param state Program state.
#'
#' @return Updated program state.
#' @export
bp_reset_costs <- function(state) {
  state$cost_log <- bp_empty_cost_log()
  state
}

bp_empty_like <- function(df) {
  if (is.null(df)) return(NULL)
  df[NULL, , drop = FALSE]
}

bp_filter_stage_stream <- function(df, stages = "all", streams = "all") {
  keep <- rep(TRUE, nrow(df))
  if (!identical(stages, "all") && !is.null(stages)) {
    keep <- keep & df$stage %in% as.character(stages)
  }
  if (!identical(streams, "all") && !is.null(streams)) {
    keep <- keep & df$stream %in% as.character(streams)
  }
  keep
}

bp_remap_source_ids <- function(x, id_map) {
  ids <- bp_parse_source_ids(x)
  ids <- ids[ids %in% names(id_map)]
  if (length(ids) == 0L) return(NA_character_)
  paste(unname(id_map[ids]), collapse = ";")
}

bp_shift_tick_col <- function(df, col, offset) {
  if (!is.null(df) && col %in% names(df)) {
    df[[col]] <- as.integer(df[[col]] - as.integer(offset))
  }
  df
}

#' Forget Historical State
#'
#' Return a state suitable for starting a new future experiment from the current
#' program while dropping historical logs, old populations, costs, and usually
#' old models.
#'
#' @param state Program state.
#' @param keep_stages Stage filter or `"all"`.
#' @param keep_streams Stream filter or `"all"`.
#' @param keep_available Keep active cohorts already available at current time.
#' @param keep_future Keep active cohorts that become available in the future.
#' @param keep_latest_by_stage Preserve the latest available cohort per kept
#'   stage/stream even if it is inactive.
#' @param reset_time If `TRUE`, shift the retained state so current time becomes
#'   tick/year zero using the existing `state$time$dt`.
#' @param reset_costs If `TRUE`, clear `state$cost_log`.
#' @param keep_models If `TRUE`, retain `state$gs_models`.
#' @param label_prefix Prefix used to rewrite retained cohort ids. Use `NULL`
#'   or `""` to keep ids unchanged.
#' @param keep_event_log Whether to keep retained event-log rows. Default
#'   `FALSE` records only an initialization/provenance event.
#'
#' @return Updated `bp_state`.
#' @export
bp_forget_history <- function(
  state,
  keep_stages = "all",
  keep_streams = "all",
  keep_available = TRUE,
  keep_future = TRUE,
  keep_latest_by_stage = TRUE,
  reset_time = FALSE,
  reset_costs = TRUE,
  keep_models = FALSE,
  label_prefix = "Initialization",
  keep_event_log = FALSE
) {
  if (is.null(state$cohorts) || nrow(state$cohorts) == 0L) {
    if (isTRUE(reset_costs)) state <- bp_reset_costs(state)
    if (!isTRUE(keep_models)) state$gs_models <- list()
    return(state)
  }

  now <- as.integer(state$time$tick)
  eligible <- bp_filter_stage_stream(state$cohorts, keep_stages, keep_streams)
  keep <- rep(FALSE, nrow(state$cohorts))
  active <- state$cohorts$active & eligible
  if (isTRUE(keep_available)) {
    keep <- keep | (active & state$cohorts$available_tick <= now)
  }
  if (isTRUE(keep_future)) {
    keep <- keep | (active & state$cohorts$available_tick > now)
  }
  if (!isTRUE(keep_available) && !isTRUE(keep_future)) {
    keep <- keep | active
  }

  if (isTRUE(keep_latest_by_stage)) {
    elig_rows <- state$cohorts[eligible & state$cohorts$available_tick <= now, , drop = FALSE]
    if (nrow(elig_rows) > 0L) {
      keys <- paste(elig_rows$stage, elig_rows$stream, sep = "\r")
      for (key in unique(keys)) {
        rows <- elig_rows[keys == key, , drop = FALSE]
        ord <- order(rows$available_tick, rows$created_tick, decreasing = TRUE)
        keep <- keep | state$cohorts$cohort_id %in% rows$cohort_id[ord[[1L]]]
      }
    }
  }

  keep_ids_old <- as.character(state$cohorts$cohort_id[keep])
  id_map <- stats::setNames(keep_ids_old, keep_ids_old)
  if (!is.null(label_prefix) && nzchar(as.character(label_prefix))) {
    id_map <- stats::setNames(paste(as.character(label_prefix), keep_ids_old, sep = "_"), keep_ids_old)
  }

  new_state <- state
  new_state$cohorts <- state$cohorts[keep, , drop = FALSE]
  rownames(new_state$cohorts) <- NULL
  new_state$cohorts$cohort_id <- unname(id_map[as.character(new_state$cohorts$cohort_id)])
  new_state$cohorts$source_cohort_id <- vapply(
    new_state$cohorts$source_cohort_id,
    bp_remap_source_ids,
    character(1),
    id_map = id_map
  )

  new_state$pops <- state$pops[keep_ids_old]
  names(new_state$pops) <- unname(id_map[keep_ids_old])

  new_state$phenotype_log <- bp_empty_phenotype_log()
  new_state$genotype_log <- bp_empty_genotype_log()
  new_state$cost_log <- if (isTRUE(reset_costs)) bp_empty_cost_log() else {
    costs <- state$cost_log[state$cost_log$cohort_id %in% keep_ids_old, , drop = FALSE]
    costs$cohort_id <- unname(id_map[as.character(costs$cohort_id)])
    costs
  }
  new_state$gs_models <- if (isTRUE(keep_models)) state$gs_models else list()

  if (isTRUE(keep_event_log)) {
    ev <- state$event_log
    if (!is.null(ev) && nrow(ev) > 0L) {
      uses_kept <- vapply(seq_len(nrow(ev)), function(i) {
        ids <- unique(c(bp_parse_source_ids(ev$source_ids[[i]]), bp_parse_source_ids(ev$output_id[[i]])))
        any(ids %in% keep_ids_old)
      }, logical(1))
      ev <- ev[uses_kept, , drop = FALSE]
      ev$source_ids <- vapply(ev$source_ids, bp_remap_source_ids, character(1), id_map = id_map)
      ev$output_id <- ifelse(ev$output_id %in% names(id_map), unname(id_map[ev$output_id]), ev$output_id)
      new_state$event_log <- ev
    } else {
      new_state$event_log <- bp_empty_event_log()
    }
  } else {
    new_state$event_log <- bp_empty_event_log()
  }

  if (isTRUE(reset_time)) {
    for (col in c("created_tick", "done_tick", "available_tick", "closed_tick")) {
      new_state$cohorts <- bp_shift_tick_col(new_state$cohorts, col, now)
    }
    new_state$cost_log <- bp_shift_tick_col(new_state$cost_log, "tick", now)
    new_state$event_log <- bp_shift_tick_col(new_state$event_log, "tick", now)
    if (!is.null(new_state$event_log) && "year" %in% names(new_state$event_log)) {
      new_state$event_log$year <- bp_tick_to_year(new_state, new_state$event_log$tick)
    }
    new_state$time$tick <- 0L
    new_state$time$t <- 0
  }

  new_state$counters$cohort <- as.integer(nrow(new_state$cohorts))
  new_state$counters$model <- as.integer(length(new_state$gs_models))
  new_state$counters$event <- as.integer(nrow(new_state$event_log))

  new_state <- bp_log_event(
    new_state,
    fn = "bp_forget_history",
    event_type = "initialization",
    stage = unique(new_state$cohorts$stage),
    source_ids = character(0),
    output_id = NA_character_,
    event_string = sprintf("Initialized state from %d retained cohort(s) after forgetting history.", nrow(new_state$cohorts)),
    template_string = "Forget history and initialize retained state",
    details = list(
      retained_cohorts = as.character(new_state$cohorts$cohort_id),
      reset_time = isTRUE(reset_time),
      reset_costs = isTRUE(reset_costs),
      keep_models = isTRUE(keep_models)
    )
  )
  bp_refresh_genotyped_flags(new_state)
}

#' Compact Historical Populations
#'
#' Downsample older stored populations while preserving recent, current, and
#' future cohorts.
#'
#' @param state Program state.
#' @param stages Optional stages to compact.
#' @param streams Optional streams to compact.
#' @param max_n Maximum stored population size after compaction.
#' @param keep_recent_nCohorts Number of most recent old cohorts to keep per
#'   stage/stream.
#' @param keep_available Keep currently available active cohorts at full size.
#' @param keep_future Keep future active cohorts at full size.
#' @param selection Sampling method: `"random"` or `"first"`.
#' @param seed Optional random seed.
#' @param log_event Whether to append a compaction event.
#'
#' @return Updated program state.
#' @export
bp_compact_history <- function(
  state,
  stages = NULL,
  streams = NULL,
  max_n = 100L,
  keep_recent_nCohorts = 0L,
  keep_available = TRUE,
  keep_future = TRUE,
  selection = c("random", "first"),
  seed = NULL,
  log_event = TRUE
) {
  selection <- match.arg(selection)
  max_n <- as.integer(max_n)
  if (is.na(max_n) || max_n < 1L) stop("max_n must be a positive integer.", call. = FALSE)
  keep_recent_nCohorts <- as.integer(keep_recent_nCohorts %||% 0L)
  if (is.na(keep_recent_nCohorts) || keep_recent_nCohorts < 0L) {
    stop("keep_recent_nCohorts must be a non-negative integer.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  if (is.null(state$cohorts) || nrow(state$cohorts) == 0L) return(state)

  now <- as.integer(state$time$tick)
  eligible <- bp_filter_stage_stream(state$cohorts, stages %||% "all", streams %||% "all")
  protected <- rep(FALSE, nrow(state$cohorts))
  if (isTRUE(keep_available)) {
    protected <- protected | (state$cohorts$active & state$cohorts$available_tick <= now)
  }
  if (isTRUE(keep_future)) {
    protected <- protected | (state$cohorts$active & state$cohorts$available_tick > now)
  }
  if (keep_recent_nCohorts > 0L) {
    keys <- paste(state$cohorts$stage, state$cohorts$stream, sep = "\r")
    for (key in unique(keys[eligible])) {
      idx <- which(eligible & keys == key)
      rows <- state$cohorts[idx, , drop = FALSE]
      ord <- order(rows$available_tick, rows$created_tick, decreasing = TRUE)
      keep_idx <- idx[ord][seq_len(min(keep_recent_nCohorts, length(idx)))]
      protected[keep_idx] <- TRUE
    }
  }

  compacted <- list()
  for (i in which(eligible & !protected)) {
    cid <- as.character(state$cohorts$cohort_id[[i]])
    pop <- state$pops[[cid]]
    if (is.null(pop)) next
    n_old <- pop_n_ind(pop)
    if (n_old <= max_n) next
    idx <- if (selection == "random") sort(sample.int(n_old, max_n)) else seq_len(max_n)
    state$pops[[cid]] <- pop_subset(pop, idx)
    n_new <- pop_n_ind(state$pops[[cid]])
    compacted[[length(compacted) + 1L]] <- data.frame(
      tick = as.integer(state$time$tick),
      cohort_id = cid,
      original_n_ind = as.integer(n_old),
      stored_n_ind = as.integer(n_new),
      selection = selection,
      stringsAsFactors = FALSE
    )
  }

  if (length(compacted) > 0L) {
    rows <- do.call(rbind, compacted)
    if (is.null(state$sim$compaction_log)) {
      state$sim$compaction_log <- rows[NULL, , drop = FALSE]
    }
    state$sim$compaction_log <- rbind(state$sim$compaction_log, rows)
    if (isTRUE(log_event)) {
      state <- bp_log_event(
        state,
        fn = "bp_compact_history",
        event_type = "history_compaction",
        stage = unique(state$cohorts$stage[match(rows$cohort_id, state$cohorts$cohort_id)]),
        source_ids = rows$cohort_id,
        output_id = NA_character_,
        event_string = sprintf("Compacted %d historical cohort population(s) to at most n=%d stored individuals.", nrow(rows), max_n),
        template_string = "Compact historical populations",
        details = list(compacted = rows)
      )
    }
  }
  state
}

bp_pop_ids <- function(pop) {
  if (methods::is(pop, "Pop") || methods::is(pop, "RawPop")) {
    return(as.character(pop@id))
  }
  rn <- tryCatch(rownames(pop), error = function(e) NULL)
  if (!is.null(rn) && length(rn) == pop_n_ind(pop) && !identical(rn, as.character(seq_len(length(rn))))) {
    return(as.character(rn))
  }
  as.character(seq_len(pop_n_ind(pop)))
}

#' Store Individual-Level Auxiliary Data on a Population
#'
#' Add or replace a named `misc` field containing one value, row, or list
#' element per individual. AlphaSimR carries these values with individuals when
#' a population is subset or selected.
#'
#' @param pop AlphaSimR `Pop` object.
#' @param name Name used in `pop@misc`. Names beginning with `bps_` are
#'   reserved for package-managed fields.
#' @param values Atomic vector, factor, list, or matrix aligned with the
#'   individuals in `pop`. Matrices must have one row per individual; other
#'   objects must have one element per individual.
#'
#' @return Updated population.
#' @export
bp_set_misc_values <- function(pop, name, values) {
  if (!methods::is(pop, "Pop")) {
    stop("pop must be an AlphaSimR Pop object.", call. = FALSE)
  }
  name <- as.character(name)
  if (length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("name must be a single non-empty character value.", call. = FALSE)
  }
  if (startsWith(name, "bps_")) {
    stop("Names beginning with 'bps_' are reserved for package-managed fields.", call. = FALSE)
  }

  n_values <- if (is.matrix(values)) nrow(values) else length(values)
  if (n_values != pop_n_ind(pop)) {
    stop("values must contain one value or row per individual.", call. = FALSE)
  }

  pop@misc[[name]] <- values
  pop
}

#' Update an Existing Stage Population
#'
#' Safely replace a stored population object after validating population size and
#' individual identity/order.
#'
#' @param state Program state.
#' @param cohort_id Cohort id to update.
#' @param pop Replacement population object.
#' @param require_same_ids Require same individual ids.
#' @param allow_reorder Allow replacement to be reordered to the original order.
#' @param log_event Whether to append an update event.
#'
#' @return Updated program state.
#' @export
bp_update_stage_pop <- function(
  state,
  cohort_id,
  pop,
  require_same_ids = TRUE,
  allow_reorder = FALSE,
  log_event = TRUE
) {
  cohort_id <- as.character(cohort_id)
  if (length(cohort_id) != 1L || is.na(cohort_id) || !nzchar(cohort_id)) {
    stop("cohort_id must be a single non-missing character value.", call. = FALSE)
  }
  if (is.null(state$pops[[cohort_id]])) {
    stop(sprintf("cohort_id '%s' not found in state$pops.", cohort_id), call. = FALSE)
  }
  old_pop <- state$pops[[cohort_id]]
  old_n <- pop_n_ind(old_pop)
  new_n <- pop_n_ind(pop)
  if (old_n != new_n) {
    stop(sprintf("Replacement pop has %d individuals; expected %d.", new_n, old_n), call. = FALSE)
  }
  if (isTRUE(require_same_ids)) {
    old_ids <- bp_pop_ids(old_pop)
    new_ids <- bp_pop_ids(pop)
    if (!setequal(old_ids, new_ids)) {
      stop("Replacement pop does not contain the same individual ids as the stored cohort.", call. = FALSE)
    }
    if (!identical(old_ids, new_ids)) {
      if (!isTRUE(allow_reorder)) {
        stop("Replacement pop has the same ids but in a different order; set allow_reorder = TRUE to reorder.", call. = FALSE)
      }
      ord <- match(old_ids, new_ids)
      pop <- pop_subset(pop, ord)
    }
  }
  state$pops[[cohort_id]] <- pop
  if (isTRUE(log_event)) {
    row <- state$cohorts[match(cohort_id, state$cohorts$cohort_id), , drop = FALSE]
    state <- bp_log_event(
      state,
      fn = "bp_update_stage_pop",
      event_type = "population_update",
      stage = if (nrow(row) == 0L) "" else as.character(row$stage[[1L]]),
      source_ids = cohort_id,
      output_id = cohort_id,
      event_string = sprintf("Updated stored population for cohort %s after validating %d individuals.", cohort_id, new_n),
      template_string = "Update stored cohort population",
      details = list(cohort_id = cohort_id, n_ind = as.integer(new_n))
    )
  }
  state
}

#' Construct a Residual Covariance Matrix
#'
#' Build a residual covariance matrix from trait heritabilities, genetic
#' variances, and a residual correlation matrix.
#'
#' @param h2 Numeric vector of trait heritabilities.
#' @param varG Numeric vector of genetic variances on the reporting scale.
#' @param corE Residual correlation matrix.
#' @param trait_names Optional trait names for rows and columns.
#'
#' @return Numeric residual covariance matrix.
#' @export
bp_make_varE <- function(h2, varG, corE, trait_names = NULL) {
  h2 <- as.numeric(h2)
  varG <- as.numeric(varG)
  if (length(h2) == 0L || length(h2) != length(varG)) {
    stop("h2 and varG must be numeric vectors of the same non-zero length.", call. = FALSE)
  }
  corE <- as.matrix(corE)
  if (!all(dim(corE) == length(h2))) {
    stop("corE dimensions must match length(h2).", call. = FALSE)
  }
  if (any(!is.finite(h2)) || any(h2 < 0) || any(h2 > 1)) {
    stop("h2 values must be finite and between 0 and 1.", call. = FALSE)
  }
  if (any(!is.finite(varG)) || any(varG < 0)) {
    stop("varG values must be finite and non-negative.", call. = FALSE)
  }
  if (any(!is.finite(corE))) {
    stop("corE must contain only finite values.", call. = FALSE)
  }
  varE_diag <- ((1 - h2) / pmax(h2, 1e-8)) * varG
  D <- diag(sqrt(varE_diag), nrow = length(varE_diag))
  out <- D %*% corE %*% D
  if (!is.null(trait_names)) {
    if (length(trait_names) != length(h2)) stop("trait_names length must match length(h2).", call. = FALSE)
    rownames(out) <- colnames(out) <- as.character(trait_names)
  }
  out
}

#' Set Trait Reporting Baseline
#'
#' Store trait-level reporting normalization in `state$sim$trait_baselines`.
#'
#' @param state Program state.
#' @param pop Optional population used to derive baseline from `@gv`.
#' @param values Optional list/vector of baseline values. Lists may contain
#'   `mean`, `sd`, and optionally `cov`.
#' @param covariance Optional trait covariance matrix used for index scaling
#'   when `values` are supplied.
#' @param traits Optional trait indices or names.
#' @param names Optional trait names.
#' @param label Baseline label.
#' @param include_index Whether to add an index baseline.
#' @param index_weights Optional index weights.
#' @param synthetic_traits Optional synthetic-trait definitions or registered
#'   names to add to the baseline.
#' @param varE Residual covariance used for nonlinear synthetic GVs.
#' @param synthetic_gv_n_trials Number of TPE trials for nonlinear synthetic
#'   GV integration.
#' @param synthetic_gv_n_plants Number of plant residual draws per trial.
#' @param seed Optional Monte Carlo seed.
#'
#' @return Updated program state.
#' @export
bp_set_trait_baseline <- function(
  state,
  pop = NULL,
  values = NULL,
  covariance = NULL,
  traits = NULL,
  names = NULL,
  label = "default",
  include_index = FALSE,
  index_weights = NULL,
  synthetic_traits = NULL,
  varE = NULL,
  synthetic_gv_n_trials = 100L,
  synthetic_gv_n_plants = 10L,
  seed = NULL
) {
  label <- as.character(label %||% "default")
  if (!nzchar(label)) stop("label must be non-empty.", call. = FALSE)
  baseline_gv <- NULL

  if (!is.null(pop)) {
    if (!methods::is(pop, "Pop") && !methods::is(pop, "RawPop")) {
      stop("pop-derived baselines require an AlphaSimR Pop-like object with @gv.", call. = FALSE)
    }
    gv <- as.matrix(pop@gv)
    if (is.null(traits)) traits <- seq_len(ncol(gv))
    if (is.numeric(traits)) {
      gv <- gv[, as.integer(traits), drop = FALSE]
    } else {
      gv <- gv[, match(as.character(traits), colnames(gv)), drop = FALSE]
    }
    trait_names <- names %||% colnames(gv) %||% paste0("trait", seq_len(ncol(gv)))
    baseline_gv <- gv
    means <- colMeans(gv, na.rm = TRUE)
    sds <- apply(gv, 2L, stats::sd, na.rm = TRUE)
    names(means) <- names(sds) <- trait_names
    covariance <- stats::cov(gv, use = "pairwise.complete.obs")
    colnames(covariance) <- rownames(covariance) <- trait_names
  } else {
    if (is.null(values)) stop("Provide either pop or values.", call. = FALSE)
    means <- if (is.list(values)) values$mean %||% values$means else values
    sds <- if (is.list(values)) values$sd %||% values$sds else rep(1, length(means))
    covariance <- covariance %||% if (is.list(values)) values$cov %||% values$covariance else NULL
    input_names <- names %||% names(means) %||% names(sds)
    means <- as.numeric(means)
    sds <- as.numeric(sds)
    if (length(sds) == 1L && length(means) > 1L) sds <- rep(sds, length(means))
    if (length(means) != length(sds)) stop("Baseline mean and sd lengths must match.", call. = FALSE)
    trait_names <- input_names %||% paste0("trait", seq_along(means))
    names(means) <- names(sds) <- trait_names
    if (!is.null(covariance)) {
      covariance <- as.matrix(covariance)
      if (!all(dim(covariance) == length(means))) {
        stop("covariance dimensions must match baseline traits.", call. = FALSE)
      }
      if (is.null(rownames(covariance))) rownames(covariance) <- trait_names
      if (is.null(colnames(covariance))) colnames(covariance) <- trait_names
      covariance <- covariance[trait_names, trait_names, drop = FALSE]
    }
  }

  bad_sd <- !is.finite(sds) | sds <= 0
  if (any(bad_sd)) {
    stop("Baseline sd values must be finite and positive.", call. = FALSE)
  }

  baseline <- list(
    label = label,
    traits = as.character(trait_names),
    mean = means,
    sd = sds,
    source = if (!is.null(pop)) "pop@gv" else "values"
  )
  if (!is.null(covariance)) baseline$cov <- covariance
  if (isTRUE(include_index)) {
    w <- as.numeric(index_weights %||% rep(1, length(means)))
    if (length(w) != length(means)) stop("index_weights length must match number of traits.", call. = FALSE)
    baseline$index_weights <- w
    baseline$index_mean <- as.numeric(sum(means * w))
    baseline$index_sd <- if (!is.null(pop)) {
      idx <- drop(baseline_gv %*% matrix(w, ncol = 1L))
      as.numeric(stats::sd(idx, na.rm = TRUE))
    } else if (!is.null(covariance)) {
      as.numeric(sqrt(drop(t(w) %*% covariance %*% w)))
    } else {
      as.numeric(sqrt(sum((sds * w)^2)))
    }
    baseline$mean <- c(baseline$mean, Index = baseline$index_mean)
    baseline$sd <- c(baseline$sd, Index = baseline$index_sd)
  }
  synthetic_defs <- bp_resolve_synthetic_traits(state, synthetic_traits)
  if (length(synthetic_defs) > 0L) {
    if (is.null(pop)) stop("Synthetic baselines require pop.", call. = FALSE)
    for (def in synthetic_defs) {
      syn_gv <- bp_get_synthetic_gv(
        pop = pop,
        synthetic_trait = def,
        state = state,
        varE = varE,
        n_trials = synthetic_gv_n_trials,
        n_plants_per_trial = synthetic_gv_n_plants,
        seed = seed
      )
      baseline$mean[[def$name]] <- mean(syn_gv, na.rm = TRUE)
      baseline$sd[[def$name]] <- stats::sd(syn_gv, na.rm = TRUE)
    }
  }
  if (is.null(state$sim$trait_baselines)) state$sim$trait_baselines <- list()
  state$sim$trait_baselines[[label]] <- baseline
  state
}

bp_strip_r_comment <- function(line) {
  chars <- strsplit(line, "", fixed = TRUE)[[1L]]
  if (length(chars) == 0L) return(line)
  quote <- NULL
  escaped <- FALSE
  out <- character(0)
  for (ch in chars) {
    if (!is.null(quote)) {
      out <- c(out, ch)
      if (escaped) {
        escaped <- FALSE
      } else if (identical(ch, "\\")) {
        escaped <- TRUE
      } else if (identical(ch, quote)) {
        quote <- NULL
      }
      next
    }
    if (ch %in% c("\"", "'")) {
      quote <- ch
      out <- c(out, ch)
      next
    }
    if (identical(ch, "#")) break
    out <- c(out, ch)
  }
  paste(out, collapse = "")
}

bp_mask_r_strings <- function(line) {
  chars <- strsplit(line, "", fixed = TRUE)[[1L]]
  if (length(chars) == 0L) return(line)
  quote <- NULL
  escaped <- FALSE
  out <- character(0)
  for (ch in chars) {
    if (!is.null(quote)) {
      if (escaped) {
        escaped <- FALSE
      } else if (identical(ch, "\\")) {
        escaped <- TRUE
      } else if (identical(ch, quote)) {
        quote <- NULL
      }
      out <- c(out, " ")
      next
    }
    if (ch %in% c("\"", "'")) {
      quote <- ch
      out <- c(out, " ")
      next
    }
    out <- c(out, ch)
  }
  paste(out, collapse = "")
}

bp_cfg_paths_in_text <- function(line) {
  line <- bp_mask_r_strings(bp_strip_r_comment(line))
  m <- gregexpr("cfg\\$[A-Za-z.][A-Za-z0-9_.]*(?:\\$[A-Za-z.][A-Za-z0-9_.]*)*", line, perl = TRUE)[[1]]
  if (m[[1L]] < 0L) return(character(0))
  out <- regmatches(line, list(m))[[1L]]
  out <- sub("^cfg\\$", "", out)
  unique(gsub("$", ".", out, fixed = TRUE))
}

bp_cfg_is_assignment_call <- function(expr) {
  is.call(expr) &&
    is.symbol(expr[[1L]]) &&
    as.character(expr[[1L]]) %in% c("<-", "=")
}

bp_cfg_assignment_in_line <- function(line) {
  line <- bp_strip_r_comment(line)
  exprs <- tryCatch(parse(text = line, keep.source = FALSE), error = function(e) expression())
  if (length(exprs) == 0L) return(NULL)
  for (expr in as.list(exprs)) {
    if (!bp_cfg_is_assignment_call(expr)) next
    lhs <- bp_cfg_lhs_path(expr[[2L]])
    if (is.na(lhs) || !nzchar(lhs)) next
    rhs_fields <- bp_cfg_paths_in_text(paste(deparse(expr[[3L]], width.cutoff = 500L), collapse = " "))
    return(list(field = lhs, self_referenced = lhs %in% rhs_fields))
  }
  NULL
}

bp_cfg_ref_for_field <- function(field) {
  paste0("cfg$", gsub(".", "$", field, fixed = TRUE))
}

bp_cfg_name_for_skeleton <- function(field) {
  gsub(".", "_", field, fixed = TRUE)
}

bp_cfg_default_for_field <- function(line, field) {
  line <- bp_strip_r_comment(line)
  ref <- bp_cfg_ref_for_field(field)
  pattern <- paste0(gsub("([$])", "\\\\\\1", ref), "\\s*%\\|\\|%\\s*")
  m <- regexpr(pattern, line, perl = TRUE)
  extract_fragment_default <- function() {
    if (m[[1L]] < 0L) return(NA_character_)
    start <- m[[1L]] + attr(m, "match.length")
    txt <- substr(line, start, nchar(line))
    chars <- strsplit(txt, "", fixed = TRUE)[[1L]]
    if (length(chars) == 0L) return(NA_character_)
    depth <- 0L
    quote <- NULL
    escaped <- FALSE
    out <- character(0)
    for (ch in chars) {
      if (!is.null(quote)) {
        out <- c(out, ch)
        if (escaped) {
          escaped <- FALSE
        } else if (identical(ch, "\\")) {
          escaped <- TRUE
        } else if (identical(ch, quote)) {
          quote <- NULL
        }
        next
      }
      if (ch %in% c("\"", "'")) {
        quote <- ch
        out <- c(out, ch)
        next
      }
      if (identical(ch, "#")) break
      if (depth == 0L && ch %in% c(",", ";", "{", "}")) break
      if (depth == 0L && identical(ch, ")")) break
      if (ch %in% c("(", "[", "{")) {
        depth <- depth + 1L
        out <- c(out, ch)
        next
      }
      if (ch %in% c(")", "]", "}")) {
        depth <- max(0L, depth - 1L)
        out <- c(out, ch)
        next
      }
      out <- c(out, ch)
    }
    val <- trimws(paste(out, collapse = ""))
    if (!nzchar(val)) NA_character_ else val
  }

  exprs <- tryCatch(parse(text = line, keep.source = FALSE), error = function(e) expression())
  if (length(exprs) > 0L) {
    find_default <- function(expr) {
      if (!is.call(expr)) return(NA_character_)
      fn <- as.character(expr[[1L]])
      if (identical(fn, "%||%") && length(expr) >= 3L) {
        lhs <- bp_cfg_lhs_path(expr[[2L]])
        if (!is.na(lhs) && identical(lhs, field)) {
          return(paste(deparse(expr[[3L]], width.cutoff = 500L), collapse = " "))
        }
      }
      for (i in seq_along(expr)[-1L]) {
        out <- find_default(expr[[i]])
        if (!is.na(out)) return(out)
      }
      NA_character_
    }

    for (expr in as.list(exprs)) {
      out <- find_default(expr)
      if (!is.na(out) && nzchar(out)) return(out)
    }
  }

  extract_fragment_default()
}

bp_cfg_leaf_info <- function(refs, field, values = NULL) {
  rows <- refs[refs$field == field, , drop = FALSE]
  script_assigned <- "assigned_in_script" %in% names(rows) && any(rows$assigned_in_script)
  def <- NA_character_
  def_row <- NA_integer_
  if (nrow(rows) > 0L) {
    defs <- vapply(seq_len(nrow(rows)), function(j) {
      bp_cfg_default_for_field(rows$text[[j]], field)
    }, character(1))
    has_def <- !is.na(defs) & nzchar(defs)
    if (any(has_def)) {
      def_row <- which(has_def)[[1L]]
      def <- defs[[def_row]]
    }
  }
  imported <- !is.null(values) && field %in% names(values) && nzchar(values[[field]])
  value <- if (isTRUE(imported)) {
    values[[field]]
  } else if (isTRUE(script_assigned)) {
    "NULL"
  } else if (is.na(def)) {
    "XX"
  } else {
    def
  }
  comment <- ""
  if (!isTRUE(imported) && isTRUE(script_assigned)) {
    comment <- "  # set inside scanned script"
  } else if (!isTRUE(imported) && !is.na(def) && nrow(rows) > 0L && !is.na(def_row)) {
    src <- sprintf("%s:%d", basename(rows$file[[def_row]]), as.integer(rows$line[[def_row]]))
    comment <- sprintf("  # default used in %s", src)
  }
  list(value = value, comment = comment)
}

bp_cfg_child_names <- function(fields, prefix = "") {
  if (nzchar(prefix)) {
    fields <- fields[startsWith(fields, paste0(prefix, "."))]
    rel <- sub(paste0("^", prefix, "\\."), "", fields)
  } else {
    rel <- fields
  }
  rel <- rel[!grepl("\\.", rel, fixed = FALSE) | nzchar(sub("^([^.]*)\\..*$", "\\1", rel))]
  unique(sort(sub("\\..*$", "", rel)))
}

bp_cfg_render_node_items <- function(refs, fields, prefix = "", indent = 2L, values = NULL) {
  names <- bp_cfg_child_names(fields, prefix = prefix)
  if (length(names) == 0L) return(character(0))
  out <- character(0)
  pad <- paste(rep(" ", indent), collapse = "")
  for (i in seq_along(names)) {
    nm <- names[[i]]
    full <- if (nzchar(prefix)) paste(prefix, nm, sep = ".") else nm
    has_children <- any(startsWith(fields, paste0(full, ".")))
    comma <- if (i < length(names)) "," else ""
    if (has_children) {
      inner <- bp_cfg_render_node_items(refs, fields, prefix = full, indent = indent + 2L, values = values)
      out <- c(
        out,
        sprintf("%s%s = list(", pad, nm),
        inner,
        sprintf("%s)%s", pad, comma)
      )
    } else {
      leaf <- bp_cfg_leaf_info(refs, full, values = values)
      out <- c(out, sprintf("%s%s = %s%s%s", pad, nm, leaf$value, comma, leaf$comment))
    }
  }
  out
}

bp_cfg_add_trailing_comma <- function(lines) {
  if (length(lines) == 0L) return(lines)
  idx <- utils::tail(which(nzchar(trimws(lines)) & !grepl("^\\s*#", lines)), 1L)
  if (length(idx) == 0L || grepl(",\\s*(#.*)?$", lines[[idx]])) return(lines)
  lines[[idx]] <- sub("(\\s*#.*)$", ",\\1", lines[[idx]])
  if (!grepl(",\\s*(#.*)?$", lines[[idx]])) {
    lines[[idx]] <- paste0(lines[[idx]], ",")
  }
  lines
}

bp_cfg_section_lines <- function(refs, fields, comment = NULL, values = NULL) {
  fields <- sort(unique(fields))
  if (length(fields) == 0L) return(character(0))
  lines <- bp_cfg_render_node_items(refs, fields, prefix = "", indent = 2L, values = values)
  if (!is.null(comment) && nzchar(comment)) {
    lines <- c(sprintf("  # %s", comment), lines)
  }
  lines
}

bp_cfg_finalize_sections <- function(sections) {
  sections <- sections[vapply(sections, length, integer(1)) > 0L]
  if (length(sections) == 0L) return(character(0))
  out <- character(0)
  for (i in seq_along(sections)) {
    lines <- sections[[i]]
    if (i < length(sections)) lines <- bp_cfg_add_trailing_comma(lines)
    out <- c(out, lines)
  }
  out
}

bp_cfg_build_skeleton <- function(refs, fields, values = NULL) {
  if (length(fields) == 0L) return("cfg <- list(\n)")
  fields <- sort(unique(fields))
  file_sets <- lapply(fields, function(field) {
    sort(unique(basename(refs$file[refs$field == field])))
  })
  names(file_sets) <- fields

  sections <- list()
  shared <- fields[vapply(file_sets, length, integer(1)) > 1L]
  if (length(shared) > 0L) {
    shared_keys <- vapply(file_sets[shared], paste, character(1), collapse = ", ")
    for (key in unique(shared_keys)) {
      section_fields <- shared[shared_keys == key]
      sections[[length(sections) + 1L]] <- bp_cfg_section_lines(
        refs,
        section_fields,
        comment = sprintf("Shared across: %s", key),
        values = values
      )
    }
  }

  unique_fields <- setdiff(fields, shared)
  if (length(unique_fields) > 0L) {
    scheme_names <- sort(unique(basename(refs$file)))
    for (scheme in scheme_names) {
      section_fields <- unique_fields[vapply(file_sets[unique_fields], function(x) identical(x, scheme), logical(1))]
      if (length(section_fields) == 0L) next
      comment <- if (length(scheme_names) > 1L) sprintf("%s only", scheme) else NULL
      sections[[length(sections) + 1L]] <- bp_cfg_section_lines(refs, section_fields, comment = comment, values = values)
    }
  }

  lines <- bp_cfg_finalize_sections(sections)
  paste(c("cfg <- list(", lines, ")"), collapse = "\n")
}

bp_cfg_build_missing_snippet <- function(refs, fields) {
  fields <- sort(unique(fields))
  if (length(fields) == 0L) return("# Missing cfg entries: none")
  lines <- bp_cfg_finalize_sections(list(bp_cfg_section_lines(refs, fields, comment = "Missing cfg entries")))
  paste(lines, collapse = "\n")
}

bp_cfg_build_update_block <- function(refs, fields, files, existing_cfg = TRUE) {
  fields <- sort(unique(fields))
  if (length(fields) == 0L) return(character(0))
  scheme_txt <- paste(sort(unique(basename(files))), collapse = ", ")
  if (isTRUE(existing_cfg)) {
    lines <- bp_cfg_finalize_sections(list(bp_cfg_section_lines(refs, fields, comment = NULL)))
    c(
      "# ---- BPS missing cfg entries ----",
      sprintf("# Schemes checked: %s", scheme_txt),
      "# Fill in XX values. Existing cfg values above are not overwritten.",
      "cfg <- utils::modifyList(cfg, list(",
      lines,
      "))",
      "# ---- end BPS missing cfg entries ----"
    )
  } else {
    unlist(strsplit(bp_cfg_build_skeleton(refs, fields), "\n", fixed = TRUE), use.names = FALSE)
  }
}

bp_cfg_fields_from_list_call <- function(expr, prefix = "") {
  if (!is.call(expr) || !identical(as.character(expr[[1L]]), "list")) return(character(0))
  args <- as.list(expr)[-1L]
  nms <- names(args)
  out <- character(0)
  for (i in seq_along(args)) {
    nm <- nms[[i]]
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) next
    field <- if (nzchar(prefix)) paste(prefix, nm, sep = ".") else nm
    if (is.call(args[[i]]) && identical(as.character(args[[i]][[1L]]), "list")) {
      child <- bp_cfg_fields_from_list_call(args[[i]], prefix = field)
      out <- c(out, if (length(child) == 0L) field else child)
    } else {
      out <- c(out, field)
    }
  }
  unique(out)
}

bp_cfg_deparse_expr <- function(expr) {
  paste(deparse(expr, width.cutoff = 500L), collapse = " ")
}

bp_cfg_values_from_list_call <- function(expr, prefix = "") {
  if (!is.call(expr) || !identical(as.character(expr[[1L]]), "list")) return(character(0))
  args <- as.list(expr)[-1L]
  nms <- names(args)
  out <- character(0)
  for (i in seq_along(args)) {
    nm <- nms[[i]]
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) next
    field <- if (nzchar(prefix)) paste(prefix, nm, sep = ".") else nm
    if (is.call(args[[i]]) && identical(as.character(args[[i]][[1L]]), "list")) {
      out <- c(out, bp_cfg_values_from_list_call(args[[i]], prefix = field))
    } else {
      out <- c(out, stats::setNames(bp_cfg_deparse_expr(args[[i]]), field))
    }
  }
  out
}

bp_cfg_values_from_rhs <- function(expr) {
  if (is.call(expr) && identical(as.character(expr[[1L]]), "list")) {
    return(bp_cfg_values_from_list_call(expr))
  }
  if (!is.call(expr)) return(character(0))
  out <- character(0)
  for (i in seq_along(expr)[-1L]) {
    if (is.call(expr[[i]]) && identical(as.character(expr[[i]][[1L]]), "list")) {
      out <- c(out, bp_cfg_values_from_list_call(expr[[i]]))
    }
  }
  out
}

bp_cfg_lhs_path <- function(expr) {
  if (is.symbol(expr)) {
    return(if (identical(as.character(expr), "cfg")) "" else NA_character_)
  }
  if (!is.call(expr) || !identical(as.character(expr[[1L]]), "$")) return(NA_character_)
  parent <- bp_cfg_lhs_path(expr[[2L]])
  if (is.na(parent)) return(NA_character_)
  child <- as.character(expr[[3L]])
  if (!nzchar(parent)) child else paste(parent, child, sep = ".")
}

bp_cfg_fields_from_rhs <- function(expr) {
  if (is.call(expr) && identical(as.character(expr[[1L]]), "list")) {
    return(bp_cfg_fields_from_list_call(expr))
  }
  if (!is.call(expr)) return(character(0))
  out <- character(0)
  for (i in seq_along(expr)[-1L]) {
    if (is.call(expr[[i]]) && identical(as.character(expr[[i]][[1L]]), "list")) {
      out <- c(out, bp_cfg_fields_from_list_call(expr[[i]]))
    }
  }
  unique(out)
}

bp_cfg_fields_from_cfg_file <- function(cfg_file) {
  if (is.null(cfg_file) || !file.exists(cfg_file)) return(character(0))
  exprs <- tryCatch(parse(cfg_file), error = function(e) expression())
  out <- character(0)
  for (expr in as.list(exprs)) {
    if (!bp_cfg_is_assignment_call(expr)) next
    lhs <- bp_cfg_lhs_path(expr[[2L]])
    if (is.na(lhs)) next
    rhs <- expr[[3L]]
    if (!nzchar(lhs)) {
      out <- c(out, bp_cfg_fields_from_rhs(rhs))
    } else {
      out <- c(out, lhs)
    }
  }
  sort(unique(out))
}

bp_cfg_values_from_cfg_file <- function(cfg_file) {
  if (is.null(cfg_file) || !file.exists(cfg_file)) return(character(0))
  exprs <- tryCatch(parse(cfg_file), error = function(e) expression())
  out <- character(0)
  for (expr in as.list(exprs)) {
    if (!bp_cfg_is_assignment_call(expr)) next
    lhs <- bp_cfg_lhs_path(expr[[2L]])
    if (is.na(lhs)) next
    rhs <- expr[[3L]]
    if (!nzchar(lhs)) {
      out <- c(out, bp_cfg_values_from_rhs(rhs))
    } else {
      out <- c(out, stats::setNames(bp_cfg_deparse_expr(rhs), lhs))
    }
  }
  bp_cfg_keep_last_named(out)
}

bp_cfg_file_has_cfg <- function(cfg_file) {
  if (is.null(cfg_file) || !file.exists(cfg_file)) return(FALSE)
  exprs <- tryCatch(parse(cfg_file), error = function(e) expression())
  any(vapply(as.list(exprs), function(expr) {
    bp_cfg_is_assignment_call(expr) &&
      !is.na(bp_cfg_lhs_path(expr[[2L]])) &&
      !nzchar(bp_cfg_lhs_path(expr[[2L]]))
  }, logical(1)))
}

bp_cfg_fields_from_list_object <- function(x, prefix = "") {
  if (!is.list(x)) return(character(0))
  out <- character(0)
  for (nm in names(x)) {
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) next
    field <- if (nzchar(prefix)) paste(prefix, nm, sep = ".") else nm
    if (is.list(x[[nm]]) && !is.data.frame(x[[nm]])) {
      child <- bp_cfg_fields_from_list_object(x[[nm]], prefix = field)
      out <- c(out, if (length(child) == 0L) field else child)
    } else {
      out <- c(out, field)
    }
  }
  sort(unique(out))
}

bp_cfg_value_from_object <- function(x) {
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}

bp_cfg_values_from_list_object <- function(x, prefix = "") {
  if (!is.list(x)) return(character(0))
  out <- character(0)
  for (nm in names(x)) {
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) next
    field <- if (nzchar(prefix)) paste(prefix, nm, sep = ".") else nm
    if (is.list(x[[nm]]) && !is.data.frame(x[[nm]])) {
      out <- c(out, bp_cfg_values_from_list_object(x[[nm]], prefix = field))
    } else {
      out <- c(out, stats::setNames(bp_cfg_value_from_object(x[[nm]]), field))
    }
  }
  out
}

bp_cfg_keep_last_named <- function(x) {
  if (length(x) == 0L) return(character(0))
  nms <- names(x)
  keep <- !is.na(nms) & nzchar(nms)
  x <- x[keep]
  if (length(x) == 0L) return(character(0))
  rev_x <- rev(x)
  rev_x <- rev_x[!duplicated(names(rev_x))]
  rev(rev_x)
}

bp_write_missing_cfg_entries <- function(cfg_file, scan, existing_fields) {
  missing <- setdiff(scan$fields, existing_fields)
  has_cfg <- bp_cfg_file_has_cfg(cfg_file)
  block <- bp_cfg_build_update_block(scan$refs, missing, scan$files, existing_cfg = has_cfg)
  if (length(block) == 0L) return(character(0))
  old <- if (file.exists(cfg_file)) readLines(cfg_file, warn = FALSE) else character(0)
  new_lines <- if (length(old) == 0L) block else c(old, "", block)
  writeLines(new_lines, cfg_file, useBytes = TRUE)
  missing
}

bp_cfg_import_values <- function(cfg = list(), cfg_file = NULL, import_cfg_file = NULL) {
  vals <- character(0)
  if (!is.null(import_cfg_file)) vals <- c(vals, bp_cfg_values_from_cfg_file(import_cfg_file))
  if (!is.null(cfg_file)) vals <- c(vals, bp_cfg_values_from_cfg_file(cfg_file))
  vals <- c(vals, bp_cfg_values_from_list_object(cfg))
  bp_cfg_keep_last_named(vals)
}

#' Scan Scheme Config Requirements
#'
#' Scan R files for direct `cfg$...` references without evaluating them.
#' Printing the returned object emits a copyable `cfg <- list(...)` skeleton.
#' Fields with inline default fallbacks are filled with the detected default and
#' a source-location comment; other fields are filled with `XX`. Nested references
#' such as `cfg$PYT$locs` are rendered as nested list blocks. When multiple
#' files are supplied, fields used by more than one scheme are grouped first
#' with a comment listing the schemes they apply to, followed by per-scheme
#' sections for fields used by only one script. Fields assigned inside a script
#' are included with value `NULL` and treated as optional/script-set fields.
#' An assignment that also reads the same cfg field remains a required input.
#'
#' @param files One or more R script paths.
#'
#' @return Structured object containing config references and a skeleton.
#' @export
bp_scan_cfg_requirements <- function(files) {
  files <- as.character(files)
  rows <- list()
  script_assigned <- character(0)
  assigned_rows <- list()
  for (file in files) {
    lines <- readLines(file, warn = FALSE)
    for (i in seq_along(lines)) {
      assignment <- bp_cfg_assignment_in_line(lines[[i]])
      if (!is.null(assignment) && !isTRUE(assignment$self_referenced)) {
        script_assigned <- c(script_assigned, assignment$field)
        assigned_rows[[length(assigned_rows) + 1L]] <- data.frame(
          file = file,
          line = as.integer(i),
          field = assignment$field,
          text = lines[[i]],
          stringsAsFactors = FALSE
        )
      }
      fields <- bp_cfg_paths_in_text(lines[[i]])
      if (length(fields) == 0L) next
      for (field in fields) {
        rows[[length(rows) + 1L]] <- data.frame(
          file = file,
          line = as.integer(i),
          field = field,
          has_inline_default = !is.na(bp_cfg_default_for_field(lines[[i]], field)),
          text = lines[[i]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  refs <- if (length(rows) == 0L) {
    data.frame(file = character(), line = integer(), field = character(), has_inline_default = logical(), text = character(), stringsAsFactors = FALSE)
  } else {
    do.call(rbind, rows)
  }
  assigned_refs <- if (length(assigned_rows) == 0L) {
    data.frame(file = character(), line = integer(), field = character(), text = character(), stringsAsFactors = FALSE)
  } else {
    do.call(rbind, assigned_rows)
  }
  script_assigned <- sort(unique(script_assigned))
  refs$assigned_in_script <- refs$field %in% script_assigned
  fields <- sort(unique(refs$field))
  required <- sort(unique(refs$field[!refs$has_inline_default & !refs$assigned_in_script]))
  defaulted <- sort(setdiff(fields, required))
  skeleton <- bp_cfg_build_skeleton(refs, fields)
  out <- list(refs = refs, fields = fields, files = files, assigned_in_scripts = script_assigned, assigned_refs = assigned_refs, required = required, defaulted = defaulted, skeleton = skeleton)
  class(out) <- "bp_cfg_requirements"
  out
}

bp_cfg_has_path <- function(cfg, path) {
  parts <- strsplit(path, ".", fixed = TRUE)[[1L]]
  cur <- cfg
  for (part in parts) {
    if (!is.list(cur) || is.null(cur[[part]])) return(FALSE)
    cur <- cur[[part]]
  }
  TRUE
}

#' Check Config Requirements
#'
#' Compare a config list against direct `cfg$...` references in one or more
#' scheme files.
#' If no `cfg_file` is supplied, printing the result emits the same copyable,
#' grouped config skeleton as [bp_scan_cfg_requirements()]. If `cfg_file`
#' exists, printing reports only missing entries as list elements that can be
#' pasted into the existing file, fields present in the file but overwritten by
#' the scanned schemes, and fields present in the file but unused by the scanned
#' schemes. If `cfg_file` does not exist, it is created with the full grouped
#' template.
#'
#' @param cfg Config list.
#' @param files One or more R script paths.
#' @param cfg_file Optional cfg `.R` file to inspect and update.
#' @param update_file If `TRUE` and `cfg_file` is supplied, append missing cfg
#'   entries to `cfg_file` without changing existing lines.
#' @param import_cfg_file Optional existing cfg `.R` file used as a source of
#'   values when creating or rewriting a grouped cfg template.
#' @param rewrite_file If `TRUE`, write a fresh grouped cfg template to
#'   `cfg_file`, preserving values that can be parsed from `cfg_file`,
#'   `import_cfg_file`, or `cfg`.
#'
#' @return Structured check object.
#' @export
bp_check_cfg_requirements <- function(
  cfg = list(),
  files,
  cfg_file = NULL,
  update_file = !is.null(cfg_file),
  import_cfg_file = NULL,
  rewrite_file = FALSE
) {
  scan <- bp_scan_cfg_requirements(files)
  cfg_fields <- bp_cfg_fields_from_list_object(cfg)
  cfg_file_existed <- !is.null(cfg_file) && file.exists(cfg_file)
  file_fields <- bp_cfg_fields_from_cfg_file(cfg_file)
  imported_values <- bp_cfg_import_values(cfg = cfg, cfg_file = cfg_file, import_cfg_file = import_cfg_file)
  grouped_skeleton <- bp_cfg_build_skeleton(scan$refs, scan$fields, values = imported_values)
  all_present_fields <- sort(unique(c(cfg_fields, file_fields)))
  present <- scan$fields %in% all_present_fields
  missing_before_update <- scan$fields[!present]
  missing_skeleton <- if (!is.null(cfg_file) && cfg_file_existed) {
    bp_cfg_build_missing_snippet(scan$refs, missing_before_update)
  } else {
    grouped_skeleton
  }
  added_to_file <- character(0)
  rewritten_file <- FALSE
  if (!is.null(cfg_file) && isTRUE(rewrite_file)) {
    writeLines(unlist(strsplit(grouped_skeleton, "\n", fixed = TRUE), use.names = FALSE), cfg_file, useBytes = TRUE)
    rewritten_file <- TRUE
    file_fields <- bp_cfg_fields_from_cfg_file(cfg_file)
    all_present_fields <- sort(unique(c(cfg_fields, file_fields)))
    present <- scan$fields %in% all_present_fields
  } else if (!is.null(cfg_file) && isTRUE(update_file)) {
    if (isTRUE(cfg_file_existed)) {
      added_to_file <- bp_write_missing_cfg_entries(cfg_file, scan, existing_fields = file_fields)
    } else {
      writeLines(unlist(strsplit(grouped_skeleton, "\n", fixed = TRUE), use.names = FALSE), cfg_file, useBytes = TRUE)
      added_to_file <- missing_before_update
    }
    file_fields <- bp_cfg_fields_from_cfg_file(cfg_file)
    all_present_fields <- sort(unique(c(cfg_fields, file_fields)))
    present <- scan$fields %in% all_present_fields
  }
  missing <- scan$fields[!present]
  overwritten_cfg_fields <- intersect(all_present_fields, scan$assigned_in_scripts)
  if (length(overwritten_cfg_fields) > 0L) {
    assigned_values <- imported_values[overwritten_cfg_fields]
    overwritten_cfg_fields <- overwritten_cfg_fields[
      !is.na(assigned_values) & trimws(unname(assigned_values)) != "NULL"
    ]
  }
  out <- list(
    scan = scan,
    skeleton = grouped_skeleton,
    missing_skeleton = missing_skeleton,
    cfg_file = cfg_file,
    cfg_file_existed = cfg_file_existed,
    present = scan$fields[present],
    missing = missing,
    missing_before_update = missing_before_update,
    missing_required = intersect(missing, scan$required),
    missing_defaulted = intersect(missing, scan$defaulted),
    added_to_file = added_to_file,
    rewritten_file = rewritten_file,
    imported_values = imported_values,
    overwritten_cfg_fields = overwritten_cfg_fields,
    unused_cfg_fields = setdiff(all_present_fields, c(scan$fields, overwritten_cfg_fields))
  )
  class(out) <- "bp_cfg_check"
  out
}

bp_cfg_assignment_locations <- function(scan, field) {
  refs <- scan$assigned_refs
  if (is.null(refs) || nrow(refs) == 0L) return(character(0))
  refs <- refs[refs$field == field, , drop = FALSE]
  if (nrow(refs) == 0L) return(character(0))
  unique(sprintf("%s:%d", basename(refs$file), as.integer(refs$line)))
}

#' @export
print.bp_cfg_requirements <- function(x, ...) {
  cat(x$skeleton, "\n", sep = "")
  invisible(x)
}

#' @export
print.bp_cfg_check <- function(x, ...) {
  if (isTRUE(x$rewritten_file)) {
    cat(x$skeleton, "\n", sep = "")
  } else {
    cat(x$missing_skeleton, "\n", sep = "")
  }
  if (length(x$overwritten_cfg_fields) > 0L) {
    cat("# cfg fields present in cfg_file but overwritten inside scanned scheme scripts:\n")
    for (field in x$overwritten_cfg_fields) {
      locs <- bp_cfg_assignment_locations(x$scan, field)
      suffix <- if (length(locs) > 0L) sprintf("  # overwritten in %s", paste(locs, collapse = ", ")) else ""
      cat(sprintf("# - %s%s\n", field, suffix))
    }
  }
  if (length(x$unused_cfg_fields) > 0L) {
    cat("# Unused cfg fields present but not referenced by scanned scheme scripts:\n")
    for (field in x$unused_cfg_fields) {
      cat(sprintf("# - %s\n", field))
    }
  }
  if (isTRUE(x$rewritten_file) && !is.null(x$cfg_file)) {
    cat(sprintf("# Rewrote grouped cfg template in %s\n", x$cfg_file))
  }
  if (length(x$added_to_file) > 0L && !is.null(x$cfg_file) && isTRUE(x$cfg_file_existed)) {
    cat(sprintf("# Added %d missing cfg field(s) to %s\n", length(x$added_to_file), x$cfg_file))
  }
  invisible(x)
}

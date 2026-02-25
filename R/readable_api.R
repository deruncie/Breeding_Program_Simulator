# Readable functional API: list-based state + explicit stage verbs.

#' Initialize Program State
#'
#' Create an empty state container for readable breeding-program simulations.
#'
#' @param SP AlphaSimR `SimParam` object.
#' @param dt Tick length in years.
#' @param start_time Starting time in years.
#' @param sim Optional list merged into `state$sim`.
#'
#' @return A `bp_state` list.
#' @export
bp_init_state <- function(SP, dt = 0.25, start_time = 0, sim = list()) {
  tick <- as.integer(round(start_time / dt))
  state <- list(
    time = list(tick = tick, dt = as.numeric(dt), t = as.numeric(start_time)),
    sim = utils::modifyList(list(SP = SP, default_chip = 1L), sim),
    pops = list(),
    cohorts = bp_empty_cohorts(),
    phenotype_log = bp_empty_phenotype_log(),
    genotype_log = bp_empty_genotype_log(),
    gs_models = list(),
    cost_log = bp_empty_cost_log(),
    outputs = list(varieties = bp_empty_variety_log()),
    counters = list(cohort = 0L, model = 0L)
  )
  class(state) <- "bp_state"
  state
}

#' Conditional Debug Breakpoint
#'
#' Enter `browser()` when debug flags in `cfg` match the current state.
#'
#' @param state Program state.
#' @param cfg Configuration list with debug flags.
#' @param where Optional label for function/stage.
#' @param year Optional current year.
#' @param tick Optional current tick.
#'
#' @details
#' When a breakpoint is triggered, this helper prints a compact summary of
#' caller-local variables and then evaluates `browser()` in the caller frame so
#' debugging starts in the caller context (e.g. `run_recurrent_gs_tick`).
#'
#' @return Invisibly returns `NULL` or enters browser mode.
#' @export
bp_debug_break <- function(state, cfg, where = NULL, year = NULL, tick = NULL) {
  if (!isTRUE(cfg$debug %||% FALSE)) return(invisible(NULL))

  if (is.null(tick)) tick <- as.integer(state$time$tick)
  if (is.null(year) && !is.null(cfg$ticks_per_year)) {
    year <- as.integer(floor(tick / as.integer(cfg$ticks_per_year)) + 1L)
  }
  if (is.null(where)) {
    where <- as.character(sys.call(-1)[[1]])
  }

  if (!is.null(cfg$debug_after_year) && !is.null(year) && year < as.integer(cfg$debug_after_year)) {
    return(invisible(NULL))
  }
  if (!is.null(cfg$debug_after_tick) && tick < as.integer(cfg$debug_after_tick)) {
    return(invisible(NULL))
  }
  if (!is.null(cfg$debug_where) && !(where %in% as.character(cfg$debug_where))) {
    return(invisible(NULL))
  }

  # Compact caller-local variable summary for fast interactive inspection.
  caller <- parent.frame()
  obj_names <- ls(envir = caller, all.names = FALSE)
  if (length(obj_names) > 0L) {
    n_show <- as.integer(cfg$debug_n_vars %||% 12L)
    n_show <- max(1L, n_show)
    show_names <- head(sort(obj_names), n_show)
    parts <- vapply(show_names, function(nm) {
      x <- get(nm, envir = caller)
      cls <- class(x)[1]
      len <- tryCatch(length(x), error = function(e) NA_integer_)
      if (!is.na(len)) sprintf("%s<%s:%d>", nm, cls, len) else sprintf("%s<%s>", nm, cls)
    }, character(1))
    cat(sprintf(
      "debug: where=%s tick=%d year=%s vars=%s\n",
      where,
      as.integer(tick),
      ifelse(is.null(year), "NA", as.character(year)),
      paste(parts, collapse = ", ")
    ))
  } else {
    cat(sprintf(
      "debug: where=%s tick=%d year=%s vars=<none>\n",
      where,
      as.integer(tick),
      ifelse(is.null(year), "NA", as.character(year))
    ))
  }

  evalq(browser(), envir = caller)
}

# Cohort metadata table.
bp_empty_cohorts <- function() {
  data.frame(
    cohort_id = character(),
    stage = character(),
    stream = character(),
    cycle_id = character(),
    source_cohort_id = character(),
    created_tick = integer(),
    done_tick = integer(),
    available_tick = integer(),
    closed_tick = integer(),
    active = logical(),
    n_ind = integer(),
    genotyped = logical(),
    chips = character(),
    stringsAsFactors = FALSE
  )
}

# Long-form phenotype history table.
bp_empty_phenotype_log <- function() {
  data.frame(
    cohort_id = character(),
    stage = character(),
    individual_id = integer(),
    environment = integer(),
    trait = character(),
    phenotype_value = numeric(),
    p_value = numeric(),
    measured_tick = integer(),
    available_tick = integer(),
    n_loc = integer(),
    reps = integer(),
    stringsAsFactors = FALSE
  )
}

# Long-form genotyping event table.
bp_empty_genotype_log <- function() {
  data.frame(
    cohort_id = character(),
    chip = character(),
    started_tick = integer(),
    done_tick = integer(),
    available_tick = integer(),
    n_ind = integer(),
    stringsAsFactors = FALSE
  )
}

# Append-only cost event table.
bp_empty_cost_log <- function() {
  data.frame(
    tick = integer(),
    stage = character(),
    cohort_id = character(),
    event = character(),
    unit = character(),
    n_units = numeric(),
    unit_cost = numeric(),
    total_cost = numeric(),
    stringsAsFactors = FALSE
  )
}

# Released variety history table.
bp_empty_variety_log <- function() {
  data.frame(
    tick = integer(),
    source_cohort_id = character(),
    variety_id = integer(),
    stringsAsFactors = FALSE
  )
}

# Pure helper: infer number of individuals in a pop-like object.
pop_n_ind <- function(pop) {
  if (methods::is(pop, "Pop") || methods::is(pop, "RawPop")) {
    return(as.integer(pop@nInd))
  }
  out <- tryCatch(nrow(pop), error = function(e) NA_integer_)
  if (length(out) == 1L && !is.na(out)) {
    return(as.integer(out))
  }
  stop("Unable to determine number of individuals for pop", call. = FALSE)
}

# Pure helper: subset a pop-like object by integer indices.
pop_subset <- function(pop, idx) {
  idx <- as.integer(idx)
  if (methods::is(pop, "Pop") || methods::is(pop, "RawPop")) {
    return(pop[idx])
  }
  if (is.data.frame(pop) || is.matrix(pop)) {
    return(pop[idx, , drop = FALSE])
  }
  out <- tryCatch(pop[idx], error = function(e) NULL)
  if (!is.null(out)) {
    return(out)
  }
  stop("Unable to subset pop with [idx]", call. = FALSE)
}

# Pure helper: convert years to integer ticks.
years_to_ticks <- function(dt, years) {
  as.integer(round(as.numeric(years) / as.numeric(dt)))
}

# Pure helper: normalize chip key for logs.
chip_key <- function(chip) {
  as.character(chip)
}

# Resolve chip config to numeric index required by AlphaSimR model functions.
chip_index <- function(state, chip) {
  if (is.numeric(chip)) return(as.integer(chip))
  if (is.character(chip) && !is.null(state$sim$chip_map[[chip]])) {
    return(as.integer(state$sim$chip_map[[chip]]))
  }
  if (is.character(chip) && grepl("^[0-9]+$", chip)) {
    return(as.integer(chip))
  }
  stop("chip must be numeric, numeric string, or present in state$sim$chip_map", call. = FALSE)
}

# Append one cost row.
bp_add_cost <- function(state, stage, cohort_id, event, unit, n_units, unit_cost) {
  row <- data.frame(
    tick = as.integer(state$time$tick),
    stage = as.character(stage),
    cohort_id = as.character(cohort_id),
    event = as.character(event),
    unit = as.character(unit),
    n_units = as.numeric(n_units),
    unit_cost = as.numeric(unit_cost),
    total_cost = as.numeric(n_units) * as.numeric(unit_cost),
    stringsAsFactors = FALSE
  )
  state$cost_log <- rbind(state$cost_log, row)
  state
}

# Add a new cohort row and store its pop.
bp_add_cohort <- function(
  state,
  pop,
  stage,
  stream = "main",
  cycle_id = "cycle_1",
  source_cohort_id = NA_character_,
  duration_years = 0,
  active = TRUE
) {
  state$counters$cohort <- state$counters$cohort + 1L
  cohort_id <- sprintf("cohort_%07d", state$counters$cohort)

  dur_ticks <- years_to_ticks(state$time$dt, duration_years)
  done_tick <- as.integer(state$time$tick + dur_ticks)

  row <- data.frame(
    cohort_id = cohort_id,
    stage = as.character(stage),
    stream = as.character(stream),
    cycle_id = as.character(cycle_id),
    source_cohort_id = as.character(source_cohort_id),
    created_tick = as.integer(state$time$tick),
    done_tick = done_tick,
    available_tick = done_tick,
    closed_tick = NA_integer_,
    active = isTRUE(active),
    n_ind = pop_n_ind(pop),
    genotyped = FALSE,
    chips = "",
    stringsAsFactors = FALSE
  )

  state$cohorts <- rbind(state$cohorts, row)
  state$pops[[cohort_id]] <- pop
  state
}

# Propagate genotype-log availability from one source cohort to a copied subset cohort.
bp_inherit_genotypes_from_source <- function(state, new_cohort_id, source_ids) {
  ids <- unique(as.character(source_ids))
  ids <- ids[!is.na(ids) & nzchar(ids) & ids != "NA"]
  if (length(ids) != 1L) {
    return(state)
  }
  src <- ids[[1L]]
  src_rows <- state$genotype_log[state$genotype_log$cohort_id == src, , drop = FALSE]
  if (nrow(src_rows) == 0L) return(state)

  new_idx <- match(new_cohort_id, state$cohorts$cohort_id)
  if (is.na(new_idx)) return(state)
  new_n <- as.integer(state$cohorts$n_ind[new_idx])

  chips <- unique(as.character(src_rows$chip))
  for (ckey in chips) {
    if (any(state$genotype_log$cohort_id == new_cohort_id & state$genotype_log$chip == ckey)) next

    s2 <- src_rows[src_rows$chip == ckey, , drop = FALSE]
    j <- which.min(s2$available_tick)
    row <- s2[j, , drop = FALSE]
    row$cohort_id <- as.character(new_cohort_id)
    row$n_ind <- new_n
    state$genotype_log <- rbind(state$genotype_log, row)
  }

  bp_refresh_genotyped_flags(state)
}

# Return the most recently created cohort id.
bp_last_cohort_id <- function(state) {
  if (nrow(state$cohorts) == 0L) return(NA_character_)
  as.character(state$cohorts$cohort_id[nrow(state$cohorts)])
}

# Mark a cohort inactive.
bp_close_cohort <- function(state, cohort_id) {
  idx <- match(cohort_id, state$cohorts$cohort_id)
  if (is.na(idx)) return(state)
  state$cohorts$active[idx] <- FALSE
  state$cohorts$closed_tick[idx] <- as.integer(state$time$tick)
  state
}

#' Get Ready Cohorts
#'
#' Return cohorts available for use at `as_of_tick`, optionally filtered by
#' stage/stream and activity status.
#'
#' @param state Program state.
#' @param stage Optional stage filter.
#' @param stream Optional stream filter.
#' @param active_only Whether to return only active cohorts.
#' @param as_of_tick Tick used for availability filtering.
#'
#' @return Cohort metadata `data.frame`.
#' @export
bp_get_ready_cohorts <- function(state, stage = NULL, stream = NULL, active_only = TRUE, as_of_tick = state$time$tick) {
  df <- state$cohorts
  if (!is.null(stage)) {
    df <- df[df$stage %in% stage, , drop = FALSE]
  }
  if (!is.null(stream)) {
    df <- df[df$stream %in% stream, , drop = FALSE]
  }
  if (isTRUE(active_only)) {
    df <- df[df$active, , drop = FALSE]
  }
  df <- df[df$available_tick <= as.integer(as_of_tick), , drop = FALSE]
  df
}

# Select which source cohorts to process from a ready cohort table.
bp_select_source_rows <- function(state, ready, cfg) {
  if (nrow(ready) == 0L) return(ready)

  policy <- as.character(cfg$input_policy %||% "latest_one")
  if (policy == "all_ready") {
    return(ready)
  }

  ready_ord <- ready[order(ready$available_tick, ready$created_tick, decreasing = TRUE), , drop = FALSE]

  if (policy == "latest_one") {
    return(ready_ord[1, , drop = FALSE])
  }

  if (policy == "latest_n") {
    n <- as.integer(cfg$input_n %||% 1L)
    n <- max(1L, min(n, nrow(ready_ord)))
    return(ready_ord[seq_len(n), , drop = FALSE])
  }

  if (policy == "by_cycle") {
    cycles <- as.character(cfg$input_cycle_ids %||% cfg$cycle_id %||% NA_character_)
    cycles <- cycles[!is.na(cycles) & nzchar(cycles)]
    if (length(cycles) == 0L) {
      return(ready_ord[0, , drop = FALSE])
    }
    return(ready_ord[ready_ord$cycle_id %in% cycles, , drop = FALSE])
  }

  if (policy == "custom") {
    fn <- cfg$select_source_cohorts_fn
    if (!is.function(fn)) {
      stop("input_policy='custom' requires cfg$select_source_cohorts_fn", call. = FALSE)
    }
    out <- fn(state, ready_ord, cfg)
    if (is.null(out) || length(out) == 0L) {
      return(ready_ord[0, , drop = FALSE])
    }
    if (is.logical(out)) {
      if (length(out) != nrow(ready_ord)) {
        stop("custom selector logical output must match nrow(ready)", call. = FALSE)
      }
      return(ready_ord[out, , drop = FALSE])
    }
    if (is.numeric(out)) {
      idx <- as.integer(out)
      idx <- idx[idx >= 1L & idx <= nrow(ready_ord)]
      idx <- unique(idx)
      return(ready_ord[idx, , drop = FALSE])
    }
    ids <- as.character(out)
    return(ready_ord[ready_ord$cohort_id %in% ids, , drop = FALSE])
  }

  stop(sprintf("Unknown input_policy: %s", policy), call. = FALSE)
}

# Standardized behavior when a stage call has no eligible source cohorts.
bp_handle_no_ready <- function(cfg, fn_name, stage_label, context = "no available source cohorts") {
  msg <- sprintf("%s: %s for stage '%s'.", fn_name, context, stage_label)
  if (isTRUE(cfg$fail_if_no_ready %||% FALSE)) {
    stop(msg, call. = FALSE)
  }
  if (!isTRUE(cfg$silent %||% FALSE)) {
    cat(msg, "\n")
  }
}

#' Get Ready Pop Bundle
#'
#' Select source cohort(s), retrieve pop(s), and optionally merge them.
#'
#' @param state Program state.
#' @param stage Input stage(s).
#' @param stream Optional stream filter.
#' @param policy Source selection policy.
#' @param combine Whether to merge selected pops into one.
#' @param input_n Number of cohorts for `policy = "latest_n"`.
#' @param cycle_id Optional cycle id.
#' @param input_cycle_ids Optional cycle id vector.
#' @param select_source_cohorts_fn Custom source selector for `policy = "custom"`.
#' @param silent Suppress no-ready messages.
#' @param fail_if_no_ready Error when no cohorts are eligible.
#'
#' @return `NULL` or a source bundle list with `pop`, `source_ids`, and metadata.
#' @export
get_ready_pop <- function(
  state,
  stage,
  stream = NULL,
  policy = "latest_one",
  combine = TRUE,
  input_n = NULL,
  cycle_id = NULL,
  input_cycle_ids = NULL,
  select_source_cohorts_fn = NULL,
  silent = FALSE,
  fail_if_no_ready = FALSE
) {
  cfg <- list(
    input_policy = policy,
    input_n = input_n,
    cycle_id = cycle_id,
    input_cycle_ids = input_cycle_ids,
    select_source_cohorts_fn = select_source_cohorts_fn,
    silent = isTRUE(silent),
    fail_if_no_ready = isTRUE(fail_if_no_ready)
  )
  ready <- bp_get_ready_cohorts(state, stage = stage, stream = stream)
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "get_ready_pop", stage)
    return(NULL)
  }
  selected <- bp_select_source_rows(state, ready, cfg)
  if (nrow(selected) == 0L) {
    bp_handle_no_ready(cfg, "get_ready_pop", stage, context = "source selection policy returned no cohorts")
    return(NULL)
  }

  pops <- lapply(selected$cohort_id, function(cid) state$pops[[cid]])
  pop <- if (isTRUE(combine)) merge_pops(pops) else pops[[1L]]
  source_ids <- as.character(selected$cohort_id)
  cycle_values <- unique(as.character(selected$cycle_id))
  cycle_out <- if (length(cycle_values) == 1L) cycle_values else paste(cycle_values, collapse = ";")

  list(
    pop = pop,
    source_ids = source_ids,
    source_rows = selected,
    stage = as.character(stage),
    stream = if (is.null(stream)) as.character(selected$stream[[1L]]) else as.character(stream),
    cycle_id = cycle_out
  )
}

#' Add Output Stage Cohort
#'
#' Create a new cohort from an output pop, carrying defaults from source metadata.
#'
#' @param state Program state.
#' @param pop Output pop.
#' @param stage Output stage.
#' @param source Optional source bundle from [get_ready_pop()].
#' @param ready_in_years Delay until cohort availability.
#' @param stream Optional output stream override.
#' @param cycle_id Optional output cycle id override.
#' @param active Whether the new cohort is active.
#' @param inherit_genotypes Whether to inherit genotype availability from source.
#'
#' @return Updated program state.
#' @export
put_stage_pop <- function(
  state,
  pop,
  stage,
  source = NULL,
  ready_in_years = 0,
  stream = NULL,
  cycle_id = NULL,
  active = TRUE,
  inherit_genotypes = FALSE
) {
  src_ids <- if (is.null(source)) NA_character_ else paste(source$source_ids, collapse = ";")
  stream_val <- if (!is.null(stream)) stream else if (!is.null(source)) source$stream else "main"
  cycle_val <- if (!is.null(cycle_id)) cycle_id else if (!is.null(source)) source$cycle_id else "cycle_1"

  state <- bp_add_cohort(
    state = state,
    pop = pop,
    stage = stage,
    stream = stream_val,
    cycle_id = cycle_val,
    source_cohort_id = src_ids,
    duration_years = ready_in_years,
    active = active
  )

  if (isTRUE(inherit_genotypes) && !is.null(source) && !is.null(source$source_ids)) {
    state <- bp_inherit_genotypes_from_source(
      state = state,
      new_cohort_id = bp_last_cohort_id(state),
      source_ids = source$source_ids
    )
  }
  state
}

#' Add Stage Cost
#'
#' Append a cost event row to `state$cost_log`.
#'
#' @param state Program state.
#' @param event Event label.
#' @param n_units Number of units.
#' @param unit_cost Cost per unit.
#' @param stage Optional stage label.
#' @param unit Unit label.
#' @param cohort_id Optional cohort id.
#' @param n Alias for `n_units`.
#'
#' @return Updated program state.
#' @export
add_stage_cost <- function(
  state,
  event,
  n_units = NULL,
  unit_cost,
  stage = NULL,
  unit = "unit",
  cohort_id = NULL,
  n = NULL
) {
  if (is.null(n_units)) {
    n_units <- n
  }
  if (is.null(n_units)) {
    stop("add_stage_cost requires n_units (or n)", call. = FALSE)
  }
  cid <- cohort_id %||% bp_last_cohort_id(state)
  stg <- stage
  if (is.null(stg) || !nzchar(stg)) {
    idx <- match(cid, state$cohorts$cohort_id)
    stg <- if (is.na(idx)) "unknown" else state$cohorts$stage[[idx]]
  }
  bp_add_cost(
    state = state,
    stage = stg,
    cohort_id = cid,
    event = event,
    unit = unit,
    n_units = n_units,
    unit_cost = unit_cost
  )
}

#' Log Trial Phenotypes
#'
#' Log trial phenotypes from `pop_trial@pheno` into `state$phenotype_log`.
#'
#' @param state Program state.
#' @param pop_trial Pop with phenotype values in `@pheno`.
#' @param stage Stage label.
#' @param source Optional source bundle.
#' @param traits Trait index/indices.
#' @param n_loc Number of locations represented.
#' @param reps Number of replicates represented.
#' @param environment Environment code (`0` for aggregate).
#' @param p_value Optional environment offset.
#'
#' @return Updated program state.
#' @export
log_trial_pheno <- function(
  state,
  pop_trial,
  stage,
  source = NULL,
  traits = 1L,
  n_loc = 1L,
  reps = 1L,
  environment = 0L,
  p_value = NA_real_
) {
  cid <- bp_last_cohort_id(state)
  if (is.na(cid)) return(state)
  avail_tick <- state$cohorts$available_tick[match(cid, state$cohorts$cohort_id)]
  ph <- pop_trial@pheno
  tr <- as.integer(traits)
  if (is.null(dim(ph))) {
    ph <- matrix(ph, ncol = length(tr))
  }
  bp_record_pheno(
    state = state,
    cohort_id = cid,
    stage = stage,
    individual_id = pop_trial@id,
    traits = tr,
    pheno_matrix = ph,
    available_tick = avail_tick,
    n_loc = as.integer(n_loc),
    reps = as.integer(reps),
    environment = as.integer(environment),
    p_value = as.numeric(p_value)
  )
}

#' Close Source Cohorts
#'
#' Mark all source cohorts in a source bundle as inactive.
#'
#' @param state Program state.
#' @param source Source bundle from [get_ready_pop()].
#'
#' @return Updated program state.
#' @export
close_sources <- function(state, source) {
  if (is.null(source) || is.null(source$source_ids)) return(state)
  for (cid in as.character(source$source_ids)) {
    state <- bp_close_cohort(state, cid)
  }
  state
}

#' Advance Time
#'
#' Advance simulation time by integer ticks.
#'
#' @param state Program state.
#' @param n_ticks Number of ticks to advance.
#'
#' @return Updated program state.
#' @export
bp_advance_time <- function(state, n_ticks = 1L) {
  state$time$tick <- as.integer(state$time$tick + as.integer(n_ticks))
  state$time$t <- as.numeric(state$time$tick * state$time$dt)
  state <- bp_refresh_genotyped_flags(state)
  state
}

# Update cohort-level genotyped/chips fields based on available genotype records.
bp_refresh_genotyped_flags <- function(state) {
  if (nrow(state$cohorts) == 0L) return(state)
  chips_done <- state$genotype_log[state$genotype_log$available_tick <= state$time$tick, , drop = FALSE]

  state$cohorts$genotyped <- FALSE
  state$cohorts$chips <- ""

  if (nrow(chips_done) == 0L) return(state)

  by_cohort <- split(chips_done$chip, chips_done$cohort_id)
  for (cid in names(by_cohort)) {
    idx <- match(cid, state$cohorts$cohort_id)
    if (is.na(idx)) next
    uniq <- unique(as.character(by_cohort[[cid]]))
    state$cohorts$genotyped[idx] <- length(uniq) > 0L
    state$cohorts$chips[idx] <- paste(uniq, collapse = ";")
  }
  state
}

# Internal helper: append phenotype rows for one environment/aggregate view.
bp_record_pheno <- function(
  state,
  cohort_id,
  stage,
  individual_id,
  traits,
  pheno_matrix,
  available_tick,
  n_loc,
  reps,
  environment = NA_integer_,
  p_value = NA_real_
) {
  tr <- as.integer(traits)
  ph <- pheno_matrix
  if (is.null(dim(ph))) {
    ph <- matrix(ph, ncol = 1)
  }

  rows <- do.call(rbind, lapply(seq_along(tr), function(k) {
    data.frame(
      cohort_id = as.character(cohort_id),
      stage = as.character(stage),
      individual_id = as.integer(individual_id),
      environment = as.integer(environment),
      trait = paste0("trait", tr[k]),
      phenotype_value = as.numeric(ph[, k]),
      p_value = as.numeric(p_value),
      measured_tick = as.integer(state$time$tick),
      available_tick = as.integer(available_tick),
      n_loc = as.integer(n_loc),
      reps = as.integer(reps),
      stringsAsFactors = FALSE
    )
  }))

  state$phenotype_log <- rbind(state$phenotype_log, rows)
  state
}

# Resolve persistent per-trial base environment means.
bp_get_env_means <- function(state, cfg) {
  trial_name <- as.character(cfg$trial_name %||% cfg$output_stage %||% "trial")
  n_loc <- as.integer(cfg$n_loc %||% 1L)

  if (is.null(state$sim$env_means)) {
    state$sim$env_means <- list()
  }

  if (!is.null(cfg$env_means)) {
    means <- as.numeric(cfg$env_means)
    if (length(means) != n_loc) {
      stop("cfg$env_means length must equal cfg$n_loc", call. = FALSE)
    }
    state$sim$env_means[[trial_name]] <- means
  } else if (!is.null(state$sim$env_means[[trial_name]])) {
    means <- as.numeric(state$sim$env_means[[trial_name]])
  } else {
    mu <- as.numeric(cfg$env_mean_mu %||% 0)
    sd <- as.numeric(cfg$env_mean_sd %||% 1)
    means <- stats::rnorm(n_loc, mean = mu, sd = sd)
    state$sim$env_means[[trial_name]] <- means
  }

  list(state = state, env_means = means)
}

# Pure helper: merge one or more AlphaSimR pops.
merge_pops <- function(pop_list) {
  if (length(pop_list) == 1L) return(pop_list[[1L]])
  if (all(vapply(pop_list, function(x) methods::is(x, "Pop"), logical(1)))) {
    return(AlphaSimR::mergePops(pop_list))
  }
  if (all(vapply(pop_list, is.data.frame, logical(1)))) {
    return(do.call(rbind, pop_list))
  }
  if (all(vapply(pop_list, is.matrix, logical(1)))) {
    return(do.call(rbind, pop_list))
  }
  stop("merge_pops requires homogeneous pop types (all Pop, all data.frame, or all matrix)", call. = FALSE)
}

#' Run Phenotyping Trial
#'
#' Generic field-trial runner with source selection, phenotype generation,
#' phenotype logging, output cohort creation, and cost logging.
#'
#' @param state Program state.
#' @param cfg Trial configuration list.
#'
#' @section Required `cfg` fields:
#' \describe{
#'   \item{`input_stage`}{Source cohort stage to phenotype.}
#'   \item{`output_stage`}{Stage name for created trial cohort(s).}
#'   \item{`varE`}{Residual variance passed to `AlphaSimR::setPheno`.}
#' }
#'
#' @section Common optional `cfg` fields:
#' \describe{
#'   \item{`traits`}{Trait index/vector. Default `1L`.}
#'   \item{`n_loc`}{Number of locations. Default `1L`.}
#'   \item{`reps`}{Replications per location. Default `1L`.}
#'   \item{`duration_years`}{Delay until output cohort is available.}
#'   \item{`cost_per_plot`}{Cost per plot for logging.}
#'   \item{`select_entries_fn`}{Function `(state, src, pop_in, cfg) -> idx` for entry selection.}
#'   \item{`input_policy`}{Source policy (`latest_one`, `latest_n`, `all_ready`, ...).}
#'   \item{`consume_input`}{Close source cohorts after creating output.}
#'   \item{`inherit_genotypes`}{Propagate genotype availability to copied subset cohorts.}
#' }
#'
#' @section Environment-control options:
#' To model explicit environment means, set one of `use_env_control = TRUE`,
#' `env_means`, `env_mean_mu`, `env_mean_sd`, or `env_year_sd`.
#' Then each environment phenotype is generated with `onlyPheno = TRUE` and
#' aggregated into `pop_trial@pheno`.
#'
#' @section Logging behavior:
#' \describe{
#'   \item{`log_per_environment`}{Log each environment separately (default `TRUE` when environment control is used).}
#'   \item{`log_aggregate`}{Log aggregate line means as environment `0` (default `TRUE`).}
#' }
#'
#' @return Updated program state.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("AlphaSimR", quietly = TRUE)) {
#'   library(AlphaSimR)
#'   h <- quickHaplo(20, 2, 50)
#'   SP <- SimParam$new(h)
#'   SP$addTraitA(10)
#'
#'   state <- bp_init_state(SP = SP, dt = 0.25)
#'   parents <- newPop(h, simParam = SP)
#'   state <- put_stage_pop(state, parents, stage = "F5", ready_in_years = 0)
#'
#'   cfg <- list(
#'     input_stage = "F5",
#'     output_stage = "PYT",
#'     traits = 1L,
#'     n_loc = 4L,
#'     reps = 2L,
#'     varE = 1.0,
#'     duration_years = 0.5,
#'     cost_per_plot = 20,
#'     env_mean_mu = 0,
#'     env_mean_sd = 1,
#'     env_year_sd = 0.2
#'   )
#'   state <- run_phenotype_trial(state, cfg)
#' }
#' }
#' @export
run_phenotype_trial <- function(state, cfg) {
  stage_label <- as.character(cfg$input_stage %||% "unknown")
  ready <- bp_get_ready_cohorts(state, stage = cfg$input_stage, stream = cfg$stream %||% NULL)
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "run_phenotype_trial", stage_label)
    return(state)
  }
  ready <- bp_select_source_rows(state, ready, cfg)
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "run_phenotype_trial", stage_label, context = "source selection policy returned no cohorts")
    return(state)
  }

  traits <- as.integer(cfg$traits %||% 1L)
  n_loc <- as.integer(cfg$n_loc %||% 1L)
  reps <- as.integer(cfg$reps %||% 1L)
  use_env_control <- isTRUE(cfg$use_env_control %||% FALSE) ||
    !is.null(cfg$env_means) || !is.null(cfg$env_year_sd) ||
    !is.null(cfg$env_mean_mu) || !is.null(cfg$env_mean_sd)

  for (i in seq_len(nrow(ready))) {
    src <- ready[i, , drop = FALSE]
    pop_in <- state$pops[[src$cohort_id]]

    idx <- if (is.function(cfg$select_entries_fn)) {
      as.integer(cfg$select_entries_fn(state, src, pop_in, cfg))
    } else {
      seq_len(pop_n_ind(pop_in))
    }
    if (length(idx) == 0L) next

    pop_trial <- pop_subset(pop_in, idx)
    if (isTRUE(use_env_control)) {
      env_out <- bp_get_env_means(state, cfg)
      state <- env_out$state
      env_means <- env_out$env_means
      env_year_sd <- as.numeric(cfg$env_year_sd %||% 0)
      env_year_delta <- stats::rnorm(n_loc, mean = 0, sd = env_year_sd)
      p_env <- env_means + env_year_delta

      env_pheno <- vector("list", n_loc)
      for (e in seq_len(n_loc)) {
        env_pheno[[e]] <- AlphaSimR::setPheno(
          pop_trial,
          varE = cfg$varE,
          reps = reps,
          traits = traits,
          p = p_env[e],
          onlyPheno = TRUE,
          simParam = state$sim$SP
        )
        if (is.null(dim(env_pheno[[e]]))) {
          env_pheno[[e]] <- matrix(env_pheno[[e]], ncol = length(traits))
        }
      }

      pheno_mean <- Reduce("+", env_pheno) / n_loc
      if (length(traits) == 1L) {
        pop_trial@pheno[] <- as.numeric(pheno_mean[, 1])
      } else {
        pop_trial@pheno[,] <- pheno_mean
      }
    } else {
      pop_trial <- AlphaSimR::setPheno(
        pop_trial,
        varE = cfg$varE,
        reps = reps,
        traits = traits,
        simParam = state$sim$SP
      )
      env_pheno <- NULL
      p_env <- rep(NA_real_, n_loc)
    }

    state <- bp_add_cohort(
      state = state,
      pop = pop_trial,
      stage = cfg$output_stage %||% cfg$trial_name,
      stream = cfg$output_stream %||% src$stream,
      cycle_id = src$cycle_id,
      source_cohort_id = src$cohort_id,
      duration_years = cfg$duration_years %||% 1
    )
    new_cohort_id <- bp_last_cohort_id(state)
    if (isTRUE(cfg$inherit_genotypes %||% TRUE)) {
      state <- bp_inherit_genotypes_from_source(
        state = state,
        new_cohort_id = new_cohort_id,
        source_ids = src$cohort_id
      )
    }

    avail_tick <- state$cohorts$available_tick[match(new_cohort_id, state$cohorts$cohort_id)]
    stage_name <- cfg$output_stage %||% cfg$trial_name

    if (isTRUE(use_env_control) && isTRUE(cfg$log_per_environment %||% TRUE)) {
      for (e in seq_len(n_loc)) {
        state <- bp_record_pheno(
          state = state,
          cohort_id = new_cohort_id,
          stage = stage_name,
          individual_id = pop_trial@id,
          traits = traits,
          pheno_matrix = env_pheno[[e]],
          available_tick = avail_tick,
          n_loc = n_loc,
          reps = reps,
          environment = e,
          p_value = p_env[e]
        )
      }
    }

    if (isTRUE(cfg$log_aggregate %||% TRUE)) {
      ph <- pop_trial@pheno
      if (is.null(dim(ph))) ph <- matrix(ph, ncol = length(traits))
      state <- bp_record_pheno(
        state = state,
        cohort_id = new_cohort_id,
        stage = stage_name,
        individual_id = pop_trial@id,
        traits = traits,
        pheno_matrix = ph,
        available_tick = avail_tick,
        n_loc = n_loc,
        reps = reps,
        environment = 0L,
        p_value = mean(p_env)
      )
    }

    n_plots <- pop_n_ind(pop_trial) * n_loc * reps
    state <- bp_add_cost(
      state,
      stage = stage_name,
      cohort_id = new_cohort_id,
      event = "phenotype_trial",
      unit = "plot",
      n_units = n_plots,
      unit_cost = cfg$cost_per_plot %||% 10
    )

    if (isTRUE(cfg$consume_input %||% TRUE)) {
      state <- bp_close_cohort(state, src$cohort_id)
    }
  }

  state
}

#' Run Genotyping
#'
#' Schedule/log genotyping events by cohort and chip, with genotyping costs.
#'
#' @param state Program state.
#' @param cfg Genotyping configuration list.
#'
#' @return Updated program state.
#' @export
run_genotyping <- function(state, cfg) {
  stage_label <- as.character(cfg$input_stage %||% "unknown")
  if (isTRUE(cfg$include_not_ready %||% FALSE)) {
    ready <- state$cohorts
    if (!is.null(cfg$input_stage)) {
      ready <- ready[ready$stage %in% cfg$input_stage, , drop = FALSE]
    }
    if (!is.null(cfg$stream)) {
      ready <- ready[ready$stream %in% cfg$stream, , drop = FALSE]
    }
    ready <- ready[ready$active, , drop = FALSE]
    ready <- ready[ready$created_tick <= as.integer(state$time$tick), , drop = FALSE]
  } else {
    ready <- bp_get_ready_cohorts(state, stage = cfg$input_stage, stream = cfg$stream %||% NULL)
  }
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "run_genotyping", stage_label)
    return(state)
  }
  ready <- bp_select_source_rows(state, ready, cfg)
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "run_genotyping", stage_label, context = "source selection policy returned no cohorts")
    return(state)
  }

  chip <- cfg$chip %||% state$sim$default_chip
  ckey <- chip_key(chip)
  dur_ticks <- years_to_ticks(state$time$dt, cfg$duration_years %||% 1)
  force <- isTRUE(cfg$force %||% FALSE)

  for (i in seq_len(nrow(ready))) {
    src <- ready[i, , drop = FALSE]
    n_ind <- src$n_ind
    already <- state$genotype_log$cohort_id == src$cohort_id & state$genotype_log$chip == ckey
    if (!force && any(already, na.rm = TRUE)) {
      next
    }

    row <- data.frame(
      cohort_id = src$cohort_id,
      chip = ckey,
      started_tick = as.integer(state$time$tick),
      done_tick = as.integer(state$time$tick + dur_ticks),
      available_tick = as.integer(state$time$tick + dur_ticks),
      n_ind = as.integer(n_ind),
      stringsAsFactors = FALSE
    )
    state$genotype_log <- rbind(state$genotype_log, row)

    state <- bp_add_cost(
      state,
      stage = src$stage,
      cohort_id = src$cohort_id,
      event = "genotyping",
      unit = "sample",
      n_units = n_ind,
      unit_cost = cfg$cost_per_sample %||% 15
    )
  }

  state <- bp_refresh_genotyped_flags(state)
  state
}

# Identify eligible training cohorts from a recent window with required genotype data.
bp_get_training_cohorts <- function(state, cfg) {
  ready <- bp_get_ready_cohorts(state, stage = cfg$from_stage, stream = cfg$stream %||% NULL)
  if (nrow(ready) == 0L) return(ready)

  lookback_ticks <- years_to_ticks(state$time$dt, cfg$lookback_years %||% 3)
  min_tick <- as.integer(state$time$tick - lookback_ticks)
  ready <- ready[ready$available_tick >= min_tick, , drop = FALSE]

  ckey <- chip_key(cfg$chip %||% state$sim$default_chip)
  if (nrow(ready) == 0L) return(ready)

  has_chip <- vapply(ready$cohort_id, function(cid) {
    any(state$genotype_log$cohort_id == cid & state$genotype_log$chip == ckey & state$genotype_log$available_tick <= state$time$tick)
  }, logical(1))
  ready <- ready[has_chip, , drop = FALSE]
  if (nrow(ready) == 0L) return(ready)

  training_policy <- as.character(cfg$training_policy %||% "all_ready")
  if (training_policy == "all_ready") {
    return(ready)
  }
  cfg2 <- cfg
  cfg2$input_policy <- training_policy
  bp_select_source_rows(state, ready, cfg2)
}

# Execute a user-provided hook with contextual error messages.
bp_call_user_fn <- function(fn, args, fn_label, stage_label) {
  tryCatch(
    do.call(fn, args),
    error = function(e) {
      stop(sprintf("%s failed for stage '%s': %s", fn_label, stage_label, conditionMessage(e)), call. = FALSE)
    }
  )
}

# Ensure one or more cohorts have available genotype records for the requested chip.
bp_assert_genotyped_cohorts <- function(state, cohort_ids, chip, stage_label = "unknown", context = "prediction") {
  ids <- unique(as.character(cohort_ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0L) {
    stop(sprintf("%s for stage '%s' requires cohort_ids for genotype validation", context, stage_label), call. = FALSE)
  }

  ckey <- chip_key(chip)
  bad <- vapply(ids, function(cid) {
    !any(
      state$genotype_log$cohort_id == cid &
        state$genotype_log$chip == ckey &
        state$genotype_log$available_tick <= state$time$tick
    )
  }, logical(1))

  if (any(bad)) {
    miss <- paste(ids[bad], collapse = ", ")
    stop(
      sprintf(
        "%s for stage '%s' requires genotyping for chip '%s'; missing cohorts: %s",
        context, stage_label, ckey, miss
      ),
      call. = FALSE
    )
  }
}

#' Predict EBV
#'
#' Predict EBVs for a target pop using either a user function or AlphaSimR.
#'
#' @param pop Target pop.
#' @param model_entry Model entry from `state$gs_models`.
#' @param state Program state.
#' @param cfg Prediction configuration list.
#' @param stage_label Stage label used in errors.
#'
#' @section Key `cfg` fields:
#' \describe{
#'   \item{`cohort_ids`}{Cohort id(s) used for genotype-chip validation. Required when `require_genotyped = TRUE`.}
#'   \item{`chip`}{Optional chip override for validation/prediction compatibility.}
#'   \item{`require_genotyped`}{Whether to enforce genotype availability checks. Default `TRUE`.}
#'   \item{`predict_ebv_fn`}{Optional custom predictor function. Signature:
#'   `(target_pop, model_obj, state, cfg, model_entry)`. Must return numeric
#'   vector (`nInd`) or numeric matrix (`nInd x nTraits`).}
#' }
#'
#' @details
#' If `predict_ebv_fn` is not provided, this function calls
#' `AlphaSimR::setEBV(pop, solution = model_entry$model, simParam = state$sim$SP)`.
#' For custom predictors, returned values are written directly into `pop@ebv`
#' (single- or multi-trait).
#'
#' @return Pop with updated `@ebv`.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("AlphaSimR", quietly = TRUE)) {
#'   library(AlphaSimR)
#'   h <- quickHaplo(16, 2, 40)
#'   SP <- SimParam$new(h)
#'   SP$addTraitA(10)
#'   pop <- newPop(h, simParam = SP)
#'
#'   state <- bp_init_state(SP = SP, dt = 1)
#'   state <- put_stage_pop(state, pop, stage = "TEST", ready_in_years = 0)
#'   cid <- tail(state$cohorts$cohort_id, 1)
#'
#'   # Record that this cohort is genotyped on chip 1.
#'   state <- run_genotyping(state, list(input_stage = "TEST", chip = 1L, duration_years = 0))
#'
#'   model_entry <- list(
#'     model = list(dummy = TRUE),
#'     chip = "1",
#'     predict_ebv_fn = function(target_pop, model_obj, state, cfg, model_entry) {
#'       seq_len(target_pop@nInd)
#'     }
#'   )
#'   pop2 <- run_predict_ebv(
#'     pop = pop,
#'     model_entry = model_entry,
#'     state = state,
#'     cfg = list(cohort_ids = cid),
#'     stage_label = "TEST"
#'   )
#'   dim(pop2@ebv)
#' }
#' }
#' @export
run_predict_ebv <- function(pop, model_entry, state, cfg, stage_label = "unknown") {
  require_genotyped <- isTRUE(cfg$require_genotyped %||% TRUE)
  if (require_genotyped) {
    chip_for_pred <- cfg$chip %||% model_entry$chip %||% state$sim$default_chip
    bp_assert_genotyped_cohorts(
      state = state,
      cohort_ids = cfg$cohort_ids %||% NULL,
      chip = chip_for_pred,
      stage_label = stage_label,
      context = "EBV prediction"
    )
  }

  predict_fn <- cfg$predict_ebv_fn %||% model_entry$predict_ebv_fn %||% NULL

  if (is.function(predict_fn)) {
    pred <- bp_call_user_fn(
      predict_fn,
      list(target_pop = pop, model_obj = model_entry$model, state = state, cfg = cfg, model_entry = model_entry),
      fn_label = "predict_ebv_fn",
      stage_label = stage_label
    )
    n <- pop_n_ind(pop)

    if (is.null(dim(pred))) {
      pred_vec <- as.numeric(pred)
      if (length(pred_vec) != n) {
        stop(sprintf("predict_ebv_fn returned vector length %d; expected %d", length(pred_vec), n), call. = FALSE)
      }
      if (anyNA(pred_vec)) {
        stop("predict_ebv_fn returned NA values", call. = FALSE)
      }
      pop@ebv <- matrix(pred_vec, ncol = 1)
    } else {
      pred_mat <- as.matrix(pred)
      if (nrow(pred_mat) != n) {
        stop(sprintf("predict_ebv_fn returned matrix with %d rows; expected %d", nrow(pred_mat), n), call. = FALSE)
      }
      if (ncol(pred_mat) < 1L) {
        stop("predict_ebv_fn returned matrix with zero columns", call. = FALSE)
      }
      if (anyNA(pred_mat)) {
        stop("predict_ebv_fn returned NA values", call. = FALSE)
      }
      storage.mode(pred_mat) <- "double"
      pop@ebv <- pred_mat
    }
    return(pop)
  }

  AlphaSimR::setEBV(pop, solution = model_entry$model, simParam = state$sim$SP)
}

#' Train GP Model
#'
#' Train and store a genomic prediction model from eligible cohorts.
#'
#' @param state Program state.
#' @param cfg Training configuration list.
#'
#' @section Required `cfg` fields:
#' \describe{
#'   \item{`from_stage`}{Stage used as training source cohorts.}
#' }
#'
#' @section Common optional `cfg` fields:
#' \describe{
#'   \item{`chip`}{Genotyping chip key/index. Default `state$sim$default_chip`.}
#'   \item{`trait`}{Trait index for default RRBLUP model.}
#'   \item{`lookback_years`}{Training cohort lookback window.}
#'   \item{`training_policy`}{Subset policy over eligible training cohorts.}
#'   \item{`model_id`}{Stored model id. Auto-generated if omitted.}
#'   \item{`train_model_fn`}{Custom trainer `(train_pop, state, cfg) -> model_object`.}
#'   \item{`predict_ebv_fn`}{Optional custom predictor stored alongside model.}
#' }
#'
#' @details
#' By default, this function trains `AlphaSimR::RRBLUP(...)` on merged training
#' cohorts filtered to the requested chip. The stored model entry includes model
#' object, trained tick, chip key, trait index, and source cohorts.
#'
#' @return Updated program state.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("AlphaSimR", quietly = TRUE)) {
#'   library(AlphaSimR)
#'   h <- quickHaplo(20, 2, 50)
#'   SP <- SimParam$new(h)
#'   SP$addTraitA(10)
#'
#'   state <- bp_init_state(SP = SP, dt = 1)
#'   pop <- newPop(h, simParam = SP)
#'   pop <- setPheno(pop, varE = 1, reps = 1, traits = 1, simParam = SP)
#'   state <- put_stage_pop(state, pop, stage = "PYT", ready_in_years = 0)
#'
#'   state <- run_genotyping(state, list(input_stage = "PYT", chip = 1L, duration_years = 0))
#'
#'   state <- run_train_gp_model(
#'     state,
#'     list(
#'       from_stage = "PYT",
#'       chip = 1L,
#'       trait = 1L,
#'       lookback_years = 3,
#'       model_id = "rrblup_pyt"
#'     )
#'   )
#'   names(state$gs_models)
#' }
#' }
#' @export
run_train_gp_model <- function(state, cfg) {
  train_cohorts <- bp_get_training_cohorts(state, cfg)
  if (nrow(train_cohorts) == 0L) {
    bp_handle_no_ready(cfg, "run_train_gp_model", as.character(cfg$from_stage %||% "unknown"), context = "no eligible training cohorts")
    return(state)
  }

  pops <- lapply(train_cohorts$cohort_id, function(cid) state$pops[[cid]])
  train_pop <- merge_pops(pops)

  chip_raw <- cfg$chip %||% state$sim$default_chip
  ckey <- chip_key(chip_raw)
  cidx <- chip_index(state, chip_raw)
  trait <- as.integer(cfg$trait %||% 1L)
  stage_label <- as.character(cfg$from_stage %||% "unknown")

  if (is.function(cfg$train_model_fn)) {
    model <- bp_call_user_fn(
      cfg$train_model_fn,
      list(train_pop = train_pop, state = state, cfg = cfg),
      fn_label = "train_model_fn",
      stage_label = stage_label
    )
  } else {
    model <- AlphaSimR::RRBLUP(train_pop, traits = trait, use = "pheno", snpChip = cidx, simParam = state$sim$SP)
  }
  if (is.null(model)) {
    stop("run_train_gp_model: model object is NULL", call. = FALSE)
  }

  state$counters$model <- state$counters$model + 1L
  model_id <- as.character(cfg$model_id %||% sprintf("gp_%03d", state$counters$model))

  state$gs_models[[model_id]] <- list(
    model = model,
    predict_ebv_fn = cfg$predict_ebv_fn %||% NULL,
    trained_tick = as.integer(state$time$tick),
    chip = ckey,
    trait = trait,
    source_cohorts = train_cohorts$cohort_id
  )

  state
}

#' Get Latest Model ID
#'
#' Return the id of the most recently trained model in `state$gs_models`.
#'
#' @param state Program state.
#'
#' @return Character scalar model id, or `NULL`.
#' @export
bp_latest_model_id <- function(state) {
  if (length(state$gs_models) == 0L) return(NULL)
  ord <- order(vapply(state$gs_models, function(x) x$trained_tick, numeric(1)), decreasing = TRUE)
  names(state$gs_models)[ord][1]
}

# Core generic wrappers used by tests and user scripts.

#' Run Generic Crossing Stage
#'
#' Generic crossing stage with optional parent selection and custom cross plan.
#'
#' @param state Program state.
#' @param cfg Crossing configuration list.
#'
#' @return Updated program state.
#' @export
run_crossing_stage <- function(state, cfg) {
  ready <- bp_get_ready_cohorts(state, stage = cfg$input_stage, stream = cfg$stream %||% NULL)
  if (nrow(ready) == 0L) return(state)
  ready <- bp_select_source_rows(state, ready, cfg)
  if (nrow(ready) == 0L) return(state)

  for (i in seq_len(nrow(ready))) {
    src <- ready[i, , drop = FALSE]
    parent_pop <- state$pops[[src$cohort_id]]

    if (is.function(cfg$select_parents_fn)) {
      sel <- cfg$select_parents_fn(state, src, parent_pop, cfg)
      if (is.numeric(sel)) {
        parent_pop <- pop_subset(parent_pop, as.integer(sel))
      } else if (!is.null(sel)) {
        parent_pop <- sel
      }
    }

    n_crosses <- as.integer(cfg$n_crosses)
    n_progeny <- as.integer(cfg$n_progeny_per_cross %||% 1L)
    n_parents <- pop_n_ind(parent_pop)
    cross_plan <- if (is.function(cfg$cross_plan_fn)) {
      as.matrix(cfg$cross_plan_fn(state, src, parent_pop, cfg))
    } else {
      matrix(sample.int(n_parents, size = n_crosses * 2L, replace = TRUE), ncol = 2)
    }
    if (ncol(cross_plan) != 2L) {
      stop("cross_plan must have 2 columns", call. = FALSE)
    }

    out_pop <- AlphaSimR::makeCross(parent_pop, crossPlan = cross_plan, nProgeny = n_progeny, simParam = state$sim$SP)

    state <- bp_add_cohort(
      state = state,
      pop = out_pop,
      stage = cfg$output_stage %||% "F1",
      stream = cfg$output_stream %||% src$stream,
      cycle_id = cfg$cycle_id %||% src$cycle_id %||% paste0("cycle_", state$time$tick),
      source_cohort_id = src$cohort_id,
      duration_years = cfg$duration_years %||% 1
    )
    new_cohort_id <- bp_last_cohort_id(state)

    yr_now <- bp_format_year(bp_tick_to_year(state, state$time$tick))
    yr_av <- bp_format_year(bp_tick_to_year(state, state$cohorts$available_tick[match(new_cohort_id, state$cohorts$cohort_id)]))
    src_label <- bp_source_labels(state, src$cohort_id, use = "created")
    cross_desc <- as.character(cfg$cross_strategy %||% cfg$cross_plan_desc %||% "cross_plan_fn/random cross plan")
    state <- bp_log_event(
      state = state,
      fn = "run_crossing_stage",
      event_type = "crossing",
      stage = cfg$output_stage %||% "F1",
      source_ids = src$cohort_id,
      output_id = new_cohort_id,
      event_string = sprintf(
        "Year %s: Created %s by crossing %s with %s (n_crosses=%d, progeny_per_cross=%d). Will be available Year %s.",
        yr_now,
        as.character(cfg$output_stage %||% "F1"),
        src_label,
        cross_desc,
        as.integer(n_crosses),
        as.integer(n_progeny),
        yr_av
      ),
      template_string = sprintf(
        "Cross to %s from %s (n_crosses=%d, n_prog=%d)",
        as.character(cfg$output_stage %||% "F1"),
        as.character(src$stage),
        as.integer(n_crosses),
        as.integer(n_progeny)
      ),
      details = list(
        n_crosses = as.integer(n_crosses),
        n_progeny_per_cross = as.integer(n_progeny),
        cross_strategy = cross_desc
      )
    )

    state <- bp_add_cost(
      state,
      stage = cfg$output_stage %||% "F1",
      cohort_id = new_cohort_id,
      event = cfg$cost_event %||% "crossing",
      unit = cfg$cost_unit %||% "cross",
      n_units = n_crosses,
      unit_cost = cfg$cost_per_cross %||% 1
    )

    if (isTRUE(cfg$consume_input %||% FALSE)) {
      state <- bp_close_cohort(state, src$cohort_id)
    }
  }
  state
}

#' Run Generic Selfing Stage
#'
#' Generic selfing/SSD-style stage with configurable generation depth.
#'
#' @param state Program state.
#' @param cfg Selfing configuration list.
#'
#' @return Updated program state.
#' @export
run_selfing_stage <- function(state, cfg) {
  ready <- bp_get_ready_cohorts(state, stage = cfg$input_stage, stream = cfg$stream %||% NULL)
  if (nrow(ready) == 0L) return(state)
  ready <- bp_select_source_rows(state, ready, cfg)
  if (nrow(ready) == 0L) return(state)

  for (i in seq_len(nrow(ready))) {
    src <- ready[i, , drop = FALSE]
    pop_curr <- state$pops[[src$cohort_id]]

    n_gen <- as.integer(cfg$n_generations %||% 4L)
    if (!is.null(cfg$n_progeny_schedule)) {
      prog <- as.integer(cfg$n_progeny_schedule)
      if (length(prog) != n_gen) {
        stop("length(cfg$n_progeny_schedule) must equal n_generations", call. = FALSE)
      }
    } else {
      first <- as.integer(cfg$n_lines_per_parent %||% cfg$n_lines_per_f1 %||% 5L)
      prog <- c(first, rep(1L, max(0L, n_gen - 1L)))
    }

    for (g in seq_len(n_gen)) {
      pop_curr <- AlphaSimR::self(pop_curr, nProgeny = prog[g], simParam = state$sim$SP)
    }

    state <- bp_add_cohort(
      state = state,
      pop = pop_curr,
      stage = cfg$output_stage %||% "F5",
      stream = cfg$output_stream %||% src$stream,
      cycle_id = src$cycle_id,
      source_cohort_id = src$cohort_id,
      duration_years = cfg$duration_years %||% 2
    )
    new_cohort_id <- bp_last_cohort_id(state)

    yr_now <- bp_format_year(bp_tick_to_year(state, state$time$tick))
    yr_av <- bp_format_year(bp_tick_to_year(state, state$cohorts$available_tick[match(new_cohort_id, state$cohorts$cohort_id)]))
    src_label <- bp_source_labels(state, src$cohort_id, use = "created")
    prog_txt <- paste(prog, collapse = ",")
    state <- bp_log_event(
      state = state,
      fn = "run_selfing_stage",
      event_type = "selfing",
      stage = cfg$output_stage %||% "F5",
      source_ids = src$cohort_id,
      output_id = new_cohort_id,
      event_string = sprintf(
        "Year %s: Advanced %s via selfing/SSD to %s (%d generations, progeny schedule=%s). Will be available Year %s.",
        yr_now,
        src_label,
        as.character(cfg$output_stage %||% "F5"),
        as.integer(n_gen),
        prog_txt,
        yr_av
      ),
      template_string = sprintf(
        "Selfing to %s from %s (%d generations)",
        as.character(cfg$output_stage %||% "F5"),
        as.character(src$stage),
        as.integer(n_gen)
      ),
      details = list(
        n_generations = as.integer(n_gen),
        n_progeny_schedule = prog
      )
    )

    state <- bp_add_cost(
      state,
      stage = cfg$output_stage %||% "F5",
      cohort_id = new_cohort_id,
      event = cfg$cost_event %||% "ssd",
      unit = "line",
      n_units = pop_n_ind(pop_curr),
      unit_cost = cfg$cost_per_line %||% 0.2
    )

    if (isTRUE(cfg$consume_input %||% TRUE)) {
      state <- bp_close_cohort(state, src$cohort_id)
    }
  }

  state
}

# Example wrappers for common scheme-specific operations.
# These are intended as editable templates.

# Generic crossing stage with optional custom parent selection and cross-plan hooks.
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

# Generic selfing stage with configurable generation depth and progeny schedule.
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

# Example streamlined crossing stage.
run_make_f1 <- function(
  state,
  input_stage = "PARENT",
  output_stage = "F1",
  n_crosses = 100L,
  n_progeny_per_cross = 1L,
  ready_in_years = 1,
  cost_per_cross = 1,
  consume_input = FALSE,
  stream = NULL,
  policy = "latest_one",
  cross_plan_fn = NULL,
  silent = FALSE,
  fail_if_no_ready = FALSE
) {
  src <- get_ready_pop(
    state = state,
    stage = input_stage,
    stream = stream,
    policy = policy,
    combine = TRUE,
    silent = isTRUE(silent),
    fail_if_no_ready = isTRUE(fail_if_no_ready)
  )
  if (is.null(src)) return(state)

  parent_pop <- src$pop
  n_crosses <- as.integer(n_crosses)
  n_progeny <- as.integer(n_progeny_per_cross)
  n_parents <- pop_n_ind(parent_pop)
  cross_plan <- if (is.function(cross_plan_fn)) {
    as.matrix(cross_plan_fn(state, src, parent_pop))
  } else {
    matrix(sample.int(n_parents, size = n_crosses * 2L, replace = TRUE), ncol = 2)
  }
  if (ncol(cross_plan) != 2L) {
    stop("cross_plan must have 2 columns", call. = FALSE)
  }

  f1 <- AlphaSimR::makeCross(parent_pop, crossPlan = cross_plan, nProgeny = n_progeny, simParam = state$sim$SP)

  state <- put_stage_pop(
    state = state,
    pop = f1,
    stage = output_stage,
    source = src,
    ready_in_years = ready_in_years,
    stream = stream %||% src$stream
  )

  state <- add_stage_cost(
    state = state,
    stage = output_stage,
    event = "crossing",
    unit = "cross",
    n = n_crosses,
    unit_cost = cost_per_cross
  )

  if (isTRUE(consume_input)) {
    state <- close_sources(state, src)
  }
  state
}

# Example selfing stage: SSD-like advancement to F5 with fixed selfing depth.
run_ssd_to_f5 <- function(state, cfg) {
  cfg2 <- utils::modifyList(
    list(
      output_stage = "F5",
      n_generations = 4L,
      n_progeny_schedule = NULL,
      duration_years = 2,
      input_policy = "all_ready",
      consume_input = TRUE,
      cost_event = "ssd"
    ),
    cfg
  )
  run_selfing_stage(state, cfg2)
}

# Example decision stage: pick variety by PYT phenotype and recycle best EBV lines as parents.
run_select_variety_and_recycle <- function(state, cfg) {
  stage_label <- as.character(cfg$input_stage %||% "unknown")
  ready <- bp_get_ready_cohorts(state, stage = cfg$input_stage, stream = cfg$stream %||% NULL)
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "run_select_variety_and_recycle", stage_label)
    return(state)
  }
  ready <- bp_select_source_rows(state, ready, cfg)
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "run_select_variety_and_recycle", stage_label, context = "source selection policy returned no cohorts")
    return(state)
  }

  model_id <- as.character(cfg$model_id %||% bp_latest_model_id(state))
  if (is.null(model_id) || !nzchar(model_id) || is.null(state$gs_models[[model_id]])) {
    return(state)
  }
  model_entry <- state$gs_models[[model_id]]

  for (i in seq_len(nrow(ready))) {
    src <- ready[i, , drop = FALSE]
    pop <- state$pops[[src$cohort_id]]

    variety <- AlphaSimR::selectInd(pop, nInd = as.integer(cfg$n_variety %||% 1L), use = "pheno", simParam = state$sim$SP)
    var_ids <- as.integer(variety@id)
    state$outputs$varieties <- rbind(
      state$outputs$varieties,
      data.frame(
        tick = as.integer(state$time$tick),
        source_cohort_id = src$cohort_id,
        variety_id = var_ids,
        stringsAsFactors = FALSE
      )
    )

    pop_ebv <- bp_predict_ebv(pop = pop, model_entry = model_entry, state = state, cfg = cfg, stage_label = stage_label)
    n_par <- as.integer(cfg$n_new_parents %||% 50L)
    n_par <- min(n_par, pop_n_ind(pop_ebv))
    ebv_trait <- as.integer(cfg$ebv_trait %||% 1L)
    n_ebv_cols <- if (is.null(dim(pop_ebv@ebv))) 1L else ncol(pop_ebv@ebv)
    if (ebv_trait < 1L || ebv_trait > n_ebv_cols) {
      stop(sprintf("cfg$ebv_trait=%d is out of range for %d EBV column(s)", ebv_trait, n_ebv_cols), call. = FALSE)
    }
    new_parents <- AlphaSimR::selectInd(pop_ebv, nInd = n_par, use = "ebv", trait = ebv_trait, simParam = state$sim$SP)

    if (isTRUE(cfg$replace_all_parents %||% TRUE)) {
      old_par <- bp_get_ready_cohorts(state, stage = cfg$parent_stage %||% "PARENT", active_only = TRUE)
      if (nrow(old_par) > 0L) {
        for (j in seq_len(nrow(old_par))) {
          state <- bp_close_cohort(state, old_par$cohort_id[j])
        }
      }
    }

    state <- bp_add_cohort(
      state = state,
      pop = new_parents,
      stage = cfg$parent_stage %||% "PARENT",
      stream = cfg$parent_stream %||% "main",
      cycle_id = src$cycle_id,
      source_cohort_id = src$cohort_id,
      duration_years = cfg$parent_duration_years %||% 1
    )
    new_cohort_id <- bp_last_cohort_id(state)

    state <- bp_add_cost(
      state,
      stage = cfg$parent_stage %||% "PARENT",
      cohort_id = new_cohort_id,
      event = "parent_recycle",
      unit = "line",
      n_units = pop_n_ind(new_parents),
      unit_cost = cfg$cost_per_parent %||% 0.5
    )

    if (isTRUE(cfg$consume_input %||% TRUE)) {
      state <- bp_close_cohort(state, src$cohort_id)
    }
  }

  state
}

# Example streamlined multi-environment trial stage.
# Designed as a direct template: edit constants in this function body.
run_multienv_trial <- function(state) {
  input_stage <- "F5"
  output_stage <- "PYT"
  policy <- "latest_one"
  stream <- NULL

  traits <- 1L
  n_loc <- 1L
  reps <- 1L
  varE <- 1.0
  ready_in_years <- 1

  env_means <- NULL
  env_mean_mu <- 0
  env_mean_sd <- 1
  env_year_sd <- 0.2

  log_per_environment <- TRUE
  log_aggregate <- TRUE
  cost_per_plot <- 20
  consume_input <- TRUE
  silent <- FALSE
  fail_if_no_ready <- FALSE

  src <- get_ready_pop(
    state = state,
    stage = input_stage,
    stream = stream,
    policy = policy,
    combine = TRUE,
    silent = silent,
    fail_if_no_ready = fail_if_no_ready
  )
  if (is.null(src)) return(state)

  traits <- as.integer(traits)
  n_loc <- as.integer(n_loc)
  reps <- as.integer(reps)
  pop_trial <- src$pop

  cfg_env <- list(
    trial_name = output_stage,
    output_stage = output_stage,
    n_loc = n_loc,
    env_means = env_means,
    env_mean_mu = env_mean_mu,
    env_mean_sd = env_mean_sd
  )
  env_out <- bp_get_env_means(state, cfg_env)
  state <- env_out$state
  p_base <- env_out$env_means
  p_year <- p_base + stats::rnorm(n_loc, mean = 0, sd = as.numeric(env_year_sd))

  env_pheno <- vector("list", n_loc)
  for (e in seq_len(n_loc)) {
    env_pheno[[e]] <- AlphaSimR::setPheno(
      pop_trial,
      varE = varE,
      reps = reps,
      traits = traits,
      p = p_year[e],
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

  state <- put_stage_pop(
    state = state,
    pop = pop_trial,
    stage = output_stage,
    source = src,
    ready_in_years = ready_in_years,
    stream = stream
  )
  new_cohort_id <- bp_last_cohort_id(state)
  avail_tick <- state$cohorts$available_tick[match(new_cohort_id, state$cohorts$cohort_id)]

  if (isTRUE(log_per_environment)) {
    for (e in seq_len(n_loc)) {
      state <- bp_record_pheno(
        state = state,
        cohort_id = new_cohort_id,
        stage = output_stage,
        individual_id = pop_trial@id,
        traits = traits,
        pheno_matrix = env_pheno[[e]],
        available_tick = avail_tick,
        n_loc = n_loc,
        reps = reps,
        environment = e,
        p_value = p_year[e]
      )
    }
  }

  if (isTRUE(log_aggregate)) {
    state <- log_trial_pheno(
      state = state,
      pop_trial = pop_trial,
      stage = output_stage,
      source = src,
      traits = traits,
      n_loc = n_loc,
      reps = reps,
      environment = 0L,
      p_value = mean(p_year)
    )
  }

  n_plots <- pop_n_ind(pop_trial) * n_loc * reps
  state <- add_stage_cost(
    state = state,
    stage = output_stage,
    cohort_id = new_cohort_id,
    event = "phenotype_trial",
    unit = "plot",
    n_units = n_plots,
    unit_cost = cost_per_plot
  )

  if (isTRUE(consume_input)) {
    state <- close_sources(state, src)
  }
  state
}

# Example streamlined single-environment PYT stage.
run_pyt <- function(
  state,
  input_stage = "F5",
  output_stage = "PYT",
  traits = 1L,
  reps = 1L,
  varE = 1.0,
  ready_in_years = 1,
  cost_per_plot = 20,
  consume_input = TRUE,
  silent = FALSE,
  fail_if_no_ready = FALSE
) {
  src <- get_ready_pop(
    state = state,
    stage = input_stage,
    policy = "latest_one",
    combine = TRUE,
    silent = silent,
    fail_if_no_ready = fail_if_no_ready
  )
  if (is.null(src)) return(state)

  pop_trial <- AlphaSimR::setPheno(
    src$pop,
    varE = varE,
    reps = as.integer(reps),
    traits = as.integer(traits),
    simParam = state$sim$SP
  )

  state <- put_stage_pop(
    state = state,
    pop = pop_trial,
    stage = output_stage,
    source = src,
    ready_in_years = ready_in_years
  )
  new_cohort_id <- bp_last_cohort_id(state)
  avail_tick <- state$cohorts$available_tick[match(new_cohort_id, state$cohorts$cohort_id)]

  ph <- pop_trial@pheno
  if (is.null(dim(ph))) ph <- matrix(ph, ncol = length(as.integer(traits)))
  state <- bp_record_pheno(
    state = state,
    cohort_id = new_cohort_id,
    stage = output_stage,
    individual_id = pop_trial@id,
    traits = as.integer(traits),
    pheno_matrix = ph,
    available_tick = avail_tick,
    n_loc = 1L,
    reps = as.integer(reps),
    environment = 0L,
    p_value = NA_real_
  )

  state <- add_stage_cost(
    state = state,
    stage = output_stage,
    cohort_id = new_cohort_id,
    event = "phenotype_trial",
    unit = "plot",
    n_units = pop_n_ind(pop_trial) * as.integer(reps),
    unit_cost = cost_per_plot
  )

  if (isTRUE(consume_input)) {
    state <- close_sources(state, src)
  }
  state
}

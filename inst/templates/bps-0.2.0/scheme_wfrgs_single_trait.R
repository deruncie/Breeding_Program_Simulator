# Setup -------------------------------------------------------------------
# BPS 0.2.0 reference template: whole-family recurrent genomic selection
# (wfRGS)
#
#   Parents --selected founders--> DH --trial/genotype--> wfPop --train GS
#      wfPop --GS select/cross--> candidates --GS select/cross--> candidates
#      candidates --GS select--> EYT (+ current variety check) --> Variety
#
# The run_simulation yearly loop advances the recurrent cycle count. Timing is
# still handled by BPS ticks and cfg$speed_breeding_cycle_years.

library(AlphaSimR)
library(BreedingProgramSimulator)
library(rrBLUP)
library(dplyr)

# Helper Utilities ---------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

cost_value <- function(state, event) {
  if (is.null(state$cost_log) || nrow(state$cost_log) == 0L) return(0)
  if (!"event" %in% names(state$cost_log) || !"total_cost" %in% names(state$cost_log)) return(NA_real_)
  sum(state$cost_log$total_cost[state$cost_log$event == event], na.rm = TRUE)
}

total_cost_value <- function(state) {
  if (is.null(state$cost_log) || nrow(state$cost_log) == 0L || !"total_cost" %in% names(state$cost_log)) return(0)
  sum(state$cost_log$total_cost, na.rm = TRUE)
}


score_with_wf_model <- function(state, src, cfg, stage_label) {
  model_id <- "wfRGS_model"
  if (is.null(state$gs_models[[model_id]])) return(src$pop)
  predict_ebv_pop(
    pop = src$pop,
    model_entry = state$gs_models[[model_id]],
    state = state,
    cfg = list(cohort_ids = src$source_ids, chip = as.integer(cfg$snpChip), require_genotyped = TRUE),
    stage_label = stage_label
  )
}

record_yearly_outputs <- function(results, results_base, state, year, cfg) {
  candidates_src <- select_latest_available(state, stage = "candidates", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  if(!is.null(candidates_src)) {
    candidates <- score_with_wf_model(state, candidates_src, cfg, "candidates")
    state <- bp_update_stage_pop(
      state = state,
      cohort_id = candidates_src$source_ids,
      pop = candidates,
      require_same_ids = TRUE,
      allow_reorder = FALSE
    )
  }

  results_year = list(
    year = year,
    bp_report_stage_metrics(state = state, stage = "starting_parents", stream = "main", metrics = c(
      meanParents = "meanG",
      maxParents = "maxG")),
    bp_report_stage_metrics(state = state, stage = "candidates", stream = "main", metrics = c(
      meanCandidates = "meanG",
      maxCandidates = "maxG",
      varCandidates = "varG",
      accGV_Candidates = "accEBV")),
    bp_report_stage_metrics(state = state, stage = "EYT", stream = "main", metrics = c(
      meanEYT = "meanG")),
    bp_report_stage_metrics(state = state, stage = "Variety", stream = "main", metrics = c(
      VarietyG = "meanG")),
    crossing_cost = cost_value(state, "crossing"),
    line_development_cost = cost_value(state, "line_development"),
    seed_increase_cost = cost_value(state, "seed_increase"),
    phenotype_cost = cost_value(state, "phenotype_trial"),
    genotype_cost = cost_value(state, "genotyping"),
    total_cost = total_cost_value(state)
  )
  results_year[lengths(results_year) == 0] = matrix(NA)

  bind_rows(results,data.frame(results_base,results_year))
}

# Event Verbs --------------------------------------------------------------
select_parents_make_wfPop <- function(state, cfg, year) {
  all_parents <- select_latest_available(state, stage = "Parents", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, all_parents, cfg, event_name = "select_parents_make_wfPop")
  if (chk$skip) return(chk$state)

  best_parents <- selectInd(all_parents$pop, as.integer(cfg$n_wfPop_Parents), use = "pheno", simParam = state$sim$SP)
  n_founders <- nInd(best_parents)

  state <- put_stage_pop(
    state = state,
    pop = best_parents,
    stage = "starting_parents",
    source = all_parents,
    ready_in_years = 0,
    cross_strategy = "full_pairwise",
    inherit_genotypes = FALSE,
    stream = "main"
  )

  if (n_founders == 2L) {
    parent_crosses <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
    f1s <- makeCross(best_parents, parent_crosses, nProgeny = 1L, simParam = state$sim$SP)
    f1_crosses <- matrix(c(1, 1), ncol = 2, byrow = TRUE)
    f2s <- makeCross(f1s, f1_crosses, nProgeny = cfg$n_wfPop, simParam = state$sim$SP)
    cross_desc <- "biparental selfed F1 to DH"
    n_crosses <- nrow(parent_crosses) + 1
  } else if (n_founders == 4L) {
    parent_crosses <- matrix(c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4), ncol = 2, byrow = TRUE)
    f1s <- makeCross(best_parents, parent_crosses, nProgeny = 1L, simParam = state$sim$SP)
    f1_crosses <- matrix(c(1, 6, 2, 5, 3, 4), ncol = 2, byrow = TRUE)
    f2s <- makeCross(f1s, f1_crosses, nProgeny = floor(cfg$n_wfPop / nrow(f1_crosses)), simParam = state$sim$SP)
    cross_desc <- "four-way cross to DH"
    n_crosses <- nrow(parent_crosses) + nrow(f1_crosses)
  } else if (n_founders > 4L) {
    parent_crosses <- t(utils::combn(seq_len(n_founders), 2L))
    f1s <- makeCross(best_parents, parent_crosses, nProgeny = 1L, simParam = state$sim$SP)
    f2s = randCross(f1s,cfg$n_wfPop, simParam = state$sim$SP)
    actual_crosses <- unique(paste(f2s@mother,f2s@father))
    cross_desc <- sprintf("half-diallel F1 random mating among %d parents to DH", n_founders)
    n_crosses <- nrow(parent_crosses) + length(actual_crosses)
  } else {
    stop("cfg$n_wfPop_Parents must be 2, 4, or greater than 4.", call. = FALSE)
  }

  state <- put_stage_pop(
    state,
    f2s,
    stage = "F2",
    source = all_parents,
    ready_in_years = 0,
    cross_strategy = cross_desc,
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_unit = cfg$candidates_cost_per_cross,
    cost_event = "crossing",
    cost_unit = "cross",
    cost_units = n_crosses
  )

  DH_pop <- makeDH(f2s, nDH = 1L, simParam = state$sim$SP)
  state <- put_stage_pop(
    state,
    DH_pop,
    stage = "DH",
    source = all_parents,
    ready_in_years = cfg$DH_duration_years,
    cross_strategy = cross_desc,
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_individual = cfg$DH_line_development_cost_per_line,
    cost_event = "line_development",
    cost_unit = "line"
  )
  state <- add_stage_cost(
    state = state,
    event = "seed_increase",
    n_units = nInd(DH_pop),
    unit_cost = cfg$DH_seed_increase_cost_per_line,
    stage = "DH",
    unit = "line",
    cohort_id = tail(as.character(state$cohorts$cohort_id), 1L)
  )
  run_genotyping(
    state,
    list(
      input_stage = "DH",
      stream = "main",
      input_policy = "latest_one",
      include_not_ready = TRUE,
      chip = as.integer(cfg$snpChip),
      duration_years = 0,
      cost_per_sample = cfg$genotype_cost,
      silent = TRUE
    )
  )
}

run_phenotype_trial_wfPop <- function(state, cfg, year) {
  input_dh <- select_latest_available(state, stage = "DH", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_dh, cfg, event_name = "run_phenotype_trial_wfPop")
  if (chk$skip) return(chk$state)

  run_phenotype_trial(
    state,
    pop = input_dh$pop,
    output_stage = "wfPop",
    input_cohorts = input_dh$source_ids,
    selection_strategy = "None; wfPop training population",
    traits = 1L,
    n_loc = cfg$wfPop_n_locs,
    reps = cfg$wfPop_reps,
    varE = cfg$varE,
    duration_years = cfg$wfPop_duration_years,
    stream = "main",
    cost_per_plot = cfg$wfPop_cost_per_plot,
    use_env_control = TRUE,
    env_means = cfg$wfPop_means,
    env_mean_sd = cfg$wfPop_means_sd,
    env_year_sd = cfg$wfPop_years_sd,
    log_per_environment = TRUE,
    log_aggregate = TRUE,
    silent = TRUE
  )
}

train_wfRGS_model <- function(state, cfg, year) {
  input_wfPop <- select_latest_available(state, stage = "wfPop", stream = "main", n = 2L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_wfPop, cfg, event_name = "train_wfRGS_model")
  if (chk$skip) return(chk$state)

  train_pop <- input_wfPop$pop
  if (isTRUE(cfg$debug_GP)) {
    train_pop <- train_pop[sample(seq_len(nInd(train_pop)), min(as.integer(cfg$debug_GP_n), nInd(train_pop)))]
  }
  model <- AlphaSimR::RRBLUP(train_pop,
                             snpChip = as.integer(cfg$snpChip),
                             simParam = state$sim$SP)
  state$gs_models[["wfRGS_model"]] <- list(
    model = model,
    model_id = "wfRGS_model",
    model_name = "AlphaSimR::RRBLUP",
    predict_ebv_fn = NULL,
    trained_tick = as.integer(state$time$tick),
    chip = as.character(as.integer(cfg$snpChip)),
    trait = 1L,
    source_cohorts = input_wfPop$source_ids
  )
  state$sim$current_model_id <- "wfRGS_model"
  bp_log_event(
    state,
    fn = "train_wfRGS_model",
    event_type = "train_model",
    stage = "wfPop",
    source_ids = input_wfPop$source_ids,
    output_id = "wfRGS_model",
    event_string = sprintf(
      "Year %.2f: Trained wfRGS model from wfPop (n_fit=%d, chip=%s).",
      state$time$t,
      nInd(train_pop),
      as.character(cfg$snpChip)
    ),
    template_string = "Train RRBLUP from latest wfPop",
    details = list(
      model_id = "wfRGS_model",
      n_fit = nInd(train_pop),
      chip = as.character(cfg$snpChip)
    )
  )
}

initialize_wfRGS <- function(state, cfg, year) {
  wfPop <- select_latest_available(state, stage = "wfPop", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, wfPop, cfg, event_name = "initialize_wfRGS")
  if (chk$skip) return(chk$state)

  wfPop_scored <- score_with_wf_model(state, wfPop, cfg, "wfPop")
  n_select <- min(cfg$n_select_RGS_start, nInd(wfPop_scored))
  selected <- selectInd(wfPop_scored, n_select, use = "ebv", simParam = state$sim$SP)
  if (nInd(selected) == 1L) selected <- c(selected, selected)
  if (nInd(selected) %% 2L == 1L) selected <- c(selected, selected[nInd(selected)])
  f1_plan <- matrix(seq_len(nInd(selected)), ncol = 2L, byrow = TRUE)
  f1s <- makeCross(selected, f1_plan, simParam = state$sim$SP)

  candidates <- randCross(f1s,cfg$nRGSpop, simParam = state$sim$SP)
  actual_crosses <- unique(paste(candidates@mother,candidates@father))

  state <- put_stage_pop(
    state,
    candidates,
    stage = "candidates",
    source = wfPop,
    ready_in_years = 2 * cfg$speed_breeding_cycle_years,
    selection_strategy = "Initial wfRGS selected by EBV from wfPop",
    cross_strategy = "cross selected wfPop entries, then random cross",
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_unit = cfg$candidates_cost_per_cross,
    cost_event = "crossing",
    cost_unit = "cross",
    cost_units = length(actual_crosses) + nInd(f1s)
  )
  state <- run_genotyping(
    state,
    list(
      input_stage = "candidates",
      stream = "main",
      input_policy = "latest_one",
      include_not_ready = TRUE,
      chip = as.integer(cfg$snpChip),
      duration_years = 0,
      cost_per_sample = cfg$genotype_cost,
      silent = TRUE
    )
  )
}

run_wfRGS_cycle <- function(state, cfg, year, cycle) {
  candidates <- select_latest_available(state, stage = "candidates", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, candidates, cfg, event_name = "run_wfRGS_cycle")
  if (chk$skip) return(chk$state)

  scored <- score_with_wf_model(state, candidates, cfg, "candidates")
  state <- bp_update_stage_pop(state = state, cohort_id = candidates$source_ids,pop = scored,require_same_ids = TRUE)

  n_select <- min(max(cfg$n_select_RGS_min, cfg$n_select_RGS_start - cycle + 1L), nInd(scored))
  selected <- selectInd(scored, n_select, use = "ebv", simParam = state$sim$SP)
  if (as.integer(cfg$add_random_candidates) > 0L) {
    others <- scored[!scored@id %in% selected@id]
    if (nInd(others) > 0L) selected <- c(selected, sample(others, min(as.integer(cfg$add_random_candidates), nInd(others))))
  }
  candidates_next <- randCross(selected,cfg$nRGSpop, simParam = state$sim$SP)
  actual_crosses <- unique(paste(candidates_next@mother,candidates_next@father))

  state <- put_stage_pop(
    state,
    candidates_next,
    stage = "candidates",
    source = candidates,
    ready_in_years = cfg$speed_breeding_cycle_years,
    selection_strategy = sprintf("wfRGS cycle %d selected by EBV", as.integer(cycle)),
    cross_strategy = "random cross among selected candidates",
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_unit = cfg$candidates_cost_per_cross,
    cost_event = "crossing",
    cost_unit = "cross",
    cost_units = length(actual_crosses)
  )
  state <- run_genotyping(
    state,
    list(
      input_stage = "candidates",
      stream = "main",
      input_policy = "latest_one",
      include_not_ready = TRUE,
      chip = as.integer(cfg$snpChip),
      duration_years = 0,
      cost_per_sample = cfg$genotype_cost,
      silent = TRUE
    )
  )
}

run_EYT <- function(state, cfg, year) {
  candidates <- select_latest_available(state, stage = "candidates", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, candidates, cfg, event_name = "run_EYT")
  if (chk$skip) return(chk$state)

  scored <- score_with_wf_model(state, candidates, cfg, "candidates")
  selected <- selectInd(scored, cfg$nEYT, use = "ebv", simParam = state$sim$SP)
  current_variety <- select_latest_available(state, stage = "Variety", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  eyt_pop <- if (is.null(current_variety)) selected else c(selected, current_variety$pop)

  run_phenotype_trial(
    state,
    pop = eyt_pop,
    output_stage = "EYT",
    input_cohorts = candidates$source_ids,
    selection_strategy = "Top EBV candidates plus current Variety check",
    traits = 1L,
    n_loc = cfg$EYT_n_locs,
    reps = cfg$EYT_reps,
    varE = cfg$varE,
    duration_years = cfg$EYT_duration_years,
    stream = "main",
    cost_per_plot = cfg$EYT_cost_per_plot,
    use_env_control = TRUE,
    env_means = cfg$EYT_means,
    env_mean_sd = cfg$EYT_means_sd,
    env_year_sd = cfg$EYT_years_sd,
    log_per_environment = FALSE,
    log_aggregate = TRUE,
    silent = TRUE
  )
}

release_variety_from_EYT <- function(state, cfg, year) {
  input_eyt <- select_latest_available(state, stage = "EYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_eyt, cfg)
  if (chk$skip) return(chk$state)

  # since current variety was in the EYT, we just take the best from the EYT
  best_variety <- selectInd(input_eyt$pop, 1L, use = "pheno", simParam = state$sim$SP)

  # state <- close_active_stage(state, "Variety")
  state <- put_stage_pop(
    state = state,
    pop = best_variety,
    stage = "Variety",
    source = input_eyt,
    ready_in_years = 0,
    selection_strategy = "Best by phenotype from Variety + latest EYT",
    inherit_genotypes = TRUE,
    stream = "main"
  )
  state
}

# Runners ------------------------------------------------------------------
run_simulation <- function(state, results_base = data.frame(), cfg, results_file = NULL, seed = NULL, start_year = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cfg$ticks_per_year <- as.integer(bp_ticks_per_year(state))
  report_start_year <- start_year %||% state$time$t
  results <- record_yearly_outputs(data.frame(), results_base, state, state$time$t - report_start_year, cfg)

  # Year 1, Month 1
  state <- select_parents_make_wfPop(state, cfg, report_start_year)
  state <- bp_advance_time_years(state = state, years = 1L)

  # Year 2, Month 1
  if(cfg$second_wfPop_trial) state <- select_parents_make_wfPop(state, cfg, state$time$t)
  state <- bp_advance_time_years(state = state, years = 1L)

  # Year 3, Month 1
  state <- run_phenotype_trial_wfPop(state, cfg, state$time$t)
  state <- bp_advance_time_years(state = state, years = 1L)

  # saving wfPop as "candidates" for logging
  wfPop <- select_latest_available(state, stage = "wfPop", stream = "main", n = 1L, combine = TRUE, silent = TRUE)$pop
  # saving wfPop as "candidates"
  state <- put_stage_pop(
    state = state,
    pop = wfPop,
    stage = "candidates",
    source = 'wfPop',
    ready_in_years = 0,
    cross_strategy = "NA",
    inherit_genotypes = TRUE,
    stream = "main"
  )
  results <- record_yearly_outputs(results, results_base, state, state$time$t - report_start_year, cfg)

  if(cfg$second_wfPop_trial) state <- run_phenotype_trial_wfPop(state, cfg, state$time$t)
  state <- train_wfRGS_model(state, cfg, state$time$t)
  state <- initialize_wfRGS(state, cfg, state$time$t)
  state <- bp_advance_time_years(state = state, years = 2*cfg$speed_breeding_cycle_years)
  # Year 3, Month 7
  results <- record_yearly_outputs(results, results_base, state, state$time$t - report_start_year, cfg)


  # Year 3, Month 7 - Year 4, Month 10
  for(cycle in 3:cfg$nCycle_wfRGS) {
    print(cycle)
    if(cycle == 5 & cfg$second_wfPop_trial) state <- train_wfRGS_model(state, cfg, state$time$t)

    state <- run_wfRGS_cycle(state, cfg, state$time$t, cycle = cycle - 2L)
    state <- bp_advance_time_years(state = state, years = cfg$speed_breeding_cycle_years)

    results <- record_yearly_outputs(results, results_base, state, state$time$t - report_start_year, cfg)
  }
  # Year 5, Month 1
  state <- run_EYT(state, cfg, state$time$t)
  state <- bp_advance_time_years(state = state, years = 2L-cfg$speed_breeding_cycle_years)

  # Year 7, Month 7
  state <- release_variety_from_EYT(state, cfg, state$time$t)
  state <- bp_advance_time_years(state = state, years = cfg$speed_breeding_cycle_years)

  # Year 7 End
  results <- record_yearly_outputs(results, results_base, state, state$time$t - report_start_year, cfg)

  if (!is.null(results_file)) write.csv(results, file = results_file, row.names = FALSE)
  invisible(list(state = state, results = results))
}

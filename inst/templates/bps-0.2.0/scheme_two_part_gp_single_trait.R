# Setup -------------------------------------------------------------------
# BPS 0.2.0 reference template: two-part line development with genomic
# recurrent parent improvement
#
# Mainstream line development:
#   F1 --> DH/Headrow --PYT--> AYT --> EYT --> Variety
#
# Parallel parent improvement within each year:
#   Parents/F1 --> PI_CAND --GS select/cross each tick--> PI_SELECTED
#      `-----------------------------------------------------> next F1
#
# The annual loop runs the visible line-development events at tick 1, runs PI
# cycles on cfg$pi_cycle_ticks, builds the next F1 at the final tick, and then
# records progress.

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

score_with_latest_model <- function(state, src, cfg, stage_label, require_genotyped = TRUE) {
  model_id <- bp_latest_model_id(state)
  if (is.null(model_id)) return(src$pop)
  predict_ebv_pop(
    pop = src$pop,
    model_entry = state$gs_models[[model_id]],
    state = state,
    cfg = list(cohort_ids = src$source_ids, chip = as.integer(cfg$snpChip), require_genotyped = require_genotyped),
    stage_label = stage_label
  )
}

record_yearly_outputs <- function(results, results_base, state, year, cfg) {
  pi_src <- select_latest_available(state, stage = "PI_CAND", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  if(!is.null(pi_src)) {
    pi_cand <- score_with_latest_model(state, pi_src, cfg, "PI_CAND", require_genotyped = TRUE)
    state <- bp_update_stage_pop(
      state = state,
      cohort_id = pi_src$source_ids,
      pop = pi_cand,
      require_same_ids = TRUE,
      allow_reorder = FALSE
    )
  }

  results_year = list(
    year = year,
    bp_report_stage_metrics(state = state, stage = "Variety", stream = "main", metrics = c(
      VarietyG = "meanG")),
    bp_report_stage_metrics(state = state, stage = "Headrow", stream = "main", metrics = c(
      meanCandidates = "meanG",
      varCandidates = "varG")),
    bp_report_stage_metrics(state = state, stage = "PI_CAND", stream = "main", metrics = c(
      meanPI_CAND = "meanG",
      varPI_CAND = "varG",
      accGV_PI_CAND = "accEBV",
      accGV_PI_CAND_wf = "wf_accEBV")),
    bp_report_stage_metrics(state = state, stage = "PYT", stream = "main", metrics = c(
      meanPYT = "meanG",
      H2_PYT = "H2")),
    bp_report_stage_metrics(state = state, stage = "AYT", stream = "main", metrics = c(
      meanAYT = "meanG",
      H2_AYT = "H2")),
    bp_report_stage_metrics(state = state, stage = "EYT", stream = "main", metrics = c(
      meanEYT = "meanG",
      H2_EYT = "H2")),
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
update_model_from_TrainPop <- function(state, cfg, year) {
  if (year - cfg$first_year_train <= 0) return(state)
  if (is.finite(cfg$last_year_train) && year > cfg$last_year_train) return(state)

  train_bundle <- select_latest_available(
    state = state,
    stage = cfg$train_stages,
    stream = "main",
    n = length(cfg$train_stages) * max(1L, min(as.integer(cfg$n_prior_years_train), year - cfg$first_year_train)),
    combine = TRUE,
    silent = TRUE
  )
  chk <- bp_skip_if_no_input(state, train_bundle, cfg, event_name = "update_model_from_TrainPop")
  if (chk$skip) return(chk$state)

  # ensure that all input cohorts are genotyped
  state <- run_genotyping(
    state,
    list(
      cohort_ids = train_bundle$source_ids,
      chip = as.integer(cfg$snpChip),
      duration_years = 0,
      cost_per_sample = cfg$genotype_cost,
      input_policy = 'latest_n',
      input_n = 1000
    )
  )

  keep_n <- min(as.integer(cfg$maxTrainPop), pop_n_ind(train_bundle$pop))
  # add a fixed effect for the trial.
  train_bundle$pop@fixEff = factor(rep(1:nrow(train_bundle$source_rows),train_bundle$source_rows$n_ind))
  train_bundle$pop <- train_bundle$pop[sample(1:pop_n_ind(train_bundle$pop),keep_n)]
  train_bundle$n_total <- keep_n

  print(sprintf('Year: %d, Training model with n=%d', as.integer(year), as.integer(train_bundle$n_total)))

  model_id <- sprintf("rrblup_year%03d_tick%03d", as.integer(floor(state$time$t) + 1L), as.integer(state$time$tick))

  debug_model <- isTRUE(cfg$debug_GP)
  model_name <- if (debug_model) "AlphaSimR::RRBLUP_debug_random_subset" else "AlphaSimR::RRBLUP"
  train_pop <- train_bundle$pop
  if (debug_model) {
    debug_n <- min(as.integer(cfg$debug_GP_n), nInd(train_pop))
    train_pop <- train_pop[sample(seq_len(nInd(train_pop)), debug_n)]
  }
  model <- AlphaSimR::RRBLUP(
    train_pop,
    snpChip = as.integer(cfg$snpChip),
    simParam = state$sim$SP
  )
  state$gs_models[[model_id]] <- list(
    model = model,
    model_id = model_id,
    model_name = model_name,
    predict_ebv_fn = NULL,
    trained_tick = as.integer(state$time$tick),
    chip = as.character(as.integer(cfg$snpChip)),
    trait = 1L,
    source_cohorts = train_bundle$source_ids
  )
  state$sim$current_model_id <- model_id
  state <- bp_log_event(
    state = state,
    fn = "update_model_from_TrainPop",
    event_type = "train_model",
    stage = "PYT",
    source_ids = train_bundle$source_ids,
    output_id = model_id,
    event_string = sprintf(
      "Year %s: %s GS model %s from latest available PYT cohorts (n_total=%d, n_fit=%d, chip=%s).",
      state$time$t,
      if (debug_model) "Trained debug random-subset" else "Trained",
      model_id,
      as.integer(train_bundle$n_total),
      as.integer(nInd(train_pop)),
      as.character(as.integer(cfg$snpChip))
    ),
    template_string = if (debug_model) "Train RRBLUP on a small debug_GP subset" else "Train RRBLUP from latest available PYT cohorts",
    details = list(
      model_id = model_id,
      model_name = model_name,
      n_total = as.integer(train_bundle$n_total),
      n_fit = as.integer(nInd(train_pop)),
      chip = as.character(as.integer(cfg$snpChip))
    )
  )
  state
}

release_variety_from_EYT <- function(state, cfg, year) {
  input_eyt <- select_latest_available(state = state, stage = "EYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state = state, input_obj = input_eyt, cfg = cfg, event_name = "release_variety_from_EYT")
  if (chk$skip) return(chk$state)

  best_variety <- selectInd(pop = input_eyt$pop, nInd = 1L, use = "pheno", simParam = state$sim$SP)

  put_stage_pop(
    state = state,
    pop = best_variety,
    stage = "Variety",
    source = input_eyt,
    ready_in_years = 0,
    selection_strategy = "Best by phenotype from latest EYT",
    inherit_genotypes = TRUE,
    stream = "main"
  )
}

select_from_AYT_and_start_EYT <- function(state, cfg, year) {
  input_ayt <- select_latest_available(state = state, stage = "AYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state = state, input_obj = input_ayt, cfg = cfg, event_name = "select_from_AYT_and_start_EYT")
  if (chk$skip) return(chk$state)

  eyt_pop <- selectInd(pop = input_ayt$pop, nInd = min(as.integer(cfg$nEYT), nInd(input_ayt$pop)), use = "pheno", simParam = state$sim$SP)
  current_variety <- select_latest_available(state = state, stage = "Variety", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  if (!is.null(current_variety)) eyt_pop <- c(eyt_pop, current_variety$pop)

  run_phenotype_trial(
    state = state,
    pop = eyt_pop,
    output_stage = "EYT",
    input_cohorts = input_ayt$source_ids,
    selection_strategy = "Top phenotype from latest AYT plus current Variety check",
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

select_from_PYT_and_run_AYT <- function(state, cfg, year) {
  input_pyt <- select_latest_available(state, stage = "PYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_pyt, cfg, event_name = "select_from_PYT_and_run_AYT")
  if (chk$skip) return(chk$state)

  pop_scored <- score_with_latest_model(state, input_pyt, cfg, "PYT", require_genotyped = TRUE)
  if (is.null(pop_scored@ebv) || ncol(pop_scored@ebv) == 0) pop_scored@ebv <- pop_scored@pheno
  ayt_pop <- selectInd(pop_scored, min(as.integer(cfg$nAYT), nInd(pop_scored)), use = "ebv", simParam = state$sim$SP)

  run_phenotype_trial(
    state = state,
    pop = ayt_pop,
    output_stage = "AYT",
    input_cohorts = input_pyt$source_ids,
    selection_strategy = "Top phenotype from latest PYT",
    traits = 1L,
    n_loc = cfg$AYT_n_locs,
    reps = cfg$AYT_reps,
    varE = cfg$varE,
    duration_years = cfg$AYT_duration_years,
    stream = "main",
    cost_per_plot = cfg$AYT_cost_per_plot,
    use_env_control = TRUE,
    env_means = cfg$AYT_means,
    env_mean_sd = cfg$AYT_means_sd,
    env_year_sd = cfg$AYT_years_sd,
    log_per_environment = TRUE,
    log_aggregate = TRUE,
    silent = TRUE
  )
}

advance_Headrow_to_PYT <- function(state, cfg, year) {
  input_dh <- select_latest_available(state, stage = "Headrow", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_dh, cfg, event_name = "advance_Headrow_to_PYT")
  if (chk$skip) return(chk$state)

  state <- run_phenotype_trial(
    state,
    pop = input_dh$pop,
    output_stage = "PYT",
    input_cohorts = input_dh$source_ids,
    selection_strategy = "Advance latest Headrow to PYT",
    traits = 1L,
    n_loc = cfg$PYT_n_locs,
    reps = cfg$PYT_reps,
    varE = cfg$varE,
    duration_years = cfg$PYT_duration_years,
    stream = "main",
    cost_per_plot = cfg$PYT_cost_per_plot,
    use_env_control = TRUE,
    env_means = cfg$PYT_mean,
    env_mean_sd = cfg$PYT_sd,
    env_year_sd = 0,
    log_per_environment = FALSE,
    log_aggregate = TRUE,
    silent = TRUE
  )
  run_genotyping(
    state,
    list(
      input_stage = "PYT",
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

advance_F1_to_Headrow <- function(state, cfg, year) {
  input_f1 <- select_latest_available(state, stage = "F1", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_f1, cfg, event_name = "advance_F1_to_Headrow")
  if (chk$skip) return(chk$state)

  n_ind_per_F1 <- sample(1:nInd(input_f1$pop),cfg$nPYT,replace = TRUE)
  dh_pop <- makeDH(input_f1$pop[n_ind_per_F1],nDH = 1L,simParam = state$sim$SP)

  state <- put_stage_pop(
    state,
    dh_pop,
    stage = "Headrow",
    source = input_f1,
    ready_in_years = cfg$Headrow_duration_years + cfg$DH_duration_years,
    cross_strategy = "makeDH(nDH=1) from latest F1",
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_individual = cfg$DH_line_development_cost_per_line,
    cost_event = "line_development",
    cost_unit = "line"
  )
  add_stage_cost(
    state = state,
    event = "seed_increase",
    n_units = nInd(dh_pop),
    unit_cost = cfg$Headrow_seed_increase_cost_per_line,
    stage = "Headrow",
    unit = "line",
    cohort_id = tail(as.character(state$cohorts$cohort_id), 1L)
  )
}

seed_PI_CAND_from_F1_or_Parents <- function(state, cfg) {
  Parents_pop <- select_latest_available(state, stage = "Parents", stream = "main", n = 1L, combine = TRUE, silent = TRUE, include_not_ready = TRUE)
  if(cfg$start_best_F1s) {
    best_parents = selectInd(Parents_pop$pop,4,use='pheno', simParam = state$sim$SP)
    parent_crosses = matrix(c(1,2,1,3,1,4,2,3,2,4,3,4),ncol=2,byrow=T)
    F1s = makeCross(best_parents,parent_crosses,nProgeny = 1, simParam = state$sim$SP)
  } else {
    F1s <- randCross(
      Parents_pop$pop,
      nCrosses = as.integer(cfg$nCrosses),
      nProgeny = 1,
      simParam = state$sim$SP
    )
  }

  state <- put_stage_pop(
    state,
    F1s,
    stage = "PI_CAND",
    source_ids = Parents_pop$source_ids,
    ready_in_years = cfg$speed_breeding_cycle_years,
    cross_strategy = "Initialize PI candidates",
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_individual = cfg$PI_CAND_crossing_cost_per_individual,
    cost_event = "crossing",
    cost_unit = "individual"
  )
  run_genotyping(
    state,
    list(
      input_stage = "PI_CAND",
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

run_PI_cycle <- function(state, cfg, cycle_index, year) {
  input_pi <- select_latest_available(state, stage = "PI_CAND", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_pi, cfg, event_name = "run_PI_cycle")
  if (chk$skip) return(seed_PI_CAND_from_F1_or_Parents(chk$state, cfg))

  parents_pi <- score_with_latest_model(state, input_pi, cfg, "PI_CAND", require_genotyped = TRUE)
  if (is.null(parents_pi@ebv) || ncol(parents_pi@ebv) == 0) return(state)

  n_cross_pi <- min(as.integer(cfg$nCrosses / cfg$nCyclesPI), floor(nInd(parents_pi) / 2L))
  n_selected_parent_pi <- max(2L, as.integer((2L * n_cross_pi) / cfg$crossesPerParentPI))

  if(cfg$start_best_F1s & nInd(parents_pi) == 6) {
    select_parents = parents_pi
    crossPlan = matrix(c(1,6,2,5,3,4),ncol=2,byrow=T)
    crossPlan = crossPlan[sample(1:nrow(crossPlan),n_selected_parent_pi/2,replace=T),]
  } else {
    select_parents <- selectInd(parents_pi, n_selected_parent_pi, use = 'ebv', simParam = state$sim$SP)
    crossPlan = matrix(select_parents@id, ncol = 2L)
  }
  next_pi <- makeCross(
    select_parents,
    crossPlan = crossPlan,
    nProgeny = as.integer(cfg$nF1PI) * as.integer(cfg$crossesPerParentPI),
    simParam = state$sim$SP
  )
  state <- put_stage_pop(
    state,
    select_parents,
    stage = "PI_SELECTED",
    source = input_pi,
    ready_in_years = 0,
    selection_strategy = sprintf("PI cycle %d selected parents by EBV", as.integer(cycle_index)),
    inherit_genotypes = TRUE,
    stream = "main"
  )
  state <- put_stage_pop(
    state,
    next_pi,
    stage = "PI_CAND",
    source = input_pi,
    ready_in_years = cfg$speed_breeding_cycle_years,
    selection_strategy = sprintf("PI cycle %d selected parents by EBV", as.integer(cycle_index)),
    cross_strategy = "makeCross among selected PI parents",
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_individual = cfg$PI_CAND_crossing_cost_per_individual,
    cost_event = "crossing",
    cost_unit = "cross",
    cost_units = nrow(crossPlan)
  )
  state <- run_genotyping(
    state,
    list(
      input_stage = "PI_CAND",
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

build_F1_from_selected_PI_parents <- function(state, cfg, year) {
  selected_pi <- select_latest_available(state, stage = "PI_SELECTED", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, selected_pi, cfg, event_name = "build_F1_from_selected_PI_parents")
  if (chk$skip) {
    input_parents <- select_latest_available(state, stage = "Parents", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
    chk <- bp_skip_if_no_input(state, input_parents, cfg)
    if (chk$skip) return(chk$state)

    f1_pop <- randCross(
      input_parents$pop,
      nCrosses = as.integer(cfg$nCrosses),
      nProgeny = 1,
      simParam = state$sim$SP
    )

    # state <- close_active_stage(state, "F1")
    state <- put_stage_pop(
      state = state,
      pop = f1_pop,
      stage = "F1",
      source = input_parents,
      ready_in_years = cfg$speed_breeding_cycle_years,
      cross_strategy = sprintf(
        "randCross nCrosses=%d nProgeny=%d",
        as.integer(cfg$nCrosses),1
      ),
      inherit_genotypes = FALSE,
      stream = "main",
      cost_per_unit = cfg$F1_cost_per_cross,
      cost_event = "crossing",
      cost_unit = "individual"
    )
    return(state)
  }

  f1_plan <- matrix(selected_pi$pop@id, ncol = 2L)
  f1_pop <- makeCross(selected_pi$pop, crossPlan = f1_plan, nProgeny = 1, simParam = state$sim$SP)
  put_stage_pop(
    state,
    f1_pop,
    stage = "F1",
    source = selected_pi,
    ready_in_years = cfg$speed_breeding_cycle_years,
    cross_strategy = "F1 from final PI-selected parents",
    inherit_genotypes = FALSE,
    stream = "main"
    # cost_per_unit = cfg$F1_cost_per_cross, # these crosses were already made in PI_cycle
    # cost_event = "crossing",
    # cost_unit = "individual"
  )
}

# Runners ------------------------------------------------------------------
run_simulation <- function(state, results_base = data.frame(), cfg, results_file = NULL, seed = NULL, start_year = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cfg$ticks_per_year <- as.integer(bp_ticks_per_year(state))
  cfg$nCyclesPI <- as.integer(round(1 / as.numeric(cfg$speed_breeding_cycle_years)))
  cfg$first_year_train <- state$time$t + cfg$first_year_train
  cfg$last_year_train <- if (is.finite(cfg$last_year_train)) state$time$t + cfg$last_year_train else Inf

  sim_start_year <- state$time$t
  report_start_year <- start_year %||% sim_start_year

  results <- record_yearly_outputs(data.frame(), results_base, state, state$time$t - report_start_year, cfg)

  for (year in sim_start_year + seq_len(as.integer(cfg$nYearsFuture))) {
    print(year)
    # Month 1
    state <- update_model_from_TrainPop(state, cfg, year)
    state <- release_variety_from_EYT(state, cfg, year)
    state <- select_from_AYT_and_start_EYT(state, cfg, year)
    state <- select_from_PYT_and_run_AYT(state, cfg, year)
    state <- advance_Headrow_to_PYT(state, cfg, year)
    state <- advance_F1_to_Headrow(state, cfg, year)

    # Months 1, 4, 7, 10
    for(pi_cycles_done in seq_len(cfg$nCyclesPI)) {
      state <- run_PI_cycle(state, cfg, cycle_index = pi_cycles_done, year = year)
      if(pi_cycles_done < cfg$nCyclesPI) state <- bp_advance_time_years(state, years = cfg$speed_breeding_cycle_years)
      if(cfg$record_results_per_cycle) results <- record_yearly_outputs(results, results_base, state, state$time$t - report_start_year, cfg)
    }
    state <- build_F1_from_selected_PI_parents(state, cfg, year)

    state <- bp_advance_time_years(state, years = cfg$speed_breeding_cycle_years)
    # Year End
    results <- record_yearly_outputs(results, results_base, state, state$time$t - report_start_year, cfg)
  }

  if (!is.null(results_file)) write.csv(results, file = results_file, row.names = FALSE)
  invisible(list(state = state, results = results))
}

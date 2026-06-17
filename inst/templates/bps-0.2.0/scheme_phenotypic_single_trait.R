# Setup -------------------------------------------------------------------
# BPS 0.2.0 reference template: traditional phenotypic line development
#
#   Parents --cross--> F1 --selfing/SSD--> Headrow --PYT--> AYT --> EYT --> Variety
#      ^                                      |        |       |
#      |                                      |        |       +-- current variety check
#      +--------------- replace oldest parents from AYT ------+
#
# One calendar year runs the full advancement cadence once. Earlier cohorts
# become available according to the ready_in_years/duration_years values below.

library(AlphaSimR)
library(BreedingProgramSimulator)
library(dplyr)

# Helper Utilities ---------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x


calibrate_varE_GxE = function(state,stage,stream = 'main',max_env_range = 5,target_trial_GxE_cor,target_trial_H2, cfg) {
  input_pop <- select_latest_available(state = state, stage = stage, stream = stream, n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state = state, input_obj = input_pop, cfg = cfg, event_name = "calibrate_varE_GxE")
  if (chk$skip) return(list(trial_z = 0, varE = 1))

  pop = input_pop$pop

  if(cfg$initVarGE > 1e-5) {
    env_grid_z = seq(0,max_env_range,length=500)
    gv_by_env = sapply(env_grid_z,function(z) bp_gxe_gv_at_z(pop,z,state))
    cor_gvz_gv = cor(gv_by_env,pop@gv)
    trial_z = env_grid_z[which.min(abs(cor_gvz_gv[,1] - target_trial_GxE_cor))]
    varE = var(bp_gxe_gv_at_z(pop,trial_z,state)) * (1-target_trial_H2)/target_trial_H2
  } else {
    trial_z = 0
    varE = var(bp_gxe_gv_at_z(pop,trial_z,state)) * (1-target_trial_H2)/target_trial_H2
  }
  list(trial_z=trial_z, varE=varE)
}


replace_oldest_parents <- function(current_parents, replacement) {
  if (is.null(current_parents) || nInd(current_parents) == 0L) return(replacement)
  if (is.null(replacement) || nInd(replacement) == 0L) return(current_parents)

  replacement <- replacement[!replacement@id %in% current_parents@id]
  if (nInd(replacement) == 0L) return(current_parents)

  drop_n <- min(nInd(replacement), nInd(current_parents))
  c(current_parents[-seq_len(drop_n)], replacement)
}

cost_value <- function(state, event) {
  if (is.null(state$cost_log) || nrow(state$cost_log) == 0L) return(0)
  if (!"event" %in% names(state$cost_log) || !"total_cost" %in% names(state$cost_log)) return(NA_real_)
  sum(state$cost_log$total_cost[state$cost_log$event == event], na.rm = TRUE)
}

total_cost_value <- function(state) {
  if (is.null(state$cost_log) || nrow(state$cost_log) == 0L || !"total_cost" %in% names(state$cost_log)) return(0)
  sum(state$cost_log$total_cost, na.rm = TRUE)
}

record_yearly_outputs <- function(results, results_base, state, year, cfg) {

  results_year = list(
    year = year,
    bp_report_stage_metrics(state = state, stage = "Variety", stream = "main", metrics =
      c(VarietyG="meanG")),
    bp_report_stage_metrics(state = state, stage = "Headrow", stream = "main", metrics = c(
      meanCandidates="meanG",
      bestCandidates="maxG",
      varCandidates="varG")),
    bp_report_stage_metrics(state = state, stage = "Parents", stream = "main", metrics = c(
      meanParents="meanG",
      varParents="varG")),
    bp_report_stage_metrics(state = state, stage = "PYT", stream = "main", metrics = c(
      meanPYT="meanG",
      H2_PYT="H2")),
    bp_report_stage_metrics(state = state, stage = "AYT", stream = "main", metrics = c(
      meanAYT="meanG",
      H2_AYT="H2")),
    bp_report_stage_metrics(state = state, stage = "EYT", stream = "main", metrics = c(
      meanEYT="meanG",
      H2_EYT="H2")),
    crossing_cost = cost_value(state, "crossing"),
    line_development_cost = cost_value(state, "line_development"),
    seed_increase_cost = cost_value(state, "seed_increase"),
    phenotype_cost = cost_value(state, "phenotype_trial"),
    genotype_cost = cost_value(state, "genotyping"),
    total_cost = total_cost_value(state)
  )
  results_year[lengths(results_year) == 0] = matrix(NA)
  # results_year = data.frame(lapply(results_year,unname))

  bind_rows(results,data.frame(results_base,results_year))
}

# Event Verbs --------------------------------------------------------------
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

record_Parent_phenotypes_from_EYT <- function(state, cfg, year) {
  input_eyt <- select_latest_available(state = state, stage = "EYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state = state, input_obj = input_eyt, cfg = cfg, event_name = "record_Parent_phenotypes_from_EYT")
  if (chk$skip) return(chk$state)

  parents_src <- select_latest_available(state = state, stage = "Parents", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state = state, input_obj = parents_src, cfg = cfg, event_name = "record_Parent_phenotypes_from_EYT")
  if (chk$skip) return(chk$state)

  eyt_ids <- input_eyt$pop@id
  parent_ids <- parents_src$pop@id
  overlap <- eyt_ids %in% parent_ids
  if (!any(overlap)) return(state)

  updated_parents <- parents_src$pop
  eyt_pop <- input_eyt$pop[overlap]
  parent_idx <- match(eyt_pop@id, parent_ids)
  updated_parents@pheno[parent_idx, ] <- eyt_pop@pheno

  bp_update_stage_pop(
    state = state,
    cohort_id = parents_src$source_ids,
    pop = updated_parents,
    require_same_ids = TRUE,
    allow_reorder = FALSE,
    log_event = TRUE
  )
}

select_from_AYT_and_update_parents <- function(state, cfg, year) {
  input_ayt <- select_latest_available(state = state, stage = "AYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state = state, input_obj = input_ayt, cfg = cfg, event_name = "select_from_AYT_and_update_parents")
  if (chk$skip) return(chk$state)

  ayt_selected <- selectInd(pop = input_ayt$pop, nInd = min(as.integer(cfg$nReplaceParents), nInd(input_ayt$pop)), use = "pheno", simParam = state$sim$SP)
  current_parents_src <- select_latest_available(state = state, stage = "Parents", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  current_parents <- if (is.null(current_parents_src)) NULL else current_parents_src$pop
  new_parents <- replace_oldest_parents(current_parents, ayt_selected)

  put_stage_pop(
    state = state,
    pop = new_parents,
    stage = "Parents",
    source = input_ayt,
    ready_in_years = 0,
    selection_strategy = "Replace oldest parents with phenotypic selections from latest AYT",
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
  input_pyt <- select_latest_available(state = state, stage = "PYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state = state, input_obj = input_pyt, cfg = cfg, event_name = "select_from_PYT_and_run_AYT")
  if (chk$skip) return(chk$state)

  ayt_pop <- selectInd(pop = input_pyt$pop, nInd = min(as.integer(cfg$nAYT), nInd(input_pyt$pop)), use = "pheno", simParam = state$sim$SP)

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

advance_headrow_to_PYT <- function(state, cfg, year) {
  input_headrow <- select_latest_available(state = state, stage = "Headrow", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state = state, input_obj = input_headrow, cfg = cfg, event_name = "advance_headrow_to_PYT")
  if (chk$skip) return(chk$state)

  headrow_pop <- input_headrow$pop
  if (isTRUE(cfg$Headrow_trial)) {
    if(sum(is.na(headrow_pop@pheno)) == 0) {
      headrow_pop <- selectInd(pop = headrow_pop, nInd = min(as.integer(cfg$nPYT), nInd(headrow_pop)), use = "pheno", simParam = state$sim$SP)
    }
  }

  if(cfg$Downsample_Headrow < nInd(input_headrow$pop)) {
    state$pops[[input_headrow$source_ids]] = selectInd(pop = input_headrow$pop, nInd = cfg$Downsample_Headrow, use = "rand", simParam = state$sim$SP)
  }

  run_phenotype_trial(
    state = state,
    pop = headrow_pop,
    output_stage = "PYT",
    input_cohorts = input_headrow$source_ids,
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
}

advance_F1_to_headrow <- function(state, cfg, year) {
  input_f1 <- select_latest_available(state = state, stage = "F1", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state = state, input_obj = input_f1, cfg = cfg, event_name = "advance_F1_to_headrow")
  if (chk$skip) return(chk$state)

  # first selfing expands to desired number of individuals
  n_ind <- cfg$nProgenyPerCross * nInd(input_f1$pop)
  n_ind_per_F1 <- sample(1:nInd(input_f1$pop),n_ind,replace = TRUE)
  headrow_pop <- self(input_f1$pop[n_ind_per_F1],nProgeny = 1,simParam = state$sim$SP)

  for (selfing_cycle in seq_len(as.integer(cfg$n_selfing_cycles_headrow))[-1]) {
    headrow_pop <- self(pop = headrow_pop, nProgeny = 1L, simParam = state$sim$SP)
  }

  if (isTRUE(cfg$Headrow_trial)) {
    state <- run_phenotype_trial(
      state = state,
      pop = headrow_pop,
      output_stage = "Headrow",
      input_cohorts = input_f1$source_ids,
      selection_strategy = sprintf(
        "Self latest F1 for %d cycles to Headrow",
        as.integer(cfg$n_selfing_cycles_headrow)
      ),
      traits = 1L,
      n_loc = cfg$Headrow_n_locs,
      reps = cfg$Headrow_reps,
      varE = cfg$varE,
      duration_years = cfg$Headrow_duration_years + cfg$SSD_duration_years,
      stream = "main",
      cost_per_plot = cfg$Headrow_cost_per_plot,
      use_env_control = TRUE,
      env_means = cfg$PYT_mean,
      env_mean_sd = cfg$PYT_sd,
      env_year_sd = 0,
      log_per_environment = FALSE,
      log_aggregate = TRUE,
      silent = TRUE
    )
    state <- add_stage_cost(
      state = state,
      event = "line_development",
      n_units = nInd(headrow_pop),
      unit_cost = cfg$SSD_line_development_cost_per_line,
      stage = "Headrow",
      unit = "line",
      cohort_id = tail(as.character(state$cohorts$cohort_id), 1L)
    )
    state <- add_stage_cost(
      state = state,
      event = "seed_increase",
      n_units = nInd(headrow_pop),
      unit_cost = cfg$Headrow_seed_increase_cost_per_line,
      stage = "Headrow",
      unit = "line",
      cohort_id = tail(as.character(state$cohorts$cohort_id), 1L)
    )
    return(state)
  }

  state <- put_stage_pop(
    state = state,
    pop = headrow_pop,
    stage = "Headrow",
    source = input_f1,
    ready_in_years = cfg$Headrow_duration_years + cfg$SSD_duration_years,
    cross_strategy = sprintf(
      "Self latest F1 for %d cycles to Headrow",
      as.integer(cfg$n_selfing_cycles_headrow)
    ),
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_individual = cfg$SSD_line_development_cost_per_line,
    cost_event = "line_development",
    cost_unit = "line"
  )
  add_stage_cost(
    state = state,
    event = "seed_increase",
    n_units = nInd(headrow_pop),
    unit_cost = cfg$Headrow_seed_increase_cost_per_line,
    stage = "Headrow",
    unit = "line",
    cohort_id = tail(as.character(state$cohorts$cohort_id), 1L)
  )
}

build_F1_from_Parents <- function(state, cfg, year) {
  input_parents <- select_latest_available(state = state, stage = "Parents", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state = state, input_obj = input_parents, cfg = cfg, event_name = "build_F1_from_Parents")
  if (chk$skip) return(chk$state)

  f1_pop <- randCross(input_parents$pop,cfg$nCrosses,nProgeny = 1,simParam = state$sim$SP)

  put_stage_pop(
    state = state,
    pop = f1_pop,
    stage = "F1",
    source = input_parents,
    ready_in_years = cfg$F1_duration_years,
    cross_strategy = sprintf(
      "makeCross random plan nCrosses=%d nProgeny=%d",
      as.integer(cfg$nCrosses),1
    ),
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_unit = cfg$F1_cost_per_cross,
    cost_event = "crossing",
    cost_unit = "individual"
  )
}

# Runners ------------------------------------------------------------------
run_simulation <- function(state, results_base = data.frame(), cfg, results_file = NULL, seed = NULL, start_year = NULL) {
  if (!is.null(seed)) set.seed(seed)

  sim_start_year <- state$time$t
  report_start_year <- start_year %||% sim_start_year

  results <- record_yearly_outputs(data.frame(), results_base, state = state, year = state$time$t - report_start_year, cfg = cfg)

  for (year in sim_start_year + seq_len(as.integer(cfg$nYearsFuture))) {
    print(year)

    # Month 1
    state <- release_variety_from_EYT(state = state, cfg = cfg, year = year)
    state <- record_Parent_phenotypes_from_EYT(state = state, cfg = cfg, year = year)
    state <- select_from_AYT_and_update_parents(state = state, cfg = cfg, year = year)
    state <- select_from_AYT_and_start_EYT(state = state, cfg = cfg, year = year)
    state <- select_from_PYT_and_run_AYT(state = state, cfg = cfg, year = year)
    state <- advance_headrow_to_PYT(state = state, cfg = cfg, year = year)
    state <- advance_F1_to_headrow(state = state, cfg = cfg, year = year)
    state <- build_F1_from_Parents(state = state, cfg = cfg, year = year)

    state <- bp_advance_time_years(state = state, years = 1)

    # Year End
    results <- record_yearly_outputs(results = results, results_base = results_base, state = state, year = state$time$t - report_start_year, cfg = cfg)
  }

  if (!is.null(results_file)) {
    write.csv(results, file = results_file, row.names = FALSE)
  }

  invisible(list(state = state, results = results))
}

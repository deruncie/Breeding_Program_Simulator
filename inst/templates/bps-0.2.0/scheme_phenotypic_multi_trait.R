# Setup -------------------------------------------------------------------
# BPS 0.2.0 reference template: multi-trait traditional phenotypic line
# development
#
#   Parents --cross--> F1 --DH--> PYT(trait subset) --> AYT(all/selected traits)
#      ^                                                     |
#      +------ replace parents by phenotypic selection index-+
#   AYT --> EYT (+ current variety check) --> Variety
#
# Selection uses cfg$traitWeights on available stage phenotypes/EBVs. Trial
# trait sets, residual covariance, costs, and stage sizes are all flat cfg$.

library(AlphaSimR)
library(BreedingProgramSimulator)
library(dplyr)

# Helper Utilities ---------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

replace_oldest_parents <- function(current, selected) {
  if (is.null(current) || nInd(current) == 0L) return(selected)
  n_replace <- min(nInd(selected), nInd(current))
  if (n_replace == nInd(current)) return(selected[seq_len(n_replace)])
  c(current[seq_len(nInd(current) - n_replace)], selected[seq_len(n_replace)])
}

cost_value <- function(state, event) {
  if (is.null(state$cost_log) || !nrow(state$cost_log)) return(0)
  vals <- state$cost_log$total_cost[state$cost_log$event == event]
  if (!length(vals)) return(0)
  sum(vals, na.rm = TRUE)
}

total_cost_value <- function(state) {
  if (is.null(state$cost_log) || !nrow(state$cost_log)) return(0)
  sum(state$cost_log$total_cost, na.rm = TRUE)
}

record_yearly_outputs <- function(results, results_base, state, year, cfg) {
  all_traits <- seq_len(cfg$nTraits)
  results_year <- c(
    year = year,
    bp_report_stage_metrics(state = state, stage = "Variety", stream = "main", metrics = c(
      VarietyG = "meanG"),
      traits = all_traits,
      append_trait = "always"),
    bp_report_stage_metrics(state = state, stage = "Parents", stream = "main", metrics = c(
      meanParents = "meanG",
      varParents = "varG"),
      traits = all_traits,
      append_trait = "always"),
    bp_report_stage_metrics(state = state, stage = "PYT", stream = "main", metrics = c(
      meanPYT = "meanG",
      H2_PYT = "H2"),
      traits = cfg$PYT_traits,
      append_trait = "always"),
    bp_report_stage_metrics(state = state, stage = "AYT", stream = "main", metrics = c(
      meanAYT = "meanG",
      H2_AYT = "H2"),
      traits = cfg$AYT_traits,
      append_trait = "always"),
    bp_report_stage_metrics(state = state, stage = "EYT", stream = "main", metrics = c(
      meanEYT = "meanG",
      H2_EYT = "H2"),
      traits = cfg$EYT_traits,
      append_trait = "always")
  )

  results_year <- c(
    results_year,
    bp_report_stage_metrics(state = state, stage = "Variety", stream = "main", metrics = c(
      VarietyG = "meanG"),
      synthetic_trait = cfg$synthetic_trait, cfg = cfg,
      append_trait = "always"),
    bp_report_stage_metrics(state = state, stage = "Parents", stream = "main", metrics = c(
      meanParents = "meanG",
      varParents = "varG"),
      synthetic_trait = cfg$synthetic_trait, cfg = cfg,
      append_trait = "always"),
    bp_report_stage_metrics(state = state, stage = "PYT", stream = "main", metrics = c(
      meanPYT = "meanG",
      H2_PYT = "H2"),
      synthetic_trait = cfg$synthetic_trait, cfg = cfg,
      append_trait = "always"),
    bp_report_stage_metrics(state = state, stage = "AYT", stream = "main", metrics = c(
      meanAYT = "meanG",
      H2_AYT = "H2"),
      synthetic_trait = cfg$synthetic_trait, cfg = cfg,
      append_trait = "always"),
    bp_report_stage_metrics(state = state, stage = "EYT", stream = "main", metrics = c(
      meanEYT = "meanG",
      H2_EYT = "H2"),
      synthetic_trait = cfg$synthetic_trait, cfg = cfg,
      append_trait = "always")
  )

  results_year <- c(
    results_year,
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
release_variety_from_EYT <- function(state, cfg, year) {
  input_eyt <- select_latest_available(state, stage = "EYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_eyt, cfg, event_name = "release_variety_from_EYT")
  if (chk$skip) return(chk$state)

  best_variety <- bp_select_synthetic(
    pop = input_eyt$pop,
    n_select = 1L,
    synthetic_trait = cfg$synthetic_trait,
    use = "pheno",
    state = state,
    varE = cfg$varE,
    simParam = state$sim$SP
  )

  put_stage_pop(
    state = state,
    pop = best_variety,
    stage = "Variety",
    source = input_eyt,
    ready_in_years = 0,
    selection_strategy = "Best phenotypic index from latest EYT",
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
  synthetic_pheno_key <- "bps_synthetic_pheno"
  eyt_synthetic <- eyt_pop@misc[[synthetic_pheno_key]]
  if (!is.null(eyt_synthetic)) {
    parent_synthetic <- updated_parents@misc[[synthetic_pheno_key]]
    if (is.null(parent_synthetic)) {
      parent_synthetic <- matrix(
        NA_real_,
        nrow = nInd(updated_parents),
        ncol = ncol(eyt_synthetic),
        dimnames = list(updated_parents@id, colnames(eyt_synthetic))
      )
    }
    missing_cols <- setdiff(colnames(eyt_synthetic), colnames(parent_synthetic))
    if (length(missing_cols) > 0L) {
      add <- matrix(
        NA_real_,
        nrow = nrow(parent_synthetic),
        ncol = length(missing_cols),
        dimnames = list(rownames(parent_synthetic), missing_cols)
      )
      parent_synthetic <- cbind(parent_synthetic, add)
    }
    parent_synthetic[parent_idx, colnames(eyt_synthetic)] <- eyt_synthetic
    updated_parents@misc[[synthetic_pheno_key]] <- parent_synthetic
  }

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
  input_ayt <- select_latest_available(state, stage = "AYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_ayt, cfg, event_name = "select_from_AYT_and_update_parents")
  if (chk$skip) return(chk$state)

  selected <- bp_select_synthetic(
    pop = input_ayt$pop,
    n_select = min(as.integer(cfg$nReplaceParents), nInd(input_ayt$pop)),
    synthetic_trait = cfg$synthetic_trait,
    use = "pheno",
    state = state,
    varE = cfg$varE,
    simParam = state$sim$SP
  )
  current <- select_latest_available(state, stage = "Parents", stream = "main", n = 1L, combine = TRUE, silent = TRUE)

  put_stage_pop(
    state = state,
    pop = replace_oldest_parents(if (is.null(current)) NULL else current$pop, selected),
    stage = "Parents",
    source = input_ayt,
    ready_in_years = 0,
    selection_strategy = "Replace oldest parents with AYT phenotypic-index selections",
    inherit_genotypes = TRUE,
    stream = "main"
  )
}

select_from_AYT_and_start_EYT <- function(state, cfg, year) {
  input_ayt <- select_latest_available(state, stage = "AYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_ayt, cfg, event_name = "select_from_AYT_and_start_EYT")
  if (chk$skip) return(chk$state)

  eyt_pop <- bp_select_synthetic(
    pop = input_ayt$pop,
    n_select = min(as.integer(cfg$nEYT), nInd(input_ayt$pop)),
    synthetic_trait = cfg$synthetic_trait,
    use = "pheno",
    state = state,
    varE = cfg$varE,
    simParam = state$sim$SP
  )
  current_variety <- select_latest_available(state, stage = "Variety", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  if (!is.null(current_variety)) eyt_pop <- c(eyt_pop, current_variety$pop)

  run_phenotype_trial(
    state,
    pop = eyt_pop,
    output_stage = "EYT",
    input_cohorts = input_ayt$source_ids,
    selection_strategy = "Top phenotypic index from latest AYT plus current Variety check",
    traits = cfg$EYT_traits,
    synthetic_traits = cfg$synthetic_trait,
    n_loc = cfg$EYT_n_locs,
    reps = cfg$EYT_reps,
    varE = cfg$varE[cfg$EYT_traits, cfg$EYT_traits],
    duration_years = cfg$EYT_duration_years,
    stream = "main",
    cost_per_plot = cfg$EYT_cost_per_plot,
    use_env_control = TRUE,
    env_means = cfg$EYT_means,
    env_mean_sd = cfg$EYT_means_sd,
    env_year_sd = cfg$EYT_years_sd,
    log_per_environment = TRUE,
    log_aggregate = TRUE,
    silent = TRUE
  )
}

select_from_PYT_and_run_AYT <- function(state, cfg, year) {
  input_pyt <- select_latest_available(state, stage = "PYT", stream = "main", n = 1L, combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_pyt, cfg, event_name = "select_from_PYT_and_run_AYT")
  if (chk$skip) return(chk$state)

  ayt_pop <- bp_select_synthetic(
    pop = input_pyt$pop,
    n_select = min(as.integer(cfg$nAYT), nInd(input_pyt$pop)),
    synthetic_trait = cfg$synthetic_trait,
    use = "pheno",
    state = state,
    varE = cfg$varE,
    simParam = state$sim$SP
  )

  run_phenotype_trial(
    state,
    pop = ayt_pop,
    output_stage = "AYT",
    input_cohorts = input_pyt$source_ids,
    selection_strategy = "Top phenotypic index from latest PYT",
    traits = cfg$AYT_traits,
    synthetic_traits = cfg$synthetic_trait,
    n_loc = cfg$AYT_n_locs,
    reps = cfg$AYT_reps,
    varE = cfg$varE[cfg$AYT_traits, cfg$AYT_traits],
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
    headrow_pop <- bp_select_synthetic(
      pop = headrow_pop,
      n_select = min(as.integer(cfg$nPYT), nInd(headrow_pop)),
      synthetic_trait = cfg$synthetic_trait,
      use = "pheno",
      state = state,
      varE = cfg$varE,
      simParam = state$sim$SP
    )
  }

  run_phenotype_trial(
    state,
    pop = headrow_pop,
    output_stage = "PYT",
    input_cohorts = input_headrow$source_ids,
    selection_strategy = "Advance latest Candidates to PYT",
    traits = cfg$PYT_traits,
    synthetic_traits = cfg$synthetic_trait,
    n_loc = cfg$PYT_n_locs,
    reps = cfg$PYT_reps,
    varE = cfg$varE[cfg$PYT_traits, cfg$PYT_traits],
    duration_years = cfg$PYT_duration_years,
    stream = "main",
    cost_per_plot = cfg$PYT_cost_per_plot,
    use_env_control = TRUE,
    env_means = cfg$PYT_mean,
    env_mean_sd = cfg$PYT_sd,
    env_year_sd = 0,
    log_per_environment = TRUE,
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
      traits = cfg$Headrow_traits,
      synthetic_traits = cfg$synthetic_trait,
      n_loc = cfg$Headrow_n_locs,
      reps = cfg$Headrow_reps,
      varE = cfg$varE[cfg$Headrow_traits, cfg$Headrow_traits],
      duration_years = cfg$Headrow_duration_years,
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
    return(add_stage_cost(
      state = state,
      event = "seed_increase",
      n_units = nInd(headrow_pop),
      unit_cost = cfg$Headrow_seed_increase_cost_per_line,
      stage = "Headrow",
      unit = "line",
      cohort_id = tail(as.character(state$cohorts$cohort_id), 1L)
    ))
  }

  state <- put_stage_pop(
    state = state,
    pop = headrow_pop,
    stage = "Headrow",
    source = input_f1,
    ready_in_years = cfg$Headrow_duration_years,
    cross_strategy = sprintf(
      "Self latest F1 for %d cycles to Headrow",
      as.integer(cfg$n_selfing_cycles_headrow)
    ),
    inherit_genotypes = FALSE,
    stream = "main",
    cost_per_individual = cfg$Headrow_line_development_cost_per_line,
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
    ready_in_years = 1,
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

    state <- bp_advance_time_years(state, years = 1)

    # Year End
    results <- record_yearly_outputs(results, results_base, state, state$time$t - report_start_year, cfg)
  }

  if (!is.null(results_file)) write.csv(results, file = results_file, row.names = FALSE)
  invisible(list(state = state, results = results))
}

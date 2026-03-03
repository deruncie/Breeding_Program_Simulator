# GSTP-like scheme with readable yearly schedule + rapid recurrent GS.
#
# This is a structured, readability-focused script with explicit sections for:
# - setup/imports
# - helper utilities
# - event verbs
# - yearly schedule
# - configuration/state bootstrap
# - simulation run and reporting

# Setup -------------------------------------------------------------------
library(AlphaSimR)

source("R/operators.R")
source("R/readable_api.R")
source("R/monitoring.R")


# Helper Utilities ---------------------------------------------------------
make_cross_plan_no_self <- function(n_parents, n_crosses) {
  if (n_parents < 2L) return(matrix(c(1L, 1L), ncol = 2L))

  p1 <- sample.int(n_parents, size = n_crosses, replace = TRUE)
  p2 <- vapply(p1, function(i) {
    j <- sample.int(n_parents - 1L, size = 1L)
    if (j >= i) j <- j + 1L
    as.integer(j)
  }, integer(1))
  cbind(p1, p2)
}


# Event Verbs: Recurrent Cycle --------------------------------------------
select_from_RC_and_cross_next_RC_tick <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_rc <- select_latest_available(state, stage = cfg$rc_stage, stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_rc, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  pop <- input_rc$pop
  model <- state$gs_models[[cfg$model_id]]
  if (is.null(model)) return(state)

  pop <- predict_ebv_pop(pop, model, state, cfg = list(ebv_trait = 1, cohort_ids = input_rc$source_ids), stage_label = "RC")
  sel <- selectInd(pop, nInd = min(cfg$rc_select_n, pop_n_ind(pop)), use = "ebv", trait = 1, simParam = state$sim$SP)

  plan <- make_cross_plan_no_self(pop_n_ind(sel), cfg$rc_crosses)
  next_pop <- makeCross(sel, crossPlan = plan, nProgeny = cfg$rc_n_progeny_per_cross, simParam = state$sim$SP)

  state <- put_stage_pop(
    state = state,
    pop = next_pop,
    stage = cfg$rc_stage,
    source = input_rc,
    selection_strategy = "Top by EBV from latest RC cohort",
    cross_strategy = "Random mating without selfing among selected RC parents",
    ready_in_years = cfg$rapid_cycle_length,
    stream = "main",
    inherit_genotypes = FALSE
  )
  state <- add_stage_cost(
    state = state,
    stage = cfg$rc_stage,
    event = "rc_crossing",
    unit = "cross",
    n = cfg$rc_crosses,
    unit_cost = cfg$cost_cross
  )
  state
}

initialize_RC_from_PYT_if_missing <- function(state, cfg) {
  bp_debug_break(state, cfg)
  already <- bp_get_ready_cohorts(state, stage = cfg$rc_stage, stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)
  if (nrow(already) > 0L) return(state)

  input_pyt <- select_latest_available(state, stage = "PYT", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_pyt, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  n <- min(cfg$rc_init_n, pop_n_ind(input_pyt$pop))
  idx <- sample.int(pop_n_ind(input_pyt$pop), size = n, replace = FALSE)
  rc0 <- pop_subset(input_pyt$pop, idx)

  state <- put_stage_pop(
    state = state,
    pop = rc0,
    stage = cfg$rc_stage,
    source = input_pyt,
    selection_strategy = "Random sample from latest PYT to initialize RC",
    ready_in_years = 0,
    stream = "main",
    inherit_genotypes = TRUE
  )
  run_genotyping(state, list(
    input_stage = cfg$rc_stage,
    stream = "main",
    input_policy = "latest_one",
    include_not_ready = TRUE,
    chip = cfg$snp_chip,
    duration_years = cfg$rapid_cycle_length,
    cost_per_sample = cfg$cost_genotype,
    silent = TRUE
  ))
}

select_from_RC_and_send_to_DH <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_rc <- select_latest_available(state, stage = cfg$rc_stage, stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_rc, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  pop <- input_rc$pop
  model <- state$gs_models[[cfg$model_id]]
  if (!is.null(model)) {
    pop <- predict_ebv_pop(pop, model, state, cfg = list(ebv_trait = 1, cohort_ids = input_rc$source_ids), stage_label = "RC_to_DH")
    send <- selectInd(pop, nInd = min(cfg$dh_n_entries_per_year, pop_n_ind(pop)), use = "ebv", trait = 1, simParam = state$sim$SP)
  } else {
    idx <- sample.int(pop_n_ind(pop), size = min(cfg$dh_n_entries_per_year, pop_n_ind(pop)), replace = FALSE)
    send <- pop_subset(pop, idx)
  }

  dh_send <- makeDH(send, nDH = 1, simParam = state$sim$SP)
  state <- put_stage_pop(
    state = state,
    pop = dh_send,
    stage = cfg$dh_stage,
    source = input_rc,
    selection_strategy = if (is.null(model)) "Random selection from RC" else "Top by EBV from RC",
    cross_strategy = "makeDH(nDH=1) per selected line",
    ready_in_years = 2,
    stream = "main",
    inherit_genotypes = FALSE
  )
  add_stage_cost(state, stage = cfg$dh_stage, event = "dh_seed_pipeline", unit = "line", n = pop_n_ind(dh_send), unit_cost = cfg$cost_line_dh)
}


# Event Verbs: Pipeline Trials --------------------------------------------
select_from_DH_and_run_PYT <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_dh <- select_latest_available(state, stage = cfg$dh_stage, stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_dh, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state
  n <- min(as.integer(cfg$n_pyt), pop_n_ind(input_dh$pop))
  idx <- sample.int(pop_n_ind(input_dh$pop), size = n, replace = FALSE)
  pyt_pop <- pop_subset(input_dh$pop, idx)
  run_phenotype_trial(
    state = state,
    pop = pyt_pop,
    output_stage = "PYT",
    input_cohorts = input_dh$source_ids,
    selection_strategy = "Random subset from latest DH_PIPE",
    traits = 1,
    n_loc = 1,
    reps = 1,
    varE = cfg$varE,
    duration_years = 0.5,
    stream = "main",
    cost_per_plot = cfg$cost_plot_pyt,
    use_env_control = FALSE,
    silent = TRUE
  )
}

genotype_latest_PYT <- function(state, cfg) {
  bp_debug_break(state, cfg)
  run_genotyping(state, list(
    input_stage = "PYT",
    stream = "main",
    input_policy = "latest_one",
    include_not_ready = TRUE,
    chip = cfg$snp_chip,
    duration_years = 0.5,
    cost_per_sample = cfg$cost_genotype,
    silent = TRUE
  ))
}

genotype_latest_RC <- function(state, cfg) {
  bp_debug_break(state, cfg)
  run_genotyping(state, list(
    input_stage = cfg$rc_stage,
    stream = "main",
    input_policy = "latest_one",
    include_not_ready = TRUE,
    chip = cfg$snp_chip,
    duration_years = cfg$rapid_cycle_length,
    cost_per_sample = cfg$cost_genotype,
    silent = TRUE
  ))
}

select_from_PYT_and_run_AYT <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_pyt <- select_latest_available(state, stage = "PYT", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_pyt, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state
  pop_in <- input_pyt$pop
  model <- state$gs_models[[cfg$model_id]]
  pop_scored <- if (is.null(model)) {
    pop_in
  } else {
    predict_ebv_pop(pop_in, model, state, cfg = list(ebv_trait = 1, cohort_ids = input_pyt$source_ids), stage_label = "AYT")
  }
  n <- min(cfg$n_ayt, pop_n_ind(pop_scored))
  ayt_pop <- if (is.null(model)) {
    selectInd(pop_scored, nInd = n, use = "pheno", simParam = state$sim$SP)
  } else {
    selectInd(pop_scored, nInd = n, use = "ebv", trait = 1, simParam = state$sim$SP)
  }

  run_phenotype_trial(
    state = state,
    pop = ayt_pop,
    output_stage = "AYT",
    input_cohorts = input_pyt$source_ids,
    selection_strategy = "Top by EBV from latest PYT (fallback: phenotype)",
    traits = 1,
    n_loc = 4,
    reps = 2,
    varE = cfg$varE,
    duration_years = 0.5,
    stream = "main",
    cost_per_plot = cfg$cost_plot_ayt,
    use_env_control = TRUE,
    env_mean_mu = 0,
    env_mean_sd = 1,
    env_year_sd = cfg$ayt_env_year_sd,
    log_per_environment = FALSE,
    log_aggregate = TRUE,
    silent = TRUE
  )
}

select_from_AYT_and_run_EYT <- function(state, cfg) {
  bp_debug_break(state, cfg)
  if (is.null(state$sim$eyt_base_means)) state$sim$eyt_base_means <- stats::rnorm(10, mean = 0, sd = 1)

  year_shift_1 <- stats::rnorm(1, mean = 0, sd = cfg$eyt_year_shift_sd)
  year_shift_2 <- stats::rnorm(1, mean = 0, sd = cfg$eyt_year_shift_sd)
  env_means_20 <- c(state$sim$eyt_base_means + year_shift_1, state$sim$eyt_base_means + year_shift_2)

  input_ayt <- select_latest_available(state, stage = "AYT", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_ayt, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state
  n <- min(cfg$n_eyt, pop_n_ind(input_ayt$pop))
  eyt_pop <- selectInd(input_ayt$pop, nInd = n, use = "pheno", simParam = state$sim$SP)

  run_phenotype_trial(
    state = state,
    pop = eyt_pop,
    output_stage = "EYT",
    input_cohorts = input_ayt$source_ids,
    selection_strategy = "Top by phenotype from latest AYT",
    traits = 1,
    n_loc = 20,
    reps = 4,
    varE = cfg$varE,
    duration_years = 1.5,
    stream = "main",
    cost_per_plot = cfg$cost_plot_eyt,
    use_env_control = TRUE,
    env_means = env_means_20,
    env_year_sd = cfg$eyt_env_within_sd,
    log_per_environment = FALSE,
    log_aggregate = TRUE,
    silent = TRUE
  )
}

select_from_EYT_and_release_Variety <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_eyt <- select_latest_available(state, stage = "EYT", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_eyt, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  variety <- selectInd(input_eyt$pop, nInd = 1, use = "pheno", simParam = state$sim$SP)
  old_var <- bp_get_ready_cohorts(state, stage = "Variety", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)
  if (nrow(old_var) > 0L) {
    for (i in seq_len(nrow(old_var))) state <- bp_close_cohort(state, old_var$cohort_id[i])
  }

  state <- put_stage_pop(
    state = state,
    pop = variety,
    stage = "Variety",
    source = input_eyt,
    selection_strategy = "Top by phenotype from latest EYT",
    ready_in_years = 0,
    stream = "main",
    inherit_genotypes = TRUE
  )
  state$outputs$varieties <- rbind(
    state$outputs$varieties,
    data.frame(
      tick = as.integer(state$time$tick),
      source_cohort_id = paste(input_eyt$source_ids, collapse = ";"),
      variety_id = as.integer(variety@id),
      stringsAsFactors = FALSE
    )
  )
  state
}

update_GP_model_from_recent_PYT_midyear <- function(state, cfg) {
  bp_debug_break(state, cfg)
  run_train_gp_model(state, list(
    from_stage = "PYT",
    stream = "main",
    chip = cfg$snp_chip,
    lookback_years = 50,
    training_policy = "latest_n",
    input_n = cfg$gp_max_pyt,
    model_id = cfg$model_id,
    silent = TRUE
  ))
}


# Yearly Schedule ----------------------------------------------------------
run_one_year <- function(state, cfg, year_index) {
  bp_debug_break(state, cfg, year = year_index)

  # Stage events at year start.
  state <- initialize_RC_from_PYT_if_missing(state, cfg)
  state <- select_from_RC_and_send_to_DH(state, cfg)
  state <- select_from_DH_and_run_PYT(state, cfg)
  state <- genotype_latest_PYT(state, cfg)
  state <- select_from_PYT_and_run_AYT(state, cfg)
  state <- select_from_AYT_and_run_EYT(state, cfg)
  state <- select_from_EYT_and_release_Variety(state, cfg)

  # Within-year recurrent cycle and mid-year GP update.
  for (q in seq_len(cfg$ticks_per_year)) {
    state <- select_from_RC_and_cross_next_RC_tick(state, cfg)
    state <- genotype_latest_RC(state, cfg)
    state <- bp_advance_time(state, n_ticks = 1L)
    if (q == 2L) state <- update_GP_model_from_recent_PYT_midyear(state, cfg)
  }

  cat(sprintf(
    "year=%d t=%.2f active=%d rc=%d dh=%d pyt=%d ayt=%d eyt=%d models=%d varieties=%d\n",
    year_index,
    state$time$t,
    sum(state$cohorts$active),
    nrow(bp_get_ready_cohorts(state, stage = cfg$rc_stage, stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = cfg$dh_stage, stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "PYT", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "AYT", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "EYT", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    length(state$gs_models),
    nrow(state$outputs$varieties)
  ))

  state
}


# Config and Initialization ------------------------------------------------
make_gstp_cfg <- function(make_plots = FALSE) {
  cfg <- list(
    rapid_cycle_length = 0.25,
    model_id = "gp_main",
    snp_chip = 1L,
    make_plots = isTRUE(make_plots),
    debug = FALSE,
    debug_after_year = NULL,
    debug_after_tick = NULL,
    debug_where = NULL,

    rc_stage = "RC",
    rc_init_n = 120L,
    rc_select_n = 60L,
    rc_crosses = 60L,
    rc_n_progeny_per_cross = 2L,

    dh_stage = "DH_PIPE",
    dh_n_entries_per_year = 120L,

    n_pyt = 120L,
    gp_max_pyt = 2L,
    n_ayt = 40L,
    n_eyt = 6L,

    varE = 1.0,
    ayt_env_year_sd = 0.20,
    eyt_year_shift_sd = 0.25,
    eyt_env_within_sd = 0.15,

    cost_cross = 1,
    cost_line_dh = 0.5,
    cost_plot_pyt = 20,
    cost_plot_ayt = 30,
    cost_plot_eyt = 35,
    cost_genotype = 25
  )
  cfg$ticks_per_year <- as.integer(round(1 / cfg$rapid_cycle_length))
  cfg
}

init_gstp_sim <- function(cfg) {
  founder_haps <- quickHaplo(nInd = 120, nChr = 3, segSites = 80)
  SP <- SimParam$new(founder_haps)
  SP$addTraitA(nQtlPerChr = 20)
  SP$addSnpChip(nSnpPerChr = 40)
  SP$setVarE(h2 = 0.35)

  founders <- newPop(founder_haps, simParam = SP)
  parents0 <- selectInd(founders, nInd = 20, use = "gv", simParam = SP)

  state <- bp_init_state(SP = SP, dt = cfg$rapid_cycle_length, start_time = 0, sim = list(default_chip = cfg$snp_chip))
  state$sim$make_plots <- cfg$make_plots

  list(state = state, SP = SP, parents0 = parents0)
}

bootstrap_with_pyt_seed <- function(state, SP, parents0, cfg) {
  f1_plan <- matrix(sample.int(pop_n_ind(parents0), size = 40L * 2L, replace = TRUE), ncol = 2)
  f1_0 <- makeCross(parents0, crossPlan = f1_plan, nProgeny = 1, simParam = SP)
  dh_0 <- makeDH(f1_0, nDH = 4, simParam = SP)
  if (pop_n_ind(dh_0) > cfg$n_pyt) dh_0 <- pop_subset(dh_0, sample.int(pop_n_ind(dh_0), size = cfg$n_pyt, replace = FALSE))

  pyt_0 <- setPheno(dh_0, varE = cfg$varE, reps = 1, traits = 1, simParam = SP)
  state <- put_stage_pop(
    state,
    pyt_0,
    stage = "PYT",
    source = NULL,
    selection_strategy = "Bootstrap PYT from initial DH seed",
    ready_in_years = 1,
    stream = "main"
  )

  run_genotyping(state, list(
    input_stage = "PYT",
    stream = "main",
    input_policy = "latest_one",
    include_not_ready = TRUE,
    chip = cfg$snp_chip,
    duration_years = 0.5,
    cost_per_sample = cfg$cost_genotype,
    silent = TRUE
  ))
}


# Top-Level Run ------------------------------------------------------------
run_gstp_loop_demo_structured <- function(n_fill_years = 4L, n_run_years = 10L, make_plots = FALSE) {
  cfg <- make_gstp_cfg(make_plots = make_plots)
  sim <- init_gstp_sim(cfg)
  state <- bootstrap_with_pyt_seed(sim$state, sim$SP, sim$parents0, cfg)

  cfg$fail_on_missing_input <- FALSE
  for (fill_yr in seq_len(as.integer(n_fill_years))) {
    state <- run_one_year(state, cfg, year_index = fill_yr)
  }
  cfg$fail_on_missing_input <- TRUE
  for (run_yr in as.integer(n_fill_years) + seq_len(as.integer(n_run_years))) {
    state <- run_one_year(state, cfg, year_index = run_yr)
  }
  invisible(state)
}


# Reporting ----------------------------------------------------------------
print_gstp_summary <- function(state, show_event_timeline = FALSE) {
  ticks_per_year <- as.integer(round(1 / state$time$dt))
  metrics <- bp_extract_cohort_metrics(
    state = state,
    stages = c("DH_PIPE", "PYT", "AYT", "EYT", "Variety"),
    trait = 1L,
    origin_stage = "DH_PIPE",
    include_inactive = TRUE,
    ticks_per_year = ticks_per_year
  )

  cat("\nVariety releases:\n")
  print(state$outputs$varieties)

  cat("\nTotal cost by event:\n")
  print(stats::aggregate(total_cost ~ event, data = state$cost_log, sum))

  cat("\nCohort metrics:\n")
  print(metrics[order(metrics$available_year, metrics$stage), c(
    "cohort_id", "stage", "origin_cohort_id", "origin_year", "available_year",
    "mean_gv", "var_gv", "max_gv", "cor_ebv_gv", "h2", "H2"
  )])

  if (isTRUE(show_event_timeline)) {
    cat("\nEvent timeline:\n")
    bp_print_event_timeline(state, collapse_year_patterns = TRUE, digits = 2)
  }
}


# Script Entry -------------------------------------------------------------
if (identical(environment(), globalenv()) && !identical(Sys.getenv("BPS_SKIP_SCRIPT_ENTRY"), "1")) {
  out <- run_gstp_loop_demo_structured(n_fill_years = 4L, n_run_years = 10L, make_plots = FALSE)
  print_gstp_summary(out)
}

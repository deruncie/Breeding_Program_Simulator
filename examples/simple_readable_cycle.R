# Simple readable breeding cycle using functional API.
# Scheme:
# inbred parents -> random F1 -> SSD to F5 -> all into PYT -> genotype PYT
# -> train GP model (RRBLUP) -> select variety from PYT phenotypes
# -> recycle best predicted PYT lines to new parents.
# Durations: SSD = 2 years, all other stages = 1 year.

library(AlphaSimR)

source("R/operators.R")
source("R/readable_api.R")

# --- local stage functions (self-contained example)
run_make_f1_local <- function(state, cfg) {
  src <- get_ready_pop(state, stage = cfg$input_stage, policy = "latest_one", combine = TRUE, silent = TRUE)
  if (is.null(src)) return(state)

  n_par <- pop_n_ind(src$pop)
  plan <- matrix(sample.int(n_par, size = as.integer(cfg$n_crosses) * 2L, replace = TRUE), ncol = 2L)
  f1 <- makeCross(src$pop, crossPlan = plan, nProgeny = as.integer(cfg$n_progeny_per_cross), simParam = state$sim$SP)

  state <- put_stage_pop(
    state = state,
    pop = f1,
    stage = cfg$output_stage,
    source = src,
    ready_in_years = cfg$duration_years,
    selection_strategy = "Random parent sampling",
    cross_strategy = "Random pairing with replacement"
  )
  state <- add_stage_cost(state, stage = cfg$output_stage, event = "crossing", unit = "cross", n = cfg$n_crosses, unit_cost = cfg$cost_per_cross)
  close_sources(state, src)
}

run_ssd_to_f5_local <- function(state, cfg) {
  ready <- bp_get_ready_cohorts(state, stage = cfg$input_stage)
  if (nrow(ready) == 0L) return(state)
  ready <- bp_select_source_rows(state, ready, list(input_policy = "all_ready"))
  if (nrow(ready) == 0L) return(state)

  for (i in seq_len(nrow(ready))) {
    src <- ready[i, , drop = FALSE]
    pop <- state$pops[[src$cohort_id]]
    pop <- self(pop, nProgeny = as.integer(cfg$n_lines_per_f1), simParam = state$sim$SP)
    pop <- self(pop, nProgeny = 1L, simParam = state$sim$SP)
    pop <- self(pop, nProgeny = 1L, simParam = state$sim$SP)
    pop <- self(pop, nProgeny = 1L, simParam = state$sim$SP)

    state <- put_stage_pop(
      state = state,
      pop = pop,
      stage = cfg$output_stage,
      source_ids = src$cohort_id,
      ready_in_years = cfg$duration_years,
      selection_strategy = "SSD to F5 (4 selfing generations)"
    )
    state <- add_stage_cost(
      state = state,
      stage = cfg$output_stage,
      cohort_id = bp_last_cohort_id(state),
      event = "ssd",
      unit = "line",
      n_units = pop_n_ind(pop),
      unit_cost = cfg$cost_per_line
    )
    if (isTRUE(cfg$consume_input)) state <- bp_close_cohort(state, src$cohort_id)
  }
  state
}

run_pyt_local <- function(state, cfg) {
  run_phenotype_trial(state, list(
    trial_name = "PYT",
    input_stage = cfg$input_stage,
    output_stage = cfg$output_stage,
    input_policy = "latest_one",
    selection_strategy = "All lines from latest F5",
    traits = cfg$traits,
    n_loc = cfg$n_loc,
    reps = cfg$reps,
    varE = cfg$varE,
    duration_years = cfg$duration_years,
    consume_input = cfg$consume_input,
    cost_per_plot = cfg$cost_per_plot,
    silent = TRUE
  ))
}

run_select_variety_and_recycle_local <- function(state, cfg) {
  src <- get_ready_pop(state, stage = cfg$input_stage, policy = "latest_one", combine = TRUE, silent = TRUE)
  if (is.null(src)) return(state)
  model <- state$gs_models[[cfg$model_id]]
  if (is.null(model)) return(state)

  pop <- src$pop
  variety <- selectInd(pop, nInd = as.integer(cfg$n_variety), use = "pheno", simParam = state$sim$SP)
  state <- put_stage_pop(
    state = state,
    pop = variety,
    stage = "Variety",
    source = src,
    ready_in_years = 0,
    selection_strategy = "Top phenotype in PYT"
  )
  state$outputs$varieties <- rbind(
    state$outputs$varieties,
    data.frame(tick = as.integer(state$time$tick), source_cohort_id = src$source_ids[[1]], variety_id = as.integer(variety@id), stringsAsFactors = FALSE)
  )

  pop_ebv <- run_predict_ebv(pop = pop, model_entry = model, state = state, cfg = list(cohort_ids = src$source_ids), stage_label = "PYT")
  n_new <- min(as.integer(cfg$n_new_parents), pop_n_ind(pop_ebv))
  new_par <- selectInd(pop_ebv, nInd = n_new, use = "ebv", trait = 1, simParam = state$sim$SP)

  old <- bp_get_ready_cohorts(state, stage = cfg$parent_stage, active_only = TRUE, as_of_tick = .Machine$integer.max)
  if (nrow(old) > 0L) for (i in seq_len(nrow(old))) state <- bp_close_cohort(state, old$cohort_id[i])

  state <- put_stage_pop(
    state = state,
    pop = new_par,
    stage = cfg$parent_stage,
    source = src,
    ready_in_years = cfg$parent_duration_years,
    selection_strategy = "Top EBV from PYT using gp_main"
  )
  state <- add_stage_cost(state, stage = cfg$parent_stage, event = "parent_recycle", unit = "line", n = pop_n_ind(new_par), unit_cost = cfg$cost_per_parent)
  close_sources(state, src)
}

# --- setup
founder_haps <- quickHaplo(nInd = 120, nChr = 6, segSites = 250)
SP <- SimParam$new(founder_haps)
SP$addTraitA(nQtlPerChr = 80)
SP$addSnpChip(nSnpPerChr = 80)
SP$setVarE(h2 = 0.35)

founders <- newPop(founder_haps, simParam = SP)
initial_parents <- selectInd(founders, nInd = 50, use = "gv", simParam = SP)

state <- bp_init_state(SP = SP, dt = 0.25, start_time = 0, sim = list(default_chip = 1L))
state <- bp_add_cohort(
  state = state,
  pop = initial_parents,
  stage = "PARENT",
  stream = "main",
  cycle_id = "cycle_0",
  duration_years = 1
)

# --- configs
cfg <- list(
  cross = list(
    input_stage = "PARENT",
    output_stage = "F1",
    n_crosses = 100,
    n_progeny_per_cross = 1,
    duration_years = 1,
    cost_per_cross = 1
  ),
  ssd = list(
    input_stage = "F1",
    output_stage = "F5",
    n_lines_per_f1 = 2,
    duration_years = 2,
    consume_input = TRUE,
    cost_per_line = 0.5
  ),
  pyt = list(
    trial_name = "PYT",
    input_stage = "F5",
    output_stage = "PYT",
    traits = 1,
    n_loc = 1,
    reps = 1,
    varE = 1.0,
    duration_years = 1,
    consume_input = TRUE,
    cost_per_plot = 20
  ),
  genotype = list(
    input_stage = "PYT",
    input_policy = "latest_one",
    include_not_ready = TRUE,
    chip = 1L,
    duration_years = 1,
    cost_per_sample = 25
  ),
  gp = list(
    from_stage = "PYT",
    chip = 1L,
    trait = 1,
    lookback_years = 3,
    model_id = "gp_main"
  ),
  recycle = list(
    input_stage = "PYT",
    model_id = "gp_main",
    n_variety = 1,
    n_new_parents = 50,
    parent_stage = "PARENT",
    parent_duration_years = 1,
    replace_all_parents = TRUE,
    consume_input = TRUE,
    cost_per_parent = 0.5
  )
)

# --- yearly loop (4 ticks/year)
for (yr in 1:10) {
  state <- run_make_f1_local(state, cfg$cross)
  state <- run_ssd_to_f5_local(state, cfg$ssd)
  state <- run_pyt_local(state, cfg$pyt)
  state <- run_genotyping(state, cfg$genotype)
  state <- run_train_gp_model(state, cfg$gp)
  state <- run_select_variety_and_recycle_local(state, cfg$recycle)

  state <- bp_advance_time(state, n_ticks = 4L)

  cat(sprintf(
    "year=%d t=%.2f cohorts=%d varieties=%d models=%d cost=%.1f\n",
    yr,
    state$time$t,
    nrow(state$cohorts),
    nrow(state$outputs$varieties),
    length(state$gs_models),
    sum(state$cost_log$total_cost)
  ))
}

cat("\nRecent cohorts:\n")
print(utils::tail(state$cohorts[, c("cohort_id", "stage", "created_tick", "available_tick", "active", "n_ind", "genotyped", "chips")], 12))

cat("\nVarieties:\n")
print(state$outputs$varieties)

cat("\nCost by event:\n")
print(stats::aggregate(total_cost ~ event, data = state$cost_log, sum))

cat("\nEvent timeline:\n")
bp_print_event_timeline(state, collapse_year_patterns = TRUE, digits = 2)

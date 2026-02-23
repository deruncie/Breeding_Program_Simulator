# Simple readable breeding cycle using functional API.
# Scheme:
# inbred parents -> random F1 -> SSD to F5 -> all into PYT -> genotype PYT
# -> train GP model (RRBLUP) -> select variety from PYT phenotypes
# -> recycle best predicted PYT lines to new parents.
# Durations: SSD = 2 years, all other stages = 1 year.

library(AlphaSimR)

source("R/operators.R")
source("R/readable_api.R")
source("R/readable_wrappers.R")

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
  state <- run_make_f1(
    state = state,
    input_stage = cfg$cross$input_stage,
    output_stage = cfg$cross$output_stage,
    n_crosses = cfg$cross$n_crosses,
    n_progeny_per_cross = cfg$cross$n_progeny_per_cross,
    ready_in_years = cfg$cross$duration_years,
    cost_per_cross = cfg$cross$cost_per_cross
  )
  state <- run_ssd_to_f5(state, cfg$ssd)
  state <- run_multienv_trial(state)
  state <- run_genotyping(state, cfg$genotype)
  state <- run_train_gp_model(state, cfg$gp)
  state <- run_select_variety_and_recycle(state, cfg$recycle)

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

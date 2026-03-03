# Simple 5-generation test simulation:
# PARENT -> F1 (crossing) -> F2 (crossing) -> PYT (phenotype) -> AYT (phenotype)

if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
  stop("AlphaSimR is required. Install with install.packages('AlphaSimR')", call. = FALSE)
}

library(AlphaSimR)

source("R/operators.R")
source("R/readable_api.R")
source("R/state_print.R")

set.seed(123)

# --- 1) Build a small starting population of 20 individuals
founder_haps <- quickHaplo(nInd = 20, nChr = 5, segSites = 120)
SP <- SimParam$new(founder_haps)
SP$addTraitA(nQtlPerChr = 60)
SP$setVarE(h2 = 0.35)

founders <- newPop(founder_haps, simParam = SP)
initial_parents <- selectInd(founders, nInd = 20, use = "gv", simParam = SP)

state <- bp_init_state(
  SP = SP,
  dt = 1,              # one tick = one generation
  start_time = 0,
  sim = list(default_chip = 1L)
)

state <- put_stage_pop(
  state = state,
  pop = initial_parents,
  stage = "PARENT",
  stream = "main",
  cycle_id = "cycle_0",
  ready_in_years = 0,
  selection_strategy = "Top 20 founders by GV"
)

# --- 2) Stage configs: two crossing stages + two phenotype stages
cfg <- list(
  cross1 = list(
    input_stage = "PARENT",
    output_stage = "F1",
    output_stream = "main",
    n_crosses = 10,
    n_progeny_per_cross = 2,
    duration_years = 0,
    input_policy = "latest_one",
    consume_input = FALSE,
    cross_strategy = "Random pairing among parents",
    cost_per_cross = 1
  ),
  cross2 = list(
    input_stage = "F1",
    output_stage = "F2",
    output_stream = "main",
    n_crosses = 10,
    n_progeny_per_cross = 2,
    duration_years = 0,
    input_policy = "latest_one",
    consume_input = TRUE,
    cross_strategy = "Random pairing among F1",
    cost_per_cross = 1
  ),
  pheno1 = list(
    trial_name = "PYT",
    input_stage = "F2",
    output_stage = "PYT",
    input_policy = "latest_one",
    traits = 1,
    n_loc = 2,
    reps = 1,
    varE = 1.0,
    duration_years = 0,
    consume_input = TRUE,
    cost_per_plot = 8
  ),
  pheno2 = list(
    trial_name = "AYT",
    input_stage = "PYT",
    output_stage = "AYT",
    input_policy = "latest_one",
    traits = 1,
    n_loc = 3,
    reps = 2,
    varE = 0.8,
    duration_years = 0,
    consume_input = TRUE,
    cost_per_plot = 12
  )
)

# --- 3) Event verbs
select_from_PARENT_and_run_F1 <- function(state, cfg) {
  input_parent <- select_latest_available(state, stage = cfg$input_stage, stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_parent, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  pop_in <- input_parent$pop
  n_par <- pop_n_ind(pop_in)
  if (n_par < 2L) return(state)

  plan <- cbind(
    sample.int(n_par, size = cfg$n_crosses, replace = TRUE),
    sample.int(n_par, size = cfg$n_crosses, replace = TRUE)
  )
  pop_f1 <- makeCross(pop_in, crossPlan = plan, nProgeny = cfg$n_progeny_per_cross, simParam = state$sim$SP)

  state <- put_stage_pop(
    state = state,
    pop = pop_f1,
    stage = cfg$output_stage,
    source = input_parent,
    cross_strategy = cfg$cross_strategy,
    ready_in_years = cfg$duration_years,
    stream = cfg$output_stream
  )
  state <- add_stage_cost(
    state = state,
    stage = cfg$output_stage,
    event = "crossing",
    unit = "cross",
    n = cfg$n_crosses,
    unit_cost = cfg$cost_per_cross
  )
  if (isTRUE(cfg$consume_input %||% FALSE)) state <- close_sources(state, input_parent)
  state
}

select_from_F1_and_run_F2 <- function(state, cfg) {
  input_f1 <- select_latest_available(state, stage = cfg$input_stage, stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_f1, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  pop_in <- input_f1$pop
  n_par <- pop_n_ind(pop_in)
  if (n_par < 2L) return(state)

  plan <- cbind(
    sample.int(n_par, size = cfg$n_crosses, replace = TRUE),
    sample.int(n_par, size = cfg$n_crosses, replace = TRUE)
  )
  pop_f2 <- makeCross(pop_in, crossPlan = plan, nProgeny = cfg$n_progeny_per_cross, simParam = state$sim$SP)

  state <- put_stage_pop(
    state = state,
    pop = pop_f2,
    stage = cfg$output_stage,
    source = input_f1,
    cross_strategy = cfg$cross_strategy,
    ready_in_years = cfg$duration_years,
    stream = cfg$output_stream
  )
  state <- add_stage_cost(
    state = state,
    stage = cfg$output_stage,
    event = "crossing",
    unit = "cross",
    n = cfg$n_crosses,
    unit_cost = cfg$cost_per_cross
  )
  if (isTRUE(cfg$consume_input %||% FALSE)) state <- close_sources(state, input_f1)
  state
}

select_from_F2_and_run_PYT <- function(state, cfg) {
  input_f2 <- select_latest_available(state, stage = cfg$input_stage, stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_f2, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state
  run_phenotype_trial(
    state = state,
    pop = input_f2$pop,
    output_stage = cfg$output_stage,
    input_cohorts = input_f2$source_ids,
    selection_strategy = "All entries from latest F2",
    traits = cfg$traits,
    n_loc = cfg$n_loc,
    reps = cfg$reps,
    varE = cfg$varE,
    duration_years = cfg$duration_years,
    stream = "main",
    cost_per_plot = cfg$cost_per_plot,
    silent = TRUE
  ) -> state
  if (isTRUE(cfg$consume_input %||% FALSE)) state <- close_sources(state, input_f2)
  state
}

select_from_PYT_and_run_AYT <- function(state, cfg) {
  input_pyt <- select_latest_available(state, stage = cfg$input_stage, stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_pyt, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state
  run_phenotype_trial(
    state = state,
    pop = input_pyt$pop,
    output_stage = cfg$output_stage,
    input_cohorts = input_pyt$source_ids,
    selection_strategy = "All entries from latest PYT",
    traits = cfg$traits,
    n_loc = cfg$n_loc,
    reps = cfg$reps,
    varE = cfg$varE,
    duration_years = cfg$duration_years,
    stream = "main",
    cost_per_plot = cfg$cost_per_plot,
    silent = TRUE
  ) -> state
  if (isTRUE(cfg$consume_input %||% FALSE)) state <- close_sources(state, input_pyt)
  state
}

# --- 4) Run test simulation in separate fill/run loops
n_fill_generations <- 2L
n_run_generations <- 3L

for (fill_gen in seq_len(n_fill_generations)) {
  state <- select_from_PARENT_and_run_F1(state, cfg$cross1)
  state <- select_from_F1_and_run_F2(state, cfg$cross2)
  state <- select_from_F2_and_run_PYT(state, cfg$pheno1)
  state <- select_from_PYT_and_run_AYT(state, cfg$pheno2)
  state <- bp_advance_time(state, n_ticks = 1L)
}

for (run_gen in seq_len(n_run_generations)) {
  state <- select_from_PARENT_and_run_F1(state, cfg$cross1)
  state <- select_from_F1_and_run_F2(state, cfg$cross2)
  state <- select_from_F2_and_run_PYT(state, cfg$pheno1)
  state <- select_from_PYT_and_run_AYT(state, cfg$pheno2)
  cat(sprintf(
    "Generation %d complete | t=%.1f | cohorts=%d | events=%d | total_cost=%.1f\n",
    run_gen + n_fill_generations,
    state$time$t,
    nrow(state$cohorts),
    nrow(state$event_log),
    sum(state$cost_log$total_cost)
  ))
  state <- bp_advance_time(state, n_ticks = 1L)
}

# --- 5) Quick summaries
cat("\nState overview:\n")
print(bp_state_overview(state))

cat("\nPhenotype records by stage:\n")
print(stats::aggregate(individual_id ~ stage, data = state$phenotype_log, FUN = length))

cat("\nRecent event timeline:\n")
bp_print_event_timeline(state, collapse_year_patterns = FALSE, digits = 1)

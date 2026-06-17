# Setup -------------------------------------------------------------------
# BPS 0.2.0 template: create single-trait warm-start states.
#
# TODO:
#   1. Set project paths and the burn-in scheme.
#   2. Generate setup_cfg.R with bp_check_cfg_requirements().
#   3. Edit the scenario grid, founder architecture, and historical process.
#   4. Run one minimal SetupID before launching production replicates.

library(AlphaSimR)
library(BreedingProgramSimulator)

scripts_dir <- "scripts"
sim_dir <- "simulation_setups"
burnin_scheme_script <- "scheme_phenotypic_single_trait.R"
setup_cfg_file <- "setup_cfg.R"

source(file.path(scripts_dir, burnin_scheme_script))
source(file.path(scripts_dir, setup_cfg_file))

# Command-Line Arguments ---------------------------------------------------
# SetupID = 100 * ScenarioID + ReplicateID.
args <- commandArgs(trailingOnly = TRUE)
SetupID <- as.integer(args[1])
nThreads <- as.integer(args[2])
ScenarioID <- SetupID %/% 100L
ReplicateID <- SetupID %% 100L

set.seed(SetupID)
dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
state_file <- file.path(sim_dir, sprintf("Sim_setup_state_%04d.rds", SetupID))

# Config Check -------------------------------------------------------------
# Run after editing this file or the burn-in scheme.
#
# bp_check_cfg_requirements(
#   files = c(
#     file.path(scripts_dir, "Create_sim_bps.R"),
#     file.path(scripts_dir, burnin_scheme_script)
#   ),
#   cfg_file = file.path(scripts_dir, setup_cfg_file),
#   rewrite_file = TRUE
# )

# Scenario Grid ------------------------------------------------------------
# TODO:
#   Start with one reference scenario. Add only genetic-architecture or
#   historical-process contrasts needed by the experiment.
make_scenario_cfg <- function(ScenarioID, base_cfg) {
  cfg <- base_cfg
  cfg$Scenario <- "Reference"
  cfg$SetupID <- SetupID
  cfg$ScenarioID <- ScenarioID
  cfg$ReplicateID <- ReplicateID
  cfg
}

ScenarioIDs <- 1L
scenario_cfgs <- setNames(lapply(ScenarioIDs, make_scenario_cfg, base_cfg = cfg), ScenarioIDs)
cfg <- scenario_cfgs[[as.character(ScenarioID)]]

# Founder Population -------------------------------------------------------
# TODO:
#   Replace runMacs() only when the user requests another founder/history
#   model. Keep every variable quantity in cfg.
founderPop <- runMacs(
  nInd = cfg$nFounders,
  nChr = cfg$nChr,
  segSites = cfg$nQtlPerChrom + cfg$nSnpPerChrom,
  inbred = TRUE,
  species = cfg$species
)

SP <- SimParam$new(founderPop)
SP$nThreads <- nThreads
SP$restrSegSites(cfg$nQtlPerChrom, cfg$nSnpPerChrom)
SP$addSnpChip(cfg$nSnpPerChrom)
SP$addTraitA(
  nQtlPerChr = cfg$nQtlPerChrom,
  mean = cfg$initMeanG,
  var = cfg$initVarG
)
SP$setTrackPed(TRUE)

# Initial Parents ----------------------------------------------------------
# TODO:
#   Match initial phenotyping and selection to the requested historical state.
base_pop <- newPop(founderPop)
base_pop <- setPheno(
  pop = base_pop,
  varE = cfg$initial_varE,
  reps = cfg$initial_reps,
  simParam = SP
)
parents <- selectInd(
  pop = base_pop,
  nInd = cfg$nParents,
  use = "pheno",
  simParam = SP
)

state <- bp_init_state(
  SP = SP,
  dt = 1 / cfg$ticks_per_year,
  start_time = 0,
  sim = list(default_chip = 1L)
)
state <- put_stage_pop(
  state = state,
  pop = parents,
  stage = "Parents",
  source_ids = "UNKNOWN",
  ready_in_years = 0,
  selection_strategy = "Initial parents",
  inherit_genotypes = TRUE,
  stream = "main"
)

# Calibration Plan ---------------------------------------------------------
# TODO:
#   1. Set cfg$calibration_stage to the most recent new-candidate stage that
#      best matches breeder experience.
#   2. Set cfg$target_candidate_reliability for its single-location trial.
#   3. Put several eligible burn-in years in cfg$calibration_burnin_years.
#   4. Replace this calculation when multi-environment or GxE evidence is more
#      relevant. Match the noise calculation to the actual trial model.
calibrate_from_recent_candidates <- function(state, cfg) {
  candidates <- select_latest_available(
    state = state,
    stage = cfg$calibration_stage,
    stream = "main",
    n = 1L,
    combine = TRUE,
    silent = TRUE
  )

  candidate_varG <- stats::var(candidates$pop@gv[, 1L])
  cfg$varE <- bp_make_varE(
    h2 = cfg$target_candidate_reliability,
    varG = candidate_varG,
    corE = matrix(1, nrow = 1L)
  )[1L, 1L]
  state <- bp_set_trait_baseline(
    state = state,
    pop = candidates$pop,
    traits = 1L,
    label = "default"
  )

  list(state = state, cfg = cfg)
}

# Burn-In ------------------------------------------------------------------
# TODO:
#   Use only enough burn-in for the intended starting state. Smoke tests should
#   use the shortest duration that fills every required stage. Production
#   burn-in should repeat the same calibration at the explicit cfg years.
if (cfg$nYearsBurnin > 0L) {
  burnin_cfg <- cfg
  burnin_cfg$nYearsFuture <- 1L
  for (burnin_year in seq_len(cfg$nYearsBurnin)) {
    if (burnin_year %in% cfg$calibration_burnin_years) {
      calibration <- calibrate_from_recent_candidates(state, burnin_cfg)
      state <- calibration$state
      burnin_cfg <- calibration$cfg
    }
    burnin <- run_simulation(
      state = state,
      results_base = data.frame(SetupID = SetupID),
      cfg = burnin_cfg,
      start_year = 0
    )
    state <- burnin$state
  }
  cfg <- burnin_cfg
}

# Save State ---------------------------------------------------------------
state$sim$setup <- cfg
saveRDS(state, state_file)
cat(sprintf("Wrote setup state to %s\n", state_file))

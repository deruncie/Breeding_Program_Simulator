# Setup -------------------------------------------------------------------
# BPS 0.2.0 template: create multi-trait warm-start states.
#
# TODO:
#   1. Set project paths, trait definitions, and the burn-in scheme.
#   2. Generate setup_cfg.R with bp_check_cfg_requirements().
#   3. Edit genetic/residual correlation scenarios.
#   4. Run one minimal SetupID before production replicates.

library(AlphaSimR)
library(BreedingProgramSimulator)

scripts_dir <- "scripts"
sim_dir <- "simulation_setups"
burnin_scheme_script <- "scheme_phenotypic_multi_trait.R"
setup_cfg_file <- "setup_cfg.R"

source(file.path(scripts_dir, burnin_scheme_script))
source(file.path(scripts_dir, setup_cfg_file))

# Command-Line Arguments ---------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
SetupID <- as.integer(args[1])
nThreads <- as.integer(args[2])
ScenarioID <- SetupID %/% 100L
ReplicateID <- SetupID %% 100L

set.seed(SetupID)
dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
state_file <- file.path(sim_dir, sprintf("Sim_setup_state_%04d.rds", SetupID))

# Config Check -------------------------------------------------------------
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
#   Replace the reference row with only the correlations or architectures the
#   user wants to compare.
scenario_design <- data.frame(
  ScenarioID = 1L,
  Scenario = "Reference",
  r_g = cfg$reference_r_g,
  r_e = cfg$reference_r_e
)
design <- scenario_design[scenario_design$ScenarioID == ScenarioID, , drop = FALSE]
cfg$Scenario <- design$Scenario
cfg$r_g <- design$r_g
cfg$r_e <- design$r_e
cfg$SetupID <- SetupID
cfg$ScenarioID <- ScenarioID
cfg$ReplicateID <- ReplicateID

cfg$traitCorA <- cfg$r_g + diag(1 - cfg$r_g, cfg$nTraits)
traitCorE <- cfg$r_e + diag(1 - cfg$r_e, cfg$nTraits)
cfg$varE <- bp_make_varE(
  h2 = cfg$traitH2,
  varG = cfg$traitVarG,
  corE = traitCorE,
  trait_names = cfg$traitNames
)

# Founder Population -------------------------------------------------------
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
  var = cfg$traitVarG,
  corA = cfg$traitCorA,
  name = cfg$traitNames
)
SP$setTrackPed(TRUE)

# Initial Parents ----------------------------------------------------------
base_pop <- newPop(founderPop)
base_pop <- setPheno(
  pop = base_pop,
  varE = cfg$varE,
  reps = cfg$initial_reps,
  traits = seq_len(cfg$nTraits),
  simParam = SP
)
parents <- selectInd(
  pop = base_pop,
  nInd = cfg$nParents,
  trait = AlphaSimR::selIndex,
  use = "pheno",
  b = cfg$initial_trait_weights,
  simParam = SP
)

state <- bp_init_state(
  SP = SP,
  dt = 1 / cfg$ticks_per_year,
  start_time = 0,
  sim = list(default_chip = 1L)
)
state <- bp_register_synthetic_traits(state, cfg$synthetic_trait)
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
#   Use recent new candidates and breeder-relevant trial reliability. The
#   reference cfg below starts from a single-location candidate trial. Replace
#   cfg$calibration_corE and the noise model when multi-environment evidence,
#   year effects, or GxE are more relevant.
calibrate_from_recent_candidates <- function(state, cfg) {
  candidates <- select_latest_available(
    state = state,
    stage = cfg$calibration_stage,
    stream = "main",
    n = 1L,
    combine = TRUE,
    silent = TRUE
  )

  cfg$traitVarG <- diag(AlphaSimR::varG(candidates$pop))
  cfg$varE <- bp_make_varE(
    h2 = cfg$target_candidate_reliability,
    varG = cfg$traitVarG,
    corE = cfg$calibration_corE,
    trait_names = cfg$traitNames
  )
  state <- bp_set_trait_baseline(
    state = state,
    pop = candidates$pop,
    traits = seq_len(cfg$nTraits),
    synthetic_traits = cfg$synthetic_trait,
    varE = cfg$varE,
    label = "default"
  )

  list(state = state, cfg = cfg)
}

# Burn-In ------------------------------------------------------------------
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

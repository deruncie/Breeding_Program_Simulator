# Setup -------------------------------------------------------------------
# BPS 0.2.0 template: compare schemes from one saved starting state.
#
# TODO:
#   1. Set paths and list every scheme used by any ExperimentID.
#   2. Generate experiment_cfg.R with bp_check_cfg_requirements().
#   3. Replace the two example method blocks with user-requested schemes.
#   4. Run SetupID 101 / ExperimentID 1 with a minimal cfg first.

library(AlphaSimR)
library(BreedingProgramSimulator)

scripts_dir <- "scripts"
sim_dir <- "simulation_setups"
results_base_dir <- "results"
experiment_cfg_file <- "experiment_cfg.R"

scheme_scripts <- c(
  "scheme_phenotypic_single_trait.R",
  "scheme_two_part_gp_single_trait.R"
)

# Experiment Config --------------------------------------------------------
# Scan this file plus every scheme. Re-run after adding a method or cfg field.
#
# bp_check_cfg_requirements(
#   files = c(
#     file.path(scripts_dir, "run_experiments.R"),
#     file.path(scripts_dir, scheme_scripts)
#   ),
#   cfg_file = file.path(scripts_dir, experiment_cfg_file),
#   rewrite_file = TRUE
# )
source(file.path(scripts_dir, experiment_cfg_file))

# Command-Line Arguments ---------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
SetupID <- as.integer(args[1])
ExperimentID <- as.integer(args[2])
nThreads <- as.integer(args[3])

ScenarioID <- SetupID %/% 100L
ReplicateID <- SetupID %% 100L

# Load Starting State ------------------------------------------------------
# Every comparison below receives this unchanged state.
state_file <- file.path(sim_dir, sprintf("Sim_setup_state_%04d.rds", SetupID))
state <- readRDS(state_file)
state <- bp_reset_costs(state)
state$sim$SP$nThreads <- nThreads

# Initial-State Calibration ------------------------------------------------
# TODO:
#   1. Select the same breeder-relevant recent candidate stage used during
#      burn-in.
#   2. Reapply that calibration procedure once at the comparison start.
#   3. Prefer the saved final burn-in noise parameters unless calibration is an
#      experiment treatment.
#   4. Save the reference cohort/time, targets, and calibrated values in
#      results metadata. Do not calibrate methods separately.
calibration_candidates <- select_latest_available(
  state = state,
  stage = cfg$calibration_stage,
  stream = "main",
  n = 1L,
  combine = TRUE,
  silent = TRUE
)
# TODO: Reapply the same single- or multi-trait bp_set_trait_baseline() call
# used in Create_sim_bps.R. Keep synthetic-trait arguments for multi-trait
# experiments. The saved state retains the final burn-in baseline until then.
cfg$varE <- state$sim$setup$varE
cfg_base <- cfg
comparison_start_year <- state$time$t

# Results Metadata ---------------------------------------------------------
# TODO:
#   Add setup/scenario fields needed to interpret every output row.
results_base <- data.frame(
  SetupID = SetupID,
  ScenarioID = ScenarioID,
  ReplicateID = ReplicateID,
  CalibrationStage = cfg$calibration_stage,
  CalibrationCohort = calibration_candidates$source_ids,
  CalibrationTime = state$time$t
)

# Base Experiment ----------------------------------------------------------
# ExperimentID 1 should always contain the primary requested comparison.
if (ExperimentID == 1L) {
  experiment_dir <- file.path(results_base_dir, "base_experiment")
  dir.create(experiment_dir, showWarnings = FALSE, recursive = TRUE)

  # Method 1: replace script and label as needed.
  source(file.path(scripts_dir, "scheme_phenotypic_single_trait.R"))
  cfg <- cfg_base
  results_base$Method <- "Phenotypic"
  phenotypic <- run_simulation(
    state = state,
    results_base = results_base,
    cfg = cfg,
    results_file = file.path(
      experiment_dir,
      sprintf("%s_results_%04d.csv", results_base$Method, SetupID)
    ),
    seed = SetupID,
    start_year = comparison_start_year
  )

  # Method 2: reset cfg and reuse the original state.
  source(file.path(scripts_dir, "scheme_two_part_gp_single_trait.R"))
  cfg <- cfg_base
  results_base$Method <- "TwoPart_GP"
  two_part <- run_simulation(
    state = state,
    results_base = results_base,
    cfg = cfg,
    results_file = file.path(
      experiment_dir,
      sprintf("%s_results_%04d.csv", results_base$Method, SetupID)
    ),
    seed = SetupID,
    start_year = comparison_start_year
  )
}

# Additional Experiment Blocks --------------------------------------------
# TODO:
#   Add one ExperimentID per coherent question. Record every varied cfg field
#   in results_base and reset cfg from cfg_base for every combination.
if (ExperimentID == 2L) {
  experiment_dir <- file.path(results_base_dir, "example_parameter_sweep")
  dir.create(experiment_dir, showWarnings = FALSE, recursive = TRUE)

  source(file.path(scripts_dir, "scheme_two_part_gp_single_trait.R"))
  for (candidate_size in cfg$experiment_candidate_sizes) {
    cfg <- cfg_base
    cfg$nPYT <- candidate_size
    results_base$Method <- "TwoPart_GP"
    results_base$nPYT <- candidate_size

    run_simulation(
      state = state,
      results_base = results_base,
      cfg = cfg,
      results_file = file.path(
        experiment_dir,
        sprintf("TwoPart_GP_nPYT_%d_%04d.csv", candidate_size, SetupID)
      ),
      seed = SetupID,
      start_year = comparison_start_year
    )
  }
}

# Output and Cost Logging --------------------------------------------------
# Add cost-log exports only when needed. Name them with method and SetupID.

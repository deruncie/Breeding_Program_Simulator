test_that("BPS 0.2.0 templates are installed and parse", {
  template_dir <- system.file(
    "templates",
    "bps-0.2.0",
    package = "BreedingProgramSimulator"
  )
  expected <- c(
    "TEMPLATE_INDEX.md",
    "Create_sim_bps_single_trait_TEMPLATE.R",
    "Create_sim_bps_multi_trait_TEMPLATE.R",
    "run_experiments_TEMPLATE.R",
    "scheme_phenotypic_single_trait.R",
    "scheme_two_part_gp_single_trait.R",
    "scheme_wfrgs_single_trait.R",
    "scheme_phenotypic_multi_trait.R"
  )

  expect_true(nzchar(template_dir))
  expect_true(all(file.exists(file.path(template_dir, expected))))

  r_files <- list.files(template_dir, pattern = "[.]R$", full.names = TRUE)
  for (file in r_files) expect_error(parse(file), NA)
})

test_that("GP templates expose the standard debug cfg", {
  template_dir <- system.file(
    "templates",
    "bps-0.2.0",
    package = "BreedingProgramSimulator"
  )
  gp_files <- file.path(
    template_dir,
    c("scheme_two_part_gp_single_trait.R", "scheme_wfrgs_single_trait.R")
  )

  for (file in gp_files) {
    text <- paste(readLines(file, warn = FALSE), collapse = "\n")
    expect_match(text, "cfg\\$debug_GP")
    expect_match(text, "cfg\\$debug_GP_n")
  }
})

test_that("orchestration templates include cfg collection guidance", {
  template_dir <- system.file(
    "templates",
    "bps-0.2.0",
    package = "BreedingProgramSimulator"
  )
  orchestration <- file.path(
    template_dir,
    c(
      "Create_sim_bps_single_trait_TEMPLATE.R",
      "Create_sim_bps_multi_trait_TEMPLATE.R",
      "run_experiments_TEMPLATE.R"
    )
  )

  for (file in orchestration) {
    text <- paste(readLines(file, warn = FALSE), collapse = "\n")
    expect_match(text, "bp_check_cfg_requirements")
  }
})

test_that("experiment templates expose calibration workflow", {
  template_dir <- system.file(
    "templates",
    "bps-0.2.0",
    package = "BreedingProgramSimulator"
  )
  create_files <- file.path(
    template_dir,
    c(
      "Create_sim_bps_single_trait_TEMPLATE.R",
      "Create_sim_bps_multi_trait_TEMPLATE.R"
    )
  )

  for (file in create_files) {
    text <- paste(readLines(file, warn = FALSE), collapse = "\n")
    expect_match(text, "calibrate_from_recent_candidates")
    expect_match(text, "cfg\\$calibration_burnin_years")
    expect_match(text, "cfg\\$target_candidate_reliability")
  }

  experiment_text <- paste(
    readLines(file.path(template_dir, "run_experiments_TEMPLATE.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(experiment_text, "CalibrationCohort")
  expect_match(experiment_text, "state\\$sim\\$setup\\$varE")
})

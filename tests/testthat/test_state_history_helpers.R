test_that("bp_years_to_ticks and bp_advance_time_years use state dt", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 0.25, start_time = 0)

  expect_equal(BreedingProgramSimulator:::bp_years_to_ticks(state, 1), 4L)
  expect_equal(BreedingProgramSimulator:::bp_years_to_ticks(state, 0.5), 2L)
  expect_error(
    BreedingProgramSimulator:::bp_years_to_ticks(state, 0.1),
    "cannot be represented exactly"
  )

  state <- BreedingProgramSimulator:::bp_advance_time_years(state, 0.5)
  expect_equal(state$time$tick, 2L)
  expect_equal(state$time$t, 0.5)
  expect_equal(BreedingProgramSimulator:::bp_tick_fraction(state), 0.5)
  expect_true(BreedingProgramSimulator:::bp_is_fraction_tick(state, 0.5))
})

test_that("run_genotyping can target explicit cohort ids", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:3), stage = "A", duration_years = 0)
  cid_a <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:4), stage = "B", duration_years = 0)
  cid_b <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  state <- BreedingProgramSimulator:::run_genotyping(
    state,
    list(cohort_ids = cid_b, chip = 1L, duration_years = 0, cost_per_sample = 2)
  )

  expect_false(any(state$genotype_log$cohort_id == cid_a))
  expect_true(any(state$genotype_log$cohort_id == cid_b))
  expect_equal(sum(state$cost_log$total_cost), 8)
})

test_that("bp_reset_costs clears cost log", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cost(state, "A", "c1", "event", "unit", 2, 3)

  state <- BreedingProgramSimulator:::bp_reset_costs(state)
  expect_equal(nrow(state$cost_log), 0L)
  expect_equal(names(state$cost_log), names(BreedingProgramSimulator:::bp_empty_cost_log()))
})

test_that("bp_forget_history keeps active/latest cohorts and drops logs", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:5), stage = "PYT", duration_years = 0)
  old_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  state <- BreedingProgramSimulator:::bp_add_cost(state, "PYT", old_id, "grow", "line", 5, 1)
  state <- BreedingProgramSimulator:::bp_close_cohort(state, old_id)
  state <- BreedingProgramSimulator:::bp_advance_time(state, 2L)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:4), stage = "PYT", duration_years = 0)
  new_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  state$phenotype_log <- rbind(
    state$phenotype_log,
    data.frame(
      cohort_id = new_id,
      stage = "PYT",
      individual_id = 1L,
      environment = "0",
      trait = "trait1",
      phenotype_value = 1,
      p_value = NA_real_,
      measured_tick = 2L,
      available_tick = 2L,
      n_loc = 1L,
      reps = 1L,
      stringsAsFactors = FALSE
    )
  )

  out <- BreedingProgramSimulator:::bp_forget_history(
    state,
    keep_stages = "PYT",
    reset_time = TRUE
  )

  expect_equal(nrow(out$cohorts), 1L)
  expect_equal(out$cohorts$cohort_id, paste("Initialization", new_id, sep = "_"))
  expect_equal(out$time$tick, 0L)
  expect_equal(out$cohorts$created_tick, 0L)
  expect_equal(nrow(out$phenotype_log), 0L)
  expect_equal(nrow(out$genotype_log), 0L)
  expect_equal(nrow(out$cost_log), 0L)
  expect_equal(length(out$gs_models), 0L)
  expect_true(any(out$event_log$event_type == "initialization"))
})

test_that("bp_compact_history downsamples only unprotected older populations", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:20), stage = "A", duration_years = 0)
  old_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  state <- BreedingProgramSimulator:::bp_close_cohort(state, old_id)
  state <- BreedingProgramSimulator:::bp_advance_time(state, 1L)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:20), stage = "A", duration_years = 0)
  current_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  state <- BreedingProgramSimulator:::bp_compact_history(
    state,
    stages = "A",
    max_n = 5L,
    selection = "first",
    log_event = TRUE
  )

  expect_equal(nrow(state$pops[[old_id]]), 5L)
  expect_equal(nrow(state$pops[[current_id]]), 20L)
  expect_equal(state$cohorts$n_ind[match(old_id, state$cohorts$cohort_id)], 20L)
  expect_true(any(state$sim$compaction_log$cohort_id == old_id))
})

test_that("bp_update_stage_pop validates ids and can reorder", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  pop <- data.frame(v = c(1, 2, 3), row.names = c("a", "b", "c"))
  state <- BreedingProgramSimulator:::bp_add_cohort(state, pop, stage = "A", duration_years = 0)
  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  reordered <- data.frame(v = c(30, 10, 20), row.names = c("c", "a", "b"))
  expect_error(
    BreedingProgramSimulator:::bp_update_stage_pop(state, cid, reordered),
    "different order"
  )

  state <- BreedingProgramSimulator:::bp_update_stage_pop(state, cid, reordered, allow_reorder = TRUE)
  expect_equal(rownames(state$pops[[cid]]), c("a", "b", "c"))
  expect_equal(state$pops[[cid]]$v, c(10, 20, 30))
})

test_that("bp_set_trait_baseline stores supplied values", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_set_trait_baseline(
    state,
    values = list(mean = c(Yield = 10, Height = 2), sd = c(Yield = 5, Height = 1)),
    include_index = TRUE,
    index_weights = c(1, -1)
  )

  baseline <- state$sim$trait_baselines$default
  expect_equal(unname(baseline$mean["Yield"]), 10)
  expect_equal(unname(baseline$sd["Height"]), 1)
  expect_equal(baseline$index_mean, 8)
})

test_that("bp_scan_cfg_requirements and bp_check_cfg_requirements report missing fields", {
  f <- tempfile(fileext = ".R")
  writeLines(
    c(
      "x <- cfg$nParents",
      "y <- cfg$PYT$locs %||% 3",
      "z <- cfg$trait_weights"
    ),
    f
  )

  scan <- BreedingProgramSimulator:::bp_scan_cfg_requirements(f)
  expect_true("nParents" %in% scan$required)
  expect_true("PYT.locs" %in% scan$defaulted)
  expect_match(scan$skeleton, "cfg <- list\\(", fixed = FALSE)
  expect_match(scan$skeleton, "  nParents = XX,", fixed = TRUE)
  expect_match(scan$skeleton, "  PYT = list(", fixed = TRUE)
  expect_match(scan$skeleton, "    locs = 3  # default used in", fixed = TRUE)
  expect_match(scan$skeleton, basename(f), fixed = TRUE)

  printed <- utils::capture.output(print(scan))
  expect_equal(paste(printed, collapse = "\n"), scan$skeleton)

  chk <- BreedingProgramSimulator:::bp_check_cfg_requirements(
    list(nParents = 10),
    f
  )
  expect_true("trait_weights" %in% chk$missing_required)
  expect_true("PYT.locs" %in% chk$missing_defaulted)
  expect_equal(paste(utils::capture.output(print(chk)), collapse = "\n"), chk$skeleton)
})

test_that("bp_scan_cfg_requirements groups overlapping fields across schemes", {
  d <- tempdir()
  f1 <- file.path(d, "scheme_a.R")
  f2 <- file.path(d, "scheme_b.R")
  writeLines(
    c(
      "x <- cfg$nParents",
      "y <- cfg$PYT$locs %||% 3",
      "a <- cfg$only_a"
    ),
    f1
  )
  writeLines(
    c(
      "x <- cfg$nParents",
      "y <- cfg$PYT$locs %||% 3",
      "b <- cfg$only_b"
    ),
    f2
  )

  scan <- BreedingProgramSimulator:::bp_scan_cfg_requirements(c(f1, f2))
  shared_pos <- regexpr("# Shared across: scheme_a.R, scheme_b.R", scan$skeleton, fixed = TRUE)[[1]]
  scheme_a_pos <- regexpr("# scheme_a.R only", scan$skeleton, fixed = TRUE)[[1]]
  scheme_b_pos <- regexpr("# scheme_b.R only", scan$skeleton, fixed = TRUE)[[1]]

  expect_true(shared_pos > 0)
  expect_true(scheme_a_pos > shared_pos)
  expect_true(scheme_b_pos > scheme_a_pos)
  expect_match(scan$skeleton, "  nParents = XX,", fixed = TRUE)
  expect_match(scan$skeleton, "  PYT = list(", fixed = TRUE)
  expect_match(scan$skeleton, "    locs = 3  # default used in scheme_a.R:2", fixed = TRUE)
  expect_match(scan$skeleton, "  only_a = XX,", fixed = TRUE)
  expect_match(scan$skeleton, "  only_b = XX", fixed = TRUE)
})

test_that("bp_check_cfg_requirements appends missing cfg entries and reports unused fields", {
  d <- tempdir()
  scheme <- file.path(d, "scheme_cfg_update.R")
  cfg_file <- file.path(d, "existing_cfg.R")
  writeLines(
    c(
      "x <- cfg$nParents",
      "y <- cfg$PYT$locs %||% 3",
      "z <- cfg$PYT$traits",
      "u <- cfg$stage_cost"
    ),
    scheme
  )
  writeLines(
    c(
      "cfg <- list(",
      "  nParents = 10,",
      "  PYT = list(",
      "    locs = 4",
      "  ),",
      "  unused_old = TRUE",
      ")"
    ),
    cfg_file
  )

  chk <- BreedingProgramSimulator:::bp_check_cfg_requirements(
    files = scheme,
    cfg_file = cfg_file
  )
  txt <- paste(readLines(cfg_file, warn = FALSE), collapse = "\n")

  expect_true("PYT.traits" %in% chk$added_to_file)
  expect_true("stage_cost" %in% chk$added_to_file)
  expect_false("nParents" %in% chk$added_to_file)
  expect_false("PYT.locs" %in% chk$added_to_file)
  expect_true("unused_old" %in% chk$unused_cfg_fields)
  expect_match(txt, "cfg <- utils::modifyList\\(cfg, list\\(", fixed = FALSE)
  expect_match(txt, "PYT = list\\(", fixed = FALSE)
  expect_match(txt, "traits = XX", fixed = TRUE)
  expect_match(txt, "stage_cost = XX", fixed = TRUE)
  expect_equal(length(gregexpr("traits = XX", txt, fixed = TRUE)[[1]]), 1L)

  chk2 <- BreedingProgramSimulator:::bp_check_cfg_requirements(
    files = scheme,
    cfg_file = cfg_file
  )
  txt2 <- paste(readLines(cfg_file, warn = FALSE), collapse = "\n")
  expect_equal(chk2$added_to_file, character(0))
  expect_equal(txt2, txt)

  printed <- paste(utils::capture.output(print(chk)), collapse = "\n")
  expect_match(printed, "# - unused_old", fixed = TRUE)
})

test_that("bp_scan_cfg_requirements extracts inline defaults without trailing calls", {
  f <- tempfile(fileext = ".R")
  writeLines(
    c(
      "if(cfg$Headrow_trial %||% FALSE) {",
      "state$pops[[id]] = selectInd(pop, cfg$Downsmple_Headrow %||% 1000, use = 'rand', simParam = state$sim$SP)",
      "n_selfing_cycles <- as.integer(cfg$n_selfing_cycles_headrow %||% 5L)",
      "reps = cfg$repHeadrow %||% (4/9),",
      "if(cfg$dropF1_pop %||% FALSE) {"
    ),
    f
  )

  scan <- BreedingProgramSimulator:::bp_scan_cfg_requirements(f)
  expect_match(scan$skeleton, "Headrow_trial = FALSE,", fixed = TRUE)
  expect_match(scan$skeleton, "Downsmple_Headrow = 1000,", fixed = TRUE)
  expect_match(scan$skeleton, "n_selfing_cycles_headrow = 5L,", fixed = TRUE)
  expect_match(scan$skeleton, "repHeadrow = (4/9)", fixed = TRUE)
  expect_match(scan$skeleton, "dropF1_pop = FALSE", fixed = TRUE)
  expect_false(grepl("simParam", scan$skeleton, fixed = TRUE))
  expect_false(grepl("FALSE) {", scan$skeleton, fixed = TRUE))
})

test_that("bp_report_stage_metric reports and scales simple AlphaSimR metrics", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(12, 2, 40)
  SP <- SimParam$new(h)
  SP$addTraitA(10)
  pop <- newPop(h, simParam = SP)
  pop@ebv <- pop@gv

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::put_stage_pop(state, pop, stage = "PYT", ready_in_years = 0)
  state <- BreedingProgramSimulator:::bp_set_trait_baseline(
    state,
    values = list(mean = c(trait1 = 0), sd = c(trait1 = 2))
  )

  expect_equal(
    BreedingProgramSimulator:::bp_report_stage_metric(state, "PYT", metric = "meanG"),
    mean(pop@gv[, 1]) / 2
  )
  expect_equal(
    BreedingProgramSimulator:::bp_report_stage_metric(state, "PYT", metric = "maxG"),
    max(pop@gv[, 1]) / 2
  )
  expect_equal(
    BreedingProgramSimulator:::bp_report_stage_metric(state, "PYT", metric = "varG"),
    stats::var(pop@gv[, 1]) / 4
  )
  expect_equal(
    BreedingProgramSimulator:::bp_report_stage_metric(state, "PYT", metric = "accEBV"),
    1
  )
  expect_error(
    BreedingProgramSimulator:::bp_report_stage_metric(state, "PYT", metric = "accGV"),
    "'arg' should be one of"
  )
})

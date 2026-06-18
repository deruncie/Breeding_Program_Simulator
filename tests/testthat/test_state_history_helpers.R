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

test_that("bp_set_misc_values stores aligned data used by AlphaSimR selection", {
  skip_if_not_installed("AlphaSimR")
  founder_haps <- AlphaSimR::quickHaplo(nInd = 6L, nChr = 1L, segSites = 10L)
  SP <- AlphaSimR::SimParam$new(founder_haps)
  SP$addTraitA(nQtlPerChr = 2L)
  pop <- AlphaSimR::newPop(founder_haps, simParam = SP)

  pop <- bp_set_misc_values(pop, "progeny_mean", seq_len(pop@nInd))
  selected <- AlphaSimR::selectInd(
    pop,
    nInd = 2L,
    use = function(pop, trait = NULL) pop@misc[["progeny_mean"]],
    simParam = SP
  )

  expect_equal(selected@misc$progeny_mean, c(6L, 5L))
  expect_equal(selected@id, pop@id[c(6L, 5L)])

  scores <- cbind(first = seq_len(pop@nInd), second = 11:16)
  pop <- bp_set_misc_values(pop, "scores", scores)
  expect_equal(pop[c(4L, 2L)]@misc$scores, scores[c(4L, 2L), , drop = FALSE])
  expect_error(
    bp_set_misc_values(pop, "bad", 1:2),
    "one value or row per individual"
  )
  expect_error(
    bp_set_misc_values(pop, "bps_synthetic_gv", seq_len(pop@nInd)),
    "reserved for package-managed fields"
  )
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

test_that("bp_make_varE constructs residual covariance", {
  h2 <- c(T1 = 0.5, T2 = 0.25)
  varG <- c(T1 = 2, T2 = 4)
  corE <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  out <- BreedingProgramSimulator:::bp_make_varE(h2, varG, corE, trait_names = names(h2))
  varE_diag <- ((1 - h2) / pmax(h2, 1e-8)) * varG
  D <- diag(sqrt(varE_diag), nrow = 2)

  expect_equal(unname(out), unname(D %*% corE %*% D))
  expect_equal(rownames(out), names(h2))
  expect_equal(colnames(out), names(h2))
})

test_that("bp_set_trait_baseline computes correlated index baselines from values", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  cov_mat <- matrix(c(4, 1.5, 1.5, 9), nrow = 2, dimnames = list(c("Yield", "Height"), c("Yield", "Height")))
  w <- c(1, -0.5)
  state <- BreedingProgramSimulator:::bp_set_trait_baseline(
    state,
    values = list(mean = c(Yield = 10, Height = 2), sd = c(Yield = 2, Height = 3), cov = cov_mat),
    include_index = TRUE,
    index_weights = w
  )

  baseline <- state$sim$trait_baselines$default
  expect_equal(unname(baseline$index_mean), sum(c(10, 2) * w))
  expect_equal(unname(baseline$index_sd), sqrt(drop(t(w) %*% cov_mat %*% w)))
  expect_equal(unname(baseline$mean["Index"]), baseline$index_mean)
  expect_equal(unname(baseline$sd["Index"]), baseline$index_sd)
})

test_that("bp_set_trait_baseline computes pop-derived index sd from direct index", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(20, 2, 40)
  SP <- SimParam$new(h)
  SP$addTraitA(10)
  SP$addTraitA(10)
  pop <- newPop(h, simParam = SP)
  w <- c(1.2, -0.4)

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_set_trait_baseline(state, pop = pop, include_index = TRUE, index_weights = w)
  baseline <- state$sim$trait_baselines$default
  idx <- drop(as.matrix(pop@gv) %*% matrix(w, ncol = 1L))

  expect_equal(unname(baseline$index_mean), mean(idx), tolerance = 1e-10)
  expect_equal(unname(baseline$index_sd), stats::sd(idx), tolerance = 1e-10)
  expect_equal(unname(baseline$sd["Index"]), stats::sd(idx), tolerance = 1e-10)
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
      "u <- cfg$stage_cost",
      "cfg$overwritten_old <- FALSE"
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
      "  overwritten_old = TRUE,",
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
  expect_false("overwritten_old" %in% chk$unused_cfg_fields)
  expect_true("overwritten_old" %in% chk$overwritten_cfg_fields)
  expect_match(txt, "cfg <- utils::modifyList\\(cfg, list\\(", fixed = FALSE)
  expect_match(txt, "PYT = list\\(", fixed = FALSE)
  expect_match(txt, "traits = XX", fixed = TRUE)
  expect_match(txt, "stage_cost = XX", fixed = TRUE)
  expect_equal(length(gregexpr("traits = XX", txt, fixed = TRUE)[[1]]), 1L)
  expect_true("PYT.traits" %in% chk$missing_before_update)
  expect_false("nParents" %in% chk$missing_before_update)

  chk2 <- BreedingProgramSimulator:::bp_check_cfg_requirements(
    files = scheme,
    cfg_file = cfg_file
  )
  txt2 <- paste(readLines(cfg_file, warn = FALSE), collapse = "\n")
  expect_equal(chk2$added_to_file, character(0))
  expect_equal(txt2, txt)

  printed <- paste(utils::capture.output(print(chk)), collapse = "\n")
  expect_match(printed, "# Missing cfg entries", fixed = TRUE)
  expect_match(printed, "PYT = list(", fixed = TRUE)
  expect_match(printed, "traits = XX", fixed = TRUE)
  expect_match(printed, "stage_cost = XX", fixed = TRUE)
  expect_false(grepl("nParents = XX", printed, fixed = TRUE))
  expect_false(grepl("locs = 3", printed, fixed = TRUE))
  expect_match(printed, "# cfg fields present in cfg_file but overwritten inside scanned scheme scripts:", fixed = TRUE)
  expect_match(printed, "# - overwritten_old  # overwritten in scheme_cfg_update.R:5", fixed = TRUE)
  expect_match(printed, "# - unused_old", fixed = TRUE)
})

test_that("bp_check_cfg_requirements creates missing cfg files in printed grouped format", {
  d <- tempdir()
  f1 <- file.path(d, "new_scheme_a.R")
  f2 <- file.path(d, "new_scheme_b.R")
  cfg_file <- file.path(d, "new_cfg_file.R")
  if (file.exists(cfg_file)) unlink(cfg_file)
  writeLines(c("x <- cfg$shared", "a <- cfg$only_a"), f1)
  writeLines(c("x <- cfg$shared", "b <- cfg$only_b"), f2)

  chk <- BreedingProgramSimulator:::bp_check_cfg_requirements(
    files = c(f1, f2),
    cfg_file = cfg_file
  )
  txt <- paste(readLines(cfg_file, warn = FALSE), collapse = "\n")
  printed <- paste(utils::capture.output(print(chk)), collapse = "\n")

  expect_false(chk$cfg_file_existed)
  expect_equal(txt, chk$skeleton)
  expect_equal(printed, chk$skeleton)
  expect_match(txt, "# Shared across: new_scheme_a.R, new_scheme_b.R", fixed = TRUE)
  expect_match(txt, "# new_scheme_a.R only", fixed = TRUE)
  expect_match(txt, "# new_scheme_b.R only", fixed = TRUE)
})

test_that("bp_check_cfg_requirements can rewrite grouped cfg files with imported values", {
  d <- tempdir()
  f1 <- file.path(d, "rewrite_scheme_a.R")
  f2 <- file.path(d, "rewrite_scheme_b.R")
  cfg_file <- file.path(d, "rewrite_cfg.R")
  writeLines(c("x <- cfg$shared", "a <- cfg$only_a", "n <- cfg$PYT$locs"), f1)
  writeLines(c("x <- cfg$shared", "b <- cfg$only_b", "n <- cfg$PYT$locs"), f2)
  writeLines(
    c(
      "cfg <- list(",
      "  only_a = 20,",
      "  unused_old = 99,",
      "  PYT = list(locs = 4),",
      "  shared = 10",
      ")"
    ),
    cfg_file
  )

  chk <- BreedingProgramSimulator:::bp_check_cfg_requirements(
    files = c(f1, f2),
    cfg_file = cfg_file,
    rewrite_file = TRUE
  )
  txt <- paste(readLines(cfg_file, warn = FALSE), collapse = "\n")
  printed <- paste(utils::capture.output(print(chk)), collapse = "\n")

  expect_true(chk$rewritten_file)
  expect_equal(txt, chk$skeleton)
  expect_match(txt, "# Shared across: rewrite_scheme_a.R, rewrite_scheme_b.R", fixed = TRUE)
  expect_match(txt, "shared = 10", fixed = TRUE)
  expect_match(txt, "locs = 4", fixed = TRUE)
  expect_match(txt, "only_a = 20", fixed = TRUE)
  expect_match(txt, "only_b = XX", fixed = TRUE)
  expect_false(grepl("unused_old", txt, fixed = TRUE))
  expect_match(printed, "# Rewrote grouped cfg template", fixed = TRUE)
})

test_that("bp_check_cfg_requirements can create a new grouped cfg from an old cfg file", {
  d <- tempdir()
  f1 <- file.path(d, "import_scheme_a.R")
  f2 <- file.path(d, "import_scheme_b.R")
  old_cfg <- file.path(d, "old_cfg.R")
  new_cfg <- file.path(d, "new_imported_cfg.R")
  if (file.exists(new_cfg)) unlink(new_cfg)
  writeLines(c("x <- cfg$shared", "a <- cfg$only_a"), f1)
  writeLines(c("x <- cfg$shared", "b <- cfg$only_b"), f2)
  writeLines(
    c(
      "cfg <- list(",
      "  shared = 10,",
      "  only_a = 20",
      ")",
      "cfg$only_a <- 21"
    ),
    old_cfg
  )

  chk <- BreedingProgramSimulator:::bp_check_cfg_requirements(
    files = c(f1, f2),
    cfg_file = new_cfg,
    import_cfg_file = old_cfg
  )
  txt <- paste(readLines(new_cfg, warn = FALSE), collapse = "\n")

  expect_false(chk$cfg_file_existed)
  expect_equal(txt, chk$skeleton)
  expect_match(txt, "shared = 10", fixed = TRUE)
  expect_match(txt, "only_a = 21", fixed = TRUE)
  expect_match(txt, "only_b = XX", fixed = TRUE)
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

test_that("bp_scan_cfg_requirements includes script-assigned cfg fields as optional NULL values", {
  f <- tempfile(fileext = ".R")
  writeLines(
    c(
      "cfg$ticks_per_year <- as.integer(bp_ticks_per_year(state))",
      "cfg$nCyclesPI <- as.integer(round(1 / as.numeric(cfg$speed_breeding_cycle_years)))",
      "cfg$first_year_train <- state$time$t + cfg$first_year_train",
      "x <- cfg$nCyclesPI + cfg$ticks_per_year + cfg$speed_breeding_cycle_years + cfg$first_year_train"
    ),
    f
  )

  scan <- BreedingProgramSimulator:::bp_scan_cfg_requirements(f)
  expect_true("ticks_per_year" %in% scan$fields)
  expect_true("nCyclesPI" %in% scan$fields)
  expect_true("speed_breeding_cycle_years" %in% scan$fields)
  expect_true("first_year_train" %in% scan$fields)
  expect_true(all(c("ticks_per_year", "nCyclesPI") %in% scan$assigned_in_scripts))
  expect_false(any(c("ticks_per_year", "nCyclesPI") %in% scan$required))
  expect_true("first_year_train" %in% scan$required)
  expect_match(scan$skeleton, "ticks_per_year = NULL", fixed = TRUE)
  expect_match(scan$skeleton, "nCyclesPI = NULL", fixed = TRUE)
  expect_match(scan$skeleton, "first_year_train = XX", fixed = TRUE)
})

test_that("bp_check_cfg_requirements writes missing script-assigned fields as NULL", {
  scheme <- tempfile(fileext = ".R")
  cfg_file <- tempfile(fileext = ".R")
  writeLines(
    c(
      "cfg$derived_value <- state$time$t",
      "x <- cfg$required_value"
    ),
    scheme
  )
  writeLines("cfg <- list(required_value = 10)", cfg_file)

  chk <- BreedingProgramSimulator:::bp_check_cfg_requirements(
    files = scheme,
    cfg_file = cfg_file
  )
  txt <- paste(readLines(cfg_file, warn = FALSE), collapse = "\n")

  expect_true("derived_value" %in% chk$added_to_file)
  expect_false("derived_value" %in% chk$missing_required)
  expect_false("derived_value" %in% chk$overwritten_cfg_fields)
  expect_match(txt, "derived_value = NULL", fixed = TRUE)
})

test_that("bp_scan_cfg_requirements ignores cfg references in comments", {
  f <- tempfile(fileext = ".R")
  writeLines(
    c(
      "# cfg$comment_only should not be required",
      "x <- cfg$real_field # cfg$trailing_comment should not be required",
      "label <- '# cfg$inside_string is not a cfg reference comment'",
      "y <- cfg$after_string"
    ),
    f
  )

  scan <- BreedingProgramSimulator:::bp_scan_cfg_requirements(f)
  expect_true("real_field" %in% scan$fields)
  expect_true("after_string" %in% scan$fields)
  expect_false("comment_only" %in% scan$fields)
  expect_false("trailing_comment" %in% scan$fields)
  expect_false("inside_string" %in% scan$fields)
})

test_that("bp_scan_cfg_requirements handles method calls with cfg arguments", {
  f <- tempfile(fileext = ".R")
  writeLines("SP$restrSegSites(cfg$nQtlPerChrom, cfg$nSnpPerChrom)", f)

  scan <- BreedingProgramSimulator:::bp_scan_cfg_requirements(f)

  expect_setequal(scan$fields, c("nQtlPerChrom", "nSnpPerChrom"))
  expect_equal(scan$assigned_in_scripts, character(0))
  expect_match(scan$skeleton, "nQtlPerChrom = XX", fixed = TRUE)
  expect_match(scan$skeleton, "nSnpPerChrom = XX", fixed = TRUE)
})

test_that("bp_report_stage_metrics reports and scales simple AlphaSimR metrics", {
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
    bp_report_stage_metrics(state, "PYT", metrics = c(value = "meanG"))$value,
    mean(pop@gv[, 1]) / 2
  )
  expect_equal(
    bp_report_stage_metrics(state, "PYT", metrics = c(value = "maxG"))$value,
    max(pop@gv[, 1]) / 2
  )
  expect_equal(
    bp_report_stage_metrics(state, "PYT", metrics = c(value = "varG"))$value,
    stats::var(pop@gv[, 1]) / 4
  )
  expect_equal(
    bp_report_stage_metrics(state, "PYT", metrics = c(value = "accEBV"))$value,
    1
  )
  expect_equal(
    bp_report_stage_metrics(state, "PYT", metrics = c(value = "wf_accEBV"))$value,
    1
  )

  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  pop_no_ebv <- state$pops[[cid]]
  pop_no_ebv@ebv[,] <- NA_real_
  pop_no_ebv@pheno <- pop_no_ebv@gv
  state$pops[[cid]] <- pop_no_ebv
  expect_equal(
    bp_report_stage_metrics(state, "PYT", metrics = c(value = "accEBV"))$value,
    1
  )
  expect_equal(
    bp_report_stage_metrics(state, "PYT", metrics = c(value = "wf_accEBV"))$value,
    1
  )

  pop_no_pred <- state$pops[[cid]]
  pop_no_pred@ebv <- matrix(numeric(0), nrow = pop_no_pred@nInd, ncol = 0L)
  pop_no_pred@pheno <- matrix(numeric(0), nrow = pop_no_pred@nInd, ncol = 0L)
  state$pops[[cid]] <- pop_no_pred
  expect_true(is.na(bp_report_stage_metrics(state, "PYT", metrics = c(value = "accEBV"))$value))
  expect_true(is.na(bp_report_stage_metrics(state, "PYT", metrics = c(value = "wf_accEBV"))$value))

  expect_error(
    bp_report_stage_metrics(state, "PYT", metrics = "accGV"),
    "metrics must contain only"
  )
})

test_that("bp_report_stage_metrics treats zero-column EBV matrices as unavailable", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  founder <- quickHaplo(nInd = 12, nChr = 1, segSites = 30)
  SP <- SimParam$new(founder)
  SP$addTraitA(10)

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)
  pop <- newPop(founder, simParam = SP)
  pop@ebv <- matrix(numeric(0), nrow = pop@nInd, ncol = 0L)
  pop@pheno <- pop@gv
  state <- BreedingProgramSimulator:::put_stage_pop(state, pop, stage = "PYT", ready_in_years = 0)

  expect_equal(
    bp_report_stage_metrics(state, "PYT", metrics = c(value = "accEBV"), traits = 1L)$value,
    1
  )

  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  state$pops[[cid]]@pheno <- matrix(numeric(0), nrow = pop@nInd, ncol = 0L)
  expect_true(is.na(bp_report_stage_metrics(state, "PYT", metrics = c(value = "accEBV"), traits = 1L)$value))
})

test_that("bp_report_stage_metrics reports biological traits and synthetic index separately", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(24, 2, 50)
  SP <- SimParam$new(h)
  SP$addTraitA(10)
  SP$addTraitA(10)
  pop <- newPop(h, simParam = SP)
  pop@ebv <- pop@gv
  w <- c(1, -0.5)
  index_def <- bp_synthetic_trait(
    "Index",
    traits = 1:2,
    fun = AlphaSimR::selIndex,
    args = list(b = w),
    linear = TRUE
  )
  idx <- as.numeric(AlphaSimR::selIndex(pop@gv[, 1:2, drop = FALSE], b = w))
  pop <- bp_set_synthetic_values(pop, "Index", idx, type = "ebv")

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)
  state <- bp_register_synthetic_traits(state, index_def)
  state <- BreedingProgramSimulator:::bp_set_trait_baseline(
    state,
    pop = pop,
    synthetic_traits = "Index"
  )
  state <- BreedingProgramSimulator:::put_stage_pop(state, pop, stage = "PYT", ready_in_years = 0)

  mean_out <- bp_report_stage_metrics(
    state, "PYT", metrics = c(value = "meanG"), traits = 1:2,
    append_trait = "always"
  )
  index_out <- bp_report_stage_metrics(
    state, "PYT", metrics = c(value = "meanG"), synthetic_trait = "Index",
    append_trait = "always"
  )
  acc_index_out <- bp_report_stage_metrics(
    state, "PYT", metrics = c(value = "accEBV"), synthetic_trait = "Index",
    append_trait = "always"
  )

  expect_equal(names(mean_out), c("value_Trait1", "value_Trait2"))
  expect_equal(names(index_out), "value_Index")
  expect_equal(unlist(mean_out, use.names = FALSE), rep(0, 2), tolerance = 1e-10)
  expect_equal(unlist(index_out, use.names = FALSE), 0, tolerance = 1e-10)
  expect_equal(unlist(acc_index_out, use.names = FALSE), 1, tolerance = 1e-10)
})

test_that("bp_report_stage_metrics does not use synthetic-only baselines for biological traits", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(24, 2, 50)
  SP <- SimParam$new(h)
  SP$addTraitA(10, name = "Trait1")
  SP$addTraitA(10, name = "Trait2")
  pop <- newPop(h, simParam = SP)
  index_def <- bp_synthetic_trait(
    "Index",
    traits = 1:2,
    fun = AlphaSimR::selIndex,
    args = list(b = c(0.5, 0.5)),
    linear = TRUE
  )

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)
  state <- bp_register_synthetic_traits(state, index_def)
  state$sim$trait_baselines$default <- list(
    label = "default",
    traits = "Index",
    mean = c(Index = 100),
    sd = c(Index = 10),
    source = "values"
  )
  state <- BreedingProgramSimulator:::put_stage_pop(state, pop, stage = "PYT", ready_in_years = 0)

  biological <- bp_report_stage_metrics(
    state,
    "PYT",
    metrics = c(value = "meanG"),
    traits = c("Trait1", "Trait2"),
    append_trait = "always"
  )
  synthetic <- bp_report_stage_metrics(
    state,
    "PYT",
    metrics = c(value = "meanG"),
    synthetic_trait = "Index",
    append_trait = "always"
  )
  raw_index <- mean(AlphaSimR::selIndex(pop@gv[, 1:2, drop = FALSE], b = c(0.5, 0.5)))

  expect_equal(
    unlist(biological, use.names = FALSE),
    c(mean(pop@gv[, 1]), mean(pop@gv[, 2]))
  )
  expect_equal(unlist(synthetic, use.names = FALSE), (raw_index - 100) / 10)
})

test_that("bp_report_stage_metrics falls back to positional baseline scaling when names differ", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(16, 2, 40)
  SP <- SimParam$new(h)
  SP$addTraitA(10, name = "Trait1")
  SP$addTraitA(10, name = "Trait2")
  pop <- newPop(h, simParam = SP)

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::put_stage_pop(state, pop, stage = "PYT", ready_in_years = 0)
  state <- BreedingProgramSimulator:::bp_set_trait_baseline(
    state,
    values = list(mean = c(trait1 = 100, trait2 = 200), sd = c(trait1 = 2, trait2 = 4))
  )

  out <- bp_report_stage_metrics(
    state,
    "PYT",
    metrics = c(value = "meanG"),
    traits = 1:2,
    append_trait = "always"
  )

  expect_equal(
    unlist(out, use.names = FALSE),
    c((mean(pop@gv[, 1]) - 100) / 2, (mean(pop@gv[, 2]) - 200) / 4)
  )

  out_trait2 <- bp_report_stage_metrics(
    state,
    "PYT",
    metrics = c(value = "meanG"),
    traits = 2,
    append_trait = "always"
  )

  expect_equal(
    unlist(out_trait2, use.names = FALSE),
    (mean(pop@gv[, 2]) - 200) / 4
  )

  out_character <- bp_report_stage_metrics(
    state,
    "PYT",
    metrics = c(value = "meanG"),
    traits = c("Trait1", "Trait2"),
    append_trait = "always"
  )

  expect_equal(
    unlist(out_character, use.names = FALSE),
    c((mean(pop@gv[, 1]) - 100) / 2, (mean(pop@gv[, 2]) - 200) / 4)
  )

  state$sim$trait_baselines$default$mean <- c(Trait1 = 100, Trait2 = NA_real_, trait2 = 200)
  state$sim$trait_baselines$default$sd <- c(Trait1 = 2, Trait2 = NA_real_, trait2 = 4)
  out_stale_name <- bp_report_stage_metrics(
    state,
    "PYT",
    metrics = c(value = "meanG"),
    traits = c("Trait1", "Trait2"),
    append_trait = "always"
  )

  expect_equal(
    unlist(out_stale_name, use.names = FALSE),
    c((mean(pop@gv[, 1]) - 100) / 2, (mean(pop@gv[, 2]) - 200) / 4)
  )
})

test_that("bp_report_stage_metrics returns stable NA columns when no cohort is available", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_set_trait_baseline(
    state,
    values = list(mean = c(Trait1 = 0, Trait2 = 0), sd = c(Trait1 = 1, Trait2 = 1)),
    covariance = diag(2)
  )

  out <- bp_report_stage_metrics(
    state,
    stage = "PYT",
    metrics = c(meanPYT = "meanG"),
    traits = c("Trait1", "Trait2"),
    append_trait = "always"
  )

  expect_equal(names(out), c("meanPYT_Trait1", "meanPYT_Trait2"))
  expect_true(all(vapply(out, is.na, logical(1))))
})

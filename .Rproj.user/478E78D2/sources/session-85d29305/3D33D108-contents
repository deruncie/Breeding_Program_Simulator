test_that("source policy defaults to latest_one and supports latest_n", {
  state <- list()
  ready <- data.frame(
    cohort_id = c("c1", "c2", "c3"),
    cycle_id = c("cy1", "cy2", "cy3"),
    available_tick = c(8L, 10L, 9L),
    created_tick = c(7L, 9L, 8L),
    stringsAsFactors = FALSE
  )

  out_default <- BreedingProgramSimulator:::bp_select_source_rows(state, ready, list())
  expect_equal(out_default$cohort_id, "c2")

  out_n <- BreedingProgramSimulator:::bp_select_source_rows(
    state,
    ready,
    list(input_policy = "latest_n", input_n = 2)
  )
  expect_equal(out_n$cohort_id, c("c2", "c3"))
})

test_that("run_genotyping is idempotent for same cohort and chip", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0, sim = list(default_chip = 1L))

  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:10),
    stage = "PYT",
    duration_years = 0
  )

  cfg <- list(input_stage = "PYT", chip = 1L, duration_years = 1, cost_per_sample = 5)
  state <- BreedingProgramSimulator:::run_genotyping(state, cfg)
  state <- BreedingProgramSimulator:::run_genotyping(state, cfg)

  expect_equal(nrow(state$genotype_log), 1L)
  expect_equal(sum(state$cost_log$event == "genotyping"), 1L)

  # Different chip should be allowed.
  state <- BreedingProgramSimulator:::run_genotyping(state, modifyList(cfg, list(chip = 2L)))
  expect_equal(nrow(state$genotype_log), 2L)
  expect_equal(sum(state$cost_log$event == "genotyping"), 2L)
})

test_that("run_train_gp_model supports custom train_model_fn", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0, sim = list(default_chip = 1L))
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:7),
    stage = "PYT",
    duration_years = 0
  )
  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  state$genotype_log <- rbind(
    state$genotype_log,
    data.frame(
      cohort_id = cid,
      chip = "1",
      started_tick = 0L,
      done_tick = 0L,
      available_tick = 0L,
      n_ind = 7L,
      stringsAsFactors = FALSE
    )
  )

  cfg <- list(
    from_stage = "PYT",
    chip = 1L,
    trait = 1L,
    lookback_years = 3,
    model_id = "custom_model",
    train_model_fn = function(train_pop, state, cfg) {
      list(n_train = nrow(train_pop), tag = "ok")
    }
  )
  state <- BreedingProgramSimulator:::run_train_gp_model(state, cfg)
  expect_true("custom_model" %in% names(state$gs_models))
  expect_equal(state$gs_models$custom_model$model$n_train, 7)
  expect_equal(state$gs_models$custom_model$model$tag, "ok")
})

test_that("run_train_gp_model wraps train_model_fn errors with context", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0, sim = list(default_chip = 1L))
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:5),
    stage = "PYT",
    duration_years = 0
  )
  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  state$genotype_log <- rbind(
    state$genotype_log,
    data.frame(
      cohort_id = cid,
      chip = "1",
      started_tick = 0L,
      done_tick = 0L,
      available_tick = 0L,
      n_ind = 5L,
      stringsAsFactors = FALSE
    )
  )

  cfg <- list(
    from_stage = "PYT",
    chip = 1L,
    model_id = "bad_model",
    train_model_fn = function(train_pop, state, cfg) {
      stop("boom")
    }
  )
  expect_error(
    BreedingProgramSimulator:::run_train_gp_model(state, cfg),
    "train_model_fn failed for stage 'PYT': boom"
  )
})

test_that("bp_predict_ebv accepts multi-trait matrix output", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(20, 2, 50)
  SP <- SimParam$new(h)
  SP$addTraitA(10)
  SP$addTraitA(10)
  pop <- newPop(h, simParam = SP)

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0, sim = list(default_chip = 1L))
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = pop,
    stage = "TEST",
    duration_years = 0
  )
  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  pop <- state$pops[[cid]]
  state$genotype_log <- rbind(
    state$genotype_log,
    data.frame(
      cohort_id = cid,
      chip = "1",
      started_tick = 0L,
      done_tick = 0L,
      available_tick = 0L,
      n_ind = pop@nInd,
      stringsAsFactors = FALSE
    )
  )
  model_entry <- list(model = list(dummy = TRUE))
  cfg <- list(
    cohort_ids = cid,
    predict_ebv_fn = function(target_pop, model_obj, state, cfg, model_entry) {
      cbind(seq_len(target_pop@nInd), rev(seq_len(target_pop@nInd)))
    }
  )

  pop2 <- BreedingProgramSimulator:::bp_predict_ebv(pop, model_entry, state, cfg, stage_label = "TEST")
  expect_equal(dim(pop2@ebv), c(pop2@nInd, 2L))
  s <- selectInd(pop2, nInd = 3, use = "ebv", trait = 2, simParam = SP)
  expect_equal(s@nInd, 3L)
})

test_that("bp_predict_ebv requires logged genotyping by default", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(12, 2, 40)
  SP <- SimParam$new(h)
  SP$addTraitA(10)
  pop <- newPop(h, simParam = SP)

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0, sim = list(default_chip = 1L))
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = pop,
    stage = "TEST",
    duration_years = 0
  )
  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  pop <- state$pops[[cid]]

  model_entry <- list(model = list(dummy = TRUE), chip = "1")
  cfg <- list(
    cohort_ids = cid,
    predict_ebv_fn = function(target_pop, model_obj, state, cfg, model_entry) {
      seq_len(target_pop@nInd)
    }
  )

  expect_error(
    BreedingProgramSimulator:::bp_predict_ebv(pop, model_entry, state, cfg, stage_label = "TEST"),
    "requires genotyping"
  )
})

test_that("cohort availability and genotyped flags are gated by tick", {
  state <- BreedingProgramSimulator:::bp_init_state(
    SP = NULL,
    dt = 1,
    start_time = 0,
    sim = list(default_chip = 1L)
  )
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:4),
    stage = "F5",
    duration_years = 2
  )
  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  ready0 <- BreedingProgramSimulator:::bp_get_ready_cohorts(state, stage = "F5")
  expect_equal(nrow(ready0), 0L)

  state$genotype_log <- rbind(
    state$genotype_log,
    data.frame(
      cohort_id = cid,
      chip = "1",
      started_tick = 0L,
      done_tick = 2L,
      available_tick = 2L,
      n_ind = 4L,
      stringsAsFactors = FALSE
    )
  )
  state <- BreedingProgramSimulator:::bp_refresh_genotyped_flags(state)
  expect_false(state$cohorts$genotyped[match(cid, state$cohorts$cohort_id)])

  state <- BreedingProgramSimulator:::bp_advance_time(state, n_ticks = 2L)
  ready2 <- BreedingProgramSimulator:::bp_get_ready_cohorts(state, stage = "F5")
  expect_equal(ready2$cohort_id, cid)
  expect_true(state$cohorts$genotyped[match(cid, state$cohorts$cohort_id)])
})

test_that("run_genotyping force allows repeated cohort-chip events", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:6),
    stage = "PYT",
    duration_years = 0
  )

  cfg <- list(input_stage = "PYT", chip = 1L, duration_years = 1, cost_per_sample = 2)
  state <- BreedingProgramSimulator:::run_genotyping(state, cfg)
  state <- BreedingProgramSimulator:::run_genotyping(state, modifyList(cfg, list(force = TRUE)))

  expect_equal(nrow(state$genotype_log), 2L)
  expect_equal(sum(state$cost_log$event == "genotyping"), 2L)
})

test_that("bp_get_training_cohorts filters by lookback window and chip availability", {
  state <- BreedingProgramSimulator:::bp_init_state(
    SP = NULL,
    dt = 1,
    start_time = 0,
    sim = list(default_chip = 1L)
  )

  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:3),
    stage = "PYT",
    duration_years = 0
  )
  old_cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  state <- BreedingProgramSimulator:::bp_advance_time(state, n_ticks = 8L)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:4),
    stage = "PYT",
    duration_years = 0
  )
  new_cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  state$genotype_log <- rbind(
    state$genotype_log,
    data.frame(
      cohort_id = old_cid,
      chip = "1",
      started_tick = 0L,
      done_tick = 0L,
      available_tick = 0L,
      n_ind = 3L,
      stringsAsFactors = FALSE
    ),
    data.frame(
      cohort_id = new_cid,
      chip = "1",
      started_tick = 8L,
      done_tick = 8L,
      available_tick = 8L,
      n_ind = 4L,
      stringsAsFactors = FALSE
    )
  )

  cfg <- list(from_stage = "PYT", chip = 1L, lookback_years = 3)
  out <- BreedingProgramSimulator:::bp_get_training_cohorts(state, cfg)
  expect_equal(out$cohort_id, new_cid)
})

test_that("run_phenotype_trial respects log_per_environment and log_aggregate", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  set.seed(11)
  h <- quickHaplo(12, 2, 50)
  SP <- SimParam$new(h)
  SP$addTraitA(10)

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = newPop(h, simParam = SP),
    stage = "F5",
    duration_years = 0
  )

  cfg_env_only <- list(
    input_stage = "F5",
    output_stage = "PYT",
    traits = 1L,
    n_loc = 3L,
    reps = 2L,
    varE = 1,
    duration_years = 0,
    use_env_control = TRUE,
    log_per_environment = TRUE,
    log_aggregate = FALSE,
    consume_input = FALSE
  )
  state1 <- BreedingProgramSimulator:::run_phenotype_trial(state, cfg_env_only)
  expect_equal(sort(unique(state1$phenotype_log$environment)), c(1L, 2L, 3L))
  expect_equal(
    nrow(state1$phenotype_log),
    BreedingProgramSimulator:::pop_n_ind(state1$pops[[state1$cohorts$cohort_id[2]]]) * 3L
  )

  cfg_agg_only <- modifyList(cfg_env_only, list(log_per_environment = FALSE, log_aggregate = TRUE))
  state2 <- BreedingProgramSimulator:::run_phenotype_trial(state, cfg_agg_only)
  expect_equal(unique(state2$phenotype_log$environment), 0L)
  expect_equal(
    nrow(state2$phenotype_log),
    BreedingProgramSimulator:::pop_n_ind(state2$pops[[state2$cohorts$cohort_id[2]]])
  )
})

test_that("subset copies can inherit genotype log and skip re-genotyping costs", {
  state <- BreedingProgramSimulator:::bp_init_state(
    SP = NULL,
    dt = 1,
    start_time = 0,
    sim = list(default_chip = 1L)
  )
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:10),
    stage = "SRC",
    duration_years = 0
  )
  src <- BreedingProgramSimulator:::get_ready_pop(
    state,
    stage = "SRC",
    combine = TRUE,
    policy = "latest_one",
    silent = TRUE
  )
  src_id <- src$source_ids[[1]]

  state <- BreedingProgramSimulator:::run_genotyping(
    state,
    list(input_stage = "SRC", chip = 1L, duration_years = 0, cost_per_sample = 2)
  )
  n_cost_before <- nrow(state$cost_log)

  pop_sub <- BreedingProgramSimulator:::pop_subset(src$pop, 1:4)
  state <- BreedingProgramSimulator:::put_stage_pop(
    state = state,
    pop = pop_sub,
    stage = "SUB",
    source = src,
    ready_in_years = 0,
    inherit_genotypes = TRUE
  )
  sub_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  expect_true(any(state$genotype_log$cohort_id == src_id & state$genotype_log$chip == "1"))
  expect_true(any(state$genotype_log$cohort_id == sub_id & state$genotype_log$chip == "1"))

  state <- BreedingProgramSimulator:::run_genotyping(
    state,
    list(input_stage = "SUB", chip = 1L, duration_years = 0, cost_per_sample = 2)
  )
  expect_equal(nrow(state$cost_log), n_cost_before)
})

test_that("cross-derived cohorts do not inherit genotype log by default", {
  state <- BreedingProgramSimulator:::bp_init_state(
    SP = NULL,
    dt = 1,
    start_time = 0,
    sim = list(default_chip = 1L)
  )
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:8),
    stage = "SRC",
    duration_years = 0
  )
  src <- BreedingProgramSimulator:::get_ready_pop(
    state,
    stage = "SRC",
    combine = TRUE,
    policy = "latest_one",
    silent = TRUE
  )

  state <- BreedingProgramSimulator:::run_genotyping(
    state,
    list(input_stage = "SRC", chip = 1L, duration_years = 0, cost_per_sample = 2)
  )
  n_cost_before <- nrow(state$cost_log)

  pop_new <- data.frame(v = 101:106)
  state <- BreedingProgramSimulator:::put_stage_pop(
    state = state,
    pop = pop_new,
    stage = "CROSS",
    source = src,
    ready_in_years = 0,
    inherit_genotypes = FALSE
  )
  cross_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  expect_false(any(state$genotype_log$cohort_id == cross_id & state$genotype_log$chip == "1"))

  state <- BreedingProgramSimulator:::run_genotyping(
    state,
    list(input_stage = "CROSS", chip = 1L, duration_years = 0, cost_per_sample = 2)
  )
  expect_equal(nrow(state$cost_log), n_cost_before + 1L)
})

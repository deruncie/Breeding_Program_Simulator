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

test_that("select_latest_available supports n and preserves source ids", {
  state <- BreedingProgramSimulator:::bp_init_state(
    SP = NULL,
    dt = 1,
    start_time = 0,
    sim = list(default_chip = 1L)
  )

  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:2), stage = "PYT", duration_years = 0)
  state <- BreedingProgramSimulator:::bp_advance_time(state, 1L)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 3:4), stage = "PYT", duration_years = 0)
  state <- BreedingProgramSimulator:::bp_advance_time(state, 1L)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 5:6), stage = "PYT", duration_years = 0)

  src1 <- BreedingProgramSimulator:::select_latest_available(state, stage = "PYT", n = 1L, combine = TRUE, silent = TRUE)
  expect_equal(length(src1$source_ids), 1L)
  expect_equal(src1$source_ids[[1]], tail(state$cohorts$cohort_id, 1L))

  src2 <- BreedingProgramSimulator:::select_latest_available(state, stage = "PYT", n = 2L, combine = TRUE, silent = TRUE)
  expect_equal(length(src2$source_ids), 2L)
  expect_equal(nrow(src2$pop), 4L)

  expect_error(
    BreedingProgramSimulator:::select_latest_available(state, stage = "PYT", n = 2L, combine = FALSE, silent = TRUE),
    "combine=FALSE is ambiguous when n > 1"
  )
})

test_that("select_latest_available can include not-ready active cohorts", {
  state <- BreedingProgramSimulator:::bp_init_state(
    SP = NULL,
    dt = 1,
    start_time = 0,
    sim = list(default_chip = 1L)
  )

  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:2), stage = "PYT", duration_years = 2)
  newest_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  expect_null(
    BreedingProgramSimulator:::select_latest_available(
      state,
      stage = "PYT",
      combine = TRUE,
      silent = TRUE
    )
  )

  src <- BreedingProgramSimulator:::select_latest_available(
    state,
    stage = "PYT",
    combine = TRUE,
    include_not_ready = TRUE,
    silent = TRUE
  )
  expect_equal(src$source_ids[[1]], newest_id)
  expect_equal(nrow(src$pop), 2L)
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

test_that("predict_ebv_pop accepts multi-trait matrix output", {
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

  pop2 <- BreedingProgramSimulator:::predict_ebv_pop(pop, model_entry, state, cfg, stage_label = "TEST")
  expect_equal(dim(pop2@ebv), c(pop2@nInd, 2L))
  s <- selectInd(pop2, nInd = 3, use = "ebv", trait = 2, simParam = SP)
  expect_equal(s@nInd, 3L)
})

test_that("predict_ebv_pop requires logged genotyping by default", {
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
    BreedingProgramSimulator:::predict_ebv_pop(pop, model_entry, state, cfg, stage_label = "TEST"),
    "requires genotyping"
  )
})

test_that("run_predict_ebv persists EBVs to matching cohorts in state", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(12, 2, 40)
  SP <- SimParam$new(h)
  SP$addTraitA(10)
  pop <- newPop(h, simParam = SP)

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0, sim = list(default_chip = 1L))
  state <- BreedingProgramSimulator:::bp_add_cohort(state, pop, stage = "TEST", duration_years = 0)
  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)

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

  state$gs_models[["m1"]] <- list(
    model = list(dummy = TRUE),
    chip = "1",
    predict_ebv_fn = function(target_pop, model_obj, state, cfg, model_entry) {
      seq_len(target_pop@nInd)
    }
  )

  state2 <- BreedingProgramSimulator:::run_predict_ebv(
    state,
    list(model_id = "m1", cohort_ids = cid)
  )
  expect_equal(ncol(state2$pops[[cid]]@ebv), 1L)
  expect_false(all(is.na(state2$pops[[cid]]@ebv[, 1])))
})

test_that("bp_log_event collapses multi-stage input to one scalar label", {
  state <- BreedingProgramSimulator:::bp_init_state(
    SP = NULL,
    dt = 1,
    start_time = 0,
    sim = list(default_chip = 1L)
  )

  state2 <- BreedingProgramSimulator:::bp_log_event(
    state = state,
    fn = "x",
    event_type = "y",
    stage = c("PYT", "AYT"),
    source_ids = character(0),
    output_id = "o1",
    event_string = "msg"
  )

  expect_equal(nrow(state2$event_log), 1L)
  expect_equal(state2$event_log$stage[[1]], "PYT;AYT")
})

test_that("bp_log_event scalarizes vector text fields", {
  state <- BreedingProgramSimulator:::bp_init_state(
    SP = NULL,
    dt = 1,
    start_time = 0,
    sim = list(default_chip = 1L)
  )

  state2 <- BreedingProgramSimulator:::bp_log_event(
    state = state,
    fn = c("f1", "f2"),
    event_type = c("e1", "e2"),
    stage = "PYT",
    source_ids = c("c1", "c2"),
    output_id = c("o1", "o2"),
    event_string = c("msg1", "msg2"),
    template_string = c("tpl1", "tpl2")
  )

  expect_equal(nrow(state2$event_log), 1L)
  expect_equal(state2$event_log$fn[[1]], "f1;f2")
  expect_equal(state2$event_log$event_type[[1]], "e1;e2")
  expect_equal(state2$event_log$output_id[[1]], "o1;o2")
  expect_equal(state2$event_log$event_string[[1]], "msg1;msg2")
  expect_equal(state2$event_log$template_string[[1]], "tpl1;tpl2")
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
  expect_equal(sort(unique(state1$phenotype_log$environment)), c("1", "2", "3"))
  expect_equal(
    nrow(state1$phenotype_log),
    BreedingProgramSimulator:::pop_n_ind(state1$pops[[state1$cohorts$cohort_id[2]]]) * 3L
  )

  cfg_agg_only <- modifyList(cfg_env_only, list(log_per_environment = FALSE, log_aggregate = TRUE))
  state2 <- BreedingProgramSimulator:::run_phenotype_trial(state, cfg_agg_only)
  expect_equal(unique(state2$phenotype_log$environment), "0")
  expect_equal(
    nrow(state2$phenotype_log),
    BreedingProgramSimulator:::pop_n_ind(state2$pops[[state2$cohorts$cohort_id[2]]])
  )

  expect_true(any(state2$event_log$event_type == "phenotyping"))
})

test_that("run_phenotype_trial uses effective reps = reps * n_loc when use_env_control is FALSE", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  set.seed(21)
  h <- quickHaplo(10, 2, 50)
  SP <- SimParam$new(h)
  SP$addTraitA(10)

  base_pop <- newPop(h, simParam = SP)
  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = base_pop,
    stage = "F5",
    duration_years = 0
  )

  cfg <- list(
    input_stage = "F5",
    output_stage = "PYT",
    traits = 1L,
    n_loc = 4L,
    reps = 2L,
    varE = 1,
    duration_years = 0,
    use_env_control = FALSE,
    consume_input = FALSE
  )

  set.seed(999)
  state_out <- BreedingProgramSimulator:::run_phenotype_trial(state, cfg)
  out_cid <- tail(state_out$cohorts$cohort_id, 1L)
  pop_out <- state_out$pops[[out_cid]]

  set.seed(999)
  pop_ref <- AlphaSimR::setPheno(
    base_pop,
    varE = 1,
    reps = 8L, # reps * n_loc
    traits = 1L,
    simParam = SP
  )

  expect_equal(as.numeric(pop_out@pheno), as.numeric(pop_ref@pheno))
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

test_that("put_stage_pop accepts source_ids and strategy text", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:5),
    stage = "PYT",
    duration_years = 0
  )
  src_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  state <- BreedingProgramSimulator:::put_stage_pop(
    state = state,
    pop = data.frame(v = 1:3),
    stage = "CROSS_BLOCK",
    source_ids = src_id,
    selection_strategy = "Top 3 by ebv",
    cross_strategy = "random mating",
    ready_in_years = 1
  )
  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  idx <- match(cid, state$cohorts$cohort_id)
  expect_equal(state$cohorts$source_cohort_id[idx], src_id)
  expect_equal(state$cohorts$selection_strategy[idx], "Top 3 by ebv")
  expect_equal(state$cohorts$cross_strategy[idx], "random mating")
  expect_true(any(state$event_log$output_id == cid & state$event_log$event_type == "stage_output"))
})

test_that("put_stage_pop can log per-individual cohort creation cost", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::put_stage_pop(
    state = state,
    pop = data.frame(v = 1:4),
    stage = "F1",
    ready_in_years = 0,
    cost_per_individual = 0.5,
    cost_event = "grow_out",
    cost_unit = "plant"
  )

  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  rows <- subset(state$cost_log, cohort_id == cid & event == "grow_out")
  expect_equal(nrow(rows), 1L)
  expect_equal(rows$stage[[1]], "F1")
  expect_equal(rows$unit[[1]], "plant")
  expect_equal(rows$n_units[[1]], 4)
  expect_equal(rows$unit_cost[[1]], 0.5)
  expect_equal(rows$total_cost[[1]], 2.0)
})

test_that("bp_record_pheno supports sparse matrices and trait names", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)

  ph <- matrix(
    c(4.1, NA, 3.9,
      5.2, 5.1, NA),
    nrow = 3,
    ncol = 2
  )

  state1 <- BreedingProgramSimulator:::bp_record_pheno(
    state = state,
    cohort_id = "cohort_0001",
    stage = "PYT",
    individual_id = c(101L, 102L, 103L),
    traits = c("yield_env1", "yield_env2"),
    pheno_matrix = ph,
    available_tick = 4L,
    n_loc = 2L,
    reps = 1L,
    drop_na = TRUE
  )
  expect_equal(nrow(state1$phenotype_log), 4L)
  expect_true(all(state1$phenotype_log$trait %in% c("yield_env1", "yield_env2")))
  expect_false(any(is.na(state1$phenotype_log$phenotype_value)))

  state2 <- BreedingProgramSimulator:::bp_record_pheno(
    state = state,
    cohort_id = "cohort_0001",
    stage = "PYT",
    individual_id = c(101L, 102L, 103L),
    traits = c("yield_env1", "yield_env2"),
    pheno_matrix = ph,
    available_tick = 4L,
    n_loc = 2L,
    reps = 1L,
    drop_na = FALSE
  )
  expect_equal(nrow(state2$phenotype_log), 6L)
  expect_true(any(is.na(state2$phenotype_log$phenotype_value)))

  ph_named <- ph
  colnames(ph_named) <- c("yield_env1", "yield_env2")
  state3 <- BreedingProgramSimulator:::bp_record_pheno(
    state = state,
    cohort_id = "cohort_0001",
    stage = "PYT",
    individual_id = c(101L, 102L, 103L),
    traits = NULL,
    pheno_matrix = ph_named,
    available_tick = 4L,
    n_loc = 2L,
    reps = 1L,
    drop_na = TRUE
  )
  expect_true(all(state3$phenotype_log$trait %in% c("yield_env1", "yield_env2")))

  ph9 <- matrix(seq_len(18), nrow = 2, ncol = 9)
  state4 <- BreedingProgramSimulator:::bp_record_pheno(
    state = state,
    cohort_id = "cohort_0001",
    stage = "PYT",
    individual_id = c(201L, 202L),
    traits = rep(c("T1", "T2", "T3"), 3),
    pheno_matrix = ph9,
    available_tick = 4L,
    n_loc = 3L,
    reps = 1L,
    environment = rep(c("E1", "E2", "E3"), each = 3),
    p_value = rep(c(0.1, 0.2, 0.3), each = 3),
    drop_na = TRUE
  )
  expect_true(all(state4$phenotype_log$trait %in% c("T1", "T2", "T3")))
  expect_true(all(state4$phenotype_log$environment %in% c("E1", "E2", "E3")))
  expect_true(all(state4$phenotype_log$p_value %in% c(0.1, 0.2, 0.3)))

  expect_error(
    BreedingProgramSimulator:::bp_record_pheno(
      state = state,
      cohort_id = "cohort_0001",
      stage = "PYT",
      individual_id = c(101L, 102L, 103L),
      traits = c("t1", "t2"),
      pheno_matrix = matrix(1:9, nrow = 3, ncol = 3),
      available_tick = 4L,
      n_loc = 1L,
      reps = 1L,
      environment = c("E1", "E2", "E3"),
      drop_na = TRUE
    ),
    "not compatible with ncol\\(pheno_matrix\\)"
  )
})

test_that("run_phenotype_trial converts latent environment values to AlphaSimR p scale", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  founder <- quickHaplo(nInd = 20, nChr = 1, segSites = 50)
  SP <- SimParam$new(founder)
  SP$addTraitAEG(10, varGxE = 1, varEnv = 4)

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)
  pop <- newPop(founder, simParam = SP)
  state <- BreedingProgramSimulator:::put_stage_pop(state, pop, stage = "F5", ready_in_years = 0)
  src <- BreedingProgramSimulator:::select_latest_available(state, stage = "F5", combine = TRUE, silent = TRUE)

  state <- BreedingProgramSimulator:::run_phenotype_trial(
    state = state,
    pop = src$pop,
    output_stage = "PYT",
    input_cohorts = src$source_ids,
    selection_strategy = "all",
    traits = 1L,
    n_loc = 3L,
    reps = 1L,
    varE = 1,
    duration_years = 1,
    use_env_control = TRUE,
    env_means = c(0, 1, -1),
    env_year_sd = 0,
    log_per_environment = TRUE,
    log_aggregate = TRUE
  )

  ph <- subset(state$phenotype_log, stage == "PYT" & environment %in% c("1", "2", "3"))
  pvals <- ph$p_value[!is.na(ph$p_value)]
  expect_true(length(pvals) > 0L)
  expect_true(all(pvals > 0 & pvals < 1))
  expect_equal(
    sort(unique(round(pvals, 8))),
    sort(round(stats::pnorm(c(0, 1, -1)), 8))
  )
})

test_that("run_phenotype_trial does not cache environments by trial_name", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  founder <- quickHaplo(nInd = 12, nChr = 1, segSites = 30)
  SP <- SimParam$new(founder)
  SP$addTraitAEG(10, varGxE = 1, varEnv = 1)

  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)
  pop <- newPop(founder, simParam = SP)

  state <- BreedingProgramSimulator:::run_phenotype_trial(
    state = state,
    pop = pop,
    output_stage = "AYT",
    trial_name = "AYT",
    traits = 1L,
    n_loc = 2L,
    reps = 1L,
    varE = 1,
    duration_years = 0,
    use_env_control = TRUE,
    env_means = c(-1, 1),
    env_mean_sd = 0,
    env_year_sd = 0,
    log_per_environment = TRUE,
    log_aggregate = FALSE
  )
  p1 <- subset(state$phenotype_log, stage == "AYT" & environment %in% c("1", "2"))$p_value

  state <- BreedingProgramSimulator:::run_phenotype_trial(
    state = state,
    pop = pop,
    output_stage = "AYT2",
    trial_name = "AYT",
    traits = 1L,
    n_loc = 2L,
    reps = 1L,
    varE = 1,
    duration_years = 0,
    use_env_control = TRUE,
    env_means = c(0, 0),
    env_mean_sd = 0,
    env_year_sd = 0,
    log_per_environment = TRUE,
    log_aggregate = FALSE
  )
  p2 <- subset(state$phenotype_log, stage == "AYT2" & environment %in% c("1", "2"))$p_value

  expect_equal(sort(unique(round(p1, 8))), sort(round(stats::pnorm(c(-1, 1)), 8)))
  expect_equal(sort(unique(round(p2, 8))), rep(round(stats::pnorm(0), 8), 1L))
})

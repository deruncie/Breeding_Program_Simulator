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
  model_entry <- list(model = list(dummy = TRUE))
  cfg <- list(
    predict_ebv_fn = function(target_pop, model_obj, state, cfg, model_entry) {
      cbind(seq_len(target_pop@nInd), rev(seq_len(target_pop@nInd)))
    }
  )

  pop2 <- BreedingProgramSimulator:::bp_predict_ebv(pop, model_entry, state, cfg, stage_label = "TEST")
  expect_equal(dim(pop2@ebv), c(pop2@nInd, 2L))
  s <- selectInd(pop2, nInd = 3, use = "ebv", trait = 2, simParam = SP)
  expect_equal(s@nInd, 3L)
})

test_that("synthetic traits evaluate AlphaSimR selIndex without changing native traits", {
  testthat::skip_if_not_installed("AlphaSimR")

  set.seed(101)
  h <- AlphaSimR::quickHaplo(24, 2, 50)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(10)
  SP$addTraitA(10)
  pop <- AlphaSimR::newPop(h, simParam = SP)

  index <- bp_synthetic_trait(
    name = "Index",
    traits = 1:2,
    fun = AlphaSimR::selIndex,
    args = list(b = c(1, -0.25)),
    linear = TRUE
  )

  expect_equal(
    bp_synthetic_values(pop, index, use = "gv"),
    as.numeric(AlphaSimR::selIndex(pop@gv[, 1:2, drop = FALSE], b = c(1, -0.25)))
  )
  expect_equal(ncol(pop@gv), 2L)
  expect_equal(ncol(pop@ebv), 0L)
})

test_that("synthetic traits can zero components that are entirely missing", {
  testthat::skip_if_not_installed("AlphaSimR")

  set.seed(111)
  h <- AlphaSimR::quickHaplo(12, 2, 40)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(10)
  SP$addTraitA(10)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  index <- bp_synthetic_trait(
    name = "Index",
    traits = 1:2,
    fun = AlphaSimR::selIndex,
    args = list(b = c(1, 1)),
    linear = TRUE,
    missing_component = "zero_if_all_missing"
  )

  pop@ebv <- matrix(numeric(0), nrow = pop@nInd, ncol = 0L)
  expect_equal(bp_synthetic_values(pop, index, use = "ebv"), rep(0, pop@nInd))

  pop@pheno <- pop@gv
  pop@pheno[, 2L] <- NA_real_
  expect_equal(
    bp_synthetic_values(pop, index, use = "pheno"),
    as.numeric(AlphaSimR::selIndex(pop@pheno[, 1L, drop = FALSE], b = 1))
  )
})

test_that("synthetic traits propagate missing components by default", {
  testthat::skip_if_not_installed("AlphaSimR")

  h <- AlphaSimR::quickHaplo(12, 2, 40)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(10)
  SP$addTraitA(10)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  index <- bp_synthetic_trait(
    name = "Index",
    traits = 1:2,
    fun = AlphaSimR::selIndex,
    args = list(b = c(1, 1)),
    linear = TRUE
  )

  pop@ebv <- matrix(numeric(0), nrow = pop@nInd, ncol = 0L)
  expect_true(all(is.na(bp_synthetic_values(pop, index, use = "ebv"))))
})

test_that("synthetic trial scoring can zero unmeasured component traits", {
  testthat::skip_if_not_installed("AlphaSimR")

  h <- AlphaSimR::quickHaplo(12, 2, 40)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(10)
  SP$addTraitA(10)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  pop@pheno <- pop@gv[, 1L, drop = FALSE]
  index <- bp_synthetic_trait(
    name = "Index",
    traits = 1:2,
    fun = AlphaSimR::selIndex,
    args = list(b = c(1, 1)),
    linear = TRUE,
    missing_component = "zero_if_all_missing"
  )

  scored <- BreedingProgramSimulator:::bp_score_synthetic_trial(
    pop,
    definitions = index,
    measured_traits = 1L
  )

  expect_equal(
    scored$aggregate[, "Index"],
    as.numeric(AlphaSimR::selIndex(pop@pheno, b = 1))
  )
})

test_that("stored synthetic values survive AlphaSimR population subsets", {
  testthat::skip_if_not_installed("AlphaSimR")

  h <- AlphaSimR::quickHaplo(12, 2, 30)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(8)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  pop <- bp_set_synthetic_values(pop, "utility", seq_len(pop@nInd), type = "ebv")

  sub <- pop[c(2L, 5L, 9L)]
  expect_equal(
    bp_get_stored_synthetic_values(sub, "utility", type = "ebv"),
    c(2, 5, 9)
  )
})

test_that("nonlinear synthetic GV uses reproducible common environmental draws", {
  testthat::skip_if_not_installed("AlphaSimR")

  set.seed(102)
  h <- AlphaSimR::quickHaplo(18, 2, 40)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(10)
  SP$addTraitA(10)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  state <- bp_init_state(SP = SP, dt = 1)
  utility <- bp_synthetic_trait(
    "utility",
    traits = 1:2,
    fun = function(x) x[, 1L]^2 + x[, 2L],
    linear = FALSE
  )

  a <- bp_get_synthetic_gv(
    pop,
    utility,
    state,
    varE = diag(c(0.5, 0.25)),
    n_trials = 8L,
    n_plants_per_trial = 4L,
    seed = 44
  )
  b <- bp_get_synthetic_gv(
    pop,
    utility,
    state,
    varE = diag(c(0.5, 0.25)),
    n_trials = 8L,
    n_plants_per_trial = 4L,
    seed = 44
  )

  expect_equal(a, b)
  expect_length(a, pop@nInd)
  expect_true(all(is.finite(a)))
})

test_that("synthetic GV can be materialized and reused by multi-metric reporting", {
  testthat::skip_if_not_installed("AlphaSimR")

  set.seed(104)
  h <- AlphaSimR::quickHaplo(20, 2, 40)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(10)
  SP$addTraitA(10)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  utility <- bp_synthetic_trait(
    "utility",
    traits = 1:2,
    fun = function(x) x[, 1L]^2 + x[, 2L]
  )
  state <- bp_init_state(SP = SP, dt = 1)
  state <- bp_register_synthetic_traits(state, utility)
  pop <- bp_materialize_synthetic_gv(
    pop,
    "utility",
    state,
    varE = diag(2),
    n_trials = 5L,
    n_plants_per_trial = 3L,
    seed = 19
  )
  cached <- bp_get_stored_synthetic_values(pop, "utility", type = "gv")
  state <- bp_add_cohort(state, pop, stage = "PYT", duration_years = 0)

  out <- bp_report_stage_metrics(
    state,
    stage = "PYT",
    metrics = c("meanG", "varG"),
    synthetic_trait = "utility",
    synthetic_varE = diag(2),
    synthetic_gv_n_trials = 1L,
    synthetic_gv_n_plants = 1L,
    synthetic_gv_seed = 999
  )

  expect_equal(unname(out[["meanG"]]), mean(cached))
  expect_equal(unname(out[["varG"]]), stats::var(cached))
})

test_that("report naming separates biological traits from a synthetic index", {
  testthat::skip_if_not_installed("AlphaSimR")

  h <- AlphaSimR::quickHaplo(12, 2, 30)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(8)
  SP$addTraitA(8)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  index <- bp_synthetic_trait(
    "Index",
    traits = 1:2,
    fun = AlphaSimR::selIndex,
    args = list(b = c(0.5, 0.5)),
    linear = TRUE
  )
  state <- bp_init_state(SP = SP, dt = 1)
  state <- bp_register_synthetic_traits(state, index)
  state <- bp_add_cohort(state, pop, stage = "Candidates", duration_years = 0)

  biological <- bp_report_stage_metrics(
    state,
    stage = "Candidates",
    metrics = c(meanCandidates = "meanG"),
    traits = 1:2,
    append_trait = "always"
  )
  synthetic <- bp_report_stage_metrics(
    state,
    stage = "Candidates",
    metrics = c(meanCandidates = "meanG"),
    synthetic_trait = "Index",
    append_trait = "always"
  )

  expect_equal(names(biological), c("meanCandidates_Trait1", "meanCandidates_Trait2"))
  expect_equal(names(synthetic), "meanCandidates_Index")
  expect_true(all(vapply(c(biological, synthetic), function(x) is.numeric(x) && length(x) == 1L && is.null(names(x)), logical(1))))
})

test_that("phenotype trials cache and log synthetic phenotypes", {
  testthat::skip_if_not_installed("AlphaSimR")

  set.seed(103)
  h <- AlphaSimR::quickHaplo(16, 2, 40)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(10)
  SP$addTraitA(10)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  state <- bp_init_state(SP = SP, dt = 1)
  index <- bp_synthetic_trait(
    "Index",
    traits = 1:2,
    fun = AlphaSimR::selIndex,
    args = list(b = c(0.7, 0.3)),
    linear = TRUE
  )

  state <- run_phenotype_trial(
    state = state,
    pop = pop,
    output_stage = "PYT",
    traits = 1:2,
    n_loc = 2L,
    reps = 1L,
    varE = diag(2),
    duration_years = 0,
    use_env_control = FALSE,
    synthetic_traits = index,
    log_aggregate = TRUE
  )
  out <- state$pops[[bp_last_cohort_id(state)]]
  cached <- bp_get_stored_synthetic_values(out, "Index", type = "pheno")

  expect_equal(cached, as.numeric(AlphaSimR::selIndex(out@pheno, b = c(0.7, 0.3))))
  expect_true(any(state$phenotype_log$trait == "synthetic:Index"))
})

test_that("selection can use direct synthetic EBVs without modifying pop ebv", {
  testthat::skip_if_not_installed("AlphaSimR")

  h <- AlphaSimR::quickHaplo(20, 2, 30)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(8)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  utility <- bp_synthetic_trait("utility", 1L, function(x) x[, 1L]^2)
  score <- rev(seq_len(pop@nInd))
  pop <- bp_set_synthetic_values(pop, "utility", score, type = "ebv")

  selected <- bp_select_synthetic(
    pop,
    n_select = 4L,
    synthetic_trait = utility,
    use = "ebv",
    simParam = SP
  )

  expect_equal(selected@id, pop@id[order(score, decreasing = TRUE)[1:4]])
  expect_equal(ncol(pop@ebv), 0L)
})

test_that("custom synthetic models store predictions in synthetic ebv", {
  testthat::skip_if_not_installed("AlphaSimR")

  h <- AlphaSimR::quickHaplo(10, 2, 30)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(8)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  state <- bp_init_state(SP = SP, dt = 1)
  model_entry <- list(
    model = list(dummy = TRUE),
    response_type = "synthetic",
    synthetic_trait = "utility",
    predict_ebv_fn = function(target_pop, model_obj, state, cfg, model_entry) {
      seq_len(target_pop@nInd)
    }
  )

  out <- predict_ebv_pop(
    pop,
    model_entry,
    state,
    cfg = list(require_genotyped = FALSE),
    stage_label = "TEST"
  )

  expect_equal(
    bp_get_stored_synthetic_values(out, "utility", type = "ebv"),
    seq_len(pop@nInd)
  )
  expect_equal(ncol(out@ebv), 0L)
})

test_that("AlphaSimR RRBLUP trains and predicts a synthetic phenotype", {
  testthat::skip_if_not_installed("AlphaSimR")

  set.seed(105)
  h <- AlphaSimR::quickHaplo(50, 2, 80)
  SP <- AlphaSimR::SimParam$new(h)
  SP$addTraitA(15)
  SP$addTraitA(15)
  SP$addSnpChip(20)
  pop <- AlphaSimR::newPop(h, simParam = SP)
  utility <- bp_synthetic_trait(
    "utility",
    traits = 1:2,
    fun = function(x) x[, 1L]^2 + 0.5 * x[, 2L]
  )
  state <- bp_init_state(SP = SP, dt = 1, sim = list(default_chip = 1L))
  state <- bp_register_synthetic_traits(state, utility)
  state <- run_phenotype_trial(
    state,
    pop = pop,
    output_stage = "PYT",
    traits = 1:2,
    n_loc = 1L,
    reps = 1L,
    varE = diag(2),
    duration_years = 0,
    synthetic_traits = "utility"
  )
  cid <- bp_last_cohort_id(state)
  state <- run_genotyping(
    state,
    list(cohort_ids = cid, chip = 1L, duration_years = 0, cost_per_sample = 0)
  )
  state <- run_train_gp_model(
    state,
    list(
      from_stage = "PYT",
      chip = 1L,
      response = "synthetic_pheno",
      synthetic_trait = "utility",
      model_id = "synthetic_rrblup"
    )
  )
  state <- run_predict_ebv(
    state,
    list(cohort_ids = cid, model_id = "synthetic_rrblup", chip = 1L)
  )
  predicted <- bp_get_stored_synthetic_values(
    state$pops[[cid]],
    "utility",
    type = "ebv"
  )

  expect_length(predicted, pop@nInd)
  expect_true(all(is.finite(predicted)))
  expect_equal(ncol(state$pops[[cid]]@ebv), 0L)
})

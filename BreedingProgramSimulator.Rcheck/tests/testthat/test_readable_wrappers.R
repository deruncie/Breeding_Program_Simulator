test_that("run_selfing_stage supports configurable generation depth", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(10, 2, 50)
  SP <- SimParam$new(h)
  SP$addTraitA(10)
  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)

  p <- newPop(h, simParam = SP)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, p, stage = "F1", duration_years = 0)

  cfg <- list(
    input_stage = "F1",
    output_stage = "F3",
    n_generations = 2,
    n_progeny_schedule = c(2, 1),
    duration_years = 0,
    consume_input = TRUE
  )
  state <- BreedingProgramSimulator:::run_selfing_stage(state, cfg)

  ready <- BreedingProgramSimulator:::bp_get_ready_cohorts(state, stage = "F3", active_only = TRUE)
  expect_true(nrow(ready) >= 1)
})

test_that("run_crossing_stage accepts custom cross_plan_fn", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(20, 2, 50)
  SP <- SimParam$new(h)
  SP$addTraitA(10)
  state <- BreedingProgramSimulator:::bp_init_state(SP = SP, dt = 1, start_time = 0)

  p <- newPop(h, simParam = SP)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, p, stage = "PARENT", duration_years = 0)

  cfg <- list(
    input_stage = "PARENT",
    output_stage = "F1",
    n_crosses = 5,
    n_progeny_per_cross = 1,
    duration_years = 0,
    cross_plan_fn = function(state, src, parent_pop, cfg) {
      matrix(rep(1:2, cfg$n_crosses), ncol = 2, byrow = TRUE)
    }
  )

  state <- BreedingProgramSimulator:::run_crossing_stage(state, cfg)
  ready <- BreedingProgramSimulator:::bp_get_ready_cohorts(state, stage = "F1", active_only = TRUE)
  expect_true(nrow(ready) >= 1)
})

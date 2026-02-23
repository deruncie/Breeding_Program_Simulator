test_that("get_ready_pop combines multiple cohorts", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:3), stage = "A", duration_years = 0, cycle_id = "c1")
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:2), stage = "A", duration_years = 0, cycle_id = "c2")

  src <- BreedingProgramSimulator:::get_ready_pop(state, stage = "A", combine = TRUE, policy = "all_ready", silent = TRUE)
  expect_false(is.null(src))
  expect_equal(nrow(src$pop), 5)
  expect_equal(length(src$source_ids), 2)
})

test_that("put_stage_pop + add_stage_cost + close_sources workflow", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(state, data.frame(v = 1:4), stage = "PARENT", duration_years = 0)
  src <- BreedingProgramSimulator:::get_ready_pop(state, stage = "PARENT", combine = TRUE, policy = "latest_one", silent = TRUE)

  state <- BreedingProgramSimulator:::put_stage_pop(state, src$pop, stage = "F1", source = src, ready_in_years = 1)
  cid <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  state <- BreedingProgramSimulator:::add_stage_cost(state, event = "crossing", n_units = 10, unit_cost = 2, stage = "F1", cohort_id = cid)
  state <- BreedingProgramSimulator:::close_sources(state, src)

  expect_true(any(state$cohorts$cohort_id == cid))
  expect_equal(sum(state$cost_log$total_cost), 20)
  expect_false(state$cohorts$active[match(src$source_ids[1], state$cohorts$cohort_id)])
})

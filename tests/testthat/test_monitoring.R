test_that("bp_monitor_cohorts traces origin lineage to requested stage", {
  state <- BreedingProgramSimulator:::bp_init_state(SP = NULL, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:2),
    stage = "PARENT",
    duration_years = 0
  )
  c1 <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  state <- BreedingProgramSimulator:::bp_advance_time(state, 1L)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:2),
    stage = "F1",
    source_cohort_id = c1,
    duration_years = 0
  )
  c2 <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  state <- BreedingProgramSimulator:::bp_advance_time(state, 1L)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state = state,
    pop = data.frame(v = 1:2),
    stage = "F5",
    source_cohort_id = c2,
    duration_years = 0
  )
  c3 <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  out <- BreedingProgramSimulator:::bp_monitor_cohorts(
    state,
    ticks_per_year = 1L,
    origin_stage = "PARENT"
  )
  row3 <- out[out$cohort_id == c3, , drop = FALSE]
  expect_equal(row3$origin_cohort_id, c1)
  expect_equal(row3$origin_stage, "PARENT")
  expect_equal(row3$origin_year, 1L)
})

test_that("bp_summarize_metric_by_year aggregates by year and stage", {
  metrics_df <- data.frame(
    available_year = c(1L, 1L, 2L),
    stage = c("F5", "F5", "PYT"),
    mean_gv = c(1, 3, 5),
    stringsAsFactors = FALSE
  )
  out <- BreedingProgramSimulator:::bp_summarize_metric_by_year(
    metrics_df,
    metric = "mean_gv",
    year_col = "available_year",
    stage_col = "stage"
  )

  out <- out[order(out$year, out$stage), , drop = FALSE]
  expect_equal(nrow(out), 2L)
  expect_equal(out$value[out$stage == "F5"], 2)
  expect_equal(out$value[out$stage == "PYT"], 5)
})

test_that("bp_resolve_trait accepts numeric or label and rejects invalid values", {
  num <- BreedingProgramSimulator:::bp_resolve_trait(2L)
  expect_equal(num$index, 2L)
  expect_equal(num$label, "trait2")

  lbl <- BreedingProgramSimulator:::bp_resolve_trait("trait3")
  expect_equal(lbl$index, 3L)
  expect_equal(lbl$label, "trait3")

  expect_error(BreedingProgramSimulator:::bp_resolve_trait(0), "trait must be >= 1")
  expect_error(
    BreedingProgramSimulator:::bp_resolve_trait("yield"),
    "trait must be numeric index"
  )
})

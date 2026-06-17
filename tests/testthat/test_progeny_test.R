test_that("run_progeny_test stores self-family means on candidate genotypes", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  set.seed(101)
  h <- quickHaplo(6, 2, 40)
  SP <- SimParam$new(h)
  SP$addTraitA(10)
  candidates <- newPop(h, simParam = SP)

  state <- bp_init_state(SP = SP, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state, candidates, stage = "F2", duration_years = 0
  )
  source_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  set.seed(202)
  out <- run_progeny_test(
    state = state,
    pop = candidates,
    output_stage = "F3",
    input_cohorts = source_id,
    mating = "self",
    n_progeny = 4L,
    traits = 1L,
    n_loc = 1L,
    reps = 2L,
    varE = 1,
    duration_years = 1,
    cost_per_plot = 2,
    chunk_size = 100L
  )

  set.seed(202)
  p_trial <- stats::runif(1L)
  progeny <- AlphaSimR::self(
    candidates, nProgeny = 4L, keepParents = FALSE, simParam = SP
  )
  progeny_genetic <- AlphaSimR::setPheno(
    progeny, varE = 0, reps = 1, traits = 1L, p = p_trial,
    onlyPheno = TRUE, simParam = SP
  )
  expected_genetic <- rowsum(
    progeny_genetic,
    group = rep(seq_len(candidates@nInd), each = 4L),
    reorder = FALSE
  ) / 4
  expected <- expected_genetic + stats::rnorm(candidates@nInd, sd = sqrt(1 / 2))

  output_id <- BreedingProgramSimulator:::bp_last_cohort_id(out)
  evaluated <- out$pops[[output_id]]
  expect_equal(evaluated@id, candidates@id)
  expect_equal(evaluated@geno, candidates@geno)
  expect_equal(evaluated@pheno[, 1L], as.numeric(expected[, 1L]))
  expect_equal(out$cohorts$stage[out$cohorts$cohort_id == output_id], "F3")
  expect_equal(out$cohorts$n_ind[out$cohorts$cohort_id == output_id], candidates@nInd)
  expect_equal(out$cost_log$n_units[out$cost_log$cohort_id == output_id], 6 * 4 * 2)
  expect_equal(out$cost_log$total_cost[out$cost_log$cohort_id == output_id], 6 * 4 * 2 * 2)
  expect_equal(nrow(subset(out$phenotype_log, cohort_id == output_id)), 6L)
  expect_true(any(out$event_log$fn == "run_progeny_test"))
  details <- out$event_log$details[[which(out$event_log$fn == "run_progeny_test")]]
  expect_equal(details$residual_scale, "candidate_family_mean")
})

test_that("run_progeny_test aggregates all tester families and records tester lineage", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  set.seed(303)
  h <- quickHaplo(7, 2, 40)
  SP <- SimParam$new(h)
  SP$addTraitA(10)
  base <- newPop(h, simParam = SP)
  candidates <- base[1:4]
  testers <- base[5:7]

  state <- bp_init_state(SP = SP, dt = 1, start_time = 0)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state, candidates, stage = "F3_candidates", duration_years = 0
  )
  candidate_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)
  state <- BreedingProgramSimulator:::bp_add_cohort(
    state, testers, stage = "Tester", duration_years = 0
  )
  tester_id <- BreedingProgramSimulator:::bp_last_cohort_id(state)

  set.seed(404)
  out <- run_progeny_test(
    state = state,
    pop = candidates,
    output_stage = "F3",
    input_cohorts = candidate_id,
    mating = "testcross",
    n_progeny = 2L,
    tester = testers,
    tester_cohorts = tester_id,
    traits = 1L,
    n_loc = 1L,
    reps = 1L,
    varE = 0,
    duration_years = 0,
    chunk_size = 100L
  )

  cross_plan <- cbind(
    rep(seq_len(candidates@nInd), each = testers@nInd),
    rep(seq_len(testers@nInd), times = candidates@nInd)
  )
  set.seed(404)
  p_trial <- stats::runif(1L)
  progeny <- AlphaSimR::makeCross2(
    females = candidates,
    males = testers,
    crossPlan = cross_plan,
    nProgeny = 2L,
    simParam = SP
  )
  progeny_pheno <- AlphaSimR::setPheno(
    progeny, varE = 0, reps = 1, traits = 1L, p = p_trial,
    onlyPheno = TRUE, simParam = SP
  )
  family_index <- rep(cross_plan[, 1L], each = 2L)
  expected <- rowsum(progeny_pheno, family_index, reorder = FALSE) /
    tabulate(family_index, nbins = candidates@nInd)

  output_id <- BreedingProgramSimulator:::bp_last_cohort_id(out)
  evaluated <- out$pops[[output_id]]
  expect_equal(evaluated@id, candidates@id)
  expect_equal(evaluated@pheno[, 1L], as.numeric(expected[, 1L]))
  expect_equal(
    out$cohorts$source_cohort_id[out$cohorts$cohort_id == output_id],
    paste(candidate_id, tester_id, sep = ";")
  )
  expect_equal(
    out$cost_log$n_units[out$cost_log$cohort_id == output_id],
    candidates@nInd * testers@nInd * 2L
  )
})

test_that("run_progeny_test validates testcross inputs", {
  testthat::skip_if_not_installed("AlphaSimR")
  library(AlphaSimR)

  h <- quickHaplo(2, 1, 20)
  SP <- SimParam$new(h)
  SP$addTraitA(5)
  candidates <- newPop(h, simParam = SP)
  state <- bp_init_state(SP = SP, dt = 1, start_time = 0)

  expect_error(
    run_progeny_test(
      state, candidates, "F3", mating = "testcross",
      n_progeny = 2L, varE = 1
    ),
    "tester must be an AlphaSimR Pop"
  )
  expect_error(
    run_progeny_test(
      state, candidates, "F3", mating = "self",
      n_progeny = 0L, varE = 1
    ),
    "n_progeny must be a positive whole number"
  )
})

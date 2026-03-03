# Two-part breeding program draft, warm-started from Wheat_Conv_noGS.
#
# Draft intent:
# - Use a warm-state from wheat_conv_noGS (default 20 years) as initialization.
# - Population Improvement (PI): two GS-driven half-year crossing cycles per year.
# - Product Development (PD): DH -> headrow -> PYT -> AYT -> EYT1 -> EYT2 -> Variety.
# - Keep assumptions explicit and unresolved items listed for user confirmation.

# Setup -------------------------------------------------------------------
library(AlphaSimR)

source("R/operators.R")
source("R/readable_api.R")
source("R/state_print.R")
source("R/monitoring.R")


# Filled-From-Input Table -------------------------------------------------
make_two_part_filled_table <- function() {
  data.frame(
    parameter = c(
      "warm_start_scheme",
      "warm_start_years",
      "pi_cycles_per_year",
      "pi_selected_males_per_cycle",
      "pi_selected_females_per_cycle",
      "seeds_per_cross",
      "seeds_to_next_pi_per_family",
      "seeds_to_dh_per_family",
      "pd_dh_lines_per_family",
      "pd_no_cross_block_recycling",
      "pd_training_data_use",
      "crossing_cycle_duration_years",
      "dh_production_frequency",
      "network_rendering"
    ),
    value = c(
      "wheat_conv_noGS state",
      "20",
      "2",
      "100",
      "100",
      "40",
      "30",
      "10",
      "31",
      "TRUE",
      "Add genotype+phenotype data to GS training",
      "0.5",
      "once per year (year start)",
      "disabled for now"
    ),
    source = c(
      "user text",
      "user text",
      "user text",
      "user text",
      "user text",
      "user text",
      "user text",
      "user text",
      "user text",
      "user text",
      "user text",
      "user text",
      "user text",
      "user text"
    ),
    stringsAsFactors = FALSE
  )
}


# Missing-Specification Table ---------------------------------------------
make_two_part_missing_table <- function() {
  data.frame(
    parameter = c(
      "PI candidate source at warm-start",
      "PI selection trait/use",
      "Pollination male multiplicity",
      "Whether one PI cycle contributes to one or two PD DH batches",
      "PD trial heritability targets (h2)",
      "PD environment structure (n_loc/reps by stage)",
      "Genotyping timing detail for PI updates",
      "Cost model parameters",
      "Stage-consumption policy for legacy warm-state cohorts"
    ),
    current_default = c(
      "latest PYT + latest CROSS_BLOCK (merged)",
      "both sexes selected as top-EBV (trait 1)",
      "random male assignment to selected females",
      "each PI cycle produces one PD DH input batch",
      "headrow=0.10, pyt=0.20, ayt=0.50, eyt=0.67",
      "headrow 1x1, pyt 1x1, ayt 4x1, eyt 8x1",
      "PI genotyping duration=0 (instant)",
      "placeholder costs copied from readable examples",
      "do not close warm-state cohorts unless consumed by this draft flow"
    ),
    why_it_matters = c(
      "Defines initial genetic base and continuity with prior program.",
      "Changes selection intensity and expected gain each half-cycle.",
      "Affects family structure and within-family variance.",
      "Controls annual DH throughput and downstream trial load.",
      "Sets residual variance and trial discrimination.",
      "Sets effective selection pressure and phenotyping cost.",
      "Controls whether model refresh can impact cycle 2 in-year.",
      "Needed for realistic economic comparisons.",
      "Prevents accidental double-consumption or cohort leakage."
    ),
    user_decision_needed = c(
      "confirm/replace",
      "confirm/replace",
      "confirm/replace",
      "confirm/replace",
      "confirm/replace",
      "confirm/replace",
      "confirm/replace",
      "confirm/replace",
      "confirm/replace"
    ),
    stringsAsFactors = FALSE
  )
}


# Design Questions ---------------------------------------------------------
make_two_part_design_questions <- function() {
  c(
    "Confirm PI warm-start merge should remain strictly PYT + CROSS_BLOCK only.",
    "Confirm PI parent selection remains top-EBV for both sexes across all years.",
    "Confirm random male-to-female assignment is sufficient for open pollination.",
    "Confirm annual DH production should consume all ready PD_DH_INPUT cohorts once per year.",
    "Confirm GS training should remain PYT-only.",
    "Confirm network rendering stays disabled until behavioral sign-off."
  )
}


# Helper Utilities ---------------------------------------------------------
varE_from_base_h2 <- function(base_pop, h2, trait = 1L) {
  h2 <- as.numeric(h2)
  if (!is.finite(h2) || h2 <= 0 || h2 >= 1) return(1.0)

  gv <- base_pop@gv
  if (is.null(dim(gv))) gv <- matrix(gv, ncol = 1L)
  vg <- stats::var(gv[, as.integer(trait)], na.rm = TRUE)
  if (is.na(vg) || vg <= 0) return(1.0)
  as.numeric(vg * (1 - h2) / h2)
}

resolve_two_part_varE <- function(cfg, base_pop) {
  cfg$varE_headrow <- varE_from_base_h2(base_pop, cfg$h2_headrow)
  cfg$varE_pyt <- varE_from_base_h2(base_pop, cfg$h2_pyt)
  cfg$varE_ayt <- varE_from_base_h2(base_pop, cfg$h2_ayt)
  cfg$varE_eyt <- varE_from_base_h2(base_pop, cfg$h2_eyt)
  cfg
}

score_pop_for_selection <- function(pop, use_ebv = FALSE, ebv_trait = 1L) {
  if (isTRUE(use_ebv)) {
    ebv <- bp_pop_trait_vector(pop@ebv, as.integer(ebv_trait))
    if (length(ebv) == pop_n_ind(pop) && any(!is.na(ebv))) {
      return(as.numeric(ebv))
    }
  }

  ph <- pop@pheno
  if (!is.null(ph)) {
    if (is.null(dim(ph))) return(as.numeric(ph))
    return(as.numeric(ph[, 1L]))
  }
  stats::runif(pop_n_ind(pop))
}

select_top_idx <- function(score, n) {
  n <- max(1L, min(as.integer(n), length(score)))
  order(score, decreasing = TRUE)[seq_len(n)]
}

merge_nonempty_pops <- function(pop_list) {
  pop_list <- pop_list[!vapply(pop_list, is.null, logical(1))]
  if (length(pop_list) == 0L) return(NULL)
  if (length(pop_list) == 1L) return(pop_list[[1L]])
  merge_pops(pop_list)
}

split_halfsib_seed_indices <- function(n_total, n_families, seeds_per_family, n_to_pi_per_family) {
  stopifnot(n_total == as.integer(n_families * seeds_per_family))

  pi_idx <- integer(0)
  dh_idx <- integer(0)

  for (fam in seq_len(n_families)) {
    i0 <- (fam - 1L) * seeds_per_family + 1L
    fam_idx <- i0:(i0 + seeds_per_family - 1L)
    keep_pi <- sample(fam_idx, size = as.integer(n_to_pi_per_family), replace = FALSE)
    keep_dh <- setdiff(fam_idx, keep_pi)
    pi_idx <- c(pi_idx, keep_pi)
    dh_idx <- c(dh_idx, keep_dh)
  }

  list(pi_idx = pi_idx, dh_idx = dh_idx)
}

seed_pi_from_warm_state <- function(state, cfg) {
  src_pyt <- select_latest_available(
    state = state,
    stage = "PYT",
    stream = "main",
    combine = TRUE,
    silent = TRUE
  )

  src_cb <- select_latest_available(
    state = state,
    stage = "CROSS_BLOCK",
    stream = "main",
    combine = TRUE,
    silent = TRUE
  )

  input_seed_source <- src_pyt
  if (!is.null(src_pyt) && !is.null(src_cb)) {
    input_seed_source <- src_pyt
    input_seed_source$pop <- merge_nonempty_pops(list(src_pyt$pop, src_cb$pop))
    input_seed_source$source_ids <- unique(c(src_pyt$source_ids, src_cb$source_ids))
  } else if (is.null(src_pyt) && !is.null(src_cb)) {
    input_seed_source <- src_cb
  }

  if (is.null(input_seed_source) || is.null(input_seed_source$pop)) {
    stop("Could not find warm-start source cohorts (expected PYT and/or CROSS_BLOCK).", call. = FALSE)
  }

  seed_n <- min(as.integer(cfg$pi_seed_n), pop_n_ind(input_seed_source$pop))
  idx <- sample.int(pop_n_ind(input_seed_source$pop), size = seed_n, replace = FALSE)
  pi_seed <- pop_subset(input_seed_source$pop, idx)

  state <- put_stage_pop(
    state = state,
    pop = pi_seed,
    stage = "PI_CAND",
    source = input_seed_source,
    ready_in_years = 0,
    stream = "main",
    selection_strategy = sprintf("Warm-start random sample n=%d from PYT + CROSS_BLOCK", seed_n)
  )

  state
}


# Event Verbs: Population Improvement -------------------------------------
select_from_PI_and_cross_next_PI_cycle <- function(state, cfg, cycle_label) {
  bp_debug_break(state, cfg)

  state <- run_genotyping(state, list(
    input_stage = "PI_CAND",
    stream = "main",
    input_policy = "latest_one",
    include_not_ready = TRUE,
    chip = cfg$snp_chip,
    duration_years = cfg$pi_genotyping_duration_years,
    cost_per_sample = cfg$cost_genotype,
    silent = TRUE
  ))

  input_pi <- select_latest_available(
    state = state,
    stage = "PI_CAND",
    stream = "main",
    combine = TRUE,
    silent = TRUE
  )
  chk <- bp_skip_if_no_input(state, input_pi, cfg)
  if (chk$skip) return(chk$state)

  pop <- input_pi$pop
  n_total <- pop_n_ind(pop)
  if (n_total < 4L) return(state)

  perm <- sample.int(n_total, size = n_total, replace = FALSE)
  n_f_pool <- floor(n_total / 2)
  idx_f_pool <- perm[seq_len(n_f_pool)]
  idx_m_pool <- perm[(n_f_pool + 1L):n_total]

  pop_f_pool <- pop_subset(pop, idx_f_pool)
  pop_m_pool <- pop_subset(pop, idx_m_pool)

  model_entry <- state$gs_models[[cfg$model_id]]
  if (!is.null(model_entry)) {
    pop <- predict_ebv_pop(
      pop = pop,
      model_entry = model_entry,
      state = state,
      cfg = list(ebv_trait = as.integer(cfg$ebv_trait), cohort_ids = input_pi$source_ids),
      stage_label = "PI_CAND"
    )
    # Persist EBVs on the source cohort so monitoring can use cor_ebv_gv.
    if (!is.null(input_pi$source_rows) && nrow(input_pi$source_rows) == 1L) {
      state$pops[[as.character(input_pi$source_rows$cohort_id[1])]] <- pop
    }
    pop_f_pool <- pop_subset(pop, idx_f_pool)
    pop_m_pool <- pop_subset(pop, idx_m_pool)
  }

  score_f <- score_pop_for_selection(pop_f_pool, use_ebv = !is.null(model_entry), ebv_trait = cfg$ebv_trait)
  score_m <- score_pop_for_selection(pop_m_pool, use_ebv = !is.null(model_entry), ebv_trait = cfg$ebv_trait)

  n_f_sel <- min(as.integer(cfg$n_female_select), pop_n_ind(pop_f_pool))
  n_m_sel <- min(as.integer(cfg$n_male_select), pop_n_ind(pop_m_pool))

  idx_f_sel <- select_top_idx(score_f, n_f_sel)
  idx_m_sel <- select_top_idx(score_m, n_m_sel)

  fem <- pop_subset(pop_f_pool, idx_f_sel)
  mal <- pop_subset(pop_m_pool, idx_m_sel)

  parent_pop <- merge_pops(list(fem, mal))
  n_f <- pop_n_ind(fem)
  n_m <- pop_n_ind(mal)
  if (n_f == 0L || n_m == 0L) return(state)

  # Open-pollination approximation:
  # each selected female contributes seeds_per_cross pollinations, each with
  # a randomly sampled selected male.
  crosses_per_female <- as.integer(cfg$seeds_per_cross)
  female_ids <- rep(seq_len(n_f), each = crosses_per_female)
  male_ids <- n_f + sample.int(n_m, size = length(female_ids), replace = TRUE)
  cross_plan <- cbind(female_ids, male_ids)

  halfsib_seed <- makeCross(
    parent_pop,
    crossPlan = cross_plan,
    nProgeny = 1L,
    simParam = state$sim$SP
  )

  split_idx <- split_halfsib_seed_indices(
    n_total = pop_n_ind(halfsib_seed),
    n_families = n_f,
    seeds_per_family = as.integer(cfg$seeds_per_cross),
    n_to_pi_per_family = as.integer(cfg$seeds_to_pi_per_family)
  )

  pi_next <- pop_subset(halfsib_seed, split_idx$pi_idx)
  dh_seed <- pop_subset(halfsib_seed, split_idx$dh_idx)

  state <- put_stage_pop(
    state = state,
    pop = pi_next,
    stage = "PI_CAND",
    source = input_pi,
    ready_in_years = cfg$pi_cycle_duration_years,
    stream = "main",
    selection_strategy = sprintf("%s: split male/female then GS top-%d per sex", cycle_label, as.integer(cfg$n_female_select)),
    cross_strategy = sprintf(
      "Open-pollination approximation: %d random male pollinations per selected female",
      crosses_per_female
    ),
    inherit_genotypes = FALSE
  )
  state <- add_stage_cost(
    state = state,
    stage = "PI_CAND",
    event = "pi_crossing",
    unit = "cross",
    n = nrow(cross_plan),
    unit_cost = cfg$cost_crossing
  )

  state <- put_stage_pop(
    state = state,
    pop = dh_seed,
    stage = "PD_DH_INPUT",
    source = input_pi,
    ready_in_years = cfg$pi_cycle_duration_years,
    stream = "main",
    selection_strategy = sprintf("%s: reserve %d seeds/family for DH", cycle_label, as.integer(cfg$seeds_to_dh_per_family)),
    inherit_genotypes = FALSE
  )

  state
}


# Event Verbs: Product Development ----------------------------------------
advance_PD_seed_to_DH <- function(state, cfg) {
  bp_debug_break(state, cfg)

  input_pd_seed <- get_ready_pop(
    state = state,
    stage = "PD_DH_INPUT",
    stream = "main",
    policy = "all_ready",
    combine = TRUE,
    silent = TRUE
  )
  chk <- bp_skip_if_no_input(state, input_pd_seed, cfg)
  if (chk$skip) return(chk$state)

  seeds <- input_pd_seed$pop
  n_per_family <- as.integer(cfg$seeds_to_dh_per_family)
  n_families <- floor(pop_n_ind(seeds) / n_per_family)
  if (n_families < 1L) return(state)

  rep_idx <- seq(1L, by = n_per_family, length.out = n_families)
  fam_rep <- pop_subset(seeds, rep_idx)

  dh <- makeDH(
    pop = fam_rep,
    nDH = as.integer(cfg$dh_lines_per_family),
    simParam = state$sim$SP
  )

  state <- put_stage_pop(
    state = state,
    pop = dh,
    stage = "DH_BULK",
    source = input_pd_seed,
    ready_in_years = cfg$dh_duration_years,
    stream = "main",
    cross_strategy = sprintf("makeDH(nDH=%d) from one representative seed/family", as.integer(cfg$dh_lines_per_family)),
    inherit_genotypes = FALSE
  )
  state <- add_stage_cost(
    state = state,
    stage = "DH_BULK",
    event = "dh_production",
    unit = "line",
    n = pop_n_ind(dh),
    unit_cost = cfg$cost_dh_line
  )

  state
}

select_from_DH_and_run_headrow <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_dh <- select_latest_available(state, stage = "DH_BULK", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_dh, cfg)
  if (chk$skip) return(chk$state)

  n <- min(as.integer(cfg$n_headrow_advance), pop_n_ind(input_dh$pop))
  idx <- sample.int(pop_n_ind(input_dh$pop), size = n, replace = FALSE)
  pop_sel <- pop_subset(input_dh$pop, idx)

  run_phenotype_trial(
    state = state,
    pop = pop_sel,
    output_stage = "HEADROW_SEL",
    input_cohorts = input_dh$source_ids,
    selection_strategy = "Advance fixed n from DH_BULK to headrows",
    traits = 1L,
    n_loc = 1L,
    reps = 1L,
    varE = cfg$varE_headrow,
    duration_years = 1,
    stream = "main",
    cost_per_plot = cfg$cost_headrow_line,
    silent = TRUE
  )
}

select_from_headrow_and_run_PYT <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_headrow <- select_latest_available(state, stage = "HEADROW_SEL", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_headrow, cfg)
  if (chk$skip) return(chk$state)

  n <- min(as.integer(cfg$n_pyt), pop_n_ind(input_headrow$pop))
  idx <- sample.int(pop_n_ind(input_headrow$pop), size = n, replace = FALSE)
  pop_sel <- pop_subset(input_headrow$pop, idx)

  run_phenotype_trial(
    state = state,
    pop = pop_sel,
    output_stage = "PYT",
    input_cohorts = input_headrow$source_ids,
    selection_strategy = "Random subset from HEADROW_SEL",
    traits = 1L,
    n_loc = 1L,
    reps = 1L,
    varE = cfg$varE_pyt,
    duration_years = 1,
    stream = "main",
    cost_per_plot = cfg$cost_plot_pyt,
    silent = TRUE
  )
}

select_from_PYT_and_run_AYT <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_pyt <- select_latest_available(state, stage = "PYT", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_pyt, cfg)
  if (chk$skip) return(chk$state)

  n <- min(as.integer(cfg$n_pyt_to_ayt), pop_n_ind(input_pyt$pop))
  ph <- input_pyt$pop@pheno
  idx <- if (is.null(ph)) {
    sample.int(pop_n_ind(input_pyt$pop), size = n, replace = FALSE)
  } else {
    if (!is.null(dim(ph))) ph <- ph[, 1L]
    order(as.numeric(ph), decreasing = TRUE)[seq_len(n)]
  }
  pop_sel <- pop_subset(input_pyt$pop, idx)

  run_phenotype_trial(
    state = state,
    pop = pop_sel,
    output_stage = "AYT",
    input_cohorts = input_pyt$source_ids,
    selection_strategy = "Top phenotype from latest PYT",
    traits = 1L,
    n_loc = cfg$ayt_effective_reps,
    reps = 1L,
    varE = cfg$varE_ayt,
    duration_years = 1,
    stream = "main",
    cost_per_plot = cfg$cost_plot_ayt,
    silent = TRUE
  )
}

select_from_AYT_and_run_EYT1 <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_ayt <- select_latest_available(state, stage = "AYT", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_ayt, cfg)
  if (chk$skip) return(chk$state)

  n <- min(as.integer(cfg$n_ayt_to_eyt), pop_n_ind(input_ayt$pop))
  ph <- input_ayt$pop@pheno
  idx <- if (is.null(ph)) {
    sample.int(pop_n_ind(input_ayt$pop), size = n, replace = FALSE)
  } else {
    if (!is.null(dim(ph))) ph <- ph[, 1L]
    order(as.numeric(ph), decreasing = TRUE)[seq_len(n)]
  }
  pop_sel <- pop_subset(input_ayt$pop, idx)

  run_phenotype_trial(
    state = state,
    pop = pop_sel,
    output_stage = "EYT1",
    input_cohorts = input_ayt$source_ids,
    selection_strategy = "Top phenotype from latest AYT",
    traits = 1L,
    n_loc = cfg$eyt_effective_reps,
    reps = 1L,
    varE = cfg$varE_eyt,
    duration_years = 1,
    stream = "main",
    cost_per_plot = cfg$cost_plot_eyt,
    silent = TRUE
  )
}

select_from_EYT1_and_run_EYT2 <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_eyt1 <- select_latest_available(state, stage = "EYT1", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_eyt1, cfg)
  if (chk$skip) return(chk$state)

  run_phenotype_trial(
    state = state,
    pop = input_eyt1$pop,
    output_stage = "EYT2",
    input_cohorts = input_eyt1$source_ids,
    selection_strategy = "Re-evaluate EYT1 lines",
    traits = 1L,
    n_loc = cfg$eyt_effective_reps,
    reps = 1L,
    varE = cfg$varE_eyt,
    duration_years = 1,
    stream = "main",
    cost_per_plot = cfg$cost_plot_eyt,
    silent = TRUE
  )
}

select_from_EYT_and_release_Variety <- function(state, cfg) {
  bp_debug_break(state, cfg)

  input_eyt2 <- select_latest_available(state, stage = "EYT2", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_eyt2, cfg)
  if (chk$skip) return(chk$state)

  pop <- input_eyt2$pop

  ph <- state$phenotype_log
  ph <- ph[
    ph$stage %in% c("EYT1", "EYT2") &
      ph$environment == 0L &
      ph$trait == "trait1" &
      ph$individual_id %in% pop@id,
    , drop = FALSE
  ]

  if (nrow(ph) > 0L) {
    avg <- stats::aggregate(phenotype_value ~ individual_id, data = ph, FUN = mean)
    idx <- match(pop@id, avg$individual_id)
    score <- rep(-Inf, pop_n_ind(pop))
    hit <- which(!is.na(idx))
    score[hit] <- avg$phenotype_value[idx[hit]]
    best <- which.max(score)
  } else {
    best <- 1L
  }

  variety <- pop_subset(pop, best)

  old <- bp_get_ready_cohorts(state, stage = "Variety", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)
  if (nrow(old) > 0L) {
    for (i in seq_len(nrow(old))) state <- bp_close_cohort(state, old$cohort_id[i])
  }

  state <- put_stage_pop(
    state = state,
    pop = variety,
    stage = "Variety",
    source = input_eyt2,
    ready_in_years = 0,
    stream = "main",
    selection_strategy = "Best 2-year EYT mean phenotype",
    inherit_genotypes = FALSE
  )

  state$outputs$varieties <- rbind(
    state$outputs$varieties,
    data.frame(
      tick = as.integer(state$time$tick),
      source_cohort_id = paste(input_eyt2$source_ids, collapse = ";"),
      variety_id = as.integer(variety@id),
      stringsAsFactors = FALSE
    )
  )

  state
}

genotype_PD_and_update_GP_model <- function(state, cfg) {
  bp_debug_break(state, cfg)

  state <- run_genotyping(state, list(
    input_stage = "HEADROW_SEL",
    stream = "main",
    input_policy = "latest_one",
    include_not_ready = TRUE,
    chip = cfg$snp_chip,
    duration_years = cfg$geno_duration_years,
    cost_per_sample = cfg$cost_genotype,
    silent = TRUE
  ))

  state <- run_genotyping(state, list(
    input_stage = "PYT",
    stream = "main",
    input_policy = "latest_one",
    include_not_ready = TRUE,
    chip = cfg$snp_chip,
    duration_years = cfg$geno_duration_years,
    cost_per_sample = cfg$cost_genotype,
    silent = TRUE
  ))

  state <- run_train_gp_model(state, list(
    from_stage = "PYT",
    stream = "main",
    chip = cfg$snp_chip,
    trait = cfg$ebv_trait,
    lookback_years = cfg$gp_lookback_years,
    training_policy = "all_ready",
    model_id = cfg$model_id
  ))

  state
}


# Yearly Schedule ----------------------------------------------------------
run_one_year <- function(state, cfg, year_index) {
  bp_debug_break(state, cfg, year = year_index)
  tick_start <- as.integer(state$time$tick)

  # tick 1: annual DH production first, then PI cycle A + PD stage starts + model refresh
  state <- advance_PD_seed_to_DH(state, cfg)
  state <- select_from_PI_and_cross_next_PI_cycle(state, cfg, cycle_label = sprintf("year_%d_cycle_A", year_index))
  state <- select_from_DH_and_run_headrow(state, cfg)
  state <- select_from_headrow_and_run_PYT(state, cfg)
  state <- select_from_PYT_and_run_AYT(state, cfg)
  state <- select_from_AYT_and_run_EYT1(state, cfg)
  state <- select_from_EYT1_and_run_EYT2(state, cfg)
  state <- select_from_EYT_and_release_Variety(state, cfg)
  state <- genotype_PD_and_update_GP_model(state, cfg)
  state <- bp_advance_time(state, n_ticks = 1)

  # tick 2: no new events
  state <- bp_advance_time(state, n_ticks = 1)

  # tick 3: PI cycle B
  state <- select_from_PI_and_cross_next_PI_cycle(state, cfg, cycle_label = sprintf("year_%d_cycle_B", year_index))
  state <- bp_advance_time(state, n_ticks = 1)

  # tick 4: no new events
  state <- bp_advance_time(state, n_ticks = 1)

  ticks_elapsed <- as.integer(state$time$tick - tick_start)
  if (!identical(ticks_elapsed, as.integer(cfg$ticks_per_year))) {
    stop(
      sprintf("run_one_year tick mismatch: elapsed=%d but ticks_per_year=%d", ticks_elapsed, as.integer(cfg$ticks_per_year)),
      call. = FALSE
    )
  }

  cat(sprintf(
    "scheme=%s year=%d t=%.2f active=%d pi=%d dh_in=%d dh=%d head=%d pyt=%d ayt=%d eyt1=%d eyt2=%d variety=%d models=%d\\n",
    cfg$scheme_id,
    year_index,
    state$time$t,
    sum(state$cohorts$active),
    nrow(bp_get_ready_cohorts(state, stage = "PI_CAND", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "PD_DH_INPUT", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "DH_BULK", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "HEADROW_SEL", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "PYT", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "AYT", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "EYT1", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "EYT2", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "Variety", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    length(state$gs_models)
  ))

  state
}


# Config ------------------------------------------------------------------
make_two_part_cfg <- function() {
  cfg <- list(
    scheme_id = "TwoPart_warm_from_WheatConv_DRAFT",
    debug = FALSE,
    debug_after_year = NULL,
    debug_after_tick = NULL,
    debug_where = NULL,

    # Warm-start
    pi_seed_n = 2000L,

    # PI
    model_id = "gp_pi_main",
    ebv_trait = 1L,
    n_female_select = 100L,
    n_male_select = 100L,
    seeds_per_cross = 40L,
    seeds_to_pi_per_family = 30L,
    seeds_to_dh_per_family = 10L,
    pi_cycle_duration_years = 0.5,

    # PD
    dh_lines_per_family = 31L,
    dh_duration_years = 1,
    n_headrow_advance = 500L,
    n_pyt = 500L,
    n_pyt_to_ayt = 50L,
    n_ayt_to_eyt = 10L,

    # Trial heritability assumptions
    h2_headrow = 0.10,
    h2_pyt = 0.2,
    h2_ayt = 0.50,
    h2_eyt = 0.67,

    # Effective replication assumptions
    ayt_effective_reps = 4L,
    eyt_effective_reps = 8L,

    # Genotyping/model
    snp_chip = 1L,
    pi_genotyping_duration_years = 0,
    geno_duration_years = 0.5,
    gp_lookback_years = 2,
    enable_network_render = FALSE,

    # Placeholder costs
    cost_crossing = 1,
    cost_dh_line = 0.05,
    cost_headrow_line = 0.05,
    cost_plot_pyt = 20,
    cost_plot_ayt = 30,
    cost_plot_eyt = 35,
    cost_genotype = 25
  )

  cfg$dt <- 0.25
  cfg$ticks_per_year <- as.integer(round(1 / cfg$dt))
  cfg
}


# Initialization -----------------------------------------------------------
init_two_part_from_wheat_conv <- function(warm_years = 20L, cfg = make_two_part_cfg()) {
  env <- new.env(parent = globalenv())
  # Keep all nested source() calls from the warm-start script inside this sandbox env.
  env$source <- function(file, ...) base::source(file, local = env, ...)
  base::source("examples/Wheat_Conv_noGS_readable_structured.R", local = env)

  if (!exists("run_wheat_conv_nogs", envir = env, inherits = FALSE)) {
    stop("Could not load run_wheat_conv_nogs() from examples/Wheat_Conv_noGS_readable_structured.R", call. = FALSE)
  }

  warm_state <- env$run_wheat_conv_nogs(
    n_fill_years = as.integer(warm_years),
    n_run_years = 0L,
    cfg = env$make_wheat_conv_nogs_cfg()
  )

  cfg <- resolve_two_part_varE(cfg, warm_state$pops[[warm_state$cohorts$cohort_id[1L]]])

  state <- warm_state
  if (length(state$sim$SP$snpChips) < as.integer(cfg$snp_chip)) {
    state$sim$SP$addSnpChip(nSnpPerChr = 20L)
  }
  state <- seed_pi_from_warm_state(state, cfg)

  state <- run_genotyping(state, list(
    input_stage = "PYT",
    stream = "main",
    input_policy = "all_ready",
    include_not_ready = TRUE,
    chip = cfg$snp_chip,
    duration_years = 0,
    cost_per_sample = cfg$cost_genotype,
    silent = TRUE
  ))

  state <- run_train_gp_model(state, list(
    from_stage = "PYT",
    stream = "main",
    chip = cfg$snp_chip,
    trait = cfg$ebv_trait,
    lookback_years = cfg$gp_lookback_years,
    training_policy = "all_ready",
    model_id = cfg$model_id
  ))

  list(state = state, cfg = cfg)
}


# Reporting ----------------------------------------------------------------
print_two_part_setup <- function() {
  cat("\\n=== Two-Part Draft: Filled from Input ===\\n")
  print(make_two_part_filled_table(), row.names = FALSE)

  cat("\\n=== Two-Part Draft: Missing Specs ===\\n")
  print(make_two_part_missing_table(), row.names = FALSE)

  cat("\\n=== Two-Part Draft: Design Questions ===\\n")
  qs <- make_two_part_design_questions()
  for (i in seq_along(qs)) cat(sprintf("%d. %s\\n", i, qs[i]))
}

print_two_part_summary <- function(state) {
  cat("\\nVariety releases (recent):\\n")
  print(utils::tail(state$outputs$varieties, 10))

  cat("\\nTotal cost by event:\\n")
  print(stats::aggregate(total_cost ~ event, data = state$cost_log, sum))

  cat("\\nRecent cohorts:\\n")
  print(utils::tail(state$cohorts[, c("cohort_id", "stage", "created_tick", "available_tick", "active", "n_ind")], 20))

  cat("\\nEvent timeline (collapsed patterns):\\n")
  bp_print_event_timeline(state, collapse_year_patterns = TRUE, digits = 2)
}

# Runner -------------------------------------------------------------------
run_two_part_warm_draft <- function(
    warm_years = 20L,
    n_fill_years = 3L,
    n_run_years = 10L,
    cfg = make_two_part_cfg()) {
  sim <- init_two_part_from_wheat_conv(warm_years = warm_years, cfg = cfg)
  state <- sim$state
  cfg <- sim$cfg

  cfg$fail_on_missing_input <- FALSE
  for (fill_yr in seq_len(as.integer(n_fill_years))) {
    state <- run_one_year(state, cfg, year_index = fill_yr)
  }
  cfg$fail_on_missing_input <- TRUE
  for (run_yr in as.integer(n_fill_years) + seq_len(as.integer(n_run_years))) {
    state <- run_one_year(state, cfg, year_index = run_yr)
  }

  state
}


# Script Entry -------------------------------------------------------------
if (identical(environment(), globalenv()) && !identical(Sys.getenv("BPS_SKIP_SCRIPT_ENTRY"), "1")) {
  print_two_part_setup()

  cfg <- make_two_part_cfg()
  out <- run_two_part_warm_draft(warm_years = 20L, n_fill_years = 3L, n_run_years = 10L, cfg = cfg)
  print_two_part_summary(out)
}

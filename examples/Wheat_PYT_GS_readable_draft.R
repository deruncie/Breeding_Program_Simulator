# Wheat PYT GS draft scheme (readable, first pass).
#
# Scope for this draft:
# - Start from conventional winter-wheat style backbone.
# - Add genomic selection at PYT stage and use EBV-based advancement in PYT/AYT.
# - Select next crossing block (50) by top EBV from:
#   last crossing block + all active PYT/AYT/EYT entries.
# - Keep unresolved timing/training assumptions explicit in tables/questions below.
#
# Agreement gate note:
# - This script is a behavior draft only.
# - No callable/grid wrapper packaging is included yet.

# Setup -------------------------------------------------------------------
library(AlphaSimR)

source("R/operators.R")
source("R/readable_api.R")
source("R/state_print.R")
source("R/monitoring.R")


# Filled-From-Input Table --------------------------------------------------
make_wheat_pyt_gs_filled_table <- function() {
  data.frame(
    parameter = c(
      "program_variant",
      "crossing_block_size",
      "crosses_per_year",
      "possible_parent_pairs",
      "dh_per_cross",
      "cross_to_parent_cycle_years",
      "cross_to_release_years",
      "headrow_selected",
      "pyt_entries",
      "pyt_to_ayt",
      "ayt_to_eyt",
      "eyt_years",
      "release_count",
      "h2_headrow",
      "h2_pyt",
      "h2_ayt",
      "h2_eyt",
      "pyt_gs_enabled",
      "pyt_gs_parent_pool",
      "burn_in_years",
      "future_phase_years",
      "dh_per_cross_pyt_gs"
    ),
    value = c(
      "PYT GS",
      "50 lines",
      "100",
      "1225",
      "100 (conventional)",
      "4 (conventional), reduced to 3 (PYT GS target)",
      "8",
      "500",
      "500",
      "50",
      "10",
      "2 (EYT1 + EYT2)",
      "1 variety/year",
      "0.10",
      "0.20",
      "0.50",
      "0.67",
      "Yes (at PYT stage)",
      "last crossing block + all PYT and later trial entries",
      "20",
      "20 (+ run-out to completion)",
      "97"
    ),
    source = c(
      rep("user text", 22)
    ),
    stringsAsFactors = FALSE
  )
}


# Missing-Specification Table ----------------------------------------------
make_wheat_pyt_gs_missing_table <- function() {
  data.frame(
    parameter = c(
      "Exact tick-level event order to enforce 3-year crossing cycle",
      "GS model training source stages",
      "GS model lookback window",
      "Whether crossing block selection should include inactive/archived trial cohorts",
      "How to handle candidate shortfall when <50 unique EBV-scored lines exist",
      "Whether varieties remain eligible as future parents",
      "Burn-in handoff strategy for this script (start directly as PYT GS vs warm-start from Conv)",
      "Network validation toggle for this step",
      "Cost schedule for genotyping and stage operations"
    ),
    current_default = c(
      "Single yearly tick bundle; PYT GS uses latest available cohorts and includes non-ready cohorts for genotyping/scoring",
      "PYT only",
      "4 years",
      "No (active cohorts only)",
      "Fill from best available crossing-block incumbents",
      "Yes unless explicitly excluded",
      "Start directly in PYT GS mode from founder initialization",
      "Disabled (per current step)",
      "Placeholder values aligned to existing readable examples"
    ),
    why_it_matters = c(
      "Determines whether the intended 3-year cycle is actually achieved in simulation time.",
      "Controls prediction stability and realized selection accuracy.",
      "Changes responsiveness vs noise in the GP model.",
      "Affects parent diversity and potential double-use of historical cohorts.",
      "Prevents crossing block underfill and changes selection pressure.",
      "Impacts genetic turnover and realism of product pipeline.",
      "Changes comparability to your stated burn-in protocol.",
      "Determines whether scheme extraction/diagram checks are run now.",
      "Required for fair cost-equivalent comparisons across program variants."
    ),
    user_decision_needed = c(
      rep("confirm/replace", 9)
    ),
    stringsAsFactors = FALSE
  )
}


# Design Questions ----------------------------------------------------------
make_wheat_pyt_gs_design_questions <- function() {
  c(
    "For strict 3-year crossing cycle behavior, should crossing use EBV-ranked PYT entries from the same simulation year when available?",
    "Should GS training remain PYT-only, or include AYT/EYT phenotypes as additional training rows?",
    "For crossing block candidate pool, do you want all active PYT+AYT+EYT cohorts or only latest cohort per stage?",
    "If top-50 parent set contains duplicates from overlapping cohorts, should dedup be by individual ID (current default) or by cohort-stage priority?",
    "Do you want network extraction/compare/render enabled after this draft run, or keep it off for now?"
  )
}


# Helper Utilities ----------------------------------------------------------
varE_from_base_h2 <- function(base_pop, h2, trait = 1L) {
  h2 <- as.numeric(h2)
  if (!is.finite(h2) || h2 <= 0 || h2 >= 1) return(1.0)

  gv <- base_pop@gv
  if (is.null(dim(gv))) gv <- matrix(gv, ncol = 1L)
  vg <- stats::var(gv[, as.integer(trait)], na.rm = TRUE)
  if (is.na(vg) || vg <= 0) return(1.0)
  as.numeric(vg * (1 - h2) / h2)
}

resolve_stage_varE <- function(cfg, base_pop) {
  cfg$varE_headrow <- varE_from_base_h2(base_pop, cfg$h2_headrow)
  cfg$varE_pyt <- varE_from_base_h2(base_pop, cfg$h2_pyt)
  cfg$varE_ayt <- varE_from_base_h2(base_pop, cfg$h2_ayt)
  cfg$varE_eyt <- varE_from_base_h2(base_pop, cfg$h2_eyt)
  cfg
}

pop_unique_by_id <- function(pop) {
  keep <- !duplicated(pop@id)
  pop_subset(pop, which(keep))
}

merge_nonempty_pops <- function(pop_list) {
  pop_list <- pop_list[!vapply(pop_list, is.null, logical(1))]
  if (length(pop_list) == 0L) return(NULL)
  if (length(pop_list) == 1L) return(pop_list[[1L]])
  pop_unique_by_id(merge_pops(pop_list))
}

score_pop_for_selection <- function(pop, prefer_ebv = TRUE, ebv_trait = 1L) {
  if (isTRUE(prefer_ebv)) {
    ebv <- bp_pop_trait_vector(pop@ebv, as.integer(ebv_trait))
    if (length(ebv) == pop_n_ind(pop) && any(!is.na(ebv))) return(as.numeric(ebv))
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

get_active_stage_rows <- function(state, stage) {
  rows <- state$cohorts
  rows <- rows[rows$active & rows$stage == stage, , drop = FALSE]
  if (nrow(rows) == 0L) return(rows)
  rows[order(rows$created_tick, decreasing = TRUE), , drop = FALSE]
}

make_cross_plan_without_replacement <- function(n_parents, n_crosses) {
  if (n_parents < 2L) return(matrix(c(1L, 1L), ncol = 2L))
  all_pairs <- utils::combn(n_parents, 2)
  n_possible <- ncol(all_pairs)
  n_take <- min(as.integer(n_crosses), n_possible)
  take <- sample.int(n_possible, size = n_take, replace = FALSE)
  t(all_pairs[, take, drop = FALSE])
}


# Event Verbs --------------------------------------------------------------
advance_F1_to_DH <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_cross_seed <- select_latest_available(state, stage = "CROSS_SEED", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_cross_seed, cfg)
  if (chk$skip) return(chk$state)

  dh <- makeDH(input_cross_seed$pop, nDH = as.integer(cfg$dh_per_family), simParam = state$sim$SP)
  state <- put_stage_pop(
    state = state,
    pop = dh,
    stage = "DH_BULK",
    source = input_cross_seed,
    cross_strategy = sprintf("makeDH(nDH=%d) per F1 family", as.integer(cfg$dh_per_family)),
    ready_in_years = cfg$dh_duration_years,
    stream = "main",
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
  ph <- input_dh$pop@pheno
  idx <- if (is.null(ph)) {
    sample.int(pop_n_ind(input_dh$pop), size = n, replace = FALSE)
  } else {
    if (!is.null(dim(ph))) ph <- ph[, 1L]
    order(as.numeric(ph), decreasing = TRUE)[seq_len(n)]
  }
  pop_sel <- pop_subset(input_dh$pop, idx)

  run_phenotype_trial(
    state = state,
    pop = pop_sel,
    output_stage = "HEADROW_SEL",
    input_cohorts = input_dh$source_ids,
    selection_strategy = "Visual selection on correlated-trait phenotype (h2 target 0.1)",
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
    selection_strategy = "Random subset from HEADROW_SEL into unreplicated PYT",
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

genotype_trials_and_update_GP_model <- function(state, cfg) {
  bp_debug_break(state, cfg)

  for (stg in c("PYT", "AYT", "EYT1", "EYT2", "CROSS_BLOCK")) {
    state <- run_genotyping(state, list(
      input_stage = stg,
      stream = "main",
      input_policy = "all_ready",
      include_not_ready = TRUE,
      chip = cfg$snp_chip,
      duration_years = cfg$geno_duration_years,
      cost_per_sample = cfg$cost_genotype,
      silent = TRUE
    ))
  }

  state <- run_train_gp_model(state, list(
    from_stage = cfg$gp_training_from_stage,
    stream = "main",
    chip = cfg$snp_chip,
    trait = cfg$ebv_trait,
    lookback_years = cfg$gp_lookback_years,
    training_policy = "all_ready",
    model_id = cfg$model_id,
    on_no_source = "skip"
  ))

  model <- state$gs_models[[cfg$model_id]]
  if (is.null(model)) return(state)

  rows <- state$cohorts
  rows <- rows[rows$active & rows$stage %in% c("PYT", "AYT", "EYT1", "EYT2", "CROSS_BLOCK"), , drop = FALSE]
  if (nrow(rows) == 0L) return(state)

  rows <- rows[order(rows$created_tick, decreasing = TRUE), , drop = FALSE]
  for (i in seq_len(nrow(rows))) {
    cid <- as.character(rows$cohort_id[i])
    pop <- state$pops[[cid]]
    pop2 <- predict_ebv_pop(
      pop = pop,
      model_entry = model,
      state = state,
      cfg = list(cohort_ids = cid, chip = cfg$snp_chip, ebv_trait = cfg$ebv_trait),
      stage_label = as.character(rows$stage[i])
    )
    state$pops[[cid]] <- pop2
  }

  state
}

select_from_PYT_and_run_AYT <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_pyt <- select_latest_available(state, stage = "PYT", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_pyt, cfg)
  if (chk$skip) return(chk$state)
  n <- min(as.integer(cfg$n_pyt_to_ayt), pop_n_ind(input_pyt$pop))
  score <- score_pop_for_selection(input_pyt$pop, prefer_ebv = TRUE, ebv_trait = cfg$ebv_trait)
  idx <- select_top_idx(score, n)
  pop_sel <- pop_subset(input_pyt$pop, idx)

  run_phenotype_trial(
    state = state,
    pop = pop_sel,
    output_stage = "AYT",
    input_cohorts = input_pyt$source_ids,
    selection_strategy = "Top EBV from latest PYT (fallback phenotype/random)",
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
  score <- score_pop_for_selection(input_ayt$pop, prefer_ebv = TRUE, ebv_trait = cfg$ebv_trait)
  idx <- select_top_idx(score, n)
  pop_sel <- pop_subset(input_ayt$pop, idx)

  run_phenotype_trial(
    state = state,
    pop = pop_sel,
    output_stage = "EYT1",
    input_cohorts = input_ayt$source_ids,
    selection_strategy = "Top EBV from latest AYT (fallback phenotype/random)",
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
    selection_strategy = "Re-evaluate same EYT1 lines",
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
    selection_strategy = "Best 2-year EYT mean phenotype",
    ready_in_years = 0,
    stream = "main",
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

update_parents_and_cross_to_F1 <- function(state, cfg) {
  bp_debug_break(state, cfg)

  src_cb <- select_latest_available(state, stage = "CROSS_BLOCK", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, src_cb, cfg)
  if (chk$skip) return(chk$state)

  rows_pyt <- get_active_stage_rows(state, "PYT")
  rows_ayt <- get_active_stage_rows(state, "AYT")
  rows_e1 <- get_active_stage_rows(state, "EYT1")
  rows_e2 <- get_active_stage_rows(state, "EYT2")

  trial_ids <- unique(c(rows_pyt$cohort_id, rows_ayt$cohort_id, rows_e1$cohort_id, rows_e2$cohort_id))
  trial_pops <- lapply(trial_ids, function(cid) state$pops[[as.character(cid)]])

  candidate <- merge_nonempty_pops(c(list(src_cb$pop), trial_pops))
  if (is.null(candidate) || pop_n_ind(candidate) == 0L) {
    candidate <- src_cb$pop
  }

  score <- score_pop_for_selection(candidate, prefer_ebv = TRUE, ebv_trait = cfg$ebv_trait)
  idx <- select_top_idx(score, cfg$n_cross_block)
  next_block <- pop_subset(candidate, idx)

  if (pop_n_ind(next_block) < cfg$n_cross_block) {
    missing <- cfg$n_cross_block - pop_n_ind(next_block)
    fill_ids <- which(!(src_cb$pop@id %in% next_block@id))
    if (length(fill_ids) > 0L) {
      fill <- pop_subset(src_cb$pop, fill_ids)
      fill_score <- score_pop_for_selection(fill, prefer_ebv = TRUE, ebv_trait = cfg$ebv_trait)
      take <- select_top_idx(fill_score, min(missing, pop_n_ind(fill)))
      next_block <- merge_nonempty_pops(list(next_block, pop_subset(fill, take)))
    }
  }

  src_ids_block <- unique(c(src_cb$source_ids, trial_ids))
  state <- put_stage_pop(
    state = state,
    pop = next_block,
    stage = "CROSS_BLOCK",
    source_ids = src_ids_block,
    selection_strategy = "Top 50 EBV from last CROSS_BLOCK + all active PYT/AYT/EYT entries",
    ready_in_years = cfg$cross_block_crossing_duration_years,
    stream = "main",
    inherit_genotypes = FALSE
  )
  state <- add_stage_cost(
    state = state,
    stage = "CROSS_BLOCK",
    event = "cross_block_update",
    unit = "line",
    n = pop_n_ind(next_block),
    unit_cost = cfg$cost_parent_recycle
  )

  plan <- make_cross_plan_without_replacement(pop_n_ind(next_block), cfg$n_crosses)
  f1 <- makeCross(next_block, crossPlan = plan, nProgeny = 1L, simParam = state$sim$SP)
  block_out_id <- bp_last_cohort_id(state)

  state <- put_stage_pop(
    state = state,
    pop = f1,
    stage = "CROSS_SEED",
    source_ids = block_out_id,
    cross_strategy = "Random pair sampling without replacement from updated CROSS_BLOCK",
    ready_in_years = cfg$cross_block_crossing_duration_years,
    stream = "main",
    inherit_genotypes = FALSE
  )
  state <- add_stage_cost(
    state = state,
    stage = "CROSS_SEED",
    event = "crossing",
    unit = "cross",
    n = cfg$n_crosses,
    unit_cost = cfg$cost_crossing
  )

  state
}


# Yearly Schedule ----------------------------------------------------------
run_one_year <- function(state, cfg, year_index) {
  bp_debug_break(state, cfg, year = year_index)
  tick_start <- as.integer(state$time$tick)

  # tick 1: progression and GS refresh
  state <- select_from_DH_and_run_headrow(state, cfg)
  state <- select_from_headrow_and_run_PYT(state, cfg)
  state <- genotype_trials_and_update_GP_model(state, cfg)
  state <- select_from_PYT_and_run_AYT(state, cfg)
  state <- select_from_AYT_and_run_EYT1(state, cfg)
  state <- select_from_EYT1_and_run_EYT2(state, cfg)
  state <- select_from_EYT_and_release_Variety(state, cfg)
  state <- update_parents_and_cross_to_F1(state, cfg)
  state <- bp_advance_time(state, n_ticks = 1)

  # tick 2: no new events
  state <- bp_advance_time(state, n_ticks = 1)

  # tick 3: start DH production for latest CROSS_SEED
  state <- advance_F1_to_DH(state, cfg)
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
    "scheme=%s year=%d t=%.2f active=%d cross=%d dh=%d head=%d pyt=%d ayt=%d eyt1=%d eyt2=%d variety=%d models=%d\\n",
    cfg$scheme_id,
    year_index,
    state$time$t,
    sum(state$cohorts$active),
    nrow(bp_get_ready_cohorts(state, stage = "CROSS_BLOCK", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
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
make_wheat_pyt_gs_cfg <- function(dh_per_family = 97L, n_pyt = 500L) {
  cfg <- list(
    scheme_id = "Wheat_PYT_GS_DRAFT",
    debug = FALSE,
    debug_after_year = NULL,
    debug_after_tick = NULL,
    debug_where = NULL,

    # Fixed by user protocol
    n_cross_block = 50L,
    n_crosses = 100L,
    dh_per_family = as.integer(dh_per_family),
    cross_block_crossing_duration_years = 0.25,
    dh_duration_years = 1.25,
    n_headrow_advance = 500L,
    n_pyt = as.integer(n_pyt),
    n_pyt_to_ayt = 50L,
    n_ayt_to_eyt = 10L,

    # Heritabilities from protocol
    h2_headrow = 0.10,
    h2_pyt = 0.20,
    h2_ayt = 0.50,
    h2_eyt = 0.67,

    # Effective replication assumptions
    ayt_effective_reps = 4L,
    eyt_effective_reps = 8L,

    # GS controls (current inferred defaults)
    model_id = "gp_pyt",
    ebv_trait = 1L,
    snp_chip = 1L,
    gp_lookback_years = 4,
    gp_training_from_stage = "PYT",
    geno_duration_years = 0,
    enable_network_render = FALSE,

    # Placeholder costs
    cost_crossing = 1,
    cost_dh_line = 0.05,
    cost_headrow_line = 0.05,
    cost_plot_pyt = 20,
    cost_plot_ayt = 30,
    cost_plot_eyt = 35,
    cost_parent_recycle = 0.5,
    cost_genotype = 25
  )

  cfg$dt <- 0.25
  cfg$ticks_per_year <- as.integer(round(1 / cfg$dt))
  cfg
}

init_wheat_pyt_gs_sim <- function(cfg = make_wheat_pyt_gs_cfg()) {
  founder_haps <- quickHaplo(nInd = 50, nChr = 3, segSites = 120)
  SP <- SimParam$new(founder_haps)
  SP$addTraitA(nQtlPerChr = 60)
  SP$addSnpChip(nSnpPerChr = 20L)

  founders <- newPop(founder_haps, simParam = SP)
  crossing_block <- selectInd(founders, nInd = cfg$n_cross_block, use = "gv", simParam = SP)
  cfg <- resolve_stage_varE(cfg, founders)

  state <- bp_init_state(SP = SP, dt = cfg$dt, start_time = 0, sim = list(default_chip = 1L))
  state <- put_stage_pop(
    state = state,
    pop = crossing_block,
    stage = "CROSS_BLOCK",
    source = NULL,
    selection_strategy = "Initialize crossing block from founder GV",
    ready_in_years = 0,
    stream = "main",
    inherit_genotypes = FALSE
  )

  list(state = state, cfg = cfg)
}


# Reporting ----------------------------------------------------------------
print_wheat_pyt_gs_setup <- function() {
  cat("\\n=== Wheat PYT GS Draft: Filled from Input ===\\n")
  print(make_wheat_pyt_gs_filled_table(), row.names = FALSE)

  cat("\\n=== Wheat PYT GS Draft: Missing Specs ===\\n")
  print(make_wheat_pyt_gs_missing_table(), row.names = FALSE)

  cat("\\n=== Wheat PYT GS Draft: Design Questions ===\\n")
  qs <- make_wheat_pyt_gs_design_questions()
  for (i in seq_along(qs)) cat(sprintf("%d. %s\\n", i, qs[i]))
}

print_wheat_pyt_gs_summary <- function(state) {
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
run_wheat_pyt_gs_draft <- function(n_fill_years = 4L, n_run_years = 12L, cfg = make_wheat_pyt_gs_cfg()) {
  sim <- init_wheat_pyt_gs_sim(cfg)
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
  print_wheat_pyt_gs_setup()

  cfg <- make_wheat_pyt_gs_cfg()
  out <- run_wheat_pyt_gs_draft(n_fill_years = 4L, n_run_years = 12L, cfg = cfg)

  print_wheat_pyt_gs_summary(out)
}

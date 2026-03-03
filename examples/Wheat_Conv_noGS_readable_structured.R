# Wheat conventional scheme (no GS), with recycling from AYT.
#
# Implemented from user-provided protocol details:
# - Year 1: 100 bi-parental crosses from 50-parent crossing block,
#           sampled without replacement from 1225 possible pairs.
# - Years 1-2: 100 DH lines per biparental family.
# - Year 3: Headrows, visual selection modeled as yield-phenotype selection (h2 = 0.1),
#           advance 500 lines.
# - Year 4: PYT, unreplicated yield trial (h2 = 0.2), select top 50 to AYT,
#           top 20 to next crossing block.
# - Year 5: AYT, multilocation replicated trial (h2 = 0.5), select top 10 to EYT,
#           and top 10 as crossing-block candidates.
#           Next crossing block (50) = 20 PYT + 10 AYT + 20 best from current block's
#           non-PYT lines.
# - Year 6: EYT year 1 (h2 = 0.67), keep all 10 for year-2 reevaluation.
# - Year 7: EYT year 2 reevaluation.
# - Year 8: Release best line by 2-year EYT mean.

# Setup -------------------------------------------------------------------
library(AlphaSimR)

source("R/operators.R")
source("R/readable_api.R")
source("R/state_print.R")
source("R/monitoring.R")


# Scheme Tables ------------------------------------------------------------
make_wheat_conv_nogs_scheme_table <- function() {
  data.frame(
    stage = c("Crossing", "DH", "Headrows", "PYT", "AYT", "EYT1", "EYT2", "Variety"),
    count = c("100 crosses", "100 DH per family", "500 lines", "500 lines", "50 lines", "10 lines", "10 lines", "1 line"),
    modeled_rule = c(
      "Random pair sampling without replacement from 50-parent block",
      "makeDH(nDH = 100) per F1 family",
      "Select on pheno with target h2 = 0.1",
      "Select on pheno with target h2 = 0.2",
      "Select on pheno with target h2 = 0.5",
      "Evaluate all 10 with target h2 = 0.67",
      "Reevaluate same 10 with target h2 = 0.67",
      "Release by 2-year EYT average"
    ),
    stringsAsFactors = FALSE
  )
}

make_required_input_table <- function() {
  data.frame(
    parameter = c(
      "n_years",
      "cross_block_crossing_duration_years",
      "dh_duration_years",
      "n_loc_ayt_effective",
      "n_loc_eyt_effective",
      "cost parameters",
      "founder architecture"
    ),
    current_default = c("12", "0.25", "1.25", "4", "8", "simple placeholders", "quickHaplo synthetic"),
    why_needed = c(
      "Controls warm-up and number of release opportunities.",
      "Combined crossing-block assembly + crossing + seed harvest at first tick.",
      "DH production duration after crossing event.",
      "Used to justify AYT h2 assumption.",
      "Used to justify EYT h2 assumption.",
      "Required for realistic economic comparison.",
      "Affects realism of long-term gain and variance."
    ),
    stringsAsFactors = FALSE
  )
}

make_design_questions <- function() {
  c(
    "Should crossing-block incumbent selection exclude both PYT and AYT-selected lines, or only PYT lines as specified?",
    "When AYT/EYT candidate pools are short (early warm-up), should crossing block be filled by incumbent best or random?",
    "Should released varieties be removed from crossing-block eligibility immediately?",
    "Should we add explicit checks/standards in PYT/AYT/EYT stages?"
  )
}


# Helper Utilities ---------------------------------------------------------
varE_from_base_h2 <- function(base_pop, h2, trait = 1L) {
  h2 <- as.numeric(h2)
  if (!is.finite(h2) || h2 <= 0 || h2 >= 1) return(1.0)

  pop <- base_pop
  gv <- pop@gv
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

select_top_or_random <- function(pop, n, SP) {
  n <- max(1L, min(as.integer(n), pop_n_ind(pop)))
  ph <- pop@pheno
  if (is.null(ph)) {
    idx <- sample.int(pop_n_ind(pop), size = n, replace = FALSE)
    return(pop_subset(pop, idx))
  }
  out <- tryCatch(selectInd(pop, nInd = n, use = "pheno", simParam = SP), error = function(e) NULL)
  if (!is.null(out)) return(out)
  idx <- sample.int(pop_n_ind(pop), size = n, replace = FALSE)
  pop_subset(pop, idx)
}

make_cross_plan_without_replacement <- function(n_parents, n_crosses) {
  if (n_parents < 2L) return(matrix(c(1L, 1L), ncol = 2L))
  all_pairs <- utils::combn(n_parents, 2)
  n_possible <- ncol(all_pairs)
  n_take <- min(as.integer(n_crosses), n_possible)
  take <- sample.int(n_possible, size = n_take, replace = FALSE)
  t(all_pairs[, take, drop = FALSE])
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

update_cross_block_pheno_from_eyt <- function(state, cross_pop) {
  if (is.null(cross_pop) || pop_n_ind(cross_pop) == 0L) return(cross_pop)

  ph <- state$phenotype_log
  if (nrow(ph) == 0L) return(cross_pop)

  ph <- ph[
    ph$stage %in% c("EYT1", "EYT2") &
      ph$environment == 0L &
      ph$trait == "trait1" &
      ph$available_tick <= as.integer(state$time$tick),
    , drop = FALSE
  ]
  if (nrow(ph) == 0L) return(cross_pop)

  means <- stats::aggregate(phenotype_value ~ individual_id, data = ph, FUN = mean)
  idx <- match(cross_pop@id, means$individual_id)
  hit <- which(!is.na(idx))
  if (length(hit) == 0L) return(cross_pop)

  if (is.null(cross_pop@pheno) || length(cross_pop@pheno) == 0L) {
    cross_pop@pheno <- rep(NA_real_, pop_n_ind(cross_pop))
  }
  cross_pop@pheno[hit] <- means$phenotype_value[idx[hit]]
  cross_pop
}


# Event Verbs --------------------------------------------------------------
update_parents_and_cross_to_F1 <- function(state, cfg) {
  bp_debug_break(state, cfg)
  src_cb <- select_latest_available(state, stage = "CROSS_BLOCK", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, src_cb, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  src_pyt <- select_latest_available(state, stage = "PYT", stream = "main", combine = TRUE, silent = TRUE)
  src_ayt <- select_latest_available(state, stage = "AYT", stream = "main", combine = TRUE, silent = TRUE)

  next_block <- src_cb$pop
  if (!is.null(src_pyt) && !is.null(src_ayt)) {
    py20 <- select_top_or_random(src_pyt$pop, cfg$n_pyt_to_cross, state$sim$SP)
    ay10 <- select_top_or_random(src_ayt$pop, cfg$n_ayt_to_cross, state$sim$SP)

    incumbent <- update_cross_block_pheno_from_eyt(state, src_cb$pop)

    # Rule from protocol: incumbent reserve comes from current block non-PYT lines.
    keep_ids_exclude <- py20@id
    idx_non_pyt <- which(!(incumbent@id %in% keep_ids_exclude))
    non_pyt <- if (length(idx_non_pyt) > 0L) pop_subset(incumbent, idx_non_pyt) else NULL
    keep20 <- if (!is.null(non_pyt)) select_top_or_random(non_pyt, cfg$n_incumbent_keep, state$sim$SP) else NULL

    next_block <- merge_nonempty_pops(list(py20, ay10, keep20))
    if (!is.null(next_block)) {
      if (pop_n_ind(next_block) < cfg$n_cross_block) {
        missing <- cfg$n_cross_block - pop_n_ind(next_block)
        idx_fill <- which(!(incumbent@id %in% next_block@id))
        fill_pop <- if (length(idx_fill) > 0L) pop_subset(incumbent, idx_fill) else NULL
        if (!is.null(fill_pop) && pop_n_ind(fill_pop) > 0L) {
          fill_sel <- select_top_or_random(fill_pop, missing, state$sim$SP)
          next_block <- merge_nonempty_pops(list(next_block, fill_sel))
        }
      }
      if (pop_n_ind(next_block) > cfg$n_cross_block) {
        next_block <- select_top_or_random(next_block, cfg$n_cross_block, state$sim$SP)
      }
    }
  }

  if (is.null(next_block) || pop_n_ind(next_block) == 0L) next_block <- src_cb$pop
  src_ids_block <- unique(c(
    src_cb$source_ids,
    if (!is.null(src_pyt)) src_pyt$source_ids else character(0),
    if (!is.null(src_ayt)) src_ayt$source_ids else character(0)
  ))
  sel_desc_block <- if (!is.null(src_pyt) && !is.null(src_ayt)) {
    "20 PYT + 10 AYT + 20 incumbent non-PYT (EYT-informed when available)"
  } else {
    "Fallback to incumbent crossing block only"
  }

  state <- put_stage_pop(
    state = state,
    pop = next_block,
    stage = "CROSS_BLOCK",
    source_ids = src_ids_block,
    selection_strategy = sel_desc_block,
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
    cross_strategy = "Random pair sampling without replacement from current CROSS_BLOCK",
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

advance_F1_to_DH <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_cross_seed <- select_latest_available(state, stage = "CROSS_SEED", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_cross_seed, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  dh <- makeDH(input_cross_seed$pop, nDH = cfg$dh_per_family, simParam = state$sim$SP)
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
  state <- chk$state

  n <- min(as.integer(cfg$n_headrow_advance), pop_n_ind(input_dh$pop))
  headrow_pop <- select_top_or_random(input_dh$pop, n, state$sim$SP)
  state <- run_phenotype_trial(
    state = state,
    pop = headrow_pop,
    output_stage = "HEADROW_SEL",
    input_cohorts = input_dh$source_ids,
    selection_strategy = "Visual/headrow selection on phenotype",
    traits = 1L,
    n_loc = 1L,
    reps = 1L,
    varE = cfg$varE_headrow,
    duration_years = 1,
    stream = "main",
    cost_per_plot = cfg$cost_headrow_line,
    silent = TRUE
  )
  # Optional debug line:
  # HEADROW_SEL <- select_current(state, "HEADROW_SEL", stream = "main")
  state
}

select_from_headrow_and_run_pyt <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_headrow <- select_latest_available(state, stage = "HEADROW_SEL", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_headrow, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  n <- min(as.integer(cfg$n_pyt), pop_n_ind(input_headrow$pop))
  idx <- sample.int(pop_n_ind(input_headrow$pop), size = n, replace = FALSE)
  pyt_pop <- pop_subset(input_headrow$pop, idx)
  state <- run_phenotype_trial(
    state = state,
    pop = pyt_pop,
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
  # Optional debug line:
  # PYT <- select_current(state, "PYT", stream = "main")
  state
}

select_from_PYT_and_run_AYT <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_pyt <- select_latest_available(state, stage = "PYT", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_pyt, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  n <- min(as.integer(cfg$n_pyt_to_ayt), pop_n_ind(input_pyt$pop))
  ayt_pop <- select_top_or_random(input_pyt$pop, n, state$sim$SP)
  state <- run_phenotype_trial(
    state = state,
    pop = ayt_pop,
    output_stage = "AYT",
    input_cohorts = input_pyt$source_ids,
    selection_strategy = "Top by phenotype from latest PYT",
    traits = 1L,
    n_loc = cfg$ayt_effective_reps,
    reps = 1L,
    varE = cfg$varE_ayt,
    duration_years = 1,
    stream = "main",
    cost_per_plot = cfg$cost_plot_ayt,
    silent = TRUE
  )
  # Optional debug line:
  # AYT <- select_current(state, "AYT", stream = "main")
  state
}

select_from_AYT_and_run_EYT1 <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_ayt <- select_latest_available(state, stage = "AYT", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_ayt, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  n <- min(as.integer(cfg$n_ayt_to_eyt), pop_n_ind(input_ayt$pop))
  eyt1_pop <- select_top_or_random(input_ayt$pop, n, state$sim$SP)
  state <- run_phenotype_trial(
    state = state,
    pop = eyt1_pop,
    output_stage = "EYT1",
    input_cohorts = input_ayt$source_ids,
    selection_strategy = "Top by phenotype from latest AYT",
    traits = 1L,
    n_loc = cfg$eyt_effective_reps,
    reps = 1L,
    varE = cfg$varE_eyt,
    duration_years = 1,
    stream = "main",
    cost_per_plot = cfg$cost_plot_eyt,
    silent = TRUE
  )
  state
}

run_EYT2_from_EYT1 <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_eyt1 <- select_latest_available(state, stage = "EYT1", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_eyt1, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  state <- run_phenotype_trial(
    state = state,
    pop = input_eyt1$pop,
    output_stage = "EYT2",
    input_cohorts = input_eyt1$source_ids,
    selection_strategy = "Re-evaluate same lines from EYT1",
    traits = 1L,
    n_loc = cfg$eyt_effective_reps,
    reps = 1L,
    varE = cfg$varE_eyt,
    duration_years = 1,
    stream = "main",
    cost_per_plot = cfg$cost_plot_eyt,
    silent = TRUE
  )
  state
}

select_from_EYT_and_release_Variety <- function(state, cfg) {
  bp_debug_break(state, cfg)
  input_eyt2 <- select_latest_available(state, stage = "EYT2", stream = "main", combine = TRUE, silent = TRUE)
  chk <- bp_skip_if_no_input(state, input_eyt2, cfg)
  if (chk$skip) return(chk$state)
  state <- chk$state

  pop <- input_eyt2$pop

  # Release by average over EYT1 + EYT2 records for each line id.
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
    selection_strategy = "Best 2-year mean phenotype from EYT1+EYT2",
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


# Yearly Schedule ----------------------------------------------------------
run_one_year <- function(state, cfg, year_index) {
  bp_debug_break(state, cfg, year = year_index)
  tick_start <- as.integer(state$time$tick)

  # tick 1
  # and all field-trial starts (headrows/PYT/AYT/EYT) are attempted here.
  state <- update_parents_and_cross_to_F1(state, cfg)
  state <- select_from_DH_and_run_headrow(state, cfg)
  state <- select_from_headrow_and_run_pyt(state, cfg)
  state <- select_from_PYT_and_run_AYT(state, cfg)
  state <- select_from_AYT_and_run_EYT1(state, cfg)
  state <- run_EYT2_from_EYT1(state, cfg)
  state <- select_from_EYT_and_release_Variety(state, cfg)

  state <- bp_advance_time(state, n_ticks = 1)

  # tick 2
  # nothing new
  state <- bp_advance_time(state, n_ticks = 1)

  # tick 3
  # start DH product, will take 5 ticks
  state <- advance_F1_to_DH(state, cfg)

  state <- bp_advance_time(state, n_ticks = 1)

  # tick 4
  # nothing new
  state <- bp_advance_time(state, n_ticks = 1)

  ticks_elapsed <- as.integer(state$time$tick - tick_start)
  if (!identical(ticks_elapsed, as.integer(cfg$ticks_per_year))) {
    stop(
      sprintf(
        "run_one_year tick mismatch: elapsed=%d but ticks_per_year=%d",
        ticks_elapsed, as.integer(cfg$ticks_per_year)
      ),
      call. = FALSE
    )
  }


  cat(sprintf(
    "scheme=%s year=%d t=%.2f active=%d cross=%d dh=%d head=%d pyt=%d ayt=%d eyt1=%d eyt2=%d variety=%d\n",
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
    nrow(bp_get_ready_cohorts(state, stage = "Variety", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max))
  ))

  state
}


# Config, Grid, and Initialization ----------------------------------------
make_wheat_conv_nogs_cfg <- function(n_pyt = 500L) {
  cfg <- list(
    scheme_id = "Wheat_Conv_noGS",
    debug = FALSE,
    debug_after_year = NULL,
    debug_after_tick = NULL,
    debug_where = NULL,

    # Fixed by protocol
    n_cross_block = 50L,
    n_crosses = 100L,
    dh_per_family = 100L,
    cross_block_crossing_duration_years = 0.25,
    dh_duration_years = 1.25,
    n_headrow_advance = 500L,
    n_pyt = as.integer(n_pyt),
    n_pyt_to_ayt = 50L,
    n_pyt_to_cross = 20L,
    n_ayt_to_eyt = 10L,
    n_ayt_to_cross = 10L,
    n_incumbent_keep = 20L,

    # Stage heritabilities from protocol
    h2_headrow = 0.10,
    h2_pyt = 0.20,
    h2_ayt = 0.50,
    h2_eyt = 0.67,

    # Effective replication assumptions used in costing comments
    ayt_effective_reps = 4L,
    eyt_effective_reps = 8L,

    # Placeholder costs
    cost_crossing = 1,
    cost_dh_line = 0.05,
    cost_headrow_line = 0.05,
    cost_plot_pyt = 20,
    cost_plot_ayt = 30,
    cost_plot_eyt = 35,
    cost_parent_recycle = 0.5
  )

  cfg$dt <- 0.25
  cfg$ticks_per_year <- as.integer(round(1 / cfg$dt))
  cfg
}

make_wheat_conv_nogs_param_grid <- function() {
  expand.grid(
    n_pyt = as.integer(c(400L, 500L, 600L)),
    stringsAsFactors = FALSE
  )
}

init_wheat_conv_nogs_sim <- function(cfg) {
  founder_haps <- quickHaplo(nInd = 50, nChr = 3, segSites = 120)
  SP <- SimParam$new(founder_haps)
  SP$addTraitA(nQtlPerChr = 60)

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


# Top-Level Runners --------------------------------------------------------
run_wheat_conv_nogs <- function(n_fill_years = 6L, n_run_years = 12L, cfg = make_wheat_conv_nogs_cfg()) {
  sim <- init_wheat_conv_nogs_sim(cfg)
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
  invisible(state)
}

run_wheat_conv_nogs_grid <- function(n_fill_years = 6L, n_run_years = 12L, grid = make_wheat_conv_nogs_param_grid()) {
  out <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    cfg <- make_wheat_conv_nogs_cfg(n_pyt = grid$n_pyt[i])
    state <- run_wheat_conv_nogs(n_fill_years = n_fill_years, n_run_years = n_run_years, cfg = cfg)

    out[[i]] <- data.frame(
      run_id = i,
      n_fill_years = n_fill_years,
      n_run_years = n_run_years,
      n_pyt = cfg$n_pyt,
      varieties = nrow(state$outputs$varieties),
      total_cost = sum(state$cost_log$total_cost),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}


# Reporting ----------------------------------------------------------------
print_wheat_conv_nogs_setup <- function() {
  cat("\n=== Wheat_Conv_noGS: Parameters encoded from scheme ===\n")
  print(make_wheat_conv_nogs_scheme_table(), row.names = FALSE)

  cat("\n=== Inputs still needed (or currently defaulted) ===\n")
  print(make_required_input_table(), row.names = FALSE)

  cat("\n=== Design questions ===\n")
  qs <- make_design_questions()
  for (i in seq_along(qs)) cat(sprintf("%d. %s\n", i, qs[i]))
}

print_wheat_conv_nogs_summary <- function(state) {
  cat("\nVariety releases:\n")
  print(state$outputs$varieties)

  cat("\nTotal cost by event:\n")
  print(stats::aggregate(total_cost ~ event, data = state$cost_log, sum))

  cat("\nRecent cohorts:\n")
  print(utils::tail(state$cohorts[, c(
    "cohort_id", "stage", "created_tick", "available_tick", "active", "n_ind"
  )], 20))
}


# Script Entry -------------------------------------------------------------
if (identical(environment(), globalenv()) && !identical(Sys.getenv("BPS_SKIP_SCRIPT_ENTRY"), "1")) {
  print_wheat_conv_nogs_setup()

  out <- run_wheat_conv_nogs(n_fill_years = 8L, n_run_years = 32L, cfg = make_wheat_conv_nogs_cfg())
  print_wheat_conv_nogs_summary(out)

  cat("\nGrid example:\n")
  # print(run_wheat_conv_nogs_grid(n_fill_years = 4L, n_run_years = 10L))
}

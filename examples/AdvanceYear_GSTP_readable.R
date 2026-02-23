# GSTP-like scheme with readable yearly schedule + rapid recurrent GS.
#
# Key timing:
# - Global rapid cycle length controls simulation resolution.
# - ticks_per_year = 1 / rapid_cycle_length.
# - PYT/AYT start at beginning of year and are available after 0.5 years.
# - EYT starts at beginning of year and is available after 1.5 years
#   (two years of locations executed as one combined trial block).

library(AlphaSimR)

source("R/operators.R")
source("R/readable_api.R")
source("R/readable_wrappers.R")
source("R/monitoring.R")

make_cross_plan_no_self <- function(n_parents, n_crosses) {
  if (n_parents < 2L) {
    return(matrix(c(1L, 1L), ncol = 2L))
  }
  p1 <- sample.int(n_parents, size = n_crosses, replace = TRUE)
  p2 <- vapply(p1, function(i) {
    sample.int(n_parents - 1L, size = 1L) -> j
    if (j >= i) j <- j + 1L
    as.integer(j)
  }, integer(1))
  cbind(p1, p2)
}

bp_debug_break <- function(state, cfg, where = NULL, year = NULL, tick = NULL) {
  if (!isTRUE(cfg$debug %||% FALSE)) return(invisible(NULL))

  if (is.null(tick)) tick <- as.integer(state$time$tick)
  if (is.null(year) && !is.null(cfg$ticks_per_year)) {
    year <- as.integer(floor(tick / as.integer(cfg$ticks_per_year)) + 1L)
  }
  if (is.null(where)) {
    where <- as.character(sys.call(-1)[[1]])
  }

  if (!is.null(cfg$debug_after_year) && !is.null(year) && year < as.integer(cfg$debug_after_year)) {
    return(invisible(NULL))
  }
  if (!is.null(cfg$debug_after_tick) && tick < as.integer(cfg$debug_after_tick)) {
    return(invisible(NULL))
  }
  if (!is.null(cfg$debug_where) && !(where %in% as.character(cfg$debug_where))) {
    return(invisible(NULL))
  }

  browser()
}

run_recurrent_gs_tick <- function(state, cfg) {
  bp_debug_break(state, cfg)
  src <- get_ready_pop(state, stage = cfg$rc_stage, stream = "main", policy = "latest_one", combine = TRUE, silent = TRUE)
  if (is.null(src)) return(state)

  pop <- src$pop
  model <- state$gs_models[[cfg$model_id]]
  if (is.null(model)) return(state)
  pop <- bp_predict_ebv(pop, model, state, cfg = list(ebv_trait = 1), stage_label = "RC")
  sel <- selectInd(pop, nInd = min(cfg$rc_select_n, pop_n_ind(pop)), use = "ebv", trait = 1, simParam = state$sim$SP)

  plan <- make_cross_plan_no_self(pop_n_ind(sel), cfg$rc_crosses)
  next_pop <- makeCross(sel, crossPlan = plan, nProgeny = cfg$rc_n_progeny_per_cross, simParam = state$sim$SP)

  state <- put_stage_pop(
    state = state,
    pop = next_pop,
    stage = cfg$rc_stage,
    source = src,
    ready_in_years = cfg$rapid_cycle_length,
    stream = "main"
  )
  state <- add_stage_cost(
    state = state,
    stage = cfg$rc_stage,
    event = "rc_crossing",
    unit = "cross",
    n = cfg$rc_crosses,
    unit_cost = cfg$cost_cross
  )
  state <- close_sources(state, src)
  state
}

run_send_rc_to_dh <- function(state, cfg) {
  bp_debug_break(state, cfg)
  src <- get_ready_pop(state, stage = cfg$rc_stage, stream = "main", policy = "latest_one", combine = TRUE, silent = TRUE)
  if (is.null(src)) return(state)

  pop <- src$pop
  model <- state$gs_models[[cfg$model_id]]
  if (!is.null(model)) {
    pop <- bp_predict_ebv(pop, model, state, cfg = list(ebv_trait = 1), stage_label = "RC")
    send <- selectInd(pop, nInd = min(cfg$dh_n_entries_per_year, pop_n_ind(pop)), use = "ebv", trait = 1, simParam = state$sim$SP)
  } else {
    idx <- sample.int(pop_n_ind(pop), size = min(cfg$dh_n_entries_per_year, pop_n_ind(pop)), replace = FALSE)
    send <- pop_subset(pop, idx)
  }
  dh_send <- makeDH(send, nDH = 1, simParam = state$sim$SP)

  state <- put_stage_pop(
    state = state,
    pop = dh_send,
    stage = cfg$dh_stage,
    source = src,
    ready_in_years = 2,
    stream = "main"
  )
  state <- add_stage_cost(
    state = state,
    stage = cfg$dh_stage,
    event = "dh_seed_pipeline",
    unit = "line",
    n = pop_n_ind(dh_send),
    unit_cost = cfg$cost_line_dh
  )
  state
}

run_pyt <- function(state, cfg) {
  bp_debug_break(state, cfg)
  sel_fn <- function(state, src, pop_in, cfg_local) {
    n <- min(as.integer(cfg_local$n_pyt), pop_n_ind(pop_in))
    sample.int(pop_n_ind(pop_in), size = n, replace = FALSE)
  }

  run_phenotype_trial(state, list(
    trial_name = "PYT",
    input_stage = cfg$dh_stage,
    output_stage = "PYT",
    stream = "main",
    input_policy = "latest_one",
    select_entries_fn = sel_fn,
    n_pyt = cfg$n_pyt,
    traits = 1,
    n_loc = 1,
    reps = 1,
    varE = cfg$varE,
    duration_years = 0.5,
    consume_input = TRUE,
    cost_per_plot = cfg$cost_plot_pyt,
    use_env_control = FALSE,
    silent = TRUE
  ))
}

run_genotype_pyt <- function(state, cfg) {
  bp_debug_break(state, cfg)
  run_genotyping(state, list(
    input_stage = "PYT",
    stream = "main",
    input_policy = "latest_one",
    include_not_ready = TRUE,
    chip = cfg$snp_chip,
    duration_years = 0.5,
    cost_per_sample = cfg$cost_genotype,
    silent = TRUE
  ))
}

run_ayt <- function(state, cfg) {
  bp_debug_break(state, cfg)
  sel_fn <- function(state, src, pop_in, cfg_local) {
    model <- state$gs_models[[cfg_local$model_id]]
    pop_scored <- if (is.null(model)) pop_in else bp_predict_ebv(pop_in, model, state, cfg = list(ebv_trait = 1), stage_label = "AYT")
    n <- min(cfg_local$n_ayt, pop_n_ind(pop_scored))
    selected <- if (is.null(model)) {
      selectInd(pop_scored, nInd = n, use = "pheno", simParam = state$sim$SP)
    } else {
      selectInd(pop_scored, nInd = n, use = "ebv", trait = 1, simParam = state$sim$SP)
    }
    match(selected@id, pop_in@id)
  }

  run_phenotype_trial(state, list(
    trial_name = "AYT",
    input_stage = "PYT",
    output_stage = "AYT",
    stream = "main",
    input_policy = "latest_one",
    select_entries_fn = sel_fn,
    model_id = cfg$model_id,
    n_ayt = cfg$n_ayt,
    traits = 1,
    n_loc = 4,
    reps = 2,
    varE = cfg$varE,
    duration_years = 0.5,
    consume_input = FALSE,
    cost_per_plot = cfg$cost_plot_ayt,
    use_env_control = TRUE,
    env_mean_mu = 0,
    env_mean_sd = 1,
    env_year_sd = cfg$ayt_env_year_sd,
    log_per_environment = FALSE,
    log_aggregate = TRUE,
    silent = TRUE
  ))
}

run_eyt <- function(state, cfg) {
  bp_debug_break(state, cfg)
  if (is.null(state$sim$eyt_base_means)) {
    state$sim$eyt_base_means <- stats::rnorm(10, mean = 0, sd = 1)
  }

  # Two-year EYT block in one call:
  # - locations 1..10: base_means + year_shift_1
  # - locations 11..20: base_means + year_shift_2
  # run_phenotype_trial adds per-location perturbation via env_year_sd.
  year_shift_1 <- stats::rnorm(1, mean = 0, sd = cfg$eyt_year_shift_sd)
  year_shift_2 <- stats::rnorm(1, mean = 0, sd = cfg$eyt_year_shift_sd)
  env_means_20 <- c(
    state$sim$eyt_base_means + year_shift_1,
    state$sim$eyt_base_means + year_shift_2
  )

  sel_fn <- function(state, src, pop_in, cfg_local) {
    n <- min(cfg_local$n_eyt, pop_n_ind(pop_in))
    selected <- selectInd(pop_in, nInd = n, use = "pheno", simParam = state$sim$SP)
    match(selected@id, pop_in@id)
  }

  run_phenotype_trial(state, list(
    trial_name = "EYT",
    input_stage = "AYT",
    output_stage = "EYT",
    stream = "main",
    input_policy = "latest_one",
    select_entries_fn = sel_fn,
    n_eyt = cfg$n_eyt,
    traits = 1,
    n_loc = 20,
    reps = 4,
    varE = cfg$varE,
    duration_years = 1.5,
    consume_input = FALSE,
    cost_per_plot = cfg$cost_plot_eyt,
    use_env_control = TRUE,
    env_means = env_means_20,
    env_year_sd = cfg$eyt_env_within_sd,
    log_per_environment = FALSE,
    log_aggregate = TRUE,
    silent = TRUE
  ))
}

run_release_variety <- function(state, cfg) {
  bp_debug_break(state, cfg)
  src <- get_ready_pop(state, stage = "EYT", stream = "main", policy = "latest_one", combine = TRUE, silent = TRUE)
  if (is.null(src)) return(state)

  variety <- selectInd(src$pop, nInd = 1, use = "pheno", simParam = state$sim$SP)
  old_var <- bp_get_ready_cohorts(
    state,
    stage = "Variety",
    stream = "main",
    active_only = TRUE,
    as_of_tick = .Machine$integer.max
  )
  if (nrow(old_var) > 0L) {
    for (i in seq_len(nrow(old_var))) {
      state <- bp_close_cohort(state, old_var$cohort_id[i])
    }
  }
  state <- put_stage_pop(
    state = state,
    pop = variety,
    stage = "Variety",
    source = src,
    ready_in_years = 0,
    stream = "main"
  )
  state$outputs$varieties <- rbind(
    state$outputs$varieties,
    data.frame(
      tick = as.integer(state$time$tick),
      source_cohort_id = paste(src$source_ids, collapse = ";"),
      variety_id = as.integer(variety@id),
      stringsAsFactors = FALSE
    )
  )

  # Prevent repeated release from same EYT cohort.
  state <- close_sources(state, src)
  state
}

run_midyear_gp_update <- function(state, cfg) {
  bp_debug_break(state, cfg)
  run_train_gp_model(state, list(
    from_stage = "PYT",
    stream = "main",
    chip = cfg$snp_chip,
    lookback_years = 50,
    training_policy = "latest_n",
    input_n = cfg$gp_max_pyt,
    model_id = cfg$model_id,
    silent = TRUE
  ))
}

run_start_rc_if_ready <- function(state, cfg) {
  bp_debug_break(state, cfg)
  already <- bp_get_ready_cohorts(state, stage = cfg$rc_stage, stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)
  if (nrow(already) > 0L) return(state)

  src <- get_ready_pop(state, stage = "PYT", stream = "main", policy = "latest_one", combine = TRUE, silent = TRUE)
  if (is.null(src)) return(state)

  n <- min(cfg$rc_init_n, pop_n_ind(src$pop))
  idx <- sample.int(pop_n_ind(src$pop), size = n, replace = FALSE)
  rc0 <- pop_subset(src$pop, idx)

  state <- put_stage_pop(
    state = state,
    pop = rc0,
    stage = cfg$rc_stage,
    source = src,
    ready_in_years = 0,
    stream = "main"
  )
  state
}

run_one_year <- function(state, cfg, year_index) {
  bp_debug_break(state, cfg, year = year_index)
  # Ordered from parents-side operations toward variety release.
  state <- run_start_rc_if_ready(state, cfg)
  state <- run_send_rc_to_dh(state, cfg)
  state <- run_pyt(state, cfg)
  state <- run_genotype_pyt(state, cfg)
  state <- run_ayt(state, cfg)
  state <- run_eyt(state, cfg)
  state <- run_release_variety(state, cfg)

  # Quarterly recurrent GS; retrain GP when PYT + genotype become available (mid-year).
  for (q in seq_len(cfg$ticks_per_year)) {
    state <- run_recurrent_gs_tick(state, cfg)
    state <- bp_advance_time(state, n_ticks = 1L)

    if (q == 2L) {
      state <- run_midyear_gp_update(state, cfg)
    }
  }

  cat(sprintf(
    "year=%d t=%.2f active=%d rc=%d dh=%d pyt=%d ayt=%d eyt=%d models=%d varieties=%d\n",
    year_index,
    state$time$t,
    sum(state$cohorts$active),
    nrow(bp_get_ready_cohorts(state, stage = cfg$rc_stage, stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = cfg$dh_stage, stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "PYT", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "AYT", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    nrow(bp_get_ready_cohorts(state, stage = "EYT", stream = "main", active_only = TRUE, as_of_tick = .Machine$integer.max)),
    length(state$gs_models),
    nrow(state$outputs$varieties)
  ))

  state
}

run_gstp_loop_demo <- function(n_years = 14L, make_plots = FALSE) {
  # debug hooks can be controlled in cfg below.
  founder_haps <- quickHaplo(nInd = 120, nChr = 3, segSites = 80)
  SP <- SimParam$new(founder_haps)
  SP$addTraitA(nQtlPerChr = 20)
  SP$addSnpChip(nSnpPerChr = 40)
  SP$setVarE(h2 = 0.35)

  founders <- newPop(founder_haps, simParam = SP)
  parents0 <- selectInd(founders, nInd = 20, use = "gv", simParam = SP)

  cfg <- list(
    rapid_cycle_length = 0.25,
    model_id = "gp_main",
    snp_chip = 1L,
    make_plots = isTRUE(make_plots),
    debug = FALSE,
    debug_after_year = NULL,
    debug_after_tick = NULL,
    debug_where = NULL,

    rc_stage = "RC",
    rc_init_n = 120L,
    rc_select_n = 60L,
    rc_crosses = 60L,
    rc_n_progeny_per_cross = 2L,

    dh_stage = "DH_PIPE",
    dh_n_entries_per_year = 120L,

    n_pyt = 120L,
    gp_max_pyt = 2L,

    n_ayt = 40L,
    n_eyt = 6L,

    varE = 1.0,
    ayt_env_year_sd = 0.20,
    eyt_year_shift_sd = 0.25,
    eyt_env_within_sd = 0.15,

    cost_cross = 1,
    cost_line_dh = 0.5,
    cost_plot_pyt = 20,
    cost_plot_ayt = 30,
    cost_plot_eyt = 35,
    cost_genotype = 25
  )
  cfg$ticks_per_year <- as.integer(round(1 / cfg$rapid_cycle_length))

  state <- bp_init_state(
    SP = SP,
    dt = cfg$rapid_cycle_length,
    start_time = 0,
    sim = list(default_chip = cfg$snp_chip)
  )
  state$sim$make_plots <- cfg$make_plots

  # Bootstrap: small fast start for debugging.
  f1_plan <- matrix(sample.int(pop_n_ind(parents0), size = 40L * 2L, replace = TRUE), ncol = 2)
  f1_0 <- makeCross(parents0, crossPlan = f1_plan, nProgeny = 1, simParam = SP)
  dh_0 <- makeDH(f1_0, nDH = 4, simParam = SP)
  if (pop_n_ind(dh_0) > cfg$n_pyt) {
    dh_0 <- pop_subset(dh_0, sample.int(pop_n_ind(dh_0), size = cfg$n_pyt, replace = FALSE))
  }
  pyt_0 <- setPheno(dh_0, varE = cfg$varE, reps = 1, traits = 1, simParam = SP)

  state <- put_stage_pop(state, pyt_0, stage = "PYT", source = NULL, ready_in_years = 1, stream = "main")
  state <- run_genotyping(state, list(
    input_stage = "PYT",
    stream = "main",
    input_policy = "latest_one",
    include_not_ready = TRUE,
    chip = cfg$snp_chip,
    duration_years = 0.5,
    cost_per_sample = cfg$cost_genotype,
    silent = TRUE
  ))

  for (yr in seq_len(n_years)) {
    state <- run_one_year(state, cfg, year_index = yr)
  }

  invisible(state)
}

if (identical(environment(), globalenv())) {
  out <- run_gstp_loop_demo(n_years = 14, make_plots = FALSE)
  ticks_per_year <- as.integer(round(1 / out$time$dt))
  metrics <- bp_extract_cohort_metrics(
    state = out,
    stages = c("DH_PIPE", "PYT", "AYT", "EYT", "Variety"),
    trait = 1L,
    origin_stage = "DH_PIPE",
    include_inactive = TRUE,
    ticks_per_year = ticks_per_year
  )

  cat("\nVariety releases:\n")
  print(out$outputs$varieties)

  cat("\nTotal cost by event:\n")
  print(stats::aggregate(total_cost ~ event, data = out$cost_log, sum))

  cat("\nCohort metrics:\n")
  print(metrics[order(metrics$available_year, metrics$stage), c(
    "cohort_id", "stage", "origin_cohort_id", "origin_year", "available_year",
    "mean_gv", "var_gv", "max_gv", "cor_ebv_gv", "h2", "H2"
  )])

  if (isTRUE(out$sim$make_plots) && requireNamespace("ggplot2", quietly = TRUE)) {
    metric_names <- c("mean_gv", "var_gv", "max_gv", "cor_ebv_gv", "h2", "H2")
    for (m in metric_names) {
      df_origin <- bp_summarize_metric_by_year(metrics, metric = m, year_col = "origin_year")
      p_origin <- ggplot2::ggplot(df_origin, ggplot2::aes(x = year, y = value, color = stage, group = stage)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::labs(x = "Origin Year", y = m, title = paste(m, "by Origin Year"))
      print(p_origin)

      df_available <- bp_summarize_metric_by_year(metrics, metric = m, year_col = "available_year")
      p_available <- ggplot2::ggplot(df_available, ggplot2::aes(x = year, y = value, color = stage, group = stage)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::labs(x = "Available Year", y = m, title = paste(m, "by Available Year"))
      print(p_available)
    }
  }
}

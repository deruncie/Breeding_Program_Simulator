# Template: generic custom event verb
#
# Copy this file and edit:
# - function name
# - input/output stage names
# - event action details
# - logging details and costs

input_ACTION_event_ACTION <- function(state, cfg) {
  bp_debug_break(state, cfg)

  # 1) Find input cohort
  input_cohort <- select_latest_available(
    state = state,
    stage = "STAGE",  # keep stage names inline in event verbs
    stream = "main",
    combine = TRUE,
    silent = TRUE
  )
  chk <- bp_skip_if_no_input(state, input_cohort, cfg)
  if (chk$skip) return(chk$state)

  # 2) Build event cohort from input
  event_pop <- input_cohort$pop
  # examples:
  # event_pop <- selectInd(event_pop,nInd=10)

  # 3) EVENT ACTION (user-edit block)
  # examples:
  event_pop <- setPheno(event_pop,h2 = cfg$STAGE_h2)

  # 4) add new pop to state
  # note: updates to event_pop after this point will no longer be recorded in state!
  state <- put_stage_pop(
    state = state,
    pop = event_pop,
    stage = "STAGE",
    source = input_cohort,
    selection_strategy = "FILL DESCRIPTION",
    ready_in_years = cfg$EVENT_duration_years,
    stream = "main",
    inherit_genotypes = TRUE
  )
  output_cohort_id <- bp_last_cohort_id(state)

  # 5) log the event:
  # note: if you use functions like run_phenotype_trial(), run_genotyping(), or run_train_gp_model(),
  # these log their actions directly, so don't double-log these actions
  # also, run_phenotype_trial() calls put_stage_pop() internally

  # Generic phenotype logging example:
  # - traits can be integer codes (1,2,...) or names ("yield_env1","yield_env2",...)
  # - pheno_matrix can include NA values
  # - set drop_na = TRUE to record only observed cells (sparse MET-ready)
  available_tick <- as.integer(state$time$tick + years_to_ticks(state$time$dt, cfg$EVENT_duration_years))
  traits <- cfg$EVENT_traits %||% 1L
  if (!is.null(event_pop@pheno) && length(event_pop@pheno) > 0L) {
    state <- bp_record_pheno(
      state = state,
      cohort_id = output_cohort_id,
      stage = "STAGE",
      individual_id = event_pop@id,
      traits = traits,
      pheno_matrix = event_pop@pheno,
      available_tick = available_tick,
      n_loc = as.integer(cfg$EVENT_n_loc %||% 1L),
      reps = as.integer(cfg$EVENT_reps %||% 1L),
      environment = 0L,
      p_value = NA_real_,
      drop_na = TRUE
    )
  }
  state <- add_stage_cost(
    state = state,
    stage = "STAGE",
    event = "EVENT",
    unit = "plot",
    n = pop_n_ind(event_pop),
    unit_cost = cfg$EVENT_cost_per_plot
  )

  state <- bp_log_event(
    state = state,
    fn = "input_ACTION_event_ACTION",
    event_type = "EVENT name",
    stage = "STAGE",
    source_ids = input_cohort$source_ids,
    output_id = output_cohort_id,
    event_string = sprintf(
      "EVENT name: n=%d",
      pop_n_ind(event_pop)
    ),
    template_string = "EVENT name",
    details = list()
  )

  # 6) Optional source consumption
  if (isTRUE(cfg$consume_input %||% FALSE)) {
    state <- close_sources(state, input_cohort)
  }

  state
}

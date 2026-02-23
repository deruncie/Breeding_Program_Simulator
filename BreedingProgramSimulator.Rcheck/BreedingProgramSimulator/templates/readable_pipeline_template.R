# Editable configuration template for the readable API.
# Copy and modify for your own scheme.

cfg <- list(
  cross = list(
    input_stage = "PARENT",
    output_stage = "F1",
    n_crosses = 100,
    n_progeny_per_cross = 1,
    duration_years = 1,
    cost_per_cross = 1
  ),
  ssd = list(
    input_stage = "F1",
    output_stage = "F5",
    n_lines_per_f1 = 10,
    duration_years = 2,
    consume_input = TRUE,
    cost_per_line = 0.5
  ),
  pyt = list(
    trial_name = "PYT",
    input_stage = "F5",
    output_stage = "PYT",
    traits = 1,
    n_loc = 1,
    reps = 1,
    varE = 1.0,
    # Optional explicit environment control:
    # use_env_control = TRUE,
    # env_means = c(0.1, -0.2, 0.3, -0.1), # one value per location
    # env_year_sd = 0.2,                    # year-to-year variation around env_means
    # log_per_environment = TRUE,
    # log_aggregate = TRUE,
    duration_years = 1,
    consume_input = TRUE,
    # Optional: function(state, cohort_row, pop, cfg) -> indices
    select_entries_fn = NULL,
    cost_per_plot = 20
  ),
  genotype = list(
    input_stage = "PYT",
    chip = 1,
    duration_years = 1,
    cost_per_sample = 25
  ),
  gp = list(
    from_stage = "PYT",
    chip = 1,
    trait = 1,
    lookback_years = 3,
    model_id = "gp_main"
  ),
  recycle = list(
    input_stage = "PYT",
    model_id = "gp_main",
    n_variety = 1,
    n_new_parents = 50,
    parent_stage = "PARENT",
    parent_duration_years = 1,
    replace_all_parents = TRUE,
    consume_input = TRUE,
    cost_per_parent = 0.5
  )
)

# Breeding Program Simulator

BreedingProgramSimulator (BPS) provides a readable, event-driven framework for
simulating overlapping breeding pipelines with AlphaSimR populations.

## Architecture

- State container: `state` (list) with `time`, `cohorts`, `pops`, logs, and models
- Continuous-time ticks: `state$time$tick`, `state$time$dt`, `state$time$t`
- Cohort snapshots: one POP per cohort with `created_tick` and `available_tick`
- External logs: phenotype, genotype, cost, and variety outputs

## Core API

- `bp_init_state()`
- `select_latest_available()`, `get_ready_pop()`, `put_stage_pop()`
- `run_phenotype_trial()`
- `run_genotyping()`
- `run_train_gp_model()`
- `run_predict_ebv()`
- `bp_advance_time()` and `bp_advance_time_years()`

Stage-producing functions preserve cohort provenance and append event, cost,
phenotype, and genotype records to the state.

## Monitoring API

- `bp_monitor_cohorts()`
- `bp_extract_cohort_metrics()`
- `bp_summarize_metric_by_year()`
- `bp_plot_metric_by_year()`
- `bp_report_stage_metrics()`
- `bp_event_timeline_df()` and `bp_print_event_timeline()`

## Multi-Trait and Synthetic Traits

BPS keeps AlphaSimR biological traits in their standard population slots and
stores derived synthetic traits separately in `pop@misc`. Define and register a
synthetic trait with `bp_synthetic_trait()` and
`bp_register_synthetic_traits()`. Trials, genomic prediction, selection,
baseline scaling, and stage reporting can then use the registered definition
without changing AlphaSimR's biological trait numbering.

```r
index <- bp_synthetic_trait(
  name = "Index",
  traits = 1:2,
  fun = AlphaSimR::selIndex,
  args = list(b = c(0.5, 0.5)),
  linear = TRUE
)

state <- bp_register_synthetic_traits(state, index)
```

Use `bp_select_synthetic()` for selection and
`bp_report_stage_metrics()` for consistently named one-row result fields.

## Examples

- `examples/simple_readable_cycle.R`
- `examples/AdvanceYear_GSTP_readable_structured.R`
- `examples/shiny_monitor_app/app.R`

Run the Shiny dashboard from repo root with:

```r
shiny::runApp("examples/shiny_monitor_app")
```

## Codex Skills

BPS includes skills for drafting breeding schemes and designing experiments.
Install them for your Codex user account with:

```r
BreedingProgramSimulator::bp_install_codex_skills()
```

To make the skills available only inside one project, run:

```r
BreedingProgramSimulator::bp_install_codex_skills(
  scope = "project",
  project = "/path/to/project"
)
```

After upgrading BPS, pass `overwrite = TRUE` to install the skill versions
bundled with the new package version. Codex normally detects skill changes
automatically; restart Codex if they do not appear.

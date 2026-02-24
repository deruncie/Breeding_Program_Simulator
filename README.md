# Breeding Program Simulator (Readable AlphaSimR Workflow)

List-based, event-driven breeding simulation built around AlphaSimR POP objects.

## Architecture

- State container: `state` (list) with `time`, `cohorts`, `pops`, logs, and models
- Continuous-time ticks: `state$time$tick`, `state$time$dt`, `state$time$t`
- Cohort snapshots: one POP per cohort with `created_tick` and `available_tick`
- External logs: phenotype, genotype, cost, and variety outputs

## Core API

- `bp_init_state()`
- `get_ready_pop()`, `put_stage_pop()`, `close_sources()`
- `run_phenotype_trial()`
- `run_genotyping()`
- `run_train_gp_model()`
- `bp_advance_time()`

## Monitoring API

- `bp_monitor_cohorts()`
- `bp_extract_cohort_metrics()`
- `bp_summarize_metric_by_year()`
- `bp_plot_metric_by_year()`

## Examples

- `/Users/deruncie/Box Sync/DER_projects/Optimized_GP/Breeding_Program_Simulator/examples/simple_readable_cycle.R`
- `/Users/deruncie/Box Sync/DER_projects/Optimized_GP/Breeding_Program_Simulator/examples/AdvanceYear_GSTP_readable.R`
- `/Users/deruncie/Box Sync/DER_projects/Optimized_GP/Breeding_Program_Simulator/examples/shiny_monitor_app/app.R`

Run the Shiny dashboard from repo root with:

```r
shiny::runApp("examples/shiny_monitor_app")
```

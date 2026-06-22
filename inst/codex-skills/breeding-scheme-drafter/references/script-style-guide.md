# Breeding Scheme Script Style Guide

This guide defines scheme-script style for BreedingProgramSimulator (BPS) 0.2.0.
The priorities are simplicity, readability, explicit configuration, and faithful
cohort flow.

## Contents

- [Template-First Workflow](#template-first-workflow)
- [1. Purpose and File Boundary](#1-purpose-and-file-boundary)
- [2. Scheme File Layout](#2-scheme-file-layout)
- [3. Configuration Is Explicit](#3-configuration-is-explicit)
- [4. Prefer BPS and AlphaSimR](#4-prefer-bps-and-alphasimr)
- [5. Temporary Populations and Auxiliary Data](#5-temporary-populations-and-auxiliary-individual-data)
- [6. Event Verb Names](#6-event-verb-names)
- [7. Event Verb Pattern](#7-event-verb-pattern)
- [8. Reporting](#8-reporting)
- [9. Scheduler Organization](#9-scheduler-organization)
- [10. Multiple Streams](#10-multiple-streams)
- [11. Genomic Prediction and Debug Runs](#11-genomic-prediction-and-debug-runs)
- [12. Minimal Smoke Tests](#12-minimal-smoke-tests)
- [13. Validation](#13-validation)

## Template-First Workflow

Before writing a scheme, inspect the versioned templates installed under
`templates/bps-0.2.0`. Copy the closest scheme when its biological flow and
scheduler are substantially similar; then rename and remove what is not needed.
Draft from scratch when no template fits. Do not force a user's design into a
template merely to save code.

Locate the templates with:

```r
system.file("templates", "bps-0.2.0", package = "BreedingProgramSimulator")
```

Read `TEMPLATE_INDEX.md` in that directory for the selection rules.

## 1) Purpose and File Boundary

A scheme file is loaded with `source()`. Loading it defines its reporting
function, event verbs, and `run_simulation()`; it does not initialize a
population or start a simulation by itself.

Keep these responsibilities outside the scheme file:

- genetic architecture and founder creation
- historical simulation and burn-in configuration
- scenario and replicate grids
- experiment comparisons
- default parameter values

Put those responsibilities in `Create_sim_bps.R`, `run_experiments.R`, and
their cfg files.

## 2) Scheme File Layout

Use this section order:

1. `Setup`: short scheme description, flow diagram, libraries, and essential
   `source()` calls.
2. `Helper Utilities`: only small helpers that remove genuine repetition.
3. `Reporting`: one `record_yearly_outputs()` function.
4. `Event Verbs`: visible biological and operational actions.
5. `Runners`: one explicit `run_simulation()` scheduler.

Do not add Config, Initialization, Script Entry, or demo-run sections to a
scheme file. Keep the main scheduler readable from top to bottom.

## 3) Configuration Is Explicit

Use `cfg$...` for every quantity or choice that a user may vary, including:

- stage sizes and selection intensities
- durations and cycle lengths
- trial traits, locations, reps, environmental variation, and costs
- crossing, selfing, DH, and seed-increase settings
- marker chips, training windows, models, and GP debug settings
- reporting cadence

Do not hide defaults inside scheme code with `%||%`, `if (is.null(...))`,
or fallback literals. A missing required cfg parameter should fail at its use
site.

Use `bp_check_cfg_requirements()` outside the scheme to collect required
parameters into a cfg file. Script-assigned derived fields may appear as
`NULL` in that file with a comment explaining where they are assigned.

Keep stage names, stream names intrinsic to the design, and event terminology
in the scheme code. Put them in cfg only when the experiment is intended to
vary them.

Do not add user-error checks merely to replace ordinary R errors. Validate
biological constraints only when an invalid combination would otherwise run
silently or fail misleadingly.

## 4) Prefer BPS and AlphaSimR

Use BPS and AlphaSimR functions directly whenever they express the operation:

- `select_latest_available()`
- `bp_skip_if_no_input()`
- `put_stage_pop()` and `bp_update_stage_pop()`
- `bp_set_misc_values()`
- `run_phenotype_trial()`
- `run_progeny_test()`
- `run_genotyping()`
- `run_train_gp_model()` and `run_predict_ebv()`
- `bp_select_synthetic()`
- `add_stage_cost()`
- `bp_advance_time_years()`
- AlphaSimR crossing, selfing, DH, phenotype, and selection functions

Avoid wrappers around one BPS or AlphaSimR call. Add a helper only when it
makes repeated code materially clearer. Keep scheme logic in event verbs, not
in callback layers or generic mini-frameworks.

This restriction applies to helper utilities, not to meaningful event verbs.
An event verb may have one principal BPS call when it also makes the input
choice, biological decision, timing, provenance, and scheduler position clear.

## 5) Temporary Populations and Auxiliary Individual Data

Store a population in `state` only when its individuals remain available to a
later stage or must be retained as a cohort. A population created only to
collect data can remain temporary: generate it, measure it, summarize its data
onto the retained population, and discard it without calling `put_stage_pop()`.

A progeny test is the standard example. The progeny are grown and phenotyped
only to evaluate the original candidates. Keep the candidates, attach the
family summaries to them, and discard the progeny. Prefer
`run_progeny_test()`, which performs this workflow without storing the progeny
as cohorts.

Store per-individual quantities that are not biological phenotypes, genetic
values, or EBVs as clearly named entries in `pop@misc`. Use:

```r
pop <- bp_set_misc_values(pop, "progeny_mean", progeny_mean)
```

Names beginning with `bps_` are reserved for BPS-managed fields, including
synthetic-trait storage.

The values must follow the current individual order. AlphaSimR then carries
them through population subsetting and selection. Select directly on a stored
quantity when appropriate:

```r
selected <- AlphaSimR::selectInd(
  pop,
  nInd = cfg$n_select,
  use = function(pop, trait = NULL) pop@misc[["progeny_mean"]],
  simParam = state$sim$SP
)
```

If `pop` is already stored in `state`, attach the data and then replace the
stored copy with `bp_update_stage_pop()`. For a new cohort, attach the data
before `put_stage_pop()`. Do not use `bp_update_stage_pop()` merely to edit a
local population.

## 6) Event Verb Names

Name functions for the action and its biological source/output. Follow the
terminology in the user's description and diagram.

Examples:

- `release_variety_from_EYT()`
- `select_from_AYT_and_start_EYT()`
- `select_from_PYT_and_run_AYT()`
- `advance_F1_to_headrow()`
- `build_F1_from_Parents()`
- `train_wfRGS_model()`
- `run_wfRGS_cycle()`

Use consistent stage capitalization within a scheme. Pass `state`, `cfg`,
and `year` to scheduled verbs. Add `cycle` and `stream` only when the
event needs them.

## 7) Event Verb Pattern

### Choosing Event Boundaries

An event verb represents one planned transition from currently available stored
inputs to the next meaningful stored output. It may represent a process with
duration: the scheduler invokes it once, and BPS records when its output becomes
available.

Identify persistent populations before defining event verbs. Store a population
when it contains independently available results or material needed for later
selection, reporting, GP training, branching, or reuse. Then minimize stored
populations between those boundaries by keeping one-use intermediates local to
the event. See Section 5 for temporary data-only populations.

In the examples below, PYT, AYT, and EYT are illustrative names for sequential
field-trial stages—often preliminary, advanced, and elite yield trials. They are
not BPS keywords; use and clarify the breeding program's own terminology.

A trial event combines planning and execution. Planning selects entries from
the preceding trial results and defines the trial size and design. Execution
runs the trial and produces results after its cfg-defined duration:

```text
stored PYT results
-> temporary selection of AYT entries
-> run AYT
-> stored AYT results
```

This leads to verbs such as `select_from_PYT_and_run_AYT()` and
`select_from_AYT_and_run_EYT()`. Do not store `AYT_selected` unless it has an
independent lifecycle; pass it directly into `run_phenotype_trial()`, which
creates the stored AYT cohort. This avoids storing the same lines first as
selected seed and again as a phenotyped trial population, or modifying a
stored cohort's population and availability.

Do not combine PYT, AYT, and EYT into one event. Each completed trial is a
decision and data-availability boundary: AYT cannot be planned until PYT
results are available, and EYT cannot be planned until AYT results are
available. Obtain selection inputs through availability-aware BPS functions.

Several consecutive operations may share one event when they can all be
planned at once, no decisions or branching occur between them, intermediate
populations are not used elsewhere, and the final output's availability can be
calculated. Single-seed descent is a typical example. Simulate every required
generation, but keep intermediate populations temporary and store only the
final generation.

Split actions when an intermediate population drives a decision, branches or
feeds another event, contributes data needed for reporting or training, or
must become independently available.

Write event verbs in this order:

1. Select the required input cohort or cohorts.
2. Skip when the pipeline has not produced input yet.
3. Perform the visible AlphaSimR/BPS operation.
4. Store or update output with complete provenance and logging metadata.
5. Return the updated state.

Use:

```r
input_pyt <- select_latest_available(
  state = state,
  stage = "PYT",
  stream = "main",
  n = 1L,
  combine = TRUE,
  silent = TRUE
)
chk <- bp_skip_if_no_input(
  state = state,
  input_obj = input_pyt,
  cfg = cfg,
  event_name = "select_from_PYT_and_run_AYT"
)
if (chk$skip) return(chk$state)
```

This skip handles an expected empty pipeline during warm-up. It is not a
substitute for missing cfg parameters.

Use `source = input_bundle` or explicit `source_ids` for every output. Use
`"UNKNOWN"` only when the true source cannot be represented. Supply accurate
`selection_strategy`, `cross_strategy`, duration, stream, and cost metadata.

## 8) Reporting

Define one `record_yearly_outputs(results, results_base, state, year, cfg)`.

Build one readable result record containing:

- year or cycle time
- stage metrics from `bp_report_stage_metrics()`
- biological-trait and synthetic-trait metrics separately when needed
- crossing, line-development, seed-increase, phenotype, genotype, and total
  costs that matter for the comparison

Use stable names across schemes so experiment outputs can be row-bound.
Request trait suffixes deliberately for multi-trait results. Return:

```r
bind_rows(results, data.frame(results_base, results_year))
```

Do not duplicate BPS reporting logic in scheme-local helpers unless BPS cannot
express the required metric.

## 9) Scheduler Organization

Keep the complete logical schedule visible inside `run_simulation()`.

For annual line-development pipelines, order same-time events from downstream
to upstream: release completed material, update/select later trials, advance
earlier stages, then create the newest cohort. Advance time only after all
events at that time have run.

For rapid cycles, keep the outer annual or phase structure visible and use a
small inner cycle loop. Advance BPS time by the cfg-controlled cycle duration
at the point the biological material becomes available.

Record the initial state when useful, then report at consistent year/cycle
boundaries. Use one comparison origin via `start_year`.

Do not put burn-in loops in the scheme. `Create_sim_bps.R` calls the same
`run_simulation()` with a burn-in cfg.

## 10) Multiple Streams

Represent parallel pipelines with BPS `stream` values, not separate state
objects.

- Pass `stream` explicitly to stream-specific event verbs.
- Loop over streams in the scheduler when the same event applies to each.
- Keep stream-specific cohorts in their own streams.
- When streams merge, select one input bundle per stream and pass all source
  IDs to the merged output.
- Give merged outputs a deliberate stream name.
- Report stream-specific metrics separately and include stream identity when
  it is not fixed by the scheme.

Keep simple stream-name construction visible. Do not copy the old continuous
wfRGS script's unrelated helpers or hidden defaults.

## 11) Genomic Prediction and Debug Runs

A GP scheme must expose:

```r
cfg$debug_GP
cfg$debug_GP_n
```

When `cfg$debug_GP` is true, fit the same model pathway on the smallest valid
subset, capped by `cfg$debug_GP_n`. Continue through model storage,
prediction, selection, and downstream events. Do not replace GP with random
scores: the smoke test must verify the real training/prediction integration.

Production behavior must use the full cfg-defined training population.
Neither debug value may have a hidden default in the scheme.

## 12) Minimal Smoke Tests

Before a scheme is considered working, propose and help the user run a tiny
cfg that exercises every event path.

Use approximately:

- 2 chromosomes
- a few QTL and markers per chromosome
- the smallest founder and parent populations valid for the crossing design
- small stage populations, while keeping candidate counts larger than selected
  counts so selection actually occurs
- one or very few locations and reps where valid
- only enough years/cycles to reach the final stage, train/predict any GP
  model, and exercise recycling or stream merging

Derive exact sizes from the scheme's operations. Check requirements such as
parent pairs, cross-plan dimensions, training sample size, and requested
selection counts. Small must still be biologically and computationally valid.

For GP smoke tests, set `cfg$debug_GP = TRUE` and choose the smallest
`cfg$debug_GP_n` that allows the selected model to fit.

Inspect runtime, event timeline, stage availability, source links, selection,
costs, genotyping/model state, results columns, and multi-stream convergence.
After the smoke test passes, restore the user's production cfg; never let test
values become hidden production defaults.

## 13) Validation

Validate in this order:

1. Source the scheme without running it.
2. Run `bp_check_cfg_requirements()` with the relevant orchestration file.
3. Run the minimal smoke test.
4. Check event chronology, source/output links, durations, and streams.
5. Check selection counts and recycling logic.
6. Check costs for omissions or double counting.
7. Check genotyping, GP training/prediction, and synthetic-trait behavior.
8. Check reporting names and result row binding.
9. Compare the event flow and scheduler with the user's diagram or protocol
   when supplied.

Report concrete mismatches by function, stage, stream, and time.

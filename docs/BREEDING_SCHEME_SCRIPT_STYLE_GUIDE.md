# Breeding Scheme Script Style Guide

This guide defines the standard style for readable breeding-scheme scripts in this package.

## 1) Script Structure

Use this section order in each scheme script:

1. `Setup` (libraries + `source()` calls)
2. `Helper Utilities` (small pure helpers)
3. `Event Verbs` (scheme actions)
4. `Yearly Schedule` (or tick schedule)
5. `Config` (all parameters)
6. `Initialization` (founders + initial cohorts)
7. `Runners` (fill loop + run loop)
8. `Reporting` (summary and plots)
9. `Script Entry` (guarded execution)

Keep each section short and readable. Avoid hidden control flow.

## 2) Event Verb Naming

Write event functions as explicit action names that encode source + action, for example:

- `select_from_PYT_and_run_AYT()`
- `advance_F1_to_DH()`
- `update_parents_and_cross_to_F1()`

## 3) Event Verb Implementation Pattern

Each event verb should follow this readable pattern:

1. Pull input cohorts with `select_latest_available()` or `select_current()`.
2. Handle missing input with:

```r
chk <- bp_skip_if_no_input(state, input_obj, cfg)
if (chk$skip) return(chk$state)
state <- chk$state
```

3. Run explicit stage logic in ordinary R/AlphaSimR code. This may include:
- direct crossing/selfing/subsetting code
- `run_phenotype_trial(...)`
- `run_genotyping(...)`
- `run_train_gp_model(...)`
- other wrapper/core API calls

4. Ensure output/state changes are recorded through core logging-capable verbs (for example `put_stage_pop()`, `add_stage_cost()`, and wrapper calls that already log).
5. Close sources only when one-time consumption is intended.

Phenotyping events are one common instance of this pattern:
- select entries
- run trial
- do next-stage advancement/selection in the following event verb

Prefer explicit event code over callback-heavy generic wrappers.
## 4) Time and Scheduling

- Use `rapid_cycle_length` (years per tick) as primary time parameter.
- Derive `ticks_per_year <- as.integer(round(1 / rapid_cycle_length))`.
- Keep year-level schedule explicit, with optional inner tick loops when needed.
- Separate warm-up (`fill`) years from measured (`run`) years in two loops.

## 5) Input Selection Rules

Use `select_latest_available(..., n = k)` for normal stage chaining.

Conventions:

- `n = 1` (default): one latest cohort.
- `n > 1`: combine selected cohorts when needed.
- Always pass source cohort IDs to downstream logging via `source` or `source_ids`.

If source IDs are unknown, use `"UNKNOWN"` instead of dropping provenance.

## 6) Logging Requirements

Log all stage outputs through core API verbs so the event timeline is complete.

Minimum required metadata for readable logs:

- source cohort IDs (`source_ids`, vector allowed)
- `selection_strategy` text when selection is used
- `cross_strategy` text when crossing is used
- output stage and availability timing

Missing-input behavior should log by default and can be made strict with config:

- `cfg$log_missing_input = TRUE` (default)
- `cfg$fail_on_missing_input = TRUE` for strict runs

## 7) Naming and Readability

- Use meaningful object names (`input_pyt`, `selected_ayt`, `next_cross_block`).
- Keep event functions short enough to read top-to-bottom.
- Keep cfg values near use-sites when that improves clarity.
- Do not expose internal state-table plumbing in high-level event code.

## 8) Script Entry Guard

Use this entry guard for all examples:

```r
if (identical(environment(), globalenv()) && !identical(Sys.getenv("BPS_SKIP_SCRIPT_ENTRY"), "1")) {
  # demo run
}
```

This keeps scripts runnable interactively while allowing safe `source()` in tests.

## 9) Validation Checklist

Before finalizing a scheme script:

1. Run with short fill/run years and confirm no runtime errors.
2. Check event log chronology and source/output links.
3. Confirm stage availability timing matches intent.
4. Confirm costs are not double-counted.
5. Confirm genotyping availability matches chip/model assumptions.

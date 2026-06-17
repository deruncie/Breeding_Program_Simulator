# BPS Experiment Workflow

## Contents

1. File responsibilities
2. Create_sim_bps.R organization
3. run_experiments.R organization
4. Baseline and trial calibration
5. Configuration generation
6. Minimal smoke tests
7. Validation

## 1) File Responsibilities

First read the package template index:

```r
template_dir <- system.file(
  "templates",
  "bps-0.2.0",
  package = "BreedingProgramSimulator"
)
```

Copy the closest architecture and experiment starters. Use the full scheme
references only when their biological flow is close to the requested design.
Remove irrelevant template blocks; do not grow a script by stacking examples.

Scheme files define event verbs, reporting, and `run_simulation()`.
`Create_sim_bps.R` creates founder architecture, historical state, initial
parents, and burn-in states. `run_experiments.R` reloads an unchanged saved
state and compares schemes or cfg variations.

Use flat cfg lists unless requested otherwise. Required values belong in cfg,
not fallback expressions in scripts.

## 2) Create_sim_bps.R Organization

Use commented sections for Setup, Command-Line Arguments, Base Config / Config
Check, Scenario Grid, Scenario Table, Founder Population, Trait Architecture,
Initial Parents, Calibration Plan, Burn-In, and Save State.

At each substantial section, add a short TODO-style block explaining what the
user should customize. Match the instructional blocks in the single-trait
reference Create_sim_bps.R.

Use explicit SetupID, ScenarioID, and ReplicateID rules. A one-scenario design
may contain only ScenarioID 1; do not invent a large grid.

## 3) run_experiments.R Organization

Use commented sections for Setup and Scheme List, Experiment Config / Config
Check, Command-Line Arguments, Load Starting State, Initial-State Calibration,
Results Metadata, Base Experiment, Additional Experiment Blocks, and Output /
Cost Logging.

Before each section, explain what the user edits. Document how to add a scheme,
create an ExperimentID, reset cfg, record varied parameters, preserve the
starting state, and name outputs.

Every method visibly sources its scheme, labels the method, resets
`cfg <- cfg_base`, and calls `run_simulation()` with the same input
`state`, seed convention, and comparison start. Never use another method's
result state as input.

## 4) Baseline and Trial Calibration

Treat calibration as part of the experimental design. Ask what observations
best match the breeder's experience and what reliability or heritability they
expect at each relevant trial stage.

As a starting reference, use the most recent cohort of new candidates measured
in a single-location trial. Use that cohort to set trait baselines with
`bp_set_trait_baseline()` and to calibrate residual variance/noise to the
breeder's expected candidate-trial reliability. Keep the reference stage,
traits, locations, reps, target reliability, and resulting noise parameters in
cfg or saved setup metadata.

When the breeder has better multi-environment evidence, use it instead. Match
the calibration to the actual analysis scale, including locations, reps,
residual correlations, year effects, and GxE/environment structure. Do not
equate single-plot heritability, entry-mean reliability, and multi-environment
reliability without accounting for the trial model.

Repeat the same calibration procedure at several explicit burn-in years when
the required candidate cohort is available. This keeps baselines and noise
compatible with changing genetic variance during burn-in. Put the years or
cadence in cfg; do not silently recalibrate every year.

At the comparison start, apply the final calibration once and freeze it across
all methods unless calibration itself is the experimental treatment. Never
calibrate methods separately from their future outcomes. Record the reference
cohort/time, targets, and final calibrated values in `state$sim$setup` and
results metadata.

For a minimal smoke test, one successful calibration is sufficient to verify
the path. Production burn-in should exercise the requested repeated schedule.

## 5) Configuration Generation

For setup, scan `Create_sim_bps.R` and the burn-in scheme. For experiments,
scan `run_experiments.R` and all included schemes. Use
`bp_check_cfg_requirements(..., rewrite_file = TRUE)`.

Organize cfg into traits/synthetic traits, genome architecture,
founder/history, scenarios, populations, timing, trials, crossing/line
development, genotyping/GP, costs, scheme-specific values, and script-assigned
derived fields. Do not fill unresolved values silently.

Include calibration stage/traits, target reliabilities, reference trial design,
residual/environment correlations, and burn-in calibration years in the cfg
blocks.

## 6) Minimal Smoke Tests

Construct a separate visible tiny cfg. Start near two chromosomes, 3–10 QTL
and 5–20 markers per chromosome, 8–20 founders when valid, at least four
parents for ordinary random crossing, small candidate stages larger than their
selected subsets, one location/rep unless environment behavior is under test,
and the shortest durations/years that traverse all paths.

These are starting points, not defaults. Derive valid values from cross plans,
family sizes, selection counts, training needs, streams, and checks.

GP schemes require `debug_GP = TRUE` and a smallest-valid `debug_GP_n`.
Debug mode must fit and use the requested model. Run one setup/replicate and
base ExperimentID first; inspect state, timeline, cohorts, costs, models, and
CSVs.

## 7) Validation

Confirm scripts parse; cfg scans are complete; state create/reload succeeds;
schemes share the same starting state; calibration uses the intended cohort and
is frozen across comparisons; repeated burn-in calibrations are logged; seeds
and IDs are reproducible; sweeps record varied values; stream merges preserve
sources; tiny tests perform real selection; GP training/prediction execute;
and output names are stable.

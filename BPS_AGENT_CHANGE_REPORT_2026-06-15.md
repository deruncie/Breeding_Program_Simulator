# BPS Agent Change Report

Date: 2026-06-15

This report describes the BPS APIs added or changed for synthetic traits,
multi-trait analysis, genomic prediction, selection, and reporting. It is
intended for an agent updating breeding-scheme and experiment scripts.

## Executive Summary

BPS now supports synthetic traits defined as functions of one or more
AlphaSimR biological traits. A synthetic trait can be:

- scored from trial phenotypes;
- used as a genomic-model response;
- predicted and stored as a synthetic EBV;
- used for selection;
- assigned an exact or Monte Carlo genetic value;
- included in trait baselines and stage reports.

Synthetic values are stored in `pop@misc`. They are not appended to
`pop@gv`, `pop@pheno`, or `pop@ebv`, so ordinary AlphaSimR calls using
`trait = 1`, `trait = 1:2`, and similar numeric trait specifications continue
to refer only to biological traits.

Stage reporting has been consolidated into one public function:
`bp_report_stage_metrics()`. It returns a flat named list of unnamed numeric
scalars designed to be combined incrementally into one-row result records.

## Synthetic Trait Definitions

Create a definition with `bp_synthetic_trait()`:

```r
index_def <- bp_synthetic_trait(
  name = "Index",
  traits = 1:2,
  fun = AlphaSimR::selIndex,
  args = list(b = cfg$traitWeights),
  linear = TRUE
)
```

The function receives a matrix containing the selected biological traits and
must return one numeric value per individual.

Set `linear = TRUE` only when the function is linear in its inputs. Linear
synthetic genetic values are calculated directly from biological GVs.
Nonlinear synthetic genetic values require Monte Carlo integration.

Register definitions in the state when scripts will refer to them by name:

```r
state <- bp_register_synthetic_traits(state, index_def)
```

`bp_init_state()` now initializes:

```r
state$sim$synthetic_traits
```

Most synthetic-trait arguments accept either a definition object or a
registered name such as `"Index"`.

## Synthetic Value Storage

Synthetic values are stored in matrices under:

```r
pop@misc$bps_synthetic_pheno
pop@misc$bps_synthetic_ebv
pop@misc$bps_synthetic_gv
```

Each matrix has one row per individual and one named column per synthetic
trait.

Use these helpers instead of accessing `@misc` directly:

```r
pop <- bp_set_synthetic_values(pop, "Index", values, type = "ebv")

values <- bp_get_stored_synthetic_values(
  pop,
  "Index",
  type = "ebv",
  missing_ok = TRUE
)
```

`bp_synthetic_values()` evaluates a definition directly from the selected
biological `gv`, `pheno`, or `ebv` matrix:

```r
index_from_ebv <- bp_synthetic_values(pop, index_def, use = "ebv")
```

`bp_trait_matrix()` is a light helper for extracting biological trait matrices.
It safely handles empty AlphaSimR EBV or phenotype matrices.

## Phenotype Trials

`run_phenotype_trial()` now accepts `synthetic_traits`:

```r
state <- run_phenotype_trial(
  state = state,
  pop = entries,
  output_stage = "PYT",
  traits = 1:2,
  synthetic_traits = "Index",
  n_loc = 4,
  reps = 2,
  varE = cfg$varE,
  duration_years = 1
)
```

The requested biological traits must all be measured by the trial.

The aggregate synthetic phenotype is stored in
`pop@misc$bps_synthetic_pheno`. Synthetic phenotype-log entries use trait
labels such as:

```r
"synthetic:Index"
```

When per-environment logging is active, BPS also logs the synthetic score for
each environment. The aggregate score is logged with environment `0`.

## Synthetic Genetic Values

For a linear definition:

```r
gv <- bp_get_synthetic_gv(pop, index_def, state)
```

returns the direct function of biological GVs.

For a nonlinear definition, BPS estimates expected synthetic performance over
the target population of environments. Defaults are:

- 100 random trials;
- 10 plants per trial;
- common trial and residual draws across all individuals.

Common random numbers make comparisons among individuals less noisy. Trial
effects use AlphaSimR trait `envVar` and individual GxE slopes. Plant residuals
use the supplied `varE`, including residual covariance when a matrix is given.

```r
gv <- bp_get_synthetic_gv(
  pop,
  nonlinear_def,
  state,
  varE = cfg$varE,
  n_trials = 100,
  n_plants_per_trial = 10,
  seed = 123
)
```

If several reports will use the same nonlinear synthetic GV, materialize it
once:

```r
pop <- bp_materialize_synthetic_gv(
  pop,
  synthetic_traits = "Utility",
  state = state,
  varE = cfg$varE
)
```

This stores the values in `pop@misc$bps_synthetic_gv`. Existing cached values
are reused unless `overwrite = TRUE`.

Cached synthetic GVs belong to that population state. Recalculate them after
changing the population or after changing the synthetic definition, `varE`,
or Monte Carlo design.

## Genomic Model Training and Prediction

`run_train_gp_model()` can train AlphaSimR RR-BLUP on a cached synthetic
phenotype:

```r
state <- run_train_gp_model(
  state,
  list(
    from_stage = "PYT",
    chip = 1L,
    response = "synthetic_pheno",
    synthetic_trait = "Index",
    model_id = "index_rrblup"
  )
)
```

Every training population must already contain the requested synthetic
phenotype. Normally it is created by `run_phenotype_trial()`.

`run_predict_ebv()` recognizes synthetic model metadata:

```r
state <- run_predict_ebv(
  state,
  list(
    cohort_ids = candidate_ids,
    model_id = "index_rrblup",
    chip = 1L
  )
)
```

Predictions are stored in `pop@misc$bps_synthetic_ebv`, not `pop@ebv`.

A custom `predict_ebv_fn` may return:

```r
list(
  trait_ebv = biological_ebv_matrix,
  synthetic_ebv = list(
    Index = index_predictions,
    Utility = utility_predictions
  )
)
```

It may also return a named synthetic-EBV matrix. For a model whose response is
synthetic, a simple vector or one-column matrix is stored as that model's
synthetic EBV.

## Selection

Use `bp_select_synthetic()` to select on a synthetic value:

```r
selected <- bp_select_synthetic(
  pop,
  n_select = cfg$nParents,
  synthetic_trait = "Index",
  use = "ebv",
  state = state
)
```

For `use = "ebv"`, selection prefers a directly stored synthetic EBV. It can
fall back to applying the synthetic function to biological EBVs.

For `use = "pheno"`, it prefers the stored synthetic phenotype.

For `use = "gv"`, it prefers a cached synthetic GV and otherwise calculates
one. Nonlinear GV selection therefore needs `state` and usually `varE`.

`bp_select_by_index()` and `bp_set_indexed_ebv()` are no longer part of the
public API. Define a linear `bp_synthetic_trait()` for weighted indices and
use `bp_select_synthetic()` for selection. This avoids adding synthetic
columns to AlphaSimR biological slots.

## Reporting API

The public reporting interface is now:

```r
bp_report_stage_metrics()
```

The following public functions were removed:

```r
bp_report_stage_metric()
bp_report_stage_columns()
```

Valid metrics are:

```r
"meanG"
"maxG"
"varG"
"H2"
"accEBV"
"wf_accEBV"
```

`accGV` is not a metric. `accEBV` falls back to phenotype when no usable EBV
is present. If both EBV and phenotype matrices have zero columns or contain no
usable values, accuracy returns `NA_real_`.

### Return contract

`bp_report_stage_metrics()` always returns a flat named list. Every element is
one unnamed numeric scalar, including missing results:

```r
list(
  meanCandidates_Trait1 = 0.42,
  meanCandidates_Trait2 = 0.18
)
```

This return type is intentional. It avoids compound vector names when reports
are accumulated before constructing a one-row data frame.

### Metric aliases

Name each metric to control its output base name:

```r
report <- bp_report_stage_metrics(
  state,
  stage = "Headrow",
  metrics = c(
    meanCandidates = "meanG",
    bestCandidates = "maxG",
    varCandidates = "varG"
  ),
  traits = 1L
)
```

With one returned trait and `append_trait = "auto"` (the default), this yields:

```r
list(
  meanCandidates = ...,
  bestCandidates = ...,
  varCandidates = ...
)
```

### Multi-trait biological reporting

Call biological traits together and force trait suffixes:

```r
biological <- bp_report_stage_metrics(
  state,
  stage = "Candidates",
  metrics = c(
    meanCandidates = "meanG",
    varCandidates = "varG"
  ),
  traits = 1:2,
  append_trait = "always"
)
```

Expected names:

```r
meanCandidates_Trait1
meanCandidates_Trait2
varCandidates_Trait1
varCandidates_Trait2
```

Actual suffixes come from the biological trait names when available.

### Synthetic index reporting

Report a synthetic index separately:

```r
index_report <- bp_report_stage_metrics(
  state,
  stage = "Candidates",
  metrics = c(
    meanCandidates = "meanG",
    varCandidates = "varG",
    accCandidates = "accEBV"
  ),
  synthetic_trait = "Index",
  cfg = cfg,
  append_trait = "always"
)
```

Expected names:

```r
meanCandidates_Index
varCandidates_Index
accCandidates_Index
```

For nonlinear synthetic traits, reporting first uses a cached synthetic GV.
If none exists, it calculates one. When multiple metrics are requested in one
call, the GV is calculated only once and reused.

### Naming controls

`append_trait` accepts:

- `"auto"`: append when a metric returns multiple traits;
- `"always"`: always append the biological or synthetic trait name;
- `"never"`: never append it.

For unusual conventions, provide:

```r
name_fn <- function(base_name, metric, stage, trait) {
  paste(base_name, stage, trait, sep = "_")
}
```

Duplicate output names cause an error instead of being silently repaired.

`use = "all_available"` is intentionally unsupported because the function
represents one flat result row.

### Recommended results collection

```r
results_year <- list(
  year = year,
  crossing_cost = crossing_cost
)

results_year <- c(
  results_year,
  bp_report_stage_metrics(
    state,
    stage = "Candidates",
    metrics = c(
      meanCandidates = "meanG",
      varCandidates = "varG"
    ),
    traits = 1:2,
    append_trait = "always"
  ),
  bp_report_stage_metrics(
    state,
    stage = "Candidates",
    metrics = c(
      meanCandidates = "meanG",
      varCandidates = "varG"
    ),
    synthetic_trait = "Index",
    cfg = cfg,
    append_trait = "always"
  )
)

results_year <- data.frame(results_year, check.names = FALSE)
results <- dplyr::bind_rows(results, data.frame(results_base, results_year))
```

Do not wrap the report with `unname()` or `lapply(..., unname)`. Report
elements are already unnamed internally.

## Reporting Baselines

`bp_set_trait_baseline()` can now add synthetic-trait baselines:

```r
state <- bp_set_trait_baseline(
  state,
  pop = founder_pop,
  synthetic_traits = "Index",
  varE = cfg$varE,
  synthetic_gv_n_trials = 100,
  synthetic_gv_n_plants = 10,
  seed = 123
)
```

Synthetic baseline means and SDs are stored under the synthetic trait name.
Mean-like metrics are standardized as means; `varG` is standardized by the
squared baseline SD. Accuracy and `H2` are not rescaled.

If the active baseline does not contain the requested synthetic trait,
reporting returns the raw synthetic metric rather than producing `NA` through
an unmatched baseline name.

## Compatibility Rules

- Numeric AlphaSimR biological traits remain valid throughout BPS.
- Synthetic values do not add columns to AlphaSimR biological slots.
- AlphaSimR `selectInd(trait = 1, ...)` and `trait = 1:2` semantics are
  preserved.
- A linear `AlphaSimR::selIndex` function is valid inside
  `bp_synthetic_trait()`.
- Population subsetting preserves synthetic values stored in `@misc`.
- Crossing or creating new progeny does not automatically update inherited
  synthetic caches. Re-score or re-predict the resulting population.
- Scripts must migrate from the removed singular/columns reporting functions
  to `bp_report_stage_metrics()`.

## Validation Status

The complete `testthat` suite passes. A source tarball built with
`R CMD build` also passes `R CMD check --no-manual` with status `OK`: no
errors, warnings, or notes.

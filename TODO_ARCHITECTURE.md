# Architecture TODO (Living Plan)

## API Consistency
- [ ] Enforce naming rule:
  - `bp_*` only for functions that operate on whole `state`.
  - non-`bp_` names for pure utilities.
- [ ] Rename pure helpers:
  - `bp_pop_n_ind` -> `pop_n_ind`
  - `bp_pop_subset` -> `pop_subset`
  - `bp_merge_pops` -> `merge_pops`
  - `bp_chip_key` -> `chip_key`
  - `bp_chip_index` -> `chip_index`
- [ ] Decide time helper form:
  - keep `bp_years_to_ticks(state, years)` as state-aware, or
  - switch to pure `years_to_ticks(dt, years)`.
- [ ] Remove/inline `bp_new_cohort_id` (avoid returning `{id, state}` list); generate IDs inside `bp_add_cohort`.

## Trial Layer (Readability + Flexibility)
- [ ] Keep `bp_record_pheno` internal-only; do not expose as user-facing API.
- [ ] Keep `run_phenotype_trial(state, cfg)` as the main user-facing phenotyping call.
- [ ] Confirm and document trial config schema (required vs optional fields).
- [ ] Add explicit source cohort selection policy for stage runners (especially `run_phenotype_trial`), with config options:
  - `input_policy = "latest_one"` (recommended default),
  - `input_policy = "all_ready"`,
  - `input_policy = "latest_n"`,
  - `input_policy = "by_cycle"`,
  - `input_policy = "custom"` + `select_source_cohorts_fn`.
- [ ] Set safe defaults so trial-to-trial transitions do not accidentally consume all historical ready cohorts.
- [ ] Add tests that verify AYT/PYT-style transitions use only intended source cohort(s) under each policy.
- [ ] Add genotyping idempotency checks in `run_genotyping`:
  - skip if cohort already has available genotype data for the same chip,
  - allow re-genotyping only for a different chip (or explicit force flag),
  - prevent duplicate cost entries for same cohort/chip/time window.
- [ ] Keep explicit environment control path in `run_phenotype_trial`:
  - base env means per location,
  - year-specific deviations,
  - per-env and aggregate logging controls.
- [ ] Add docs/examples for delayed data availability using `done_tick` and `available_tick`.

## Stage Function Generalization
- [ ] Treat current `run_make_f1` and `run_ssd_to_f5` as example implementations, not universal API endpoints.
- [ ] Add general crossing primitive:
  - `run_crossing_stage(state, cfg)`
  - with hooks such as `select_parents_fn` and `cross_plan_fn`.
- [ ] Add general selfing/advance primitive:
  - `run_selfing_stage(state, cfg)`
  - configurable generation depth (e.g., F3, F5, F6) and progeny schedule.
- [ ] Keep thin wrappers/examples:
  - `run_make_f1()` calling the generic crossing stage,
  - `run_ssd_to_f5()` calling the generic selfing stage.
- [ ] Add extension examples in template/docs:
  - OCS/usefulness/diallel parent selection patterns.

## Cohort/Flow Semantics
- [ ] Keep `stream` as routing metadata for parallel flows (e.g., main/testing, poolA/poolB).
- [ ] Document recommended stream naming and usage patterns.
- [ ] Confirm convention for `consume_input` across all stage functions.

## Logging + Costs
- [ ] Keep append-only logs (`phenotype_log`, `genotype_log`, `cost_log`) as source of truth.
- [ ] Add cost presets/tiers (e.g., genotyping and PYT high relative cost) to config template.
- [ ] Add summary helper(s) for cost breakdown by event/stage/time.
- [ ] Revisit standalone `run_predict_ebv` event logging:
  - for now, infer model usage from downstream selection events (`use = "ebv"` with model context),
  - later decide whether to add explicit prediction events from non-state-mutating prediction calls.

## User-Facing Layer
- [ ] Ensure high-level scheme scripts stay readable and explicit (stage-by-stage calls).
- [ ] Avoid over-abstraction in user-facing API; keep generic internals behind simple calls.
- [ ] Add a dedicated "how to add a new stage" template in `inst/templates/`.

## Validation
- [ ] Add tests for env-control phenotyping path (`use_env_control=TRUE`).
- [ ] Add tests for per-env vs aggregate phenotype logging modes.
- [ ] Add tests for delayed `available_tick` gating behavior.
- [ ] Re-run `R CMD check` after each refactor group.

## Performance (Without Hurting Readability)
- [ ] Keep performance improvements constrained to changes that preserve readability and robustness.
- [ ] Add scalable GP training controls in `run_train_gp_model`:
  - optional max training set size / subsampling strategy,
  - optional retrain cadence (e.g., every N ticks/years),
  - optional warm-start/incremental update hook for custom models where supported.
- [ ] Add model cache metadata so expensive models (e.g., MegaLMM) are not retrained unnecessarily when training inputs have not changed.
- [ ] Add optional async/offline model-training pattern in docs (train externally, load model object at next tick).
- [ ] Reduce avoidable overhead in logging/indexing without changing user-facing flow:
  - buffered log appends,
  - active cohort index by stage/stream,
  - incremental genotype flag updates.

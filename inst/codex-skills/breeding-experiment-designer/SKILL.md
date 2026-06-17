---
name: breeding-experiment-designer
description: "Create or revise BPS Create_sim_bps.R, run_experiments.R, setup/burn-in cfg, and experiment cfg files from requested architectures, historical processes, breeding schemes, and parameter comparisons. Use for reproducible single-trait or multi-trait setups, ExperimentID blocks, scenario grids, scheme comparisons, cfg sweeps, or minimal validation runs."
---

# Breeding Experiment Designer

## Required Reference

Read `references/workflow.md` completely. It is the source of
truth for BPS 0.2.0 experiment files, template selection, cfg generation,
comment blocks, ExperimentID design, and smoke tests.

Locate package templates with:

```r
system.file("templates", "bps-0.2.0", package = "BreedingProgramSimulator")
```

If empty in a source checkout, find `inst/templates/bps-0.2.0`.

## Workflow

1. Collect scheme paths/labels, architecture, historical process, scenarios,
   replicates, burn-in, comparisons, sweeps, metadata, production scale, and
   the breeder evidence used to calibrate baselines and trial reliability.
2. Use the scheme skill first when scheme logic must be created or changed.
3. Read `TEMPLATE_INDEX.md` and copy the nearest architecture and experiment
   starters; replace rather than accumulate irrelevant template code.
4. Produce `Create_sim_bps.R`, setup cfg, `run_experiments.R`, and
   experiment cfg with instructional comment blocks.
5. Design one explicit calibration procedure, repeat it at selected burn-in
   times, and apply its final values consistently at the comparison start.
6. Run the two required `bp_check_cfg_requirements()` scans and organize cfg
   values without hidden defaults.
7. Create a base ExperimentID comparison, then only the additional blocks the
   user requests.
8. Propose and help run the minimal setup/experiment smoke test from the
   reference; validate state reuse, metadata, outputs, streams, costs, and GP.

## Deliverables

Provide all four files, calibration assumptions/schedule, scenario/experiment
map, unresolved cfg table, smoke-test commands/results, and production-run
commands.

# BPS 0.2.0 Template Index

Copy the nearest template into the user's project, rename it, and edit it in
place. Do not combine templates mechanically when the biological flow differs.

## Scheme references

- `scheme_phenotypic_single_trait.R`: annual traditional phenotypic pipeline.
- `scheme_two_part_gp_single_trait.R`: line development plus rapid recurrent
  genomic parent improvement.
- `scheme_wfrgs_single_trait.R`: whole-family recurrent genomic selection.
- `scheme_phenotypic_multi_trait.R`: traditional multi-trait selection using a
  registered synthetic trait.

These are full reference implementations because trimming event logic tends to
remove the provenance, cost, timing, or reporting details the agent needs.
Delete irrelevant events only after mapping the user's flow.

## Orchestration starters

- `Create_sim_bps_single_trait_TEMPLATE.R`
- `Create_sim_bps_multi_trait_TEMPLATE.R`
- `run_experiments_TEMPLATE.R`

These are intentionally shorter and contain user-editing comment blocks.

## Selection rule

1. Read the user's diagram and protocol.
2. Choose the closest scheme and architecture template.
3. Copy rather than rewrite when stage flow and scheduler are substantially
   similar.
4. Rename stages/events to match the user's terminology.
5. Move all variable quantities into cfg.
6. Run `bp_check_cfg_requirements()`.
7. Build and run a minimal valid smoke test.
8. Draft from scratch when no template matches the flow; do not force a design
   into an ill-fitting template.

Locate installed templates with:

```r
system.file("templates", "bps-0.2.0", package = "BreedingProgramSimulator")
```


---
name: breeding-scheme-drafter
description: "Draft or refactor readable BPS breeding-scheme R scripts from diagrams, protocols, and user descriptions; identify missing specifications; validate timing, cohort flow, streams, reporting, costs, and genomic prediction; and help run minimal smoke tests before packaging. Use for new or revised single-trait, multi-trait, genomic-selection, or multi-stream breeding schemes."
---

# Breeding Scheme Drafter

## Required References

Read `references/script-style-guide.md` completely. It is the source of truth
for BPS 0.2.0 scheme structure, cfg rules, templates, event verbs, reporting,
schedulers, streams, GP debugging, and smoke tests.

Read `references/checklist.md` when validating a draft or completed scheme.

Locate package templates with:

```r
system.file("templates", "bps-0.2.0", package = "BreedingProgramSimulator")
```

If that returns an empty path in a source checkout, find
`inst/templates/bps-0.2.0` in the BPS repository.

## Workflow

1. Read the diagram, protocol, existing scripts, and cfg/experiment files.
2. Restate stages, streams, transitions, timing, selection, recycling,
   training, crossing, and reporting.
3. Provide filled-from-input and missing-specification tables plus a short
   design-question list. Never silently invent cfg values.
4. Read `TEMPLATE_INDEX.md`; copy the closest scheme only when its biological
   flow and scheduler genuinely match. Otherwise draft from scratch.
5. Keep the source-loaded scheme simple: minimal setup/helpers, one reporting
   function, visible event verbs, and one explicit scheduler.
6. Run `bp_check_cfg_requirements()` with the calling orchestration file and
   organize the cfg into readable blocks.
7. Propose and help run the style guide's minimal valid smoke test. Iterate
   after each substantial correction.
8. Validate runtime, flow, timing, provenance, streams, selection/recycling,
   costs/trials, GP, reporting, and diagram/network agreement when requested.

## Agreement Gate

Before grid packaging or final status, confirm flow, timing, selection,
recycling, streams, open questions, and a passing smoke test. Only then add
experiment integration.

## Deliverables

Provide the scheme, cfg requirements, clarification artifacts, smoke-test cfg
and command, validation summary/timeline, agreement status, and any
post-agreement experiment packaging.

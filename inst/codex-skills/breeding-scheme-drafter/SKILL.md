---
name: breeding-scheme-drafter
description: "Draft or refactor readable BPS breeding-scheme R scripts from diagrams, protocols, and user descriptions; identify missing specifications; validate timing, cohort flow, streams, reporting, costs, and genomic prediction; and help run minimal smoke tests before packaging. For completely new schemes, use staged clarification, skeleton, and implementation-plan review gates before coding. Use for new or revised single-trait, multi-trait, genomic-selection, or multi-stream breeding schemes."
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

## Completely New Schemes: Required Design Gates

Move slowly and keep the user in the design loop. Do not begin the full
implementation merely because a template looks close.

### 1. Understand the Intended Scheme

Have a focused conversation before writing implementation code. Ask short
rounds of clarification questions about unfamiliar terms and unresolved
biology. Establish:

- the meaning and purpose of every stage, stream, and named operation
- material entering and leaving each event, including temporary data-only
  populations
- timing, overlap, selection, crossing, recycling, and stopping rules
- phenotype, genotype, GP training, auxiliary-data, cost, and reporting needs
- what is fixed scheme logic versus experiment-controlled cfg

Restate the scheme in the user's terminology and maintain a visible list of
assumptions and unresolved points. Look deliberately for unnecessary stages,
stored cohorts, helper layers, duplicated operations, and stream complexity.
Suggest simpler alternatives, but do not simplify away biological intent.

Proceed only when the user confirms the restatement and no material ambiguity
remains about flow, timing, selection, or data use.

### 2. Agree on Event Verbs and Scheduler

Create a skeleton scheme script before implementing events. Include the normal
file sections, proposed event-verb names and signatures, and the complete
`run_simulation()` scheduler. For each event stub, write a two- or
three-sentence comment describing its inputs, biological action, retained or
discarded material, and output. Leave the event unimplemented and make that
status unmistakable.

Write the scheduler in full enough to show event order, loops, conditions,
streams, time advances, and reporting points. Ask the user to review the
skeleton line by line and confirm or correct the verbs and schedule. Do not
advance to implementation planning until the user approves this flow contract.

### 3. Agree on Each Event's Implementation Plan

Add an implementation plan to every event stub. Explain:

- the principal AlphaSimR and BPS calls
- input-cohort selection and output provenance
- selection criteria and temporary versus stored populations
- phenotypes, genotypes, models, and `pop@misc` data created or updated
- cfg parameters, durations, costs, logging, and reporting consequences
- any additional calculations or bookkeeping not covered by BPS

Call out uncertainty and opportunities to remove unnecessary custom code.
Ask the user to review and approve these plans. Do not implement the event
bodies until this second approval.

### 4. Implement and Validate

Implement the approved skeleton incrementally. If implementation reveals a
new biological choice or changes the agreed flow, pause and return to the
appropriate review gate instead of silently deciding.

## Implementation Workflow

For a new scheme, begin this workflow only after both design gates above are
approved. For a focused refactor of an understood scheme, summarize any flow
impact and use this workflow directly unless the requested change introduces
new biological ambiguity.

1. Read the diagram, protocol, existing scripts, and cfg/experiment files.
2. Read `TEMPLATE_INDEX.md`; copy the closest scheme only when its biological
   flow and scheduler genuinely match. Otherwise draft from scratch.
3. Keep the source-loaded scheme simple: minimal setup/helpers, one reporting
   function, visible event verbs, and one explicit scheduler.
4. Run `bp_check_cfg_requirements()` with the calling orchestration file and
   organize the cfg into readable blocks.
5. Propose and help run the style guide's minimal valid smoke test. Iterate
   after each substantial correction.
6. Validate runtime, flow, timing, provenance, streams, selection/recycling,
   costs/trials, GP, reporting, and diagram/network agreement when requested.

## Agreement Gate

For a new scheme, record separate user approval of the flow skeleton and event
implementation plans. Before grid packaging or final status, also confirm
flow, timing, selection, recycling, streams, open questions, and a passing
smoke test. Only then add experiment integration.

## Deliverables

For a new scheme, provide the terminology/assumption record, reviewed skeleton,
reviewed event plans, implemented scheme, cfg requirements, smoke-test cfg and
command, validation summary/timeline, agreement status, and any post-agreement
experiment packaging.

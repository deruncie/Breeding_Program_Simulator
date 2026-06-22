---
name: breeding-scheme-drafter
description: "Design, draft, or refactor readable BPS breeding-scheme R scripts from diagrams, protocols, and breeder descriptions. For genuinely new schemes, default to a patient stage-by-stage conversation that uncovers terminology, tacit breeding decisions, trait needs, and stored versus temporary populations before any code; then require approval of the interpretation, skeleton, and event plans. Support a design-only mode that requires no R installation and defers implementation and validation. Use for new or extended single-trait, multi-trait, genomic-selection, or multi-stream breeding schemes and their minimal smoke tests."
---

# Breeding Scheme Drafter

## Required References

Read `references/script-style-guide.md` completely. It is the source of truth
for BPS 0.2.0 scheme structure, cfg rules, templates, event verbs, reporting,
schedulers, streams, GP debugging, and smoke tests.

For a genuinely new scheme, read `references/discovery-guide.md` before asking
design questions.

Read `references/checklist.md` when validating a draft or completed scheme.

Locate package templates with:

```r
system.file("templates", "bps-0.2.0", package = "BreedingProgramSimulator")
```

If that returns an empty path in a source checkout, find
`inst/templates/bps-0.2.0` in the BPS repository.

When R is unavailable in design-only mode, look for the source-checkout path
directly. If neither installed nor source templates are available, continue
through discovery and the skeleton without them and disclose that limitation
before giving implementation plans.

## Classify the Request

Treat a scheme as an **extension** only when the user supplies a working scheme
and requests bounded changes that leave most stages and the scheduler intact.
A packaged template, a familiar crop, or a diagram resembling an old scheme
does not make the request an extension. When uncertain, use the new-scheme
workflow.

## Design-Only Mode

Use design-only mode when the user requests drafting, planning, workshop use,
or explicitly defers coding or validation. An R installation is not required
and its absence is not a blocker.

In this mode:

- do not install R or packages, run R commands, source scripts, scan cfg files,
  implement event bodies, or perform smoke tests or runtime validation
- conduct the full stage-by-stage discovery and obtain Gate 1 approval
- continue through only the gates the user requests; by default, prepare and
  review the Gate 2 skeleton and Gate 3 implementation plans, then stop
- make the Gate 2 artifact pseudocode or an unmistakably non-executable R
  skeleton according to the user's preference
- use the bundled style guide and discovery guide for Gates 1 and 2
- for Gate 3, inspect only relevant package templates, documentation, or source
  when they are available as ordinary files; do not ingest the whole package
  indiscriminately
- distinguish verified BPS operations from proposed operations that must be
  checked during a later implementation session

Finish with a handoff containing the approved stage map, persistent and
temporary populations, trait map, event/scheduler skeleton, event plans,
unresolved questions, and the deferred implementation and validation work.

## Genuinely New Schemes: Conversation First

Do not write R code, a skeleton, or detailed implementation plans at the start.
Move through the scheme stage by stage in short conversational rounds, using
the discovery guide. Adapt terminology and question depth to the user's
experience; explain BPS concepts rather than expecting a new user to know them.

For every stage, understand its purpose, inputs, activities, decisions, timing,
data, and outputs. Ask what must persist in `state` and what can remain a
temporary population inside one event. Actively uncover operations breeders may
take for granted, especially parent choice and crossing design, seed increase,
trial entry decisions, advancement rules, recycling, and the timing of data
availability.

When many traits are requested, determine how each trait changes selection,
timing, or reporting. Suggest a smaller set of representative biological traits,
stage-specific measurements, or synthetic indices when that preserves the
scientific intent. Never simplify traits or biology without user approval.

Look for unnecessary stored cohorts, stages, streams, helper layers, and
duplicated operations. Suggest simpler BPS representations while preserving the
breeder's intent.

### Gate 1: Confirm the Interpretation

Provide a stage-by-stage map, terminology notes, stored-versus-temporary
population decisions, trait/measurement map, and unresolved assumptions. Keep
discussing until the user confirms the interpretation and no material ambiguity
remains about biological flow, timing, selection, crossing, or data use.

After Gate 1, inspect `TEMPLATE_INDEX.md` and the closest relevant templates.
Use them to inform the outline, but do not let them override the agreed design.

### Gate 2: Confirm Event Verbs and Scheduler

Create a skeleton scheme script with the normal sections, proposed event-verb
signatures, and complete `run_simulation()` scheduler. Each event stub should
contain a short comment describing its available inputs, decisions, temporary
populations, stored output, and duration; leave it unmistakably unimplemented.
Show event order, loops, conditions, streams, time advances, and reporting
points. Ask the user to review and approve this flow contract.

### Gate 3: Confirm Event Implementation Plans

For every event, explain the principal AlphaSimR and BPS calls, input selection,
provenance, temporary and stored populations, cfg fields, phenotypes,
genotypes, models, costs, logging, reporting, and any custom calculations. Call
out uncertainty and opportunities to remove custom code. Do not implement event
bodies until the user approves these plans.

## Extensions of Existing Schemes

For a bounded extension, inspect the supplied scheme and describe the proposed
changes to stages, event verbs, scheduler, cfg, and reporting. Ask only the
questions needed to resolve changed or ambiguous behavior. After the user
confirms this outline, implement directly unless the extension introduces a new
biological decision; if it does, use the relevant new-scheme gate.

## Implement and Validate

Implement the approved design incrementally. Keep minimal helpers, one reporting
function, visible event verbs, and one explicit scheduler. Run
`bp_check_cfg_requirements()`, organize cfg fields, propose and run the minimal
valid smoke test, and validate timing, provenance, populations, streams,
selection, crossing, costs, GP, reporting, and agreement with any supplied
diagram or protocol.

If coding exposes a new biological choice or changes the approved flow, pause
and return to the relevant gate. Replace long planning notes in the final scheme
with only the concise comments needed to understand the implemented code.

Before experiment integration, record the three approvals for a new scheme and
confirm that the smoke test passes with no unresolved design questions.

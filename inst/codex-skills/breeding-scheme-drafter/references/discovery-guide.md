# New-Scheme Discovery Guide

Use this guide as a conversation route, not as a questionnaire. Ask a few
questions at a time, summarize the answer, and follow the breeder's terminology.
Finish one stage or uncertainty before moving to the next.

## Start with the Whole Pipeline

Ask for the program's objective, starting material, released output, cycle
length, overlapping cohorts, streams, and any diagram or protocol. Restate the
overall flow before walking through individual stages.

## Walk Through Every Stage

For each stage, establish:

- what the stage name means to this breeder
- which population or information enters and how it becomes eligible
- what physical and analytical activities occur
- who or what makes selection, mating, or advancement decisions
- population sizes before and after each decision
- what data are collected, when they become available, and how they are used
- elapsed time, season, locations, replications, and important costs
- which output is stored, which intermediates are temporary, and what consumes
  the output next

Do not assume that a familiar label implies a familiar protocol. Ask what terms
such as nursery, family, line, cycle, trial, tester, candidate, parent, or
recycling mean in this program.

## Surface Tacit Breeding Decisions

Pay special attention to activities a breeder may omit because they feel
obvious:

- how parents are eligible, chosen, paired, and replaced
- whether crossing is random, planned, assortative, factorial, reciprocal, or
  constrained by family contribution
- numbers of crosses, parents per cross, progeny, families, and retained lines
- whether seed increase, generation advance, or bulking changes timing or
  population identity
- whether selection is among individuals, families, lines, crosses, or an index
- whether failed, missing, or unbalanced entries alter decisions
- where material branches, merges, recycles, or leaves the program
- which historical populations train GP models and when models are updated

If a missing decision would require inventing code behavior, keep asking rather
than choosing a plausible default.

## Decide What Persists in State

Identify meaningful decision and data-availability boundaries first. Store a
population when it is independently available, reused, reported, used for
training, or consumed by more than one later event. Keep selected subsets,
discarded progeny, and other one-use intermediates temporary inside an event.

For each proposed stored cohort, ask what future operation needs it. For each
temporary population, confirm that no later event, report, or training process
needs the population itself or unsummarized data from it.

## Rationalize Traits and Measurements

For every requested trait, ask:

- whether it represents distinct genetic biology or another measurement of an
  existing trait
- where and when it is measurable
- whether it affects selection, GP, costs, or only reporting
- whether correlations among traits matter to the experiment
- whether several observed traits can be represented by a smaller set of
  biological traits plus stage-specific measurement error or a synthetic index

Suggest abstraction when it reduces simulation complexity without changing the
decisions being studied. Preserve separate traits when their genetic
architecture, correlations, measurement timing, or selection consequences are
scientifically important. Obtain approval before simplifying.

## Finish Discovery

Return a concise stage map, decision map, trait/measurement map, persistent and
temporary population list, cfg-versus-fixed-design split, assumptions, and open
questions. Do not proceed to the skeleton until the user confirms it.

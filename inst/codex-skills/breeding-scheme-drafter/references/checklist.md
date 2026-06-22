# Validation Checklist

Use this checklist during validation.

## Design Agreement

- Stage terminology and breeder-specific meanings are documented
- Parent choice, crossing, selection, advancement, and recycling rules are explicit
- Stored populations coincide with meaningful reuse, decision, or data boundaries
- One-use intermediate populations remain temporary
- Trait abstractions preserve measurement timing and scientific intent
- Event verbs and scheduler match the approved diagram or protocol

## Runtime

- Script runs without errors
- Stage availability timing matches described durations
- Yearly/tick summaries are internally consistent

## Protocol Alignment

- Crossing counts and parent-block counts match spec
- Advancement counts by stage match spec
- Data-only temporary populations are summarized and discarded
- Auxiliary selection data stay aligned in named `pop@misc` fields
- Recycling logic matches spec
- Variety release rule matches spec

## Logs

- Verbose logs show expected cadence
- No impossible transitions appear
- Empty/idle ticks are expected where schedule calls for them

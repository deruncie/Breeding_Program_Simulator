# Validation Checklist

Use this checklist during Step 5.

## Runtime

- Script runs without errors
- Stage availability timing matches described durations
- Yearly/tick summaries are internally consistent

## Protocol Alignment

- Crossing counts and parent-block counts match spec
- Advancement counts by stage match spec
- Recycling logic matches spec
- Variety release rule matches spec

## Logs

- Verbose logs show expected cadence
- No impossible transitions appear
- Empty/idle ticks are expected where schedule calls for them

## Network (if enabled)

- Extracted topology matches input diagram
- Method A vs Method B mismatches explained
- Remaining mismatches are documented as intentional or to-fix

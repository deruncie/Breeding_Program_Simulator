# Build a static stage-flow network from a readable scheme script.

source("R/operators.R")
source("R/readable_api.R")
source("R/scheme_network.R")

script_path <- "examples/AdvanceYear_GSTP_readable_structured.R"
orchestrator <- "run_one_year"

net <- bp_extract_scheme_network(script_path, orchestrator = orchestrator)

cat("Handlers in orchestrator order:\n")
print(net$handlers)

cat("\nExtracted edges:\n")
print(net$edges[, c("handler", "from_stage", "to_stage", "event", "details")])

dot_path <- "examples/AdvanceYear_GSTP_readable.dot"
writeLines(net$dot, dot_path)
cat("\nWrote:", dot_path, "\n")

cat("\nGraphviz rendering (optional):\n")
cat("  dot -Tpng examples/AdvanceYear_GSTP_readable.dot -o examples/AdvanceYear_GSTP_readable.png\n")

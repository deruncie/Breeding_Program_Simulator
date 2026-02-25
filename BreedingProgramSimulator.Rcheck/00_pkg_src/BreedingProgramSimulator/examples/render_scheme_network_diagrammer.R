# Render a readable scheme network with DiagrammeR.

source("R/operators.R")
source("R/readable_api.R")
source("R/scheme_network.R")

if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
  stop("Please install DiagrammeR", call. = FALSE)
}
if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  stop("Please install htmlwidgets", call. = FALSE)
}

script_path <- "examples/AdvanceYear_GSTP_readable_structured.R"
orchestrator <- "run_one_year"

net <- bp_extract_scheme_network(script_path, orchestrator = orchestrator)

g <- DiagrammeR::grViz(net$dot)

out_html <- "examples/AdvanceYear_GSTP_readable_diagrammer.html"
htmlwidgets::saveWidget(g, file = out_html, selfcontained = FALSE)

cat("Rendered DiagrammeR graph to:\n")
cat("  ", out_html, "\n", sep = "")

cat("\nNodes and edges inferred:\n")
print(net$edges[, c("handler", "from_stage", "to_stage", "event", "details")])

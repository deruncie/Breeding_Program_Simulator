# Compare two independent scheme-network inference methods.
# Method A: careful code-reading interpretation (explicit table below).
# Method B: automated parser extraction from R/scheme_network.R.

source("R/operators.R")
source("R/readable_api.R")
source("R/scheme_network.R")

script_path <- "examples/AdvanceYear_GSTP_readable_structured.R"
orchestrator <- "run_one_year"

extract_selection <- function(details) {
  m <- regexec("selection=([a-zA-Z_]+)", details)
  r <- regmatches(details, m)[[1L]]
  if (length(r) >= 2L) return(r[2L])
  NA_character_
}

# Method A: explicit human-like reading of scheme code.
method_a <- data.frame(
  from = c(
    "PYT",
    "RC",
    "DH_PIPE",
    "PYT",
    "AYT",
    "EYT",
    "RC",
    "MODEL:gp_main",
    "MODEL:gp_main",
    "MODEL:gp_main",
    "PYT",
    "GENO:PYT"
  ),
  to = c(
    "RC",
    "DH_PIPE",
    "PYT",
    "AYT",
    "EYT",
    "Variety",
    "RC",
    "DH_PIPE",
    "AYT",
    "RC",
    "MODEL:gp_main",
    "PYT"
  ),
  edge_type = c(
    "cohort",
    "cohort",
    "cohort",
    "cohort",
    "cohort",
    "cohort",
    "cohort",
    "model_dependency",
    "model_dependency",
    "model_dependency",
    "model_training",
    "genotype_info"
  ),
  selection = c(
    "random",
    "ebv",
    "random",
    "ebv",
    "pheno",
    "pheno",
    "ebv",
    NA,
    NA,
    NA,
    NA,
    NA
  ),
  stringsAsFactors = FALSE
)

# Method B: automated extractor output normalized to the same schema.
b <- bp_extract_scheme_network(script_path, orchestrator = orchestrator)$edges

method_b <- data.frame(
  from = character(),
  to = character(),
  edge_type = character(),
  selection = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(b))) {
  r <- b[i, , drop = FALSE]
  ev <- as.character(r$event)
  from <- as.character(r$from_stage)
  to <- as.character(r$to_stage)
  sel <- extract_selection(as.character(r$details))

  if (ev %in% c("stage_transition", "phenotype_trial")) {
    method_b <- rbind(method_b, data.frame(
      from = from, to = to, edge_type = "cohort", selection = sel, stringsAsFactors = FALSE
    ))
  } else if (ev == "model_dependency") {
    method_b <- rbind(method_b, data.frame(
      from = from, to = to, edge_type = "model_dependency", selection = NA_character_, stringsAsFactors = FALSE
    ))
  } else if (ev == "gp_train") {
    method_b <- rbind(method_b, data.frame(
      from = from, to = to, edge_type = "model_training", selection = NA_character_, stringsAsFactors = FALSE
    ))
  } else if (ev == "genotyping") {
    method_b <- rbind(method_b, data.frame(
      from = paste0("GENO:", from), to = from, edge_type = "genotype_info", selection = NA_character_, stringsAsFactors = FALSE
    ))
  }
}

method_a <- method_a[order(method_a$edge_type, method_a$from, method_a$to), , drop = FALSE]
method_b <- method_b[order(method_b$edge_type, method_b$from, method_b$to), , drop = FALSE]

key <- c("from", "to", "edge_type", "selection")
a_keys <- apply(method_a[, key], 1, paste, collapse = " | ")
b_keys <- apply(method_b[, key], 1, paste, collapse = " | ")

only_a <- method_a[!(a_keys %in% b_keys), , drop = FALSE]
only_b <- method_b[!(b_keys %in% a_keys), , drop = FALSE]

dir.create("examples/network_comparison", showWarnings = FALSE, recursive = TRUE)
write.csv(method_a, "examples/network_comparison/method_a_edges.csv", row.names = FALSE)
write.csv(method_b, "examples/network_comparison/method_b_edges.csv", row.names = FALSE)
write.csv(only_a, "examples/network_comparison/only_in_method_a.csv", row.names = FALSE)
write.csv(only_b, "examples/network_comparison/only_in_method_b.csv", row.names = FALSE)

cat("Method A edges:", nrow(method_a), "\n")
cat("Method B edges:", nrow(method_b), "\n")
cat("Only in A:", nrow(only_a), "\n")
cat("Only in B:", nrow(only_b), "\n")

if (nrow(only_a) > 0L) {
  cat("\nEdges only in Method A:\n")
  print(only_a)
}
if (nrow(only_b) > 0L) {
  cat("\nEdges only in Method B:\n")
  print(only_b)
}

cat("\nWrote comparison files in examples/network_comparison/\n")

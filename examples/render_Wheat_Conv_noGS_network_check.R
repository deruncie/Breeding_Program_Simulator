# Render + validate network extraction for Wheat_Conv_noGS.
#
# This script follows the repository guidance in:
# docs/CODEX_PROMPT_scheme_to_network_ggraph.txt
#
# It performs:
# 1) Method A (explicit manual/code-reading edge table),
# 2) Method B (parser extraction with bp_extract_scheme_network),
# 3) A-vs-B comparison outputs,
# 4) ggraph rendering (PNG + PDF).

source("R/operators.R")
source("R/readable_api.R")
source("R/scheme_network.R")

if (!requireNamespace("tidygraph", quietly = TRUE)) stop("Please install tidygraph", call. = FALSE)
if (!requireNamespace("ggraph", quietly = TRUE)) stop("Please install ggraph", call. = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2", call. = FALSE)

script_path <- "examples/Wheat_Conv_noGS_readable_structured.R"
orchestrator <- "run_one_year"
scheme_id <- "Wheat_Conv_noGS"

net <- bp_extract_scheme_network(script_path, orchestrator = orchestrator)
raw_edges <- net$edges

# Try to resolve cfg$... labels to default numeric values from the target script.
cfg_value_map <- list()
stage_n_map <- list()
env_cfg <- new.env(parent = baseenv())
Sys.setenv(BPS_SKIP_SCRIPT_ENTRY = "1")
try(suppressWarnings(sys.source(script_path, envir = env_cfg)), silent = TRUE)
if (exists("make_wheat_conv_nogs_cfg", envir = env_cfg, mode = "function")) {
  cfg0 <- try(get("make_wheat_conv_nogs_cfg", envir = env_cfg)(), silent = TRUE)
  if (!inherits(cfg0, "try-error") && is.list(cfg0)) {
    for (nm in names(cfg0)) {
      v <- cfg0[[nm]]
      if (length(v) == 1L && (is.numeric(v) || is.character(v))) {
        cfg_value_map[[paste0("cfg$", nm)]] <- as.character(v)
      }
    }
  }
}
if (exists("run_wheat_conv_nogs", envir = env_cfg, mode = "function") &&
    exists("make_wheat_conv_nogs_cfg", envir = env_cfg, mode = "function")) {
  cfg0 <- try(get("make_wheat_conv_nogs_cfg", envir = env_cfg)(), silent = TRUE)
  if (!inherits(cfg0, "try-error")) {
    st <- NULL
    sim_ok <- try(
      suppressMessages(
        capture.output({
          st <- get("run_wheat_conv_nogs", envir = env_cfg)(n_fill_years = 10L, n_run_years = 0L, cfg = cfg0)
        })
      ),
      silent = TRUE
    )
    if (!inherits(sim_ok, "try-error")) {
      if (is.list(st) && !is.null(st$cohorts) && nrow(st$cohorts) > 0L) {
        co <- st$cohorts[order(st$cohorts$created_tick, st$cohorts$cohort_id), , drop = FALSE]
        latest <- co[!duplicated(co$stage, fromLast = TRUE), c("stage", "n_ind"), drop = FALSE]
        for (i in seq_len(nrow(latest))) {
          stage_n_map[[as.character(latest$stage[i])]] <- as.character(latest$n_ind[i])
        }
      }
    }
  }
}

if (nrow(raw_edges) == 0L) stop("No edges extracted from scheme.", call. = FALSE)

cat("Handlers in orchestrator order:\n")
print(net$handlers)
cat("\nMethod B extracted edges:\n")
print(raw_edges[, c("handler", "from_stage", "to_stage", "event", "details")])

# Method A: explicit reproducible interpretation from readable scheme.
method_a <- data.frame(
  from = c(
    "CROSS_BLOCK",
    "DH_BULK",
    "HEADROW_EVAL",
    "HEADROW_SEL",
    "PYT",
    "PYT",
    "AYT_SEED",
    "AYT",
    "AYT",
    "EYT1_SEED",
    "EYT1",
    "EYT2_SEED",
    "PARENT_AYT",
    "EYT2"
  ),
  to = c(
    "DH_BULK",
    "HEADROW_EVAL",
    "HEADROW_SEL",
    "PYT",
    "AYT_SEED",
    "PARENT_PYT",
    "AYT",
    "EYT1_SEED",
    "PARENT_AYT",
    "EYT1",
    "EYT2_SEED",
    "EYT2",
    "CROSS_BLOCK",
    "Variety"
  ),
  edge_type = c(
    rep("cohort", 14L)
  ),
  selection = c(
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_
  ),
  stringsAsFactors = FALSE
)

extract_selection <- function(details) {
  m <- regexec("selection=([a-zA-Z_]+)", details)
  r <- regmatches(details, m)[[1L]]
  if (length(r) >= 2L) return(r[2L])
  NA_character_
}

# Normalize Method B into common schema.
method_b <- data.frame(
  from = character(),
  to = character(),
  edge_type = character(),
  selection = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(raw_edges))) {
  r <- raw_edges[i, , drop = FALSE]
  ev <- as.character(r$event)
  from <- as.character(r$from_stage)
  to <- as.character(r$to_stage)
  sel <- extract_selection(as.character(r$details))

  if (ev %in% c("stage_transition", "phenotype_trial")) {
    method_b <- rbind(method_b, data.frame(
      from = from, to = to, edge_type = "cohort", selection = sel, stringsAsFactors = FALSE
    ))
  }
}

method_a <- method_a[order(method_a$edge_type, method_a$from, method_a$to), , drop = FALSE]
method_b <- method_b[order(method_b$edge_type, method_b$from, method_b$to), , drop = FALSE]

key <- c("from", "to", "edge_type")
a_keys <- apply(method_a[, key], 1, paste, collapse = " | ")
b_keys <- apply(method_b[, key], 1, paste, collapse = " | ")

only_a <- method_a[!(a_keys %in% b_keys), , drop = FALSE]
only_b <- method_b[!(b_keys %in% a_keys), , drop = FALSE]

cmp_dir <- "examples/network_comparison"
dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(method_a, file.path(cmp_dir, paste0(scheme_id, "_method_a_edges.csv")), row.names = FALSE)
write.csv(method_b, file.path(cmp_dir, paste0(scheme_id, "_method_b_edges.csv")), row.names = FALSE)
write.csv(only_a, file.path(cmp_dir, paste0(scheme_id, "_only_in_method_a.csv")), row.names = FALSE)
write.csv(only_b, file.path(cmp_dir, paste0(scheme_id, "_only_in_method_b.csv")), row.names = FALSE)

cat("\nMethod A edges:", nrow(method_a), "\n")
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

# Render ggraph from Method B extractor output.
vis_edges <- data.frame(
  from_stage = as.character(raw_edges$from_stage),
  to_stage = as.character(raw_edges$to_stage),
  edge_class = ifelse(raw_edges$event %in% c("stage_transition", "phenotype_trial"), "cohort", "info"),
  stringsAsFactors = FALSE
)
vis_edges <- vis_edges[vis_edges$edge_class == "cohort", , drop = FALSE]
vis_edges <- vis_edges[!duplicated(vis_edges), , drop = FALSE]

nodes <- sort(unique(c(vis_edges$from_stage, vis_edges$to_stage)))
nodes_df <- data.frame(
  name = nodes,
  node_class = "cohort",
  stringsAsFactors = FALSE
)

trial_in <- raw_edges[raw_edges$event == "phenotype_trial", c("to_stage", "details"), drop = FALSE]
extract_field <- function(details, key) {
  m <- regexec(paste0("(^|;\\s*)", key, "=([^;]+)"), details)
  r <- regmatches(details, m)[[1L]]
  if (length(r) >= 3L) return(trimws(r[3L]))
  NA_character_
}

normalize_n_value <- function(n_txt) {
  if (is.na(n_txt) || !nzchar(n_txt)) return(n_txt)
  s <- gsub("\\s+", "", n_txt)
  m <- regexec("^min\\(as\\.integer\\(([0-9]+)\\),", s)
  r <- regmatches(s, m)[[1L]]
  if (length(r) >= 2L) return(r[2L])
  m <- regexec("^min\\(as\\.integer\\((cfg\\$[A-Za-z0-9_]+)\\),", s)
  r <- regmatches(s, m)[[1L]]
  if (length(r) >= 2L) {
    k <- r[2L]
    if (!is.null(cfg_value_map[[k]])) return(cfg_value_map[[k]])
    return(k)
  }
  m <- regexec("^min\\(([0-9]+),", s)
  r <- regmatches(s, m)[[1L]]
  if (length(r) >= 2L) return(r[2L])
  m <- regexec("size=([0-9]+)", s)
  r <- regmatches(s, m)[[1L]]
  if (length(r) >= 2L) return(r[2L])
  m <- regexec("size=min\\(as\\.integer\\((cfg\\$[A-Za-z0-9_]+)\\),", s)
  r <- regmatches(s, m)[[1L]]
  if (length(r) >= 2L) {
    k <- r[2L]
    if (!is.null(cfg_value_map[[k]])) return(cfg_value_map[[k]])
    return(k)
  }
  m <- regexec("size=min\\(as\\.integer\\(([0-9]+)\\),", s)
  r <- regmatches(s, m)[[1L]]
  if (length(r) >= 2L) return(r[2L])
  n_txt
}

nodes_df$annotation <- ""
for (i in seq_len(nrow(nodes_df))) {
  nm <- nodes_df$name[i]
  rows <- trial_in[trial_in$to_stage == nm, , drop = FALSE]
  if (nrow(rows) == 0L) {
    if (!is.null(stage_n_map[[nm]])) {
      nodes_df$annotation[i] <- paste0("n=", stage_n_map[[nm]])
    }
    next
  }

  ann <- unique(vapply(seq_len(nrow(rows)), function(k) {
    d <- as.character(rows$details[k])
    sel <- extract_selection(d)
    n <- extract_field(d, "n")
    n <- normalize_n_value(n)
    loc <- extract_field(d, "loc")
    rep <- extract_field(d, "reps")
    parts <- c()
    if (!is.na(sel)) parts <- c(parts, paste0("sel=", sel))
    if (!is.na(n)) parts <- c(parts, paste0("n=", n))
    if (!is.na(loc)) parts <- c(parts, paste0("loc=", loc))
    if (!is.na(rep)) parts <- c(parts, paste0("rep=", rep))
    paste(parts, collapse = "\n")
  }, character(1)))
  ann <- ann[nzchar(ann)]
  if (length(ann) > 0L) nodes_df$annotation[i] <- paste(ann, collapse = "\n")
  if (!is.null(stage_n_map[[nm]]) && !grepl("(^|\\n)n=", nodes_df$annotation[i])) {
    nodes_df$annotation[i] <- paste(c(nodes_df$annotation[i], paste0("n=", stage_n_map[[nm]])), collapse = "\n")
  }
}

nodes_df$label <- ifelse(
  nzchar(nodes_df$annotation),
  paste(nodes_df$name, nodes_df$annotation, sep = "\n"),
  nodes_df$name
)

edges_df <- data.frame(
  from = match(vis_edges$from_stage, nodes_df$name),
  to = match(vis_edges$to_stage, nodes_df$name),
  edge_class = vis_edges$edge_class,
  stringsAsFactors = FALSE
)

graph <- tidygraph::tbl_graph(nodes = nodes_df, edges = edges_df, directed = TRUE)
layout_tbl <- ggraph::create_layout(graph, layout = "sugiyama")
xr <- range(layout_tbl$x, na.rm = TRUE)
yr <- range(layout_tbl$y, na.rm = TRUE)
xpad <- 0.5
ypad <- 0.5
max_lab_chars <- max(nchar(nodes_df$label), na.rm = TRUE)
plot_w <- max(8, 7 + 0.035 * max_lab_chars)
plot_h <- max(9, diff(yr) + 3)

p <- ggraph::ggraph(layout_tbl) +
  ggraph::geom_edge_fan(
    ggplot2::aes(color = edge_class),
    arrow = grid::arrow(type = "closed", length = grid::unit(2.5, "mm")),
    end_cap = ggraph::circle(5.5, "mm"),
    start_cap = ggraph::circle(3.5, "mm"),
    alpha = 0.95,
    width = 1.0
  ) +
  ggraph::geom_node_label(
    ggplot2::aes(label = label, fill = node_class),
    size = 3.8,
    color = "black",
    lineheight = 0.95
  ) +
  ggraph::scale_edge_color_manual(values = c(cohort = "#2e7d32")) +
  ggplot2::scale_fill_manual(values = c(cohort = "white")) +
  ggplot2::theme_void() +
  ggplot2::coord_cartesian(
    xlim = c(xr[1] - xpad, xr[2] + xpad),
    ylim = c(yr[1] - ypad, yr[2] + ypad),
    clip = "off"
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.margin = grid::unit(c(8, 8, 8, 8), "pt")
  ) +
  ggplot2::labs(title = paste0(scheme_id, ": Cohort Flow (no GS)"))

png_out <- file.path("examples", paste0(scheme_id, "_network_ggraph.png"))
pdf_out <- file.path("examples", paste0(scheme_id, "_network_ggraph.pdf"))
ggplot2::ggsave(filename = png_out, plot = p, width = plot_w, height = plot_h, dpi = 200, bg = "white", limitsize = FALSE)
ggplot2::ggsave(filename = pdf_out, plot = p, width = plot_w, height = plot_h, device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)

cat("\nWrote network plots:\n")
cat(" -", png_out, "\n")
cat(" -", pdf_out, "\n")

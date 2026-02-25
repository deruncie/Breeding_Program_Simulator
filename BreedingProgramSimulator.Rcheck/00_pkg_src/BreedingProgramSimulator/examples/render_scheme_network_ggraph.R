# Render a readable scheme network with tidygraph + ggraph.

source("R/operators.R")
source("R/readable_api.R")
source("R/scheme_network.R")

if (!requireNamespace("tidygraph", quietly = TRUE)) stop("Please install tidygraph", call. = FALSE)
if (!requireNamespace("ggraph", quietly = TRUE)) stop("Please install ggraph", call. = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2", call. = FALSE)

script_path <- "examples/AdvanceYear_GSTP_readable_structured.R"
orchestrator <- "run_one_year"
net <- bp_extract_scheme_network(script_path, orchestrator = orchestrator)
raw_edges <- net$edges
cfg_vals <- bp_extract_cfg_values(script_path)

if (nrow(raw_edges) == 0L) stop("No edges extracted from scheme.", call. = FALSE)

extract_selection <- function(details) {
  m <- regexec("selection=([a-zA-Z_]+)", details)
  r <- regmatches(details, m)[[1L]]
  if (length(r) >= 2L) return(r[2L])
  NA_character_
}

extract_detail <- function(details, key) {
  m <- regexec(paste0("(^|;\\s*)", key, "=([^; ]+)"), details)
  r <- regmatches(details, m)[[1L]]
  if (length(r) >= 3L) return(r[3L])
  NA_character_
}

make_vis_edges <- function(edge_df) {
  out <- list()
  add <- function(from, to, class, label = "") {
    out[[length(out) + 1L]] <<- data.frame(
      from_stage = as.character(from),
      to_stage = as.character(to),
      edge_class = as.character(class),
      label = as.character(label),
      stringsAsFactors = FALSE
    )
  }

  for (i in seq_len(nrow(edge_df))) {
    r <- edge_df[i, , drop = FALSE]
    from <- as.character(r$from_stage)
    to <- as.character(r$to_stage)
    event <- as.character(r$event)

    if (event %in% c("stage_transition", "phenotype_trial")) {
      add(from, to, "cohort", "")
    } else if (event == "model_dependency") {
      add(from, to, "model_dep", "model used")
    } else if (event == "gp_train") {
      add(from, to, "model_train", "train model")
    } else if (event == "genotyping") {
      add(paste0("GENO:", from), from, "genotype_info", "genotype")
    }
  }

  out_df <- do.call(rbind, out)
  out_df[!duplicated(out_df), , drop = FALSE]
}

# Replace cohort self-loops by duplicated cycle nodes:
# RC -> RC_cycle and RC_cycle -> RC.
expand_cohort_self_loops <- function(vis_df) {
  out <- list()
  cycle_nodes <- character()

  for (i in seq_len(nrow(vis_df))) {
    r <- vis_df[i, , drop = FALSE]
    is_loop <- r$edge_class == "cohort" && r$from_stage == r$to_stage
    if (is_loop) {
      cyc <- paste0(r$from_stage, "_cycle")
      cycle_nodes <- c(cycle_nodes, cyc)
      out[[length(out) + 1L]] <- data.frame(
        from_stage = r$from_stage, to_stage = cyc, edge_class = r$edge_class, label = r$label, stringsAsFactors = FALSE
      )
      out[[length(out) + 1L]] <- data.frame(
        from_stage = cyc, to_stage = r$to_stage, edge_class = r$edge_class, label = r$label, stringsAsFactors = FALSE
      )
    } else {
      out[[length(out) + 1L]] <- r
    }
  }

  list(
    edges = do.call(rbind, out),
    cycle_nodes = unique(cycle_nodes)
  )
}

vis_edges <- make_vis_edges(raw_edges)
expanded <- expand_cohort_self_loops(vis_edges)
vis_edges <- expanded$edges
cycle_nodes <- expanded$cycle_nodes

nodes <- sort(unique(c(vis_edges$from_stage, vis_edges$to_stage)))
nodes_df <- data.frame(name = nodes, stringsAsFactors = FALSE)
nodes_df$node_class <- "cohort"
nodes_df$node_class[grepl("^MODEL:", nodes_df$name)] <- "model"
nodes_df$node_class[grepl("^GENO:", nodes_df$name)] <- "genotype"
nodes_df$node_class[nodes_df$name %in% cycle_nodes] <- "cohort_cycle"

# Stage annotations from phenotype events.
trial_in <- raw_edges[raw_edges$event == "phenotype_trial", c("to_stage", "details"), drop = FALSE]
trial_in$ann <- ""
if (nrow(trial_in) > 0L) {
  trial_in$ann <- vapply(seq_len(nrow(trial_in)), function(i) {
    d <- as.character(trial_in$details[i])
    sel <- extract_selection(d)
    n <- extract_detail(d, "n")
    loc <- extract_detail(d, "loc")
    reps <- extract_detail(d, "reps")
    parts <- c()
    if (!is.na(sel)) parts <- c(parts, paste0("sel=", sel))
    if (!is.na(n)) parts <- c(parts, paste0("n=", n))
    if (!is.na(loc)) parts <- c(parts, paste0("loc=", loc))
    if (!is.na(reps)) parts <- c(parts, paste0("rep=", reps))
    paste(parts, collapse = "\n")
  }, character(1))
}

nodes_df$annotation <- ""
for (i in seq_len(nrow(nodes_df))) {
  nm <- nodes_df$name[i]
  if (nm %in% cycle_nodes) next
  anns <- unique(trial_in$ann[trial_in$to_stage == nm])
  anns <- anns[!is.na(anns) & nzchar(anns)]
  if (length(anns) > 0L) nodes_df$annotation[i] <- paste(anns, collapse = "\n")
}

# Add RC annotation from recurrent-cycle selection.
rc_self <- raw_edges[raw_edges$from_stage == "RC" & raw_edges$to_stage == "RC" & raw_edges$event == "stage_transition", , drop = FALSE]
if (nrow(rc_self) > 0L) {
  sel <- extract_selection(as.character(rc_self$details[1L]))
  n_rc <- cfg_vals[["rc_select_n"]]
  rc_parts <- c()
  if (!is.na(sel)) rc_parts <- c(rc_parts, paste0("sel=", sel))
  if (!is.null(n_rc) && nzchar(as.character(n_rc))) rc_parts <- c(rc_parts, paste0("n=", as.character(n_rc)))
  rc_ann <- paste(rc_parts, collapse = "\n")
  idx_rc <- match("RC", nodes_df$name)
  if (!is.na(idx_rc) && nzchar(rc_ann)) nodes_df$annotation[idx_rc] <- rc_ann
}

nodes_df$title <- nodes_df$name
nodes_df$title <- sub("^MODEL:(.+)$", "Model\n\\1", nodes_df$title)
nodes_df$title <- sub("^GENO:(.+)$", "Genotypes\n\\1", nodes_df$title)
nodes_df$title[nodes_df$name %in% cycle_nodes] <- sub("_cycle$", "\ncycle", nodes_df$name[nodes_df$name %in% cycle_nodes])

nodes_df$label <- ifelse(
  nzchar(nodes_df$annotation),
  paste(nodes_df$title, nodes_df$annotation, sep = "\n"),
  nodes_df$title
)

edges_df <- data.frame(
  from = match(vis_edges$from_stage, nodes_df$name),
  to = match(vis_edges$to_stage, nodes_df$name),
  edge_class = vis_edges$edge_class,
  label = vis_edges$label,
  stringsAsFactors = FALSE
)

graph <- tidygraph::tbl_graph(nodes = nodes_df, edges = edges_df, directed = TRUE)

manual_pos <- data.frame(
  name = c("DH_PIPE", "PYT", "AYT", "EYT", "Variety", "RC", "RC_cycle", "MODEL:gp_main", "GENO:PYT"),
  x = c(0.46, 0.46, 0.46, 0.46, 0.46, 0.00, 0.08, 0.00, 0.00),
  y = c(1.10, 0.45, -0.25, -0.95, -1.65, -1.05, -1.34, -0.15, 1.85),
  stringsAsFactors = FALSE
)

pos_idx <- match(nodes_df$name, manual_pos$name)
nodes_df$x <- ifelse(is.na(pos_idx), 0, manual_pos$x[pos_idx])
nodes_df$y <- ifelse(is.na(pos_idx), 0, manual_pos$y[pos_idx])
layout_tbl <- ggraph::create_layout(graph, layout = "manual", x = nodes_df$x, y = nodes_df$y)

edge_data_all <- function(layout) ggraph::get_edges()(layout)
edge_data_labeled <- function(layout) {
  e <- edge_data_all(layout)
  lbl <- as.character(e$label)
  e[!is.na(lbl) & nzchar(lbl) & lbl != "NA", , drop = FALSE]
}

p <- ggraph::ggraph(layout_tbl) +
  ggraph::geom_edge_fan(
    data = edge_data_all,
    ggplot2::aes(color = edge_class, linetype = edge_class),
    arrow = grid::arrow(type = "closed", length = grid::unit(2.4, "mm")),
    end_cap = ggraph::circle(7.0, "mm"),
    start_cap = ggraph::circle(5.0, "mm"),
    angle_calc = "along",
    alpha = 0.95,
    width = 1.02,
    label_alpha = 0
  ) +
  ggraph::geom_edge_fan(
    data = edge_data_labeled,
    ggplot2::aes(label = label),
    arrow = grid::arrow(type = "closed", length = grid::unit(2.4, "mm")),
    end_cap = ggraph::circle(7.0, "mm"),
    start_cap = ggraph::circle(5.0, "mm"),
    angle_calc = "along",
    alpha = 0,
    width = 0.01,
    label_colour = "black",
    label_alpha = 1,
    label_size = 4.0,
    check_overlap = FALSE,
    show.legend = FALSE
  ) +
  ggraph::geom_node_label(
    ggplot2::aes(label = label, fill = node_class),
    size = 3.6,
    color = "black",
    lineheight = 0.92
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      cohort = "white",
      cohort_cycle = "#eef7ef",
      model = "#e8f0ff",
      genotype = "#f3f4f6"
    )
  ) +
  ggraph::scale_edge_color_manual(
    values = c(
      cohort = "#2e7d32",
      model_dep = "#2563eb",
      model_train = "#6b7280",
      genotype_info = "#6b7280"
    )
  ) +
  ggraph::scale_edge_linetype_manual(
    values = c(
      cohort = "solid",
      model_dep = "dotdash",
      model_train = "dashed",
      genotype_info = "dashed"
    )
  ) +
  ggplot2::coord_cartesian(xlim = c(-0.04, 0.52), ylim = c(-1.92, 2.02), clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.margin = grid::unit(c(10, 16, 10, 10), "pt")
  ) +
  ggplot2::labs(title = "AdvanceYear_GSTP_readable: Cohort and Info Flow")

out_png <- "examples/AdvanceYear_GSTP_readable_ggraph.png"
out_pdf <- "examples/AdvanceYear_GSTP_readable_ggraph.pdf"
ggplot2::ggsave(out_png, plot = p, width = 6.3, height = 10.0, dpi = 220, bg = "white")
ggplot2::ggsave(out_pdf, plot = p, width = 6.3, height = 10.0, bg = "white")

cat("Rendered ggraph outputs:\n")
cat("  ", out_png, "\n", sep = "")
cat("  ", out_pdf, "\n", sep = "")

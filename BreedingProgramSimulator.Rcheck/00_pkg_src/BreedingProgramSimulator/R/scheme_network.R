# Static scheme-to-network extraction helpers for readable scripts.

# Convert an expression to compact one-line text for labels/debugging.
bp_expr_to_text <- function(expr) {
  paste(trimws(gsub("\\s+", " ", deparse(expr, width.cutoff = 500L))), collapse = " ")
}

# Normalize stage/model labels for graph output.
bp_norm_label <- function(x) {
  s <- trimws(as.character(x))
  if (grepl('^".*"$', s)) {
    s <- substr(s, 2L, nchar(s) - 1L)
  }
  s <- sub("^cfg_local\\$", "cfg$", s)
  s
}

# Return function/operator name for simple calls, else "".
bp_call_name <- function(call_expr) {
  if (!is.call(call_expr)) return("")
  head <- call_expr[[1L]]
  if (is.symbol(head)) return(as.character(head))
  ""
}

# Return the named argument expression from a call, or NULL if absent.
bp_call_arg <- function(call_expr, arg_name, position = NULL) {
  if (!is.call(call_expr)) return(NULL)
  args <- as.list(call_expr)[-1]
  nms <- names(args)

  if (!is.null(nms) && arg_name %in% nms) {
    return(args[[which(nms == arg_name)[1L]]])
  }
  if (!is.null(position) && length(args) >= position) {
    return(args[[position]])
  }
  NULL
}

# Parse top-level function definitions from an R script.
bp_parse_scheme_functions <- function(script_path) {
  exprs <- parse(file = script_path, keep.source = TRUE)
  out <- list()

  for (e in as.list(exprs)) {
    if (!is.call(e)) next
    op <- bp_call_name(e)
    if (!(op %in% c("<-", "="))) next
    lhs <- e[[2L]]
    rhs <- e[[3L]]
    if (!is.symbol(lhs) || !is.call(rhs) || as.character(rhs[[1L]]) != "function") next
    out[[as.character(lhs)]] <- rhs
  }
  out
}

# Walk expression tree in lexical order and collect state <- fn(...) handler names.
bp_collect_state_call_sequence <- function(fn_expr) {
  calls <- character()

  walk <- function(x) {
    if (!is.call(x)) return(invisible(NULL))

    op <- bp_call_name(x)
    if (op %in% c("<-", "=") && length(x) >= 3L) {
      lhs <- x[[2L]]
      rhs <- x[[3L]]
      if (is.symbol(lhs) && as.character(lhs) == "state" && is.call(rhs) && is.symbol(rhs[[1L]])) {
        calls <<- c(calls, as.character(rhs[[1L]]))
      }
    }

    for (i in seq_along(as.list(x))[-1L]) {
      walk(x[[i]])
    }
    invisible(NULL)
  }

  walk(fn_expr[[3L]])
  unique(calls)
}

# Infer selection hint from function body text.
bp_selection_hint <- function(body_txt) {
  if (grepl('use\\s*=\\s*"ebv"', body_txt)) return("selection=ebv")
  if (grepl('use\\s*=\\s*"pheno"', body_txt)) return("selection=pheno")
  if (grepl("sample\\.int\\(", body_txt)) return("selection=random")
  ""
}

# Infer operation hint from function body text.
bp_operation_hint <- function(body_txt) {
  ops <- character()
  if (grepl("run_phenotype_trial\\(", body_txt) || grepl("setPheno\\(", body_txt)) ops <- c(ops, "phenotype")
  if (grepl("run_genotyping\\(", body_txt)) ops <- c(ops, "genotype")
  if (grepl("run_train_gp_model\\(", body_txt)) ops <- c(ops, "gp_train")
  if (grepl("makeCross\\(", body_txt)) ops <- c(ops, "cross")
  if (grepl("makeDH\\(", body_txt)) ops <- c(ops, "dh")
  if (grepl("selectInd\\(", body_txt)) ops <- c(ops, "select")
  paste(ops, collapse = ",")
}

# Extract stage-flow edges from one handler function.
bp_extract_edges_from_handler <- function(handler_name, fn_expr) {
  body_txt <- paste(deparse(fn_expr[[3L]], width.cutoff = 500L), collapse = "\n")
  sel_hint <- bp_selection_hint(body_txt)
  op_hint <- bp_operation_hint(body_txt)

  input_stages <- character()
  rows <- list()
  model_refs <- character()

  add_row <- function(from, to, event, details = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      handler = handler_name,
      from_stage = as.character(from),
      to_stage = as.character(to),
      event = as.character(event),
      details = as.character(details),
      stringsAsFactors = FALSE
    )
  }

  walk <- function(x) {
    if (!is.call(x)) return(invisible(NULL))
    fn <- bp_call_name(x)

    if (fn == "get_ready_pop") {
      stage_expr <- bp_call_arg(x, "stage", position = 2L)
      if (!is.null(stage_expr)) {
        input_stages <<- c(input_stages, bp_norm_label(bp_expr_to_text(stage_expr)))
      }
    }

    # Track model ids referenced via state$gs_models[[...]].
    if (fn == "[[") {
      arg1 <- bp_call_arg(x, "", position = 1L)
      arg2 <- bp_call_arg(x, "", position = 2L)
      if (is.call(arg1) && bp_call_name(arg1) == "$") {
        base_obj <- bp_call_arg(arg1, "", position = 1L)
        field <- bp_call_arg(arg1, "", position = 2L)
        if (is.symbol(base_obj) && as.character(base_obj) == "state" &&
            bp_expr_to_text(field) == "gs_models" && !is.null(arg2)) {
          model_refs <<- c(model_refs, bp_norm_label(bp_expr_to_text(arg2)))
        }
      }
    }

    if (fn == "run_phenotype_trial") {
      cfg_expr <- bp_call_arg(x, "", position = 2L)
      if (is.call(cfg_expr) && as.character(cfg_expr[[1L]]) == "list") {
        stage_in <- bp_norm_label(bp_expr_to_text(bp_call_arg(cfg_expr, "input_stage")))
        stage_out_expr <- bp_call_arg(cfg_expr, "output_stage")
        trial_name_expr <- bp_call_arg(cfg_expr, "trial_name")
        stage_out <- if (!is.null(stage_out_expr)) bp_norm_label(bp_expr_to_text(stage_out_expr)) else bp_norm_label(bp_expr_to_text(trial_name_expr))
        n_loc <- bp_call_arg(cfg_expr, "n_loc")
        reps <- bp_call_arg(cfg_expr, "reps")
        n_label <- bp_call_arg(cfg_expr, "n_pyt")
        if (is.null(n_label)) n_label <- bp_call_arg(cfg_expr, "n_ayt")
        if (is.null(n_label)) n_label <- bp_call_arg(cfg_expr, "n_eyt")
        details <- paste(
          c(
            op_hint,
            sel_hint,
            if (!is.null(n_label)) paste0("n=", bp_expr_to_text(n_label)) else "",
            if (!is.null(n_loc)) paste0("loc=", bp_expr_to_text(n_loc)) else "",
            if (!is.null(reps)) paste0("reps=", bp_expr_to_text(reps)) else ""
          ),
          collapse = "; "
        )
        details <- gsub("^;\\s*|;\\s*$", "", gsub(";\\s*;+", "; ", details))
        add_row(stage_in, stage_out, "phenotype_trial", details)
      }
    }

    if (fn == "run_genotyping") {
      cfg_expr <- bp_call_arg(x, "", position = 2L)
      if (is.call(cfg_expr) && as.character(cfg_expr[[1L]]) == "list") {
        stage_in <- bp_norm_label(bp_expr_to_text(bp_call_arg(cfg_expr, "input_stage")))
        chip <- bp_call_arg(cfg_expr, "chip")
        details <- paste(
          c(op_hint, if (!is.null(chip)) paste0("chip=", bp_expr_to_text(chip)) else ""),
          collapse = "; "
        )
        details <- gsub("^;\\s*|;\\s*$", "", gsub(";\\s*;+", "; ", details))
        add_row(stage_in, stage_in, "genotyping", details)
      }
    }

    if (fn == "run_train_gp_model") {
      cfg_expr <- bp_call_arg(x, "", position = 2L)
      if (is.call(cfg_expr) && as.character(cfg_expr[[1L]]) == "list") {
        stage_in <- bp_norm_label(bp_expr_to_text(bp_call_arg(cfg_expr, "from_stage")))
        model_id <- bp_call_arg(cfg_expr, "model_id")
        to_stage <- if (!is.null(model_id)) paste0("MODEL:", bp_norm_label(bp_expr_to_text(model_id))) else "MODEL"
        add_row(stage_in, to_stage, "gp_train", op_hint)
      }
    }

    if (fn == "put_stage_pop") {
      stage_out <- bp_call_arg(x, "stage", position = 3L)
      if (!is.null(stage_out)) {
        from <- if (length(input_stages) > 0L) input_stages[length(input_stages)] else "UNKNOWN"
        details <- paste(c(op_hint, sel_hint), collapse = "; ")
        details <- gsub("^;\\s*|;\\s*$", "", gsub(";\\s*;+", "; ", details))
        add_row(from, bp_norm_label(bp_expr_to_text(stage_out)), "stage_transition", details)
      }
    }

    for (i in seq_along(as.list(x))[-1L]) {
      walk(x[[i]])
    }
    invisible(NULL)
  }

  walk(fn_expr[[3L]])
  if (length(model_refs) > 0L) {
    model_refs <- unique(model_refs)
    # Attach model dependencies to produced stages when available.
    stage_targets <- character()
    if (length(rows) > 0L) {
      stage_targets <- unique(vapply(rows, function(r) as.character(r$to_stage), character(1)))
    } else {
      stage_targets <- unique(input_stages)
    }
    stage_targets <- stage_targets[nzchar(stage_targets) & stage_targets != "UNKNOWN"]
    for (m in model_refs) {
      for (stg in stage_targets) {
        add_row(
          from = paste0("MODEL:", m),
          to = stg,
          event = "model_dependency",
          details = paste(c("model_used", sel_hint), collapse = "; ")
        )
      }
    }
  }
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  out[!duplicated(out), , drop = FALSE]
}

# Extract literal cfg values (e.g. cfg$rc_stage = "RC") from a script.
bp_extract_cfg_values <- function(script_path) {
  exprs <- parse(file = script_path, keep.source = TRUE)
  cfg <- list()

  walk <- function(x) {
    if (!is.call(x)) return(invisible(NULL))
    op <- bp_call_name(x)
    if (op %in% c("<-", "=") && length(x) >= 3L) {
      lhs <- x[[2L]]
      rhs <- x[[3L]]
      if (is.symbol(lhs) && as.character(lhs) == "cfg" && is.call(rhs) && bp_call_name(rhs) == "list") {
        args <- as.list(rhs)[-1L]
        nms <- names(args)
        for (i in seq_along(args)) {
          nm <- nms[[i]]
          if (is.null(nm) || !nzchar(nm)) next
          val <- args[[i]]
          # Keep only simple scalar constants for substitution.
          if (is.character(val) && length(val) == 1L) cfg[[nm]] <<- val
          if (is.numeric(val) && length(val) == 1L) cfg[[nm]] <<- as.character(val)
        }
      }
    }
    for (i in seq_along(as.list(x))[-1L]) walk(x[[i]])
    invisible(NULL)
  }

  for (e in as.list(exprs)) walk(e)
  cfg
}

# Apply cfg substitutions to stage/model labels.
bp_apply_cfg_labels <- function(edge_df, cfg_values) {
  if (nrow(edge_df) == 0L || length(cfg_values) == 0L) return(edge_df)

  subst_one <- function(s) {
    out <- as.character(s)
    for (nm in names(cfg_values)) {
      pat <- paste0("cfg\\$", nm)
      out <- gsub(pat, as.character(cfg_values[[nm]]), out, fixed = FALSE)
    }
    out
  }

  edge_df$from_stage <- vapply(edge_df$from_stage, subst_one, character(1))
  edge_df$to_stage <- vapply(edge_df$to_stage, subst_one, character(1))
  if ("details" %in% names(edge_df)) {
    edge_df$details <- vapply(edge_df$details, subst_one, character(1))
  }
  if ("handler" %in% names(edge_df)) {
    edge_df$handler <- vapply(edge_df$handler, subst_one, character(1))
  }
  edge_df
}

# Convert edge table to graphviz DOT text.
bp_network_to_dot <- function(edge_df, title = "Scheme Network") {
  if (nrow(edge_df) == 0L) {
    return("digraph Scheme {\n  label=\"Scheme Network\";\n}\n")
  }

  nodes <- sort(unique(c(edge_df$from_stage, edge_df$to_stage)))
  esc <- function(x) gsub('"', '\\"', x, fixed = TRUE)
  fmt_node <- function(n) sprintf('  "%s";', esc(n))

  fmt_edge <- function(r) {
    style <- "solid"
    color <- "black"
    if (identical(r[["event"]], "genotyping") || grepl("^MODEL", r[["to_stage"]])) {
      style <- "dashed"
      color <- "gray40"
    } else if (identical(r[["event"]], "phenotype_trial")) {
      color <- "forestgreen"
    } else if (identical(r[["event"]], "model_dependency")) {
      style <- "dotted"
      color <- "royalblue4"
    }
    lbl <- paste(c(r[["handler"]], r[["event"]], r[["details"]]), collapse = " | ")
    lbl <- gsub("\\s*\\|\\s*$", "", gsub("\\|\\s*\\|", "|", lbl))
    sprintf(
      '  "%s" -> "%s" [label="%s", color="%s", style="%s"];',
      esc(r[["from_stage"]]), esc(r[["to_stage"]]), esc(lbl), color, style
    )
  }

  lines <- c(
    "digraph Scheme {",
    sprintf('  label="%s";', esc(title)),
    '  labelloc="t";',
    '  rankdir="LR";',
    '  node [shape=box, style="rounded"];',
    vapply(nodes, fmt_node, character(1)),
    apply(edge_df, 1, fmt_edge),
    "}"
  )
  paste(lines, collapse = "\n")
}

# End-to-end extraction from a scheme script.
bp_extract_scheme_network <- function(script_path, orchestrator = "run_one_year") {
  fns <- bp_parse_scheme_functions(script_path)
  if (!orchestrator %in% names(fns)) {
    stop(sprintf("Could not find orchestrator function '%s' in %s", orchestrator, script_path), call. = FALSE)
  }

  handlers <- bp_collect_state_call_sequence(fns[[orchestrator]])
  handlers <- handlers[handlers %in% names(fns)]

  edge_rows <- lapply(handlers, function(h) {
    e <- bp_extract_edges_from_handler(h, fns[[h]])
    if (nrow(e) == 0L) return(e)
    e$order_index <- match(h, handlers)
    e
  })
  edge_df <- do.call(rbind, edge_rows)
  if (!is.null(edge_df) && nrow(edge_df) > 0L) {
    edge_df <- edge_df[order(edge_df$order_index, edge_df$handler, edge_df$from_stage, edge_df$to_stage), , drop = FALSE]
    edge_df <- bp_apply_cfg_labels(edge_df, bp_extract_cfg_values(script_path))
    edge_df <- edge_df[!duplicated(edge_df[, c("handler", "from_stage", "to_stage", "event", "details")]), , drop = FALSE]
  } else {
    edge_df <- data.frame()
  }

  list(
    script_path = script_path,
    orchestrator = orchestrator,
    handlers = handlers,
    edges = edge_df,
    dot = bp_network_to_dot(edge_df, title = basename(script_path))
  )
}

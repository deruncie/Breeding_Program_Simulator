# Readable functional API: list-based state + explicit stage verbs.

#' Initialize Program State
#'
#' Create an empty state container for readable breeding-program simulations.
#'
#' @param SP AlphaSimR `SimParam` object.
#' @param dt Tick length in years.
#' @param start_time Starting time in years.
#' @param sim Optional list merged into `state$sim`.
#'
#' @return A `bp_state` list.
#' @export
bp_init_state <- function(SP, dt = 0.25, start_time = 0, sim = list()) {
  tick <- as.integer(round(start_time / dt))
  state <- list(
    time = list(tick = tick, dt = as.numeric(dt), t = as.numeric(start_time)),
    sim = utils::modifyList(list(SP = SP, default_chip = 1L), sim),
    pops = list(),
    cohorts = bp_empty_cohorts(),
    event_log = bp_empty_event_log(),
    phenotype_log = bp_empty_phenotype_log(),
    genotype_log = bp_empty_genotype_log(),
    gs_models = list(),
    cost_log = bp_empty_cost_log(),
    outputs = list(varieties = bp_empty_variety_log()),
    counters = list(cohort = 0L, model = 0L, event = 0L)
  )
  class(state) <- "bp_state"
  state
}

# Event log table (append-only).
bp_empty_event_log <- function() {
  data.frame(
    event_id = integer(),
    tick = integer(),
    year = numeric(),
    fn = character(),
    event_type = character(),
    stage = character(),
    source_ids = character(),
    output_id = character(),
    event_string = character(),
    template_string = character(),
    details = I(vector("list", 0L)),
    stringsAsFactors = FALSE
  )
}

#' Conditional Debug Breakpoint
#'
#' Enter `browser()` when debug flags in `cfg` match the current state.
#'
#' @param state Program state.
#' @param cfg Configuration list with debug flags.
#' @param where Optional label for function/stage.
#' @param year Optional current year.
#' @param tick Optional current tick.
#'
#' @details
#' When a breakpoint is triggered, this helper prints a compact summary of
#' caller-local variables and then evaluates `browser()` in the caller frame so
#' debugging starts in the caller context (e.g. `run_recurrent_gs_tick`).
#'
#' @return Invisibly returns `NULL` or enters browser mode.
#' @export
bp_debug_break <- function(state, cfg, where = NULL, year = NULL, tick = NULL) {
  if (!isTRUE(cfg$debug %||% FALSE)) return(invisible(NULL))

  if (is.null(tick)) tick <- as.integer(state$time$tick)
  if (is.null(year) && !is.null(cfg$ticks_per_year)) {
    year <- as.integer(floor(tick / as.integer(cfg$ticks_per_year)) + 1L)
  }
  if (is.null(where)) {
    where <- as.character(sys.call(-1)[[1]])
  }

  if (!is.null(cfg$debug_after_year) && !is.null(year) && year < as.integer(cfg$debug_after_year)) {
    return(invisible(NULL))
  }
  if (!is.null(cfg$debug_after_tick) && tick < as.integer(cfg$debug_after_tick)) {
    return(invisible(NULL))
  }
  if (!is.null(cfg$debug_where) && !(where %in% as.character(cfg$debug_where))) {
    return(invisible(NULL))
  }

  # Compact caller-local variable summary for fast interactive inspection.
  caller <- parent.frame()
  obj_names <- ls(envir = caller, all.names = FALSE)
  if (length(obj_names) > 0L) {
    n_show <- as.integer(cfg$debug_n_vars %||% 12L)
    n_show <- max(1L, n_show)
    show_names <- head(sort(obj_names), n_show)
    parts <- vapply(show_names, function(nm) {
      x <- get(nm, envir = caller)
      cls <- class(x)[1]
      len <- tryCatch(length(x), error = function(e) NA_integer_)
      if (!is.na(len)) sprintf("%s<%s:%d>", nm, cls, len) else sprintf("%s<%s>", nm, cls)
    }, character(1))
    cat(sprintf(
      "debug: where=%s tick=%d year=%s vars=%s\n",
      where,
      as.integer(tick),
      ifelse(is.null(year), "NA", as.character(year)),
      paste(parts, collapse = ", ")
    ))
  } else {
    cat(sprintf(
      "debug: where=%s tick=%d year=%s vars=<none>\n",
      where,
      as.integer(tick),
      ifelse(is.null(year), "NA", as.character(year))
    ))
  }

  evalq(browser(), envir = caller)
}

# Cohort metadata table.
bp_empty_cohorts <- function() {
  data.frame(
    cohort_id = character(),
    stage = character(),
    stream = character(),
    cycle_id = character(),
    source_cohort_id = character(),
    selection_strategy = character(),
    cross_strategy = character(),
    created_tick = integer(),
    done_tick = integer(),
    available_tick = integer(),
    closed_tick = integer(),
    active = logical(),
    n_ind = integer(),
    genotyped = logical(),
    chips = character(),
    stringsAsFactors = FALSE
  )
}

# Long-form phenotype history table.
bp_empty_phenotype_log <- function() {
  data.frame(
    cohort_id = character(),
    stage = character(),
    individual_id = integer(),
    environment = character(),
    trait = character(),
    phenotype_value = numeric(),
    p_value = numeric(),
    measured_tick = integer(),
    available_tick = integer(),
    n_loc = integer(),
    reps = integer(),
    stringsAsFactors = FALSE
  )
}

# Long-form genotyping event table.
bp_empty_genotype_log <- function() {
  data.frame(
    cohort_id = character(),
    chip = character(),
    started_tick = integer(),
    done_tick = integer(),
    available_tick = integer(),
    n_ind = integer(),
    stringsAsFactors = FALSE
  )
}

# Append-only cost event table.
bp_empty_cost_log <- function() {
  data.frame(
    tick = integer(),
    stage = character(),
    cohort_id = character(),
    event = character(),
    unit = character(),
    n_units = numeric(),
    unit_cost = numeric(),
    total_cost = numeric(),
    stringsAsFactors = FALSE
  )
}

# Released variety history table.
bp_empty_variety_log <- function() {
  data.frame(
    tick = integer(),
    source_cohort_id = character(),
    variety_id = integer(),
    stringsAsFactors = FALSE
  )
}

# Pure helper: infer number of individuals in a pop-like object.
pop_n_ind <- function(pop) {
  if (methods::is(pop, "Pop") || methods::is(pop, "RawPop")) {
    return(as.integer(pop@nInd))
  }
  out <- tryCatch(nrow(pop), error = function(e) NA_integer_)
  if (length(out) == 1L && !is.na(out)) {
    return(as.integer(out))
  }
  stop("Unable to determine number of individuals for pop", call. = FALSE)
}

# Pure helper: subset a pop-like object by integer indices.
pop_subset <- function(pop, idx) {
  idx <- as.integer(idx)
  if (methods::is(pop, "Pop") || methods::is(pop, "RawPop")) {
    return(pop[idx])
  }
  if (is.data.frame(pop) || is.matrix(pop)) {
    return(pop[idx, , drop = FALSE])
  }
  out <- tryCatch(pop[idx], error = function(e) NULL)
  if (!is.null(out)) {
    return(out)
  }
  stop("Unable to subset pop with [idx]", call. = FALSE)
}

# Pure helper: convert years to integer ticks.
years_to_ticks <- function(dt, years) {
  as.integer(round(as.numeric(years) / as.numeric(dt)))
}

# Pure helper: normalize chip key for logs.
chip_key <- function(chip) {
  as.character(chip)
}

# Resolve chip config to numeric index required by AlphaSimR model functions.
chip_index <- function(state, chip) {
  if (is.numeric(chip)) return(as.integer(chip))
  if (is.character(chip) && !is.null(state$sim$chip_map[[chip]])) {
    return(as.integer(state$sim$chip_map[[chip]]))
  }
  if (is.character(chip) && grepl("^[0-9]+$", chip)) {
    return(as.integer(chip))
  }
  stop("chip must be numeric, numeric string, or present in state$sim$chip_map", call. = FALSE)
}

# Convert tick to decimal year where tick 0 is Year 1.0.
bp_tick_to_year <- function(state, tick) {
  1 + as.numeric(tick) * as.numeric(state$time$dt)
}

# Format decimal year for human-readable output.
bp_format_year <- function(x, digits = 2L) {
  formatC(as.numeric(x), format = "f", digits = as.integer(digits))
}

# Parse semicolon-delimited source ids into unique vector.
bp_parse_source_ids <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character(0))
  ch <- as.character(x)
  ch <- ch[!is.na(ch) & nzchar(ch) & ch != "NA"]
  if (length(ch) == 0L) return(character(0))
  out <- unique(unlist(strsplit(ch, ";", fixed = TRUE), use.names = FALSE))
  out[!is.na(out) & nzchar(out) & out != "NA"]
}

#' Build Cohort Label
#'
#' Return a readable cohort label such as `PYT:Year3.00`.
#'
#' @param state Program state.
#' @param cohort_id Cohort id.
#' @param use Whether to use cohort `created` or `available` time.
#' @param digits Number of decimal places in the year label.
#'
#' @return Character scalar label.
#' @export
bp_cohort_label <- function(state, cohort_id, use = c("created", "available"), digits = 2L) {
  use <- match.arg(use)
  idx <- match(as.character(cohort_id), state$cohorts$cohort_id)
  if (is.na(idx)) return(as.character(cohort_id))
  tick <- if (use == "created") state$cohorts$created_tick[[idx]] else state$cohorts$available_tick[[idx]]
  stage <- as.character(state$cohorts$stage[[idx]])
  paste0(stage, ":Year", bp_format_year(bp_tick_to_year(state, tick), digits = digits))
}

# Build joined labels for one or more source cohorts.
bp_source_labels <- function(state, source_ids, use = "created", digits = 2L) {
  ids <- bp_parse_source_ids(source_ids)
  if (length(ids) == 0L) return("none")
  paste(vapply(ids, function(cid) bp_cohort_label(state, cid, use = use, digits = digits), character(1)), collapse = ", ")
}

# Append one cost row.
bp_add_cost <- function(state, stage, cohort_id, event, unit, n_units, unit_cost) {
  row <- data.frame(
    tick = as.integer(state$time$tick),
    stage = as.character(stage),
    cohort_id = as.character(cohort_id),
    event = as.character(event),
    unit = as.character(unit),
    n_units = as.numeric(n_units),
    unit_cost = as.numeric(unit_cost),
    total_cost = as.numeric(n_units) * as.numeric(unit_cost),
    stringsAsFactors = FALSE
  )
  state$cost_log <- rbind(state$cost_log, row)
  state
}

#' Log Event
#'
#' Append one row to `state$event_log`.
#'
#' @param state Program state.
#' @param fn Function name that emitted this event.
#' @param event_type Event type tag.
#' @param stage Stage label.
#' @param source_ids Source cohort id(s).
#' @param output_id Output id (cohort id or model id).
#' @param event_string Human-readable event text.
#' @param template_string Normalized template text used for pattern matching.
#' @param details Optional named list of extra details.
#'
#' @return Updated program state.
#' @export
bp_log_event <- function(
  state,
  fn,
  event_type,
  stage = "",
  source_ids = character(0),
  output_id = NA_character_,
  event_string,
  template_string = NULL,
  details = list()
) {
  if (is.null(state$event_log)) {
    state$event_log <- bp_empty_event_log()
  }
  if (is.null(state$counters$event)) {
    state$counters$event <- 0L
  }

  state$counters$event <- state$counters$event + 1L
  tick <- as.integer(state$time$tick)
  yr <- bp_tick_to_year(state, tick)
  ids <- bp_parse_source_ids(source_ids)
  src_txt <- if (length(ids) == 0L) "" else paste(ids, collapse = ";")
  stg_vals <- as.character(stage %||% "")
  stg_vals <- stg_vals[!is.na(stg_vals) & nzchar(stg_vals)]
  stg_txt <- if (length(stg_vals) == 0L) "" else paste(unique(stg_vals), collapse = ";")
  one_text <- function(x, default = "") {
    if (is.null(x)) return(default)
    v <- as.character(x)
    v <- v[!is.na(v) & nzchar(v)]
    if (length(v) == 0L) return(default)
    paste(v, collapse = ";")
  }

  fn_txt <- one_text(fn, default = "unknown_fn")
  type_txt <- one_text(event_type, default = "unknown_event")
  out_txt <- one_text(output_id %||% NA_character_, default = NA_character_)
  msg_txt <- one_text(event_string, default = "")
  tpl <- one_text(template_string %||% msg_txt, default = msg_txt)

  row <- data.frame(
    event_id = as.integer(state$counters$event),
    tick = tick,
    year = as.numeric(yr),
    fn = fn_txt,
    event_type = type_txt,
    stage = stg_txt,
    source_ids = as.character(src_txt),
    output_id = out_txt,
    event_string = msg_txt,
    template_string = tpl,
    details = I(list(details)),
    stringsAsFactors = FALSE
  )
  state$event_log <- rbind(state$event_log, row)
  state
}

#' Build Event Timeline Data Frame
#'
#' @param state Program state.
#' @param digits Number of decimal places for year labels.
#'
#' @return Event log with derived timeline columns.
#' @export
bp_event_timeline_df <- function(state, digits = 2L) {
  ev <- state$event_log
  if (is.null(ev) || nrow(ev) == 0L) return(ev)
  ev <- ev[order(ev$tick, ev$event_id), , drop = FALSE]
  ev$year_label <- paste0("Year ", bp_format_year(ev$year, digits = digits))
  ev$year_int <- as.integer(floor(ev$year))
  ev$within_year <- round(ev$year - ev$year_int, digits = digits)
  ev
}

#' Print Event Timeline
#'
#' Print a readable event timeline grouped by time, with optional collapse of
#' repeated year-level patterns.
#'
#' @param state Program state.
#' @param collapse_year_patterns Whether to collapse consecutive years with the
#'   same event pattern.
#' @param digits Number of decimal places for printed years.
#' @param md_file Optional path to save the rendered timeline as a markdown
#'   text file.
#'
#' @return Invisibly returns `state`.
#' @export
bp_print_event_timeline <- function(state, collapse_year_patterns = TRUE, digits = 2L, md_file = NULL) {
  ev <- bp_event_timeline_df(state, digits = digits)
  if (is.null(ev) || nrow(ev) == 0L) {
    msg <- "No events logged."
    cat(msg, "\n", sep = "")
    if (!is.null(md_file)) writeLines(msg, con = md_file, useBytes = TRUE)
    return(invisible(state))
  }

  fmt_sentences <- function(msg) {
    txt <- trimws(as.character(msg))
    txt <- sub("^Year[[:space:]]+[0-9]+(?:\\.[0-9]+)?:[[:space:]]*", "", txt, perl = TRUE)
    if (!nzchar(txt)) return(character(0))

    parts <- trimws(unlist(strsplit(txt, "\\.\\s+", perl = TRUE), use.names = FALSE))
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0L) return(character(0))
    parts <- vapply(parts, function(x) {
      if (grepl("[.!?]$", x)) x else paste0(x, ".")
    }, character(1))
    parts
  }

  cat_wrapped <- function(prefix, text, cont_prefix) {
    w <- getOption("width", 100L)
    lines <- strwrap(text, width = max(40L, as.integer(w) - nchar(prefix)))
    if (length(lines) == 0L) return(invisible(NULL))
    cat(prefix, lines[[1L]], "\n", sep = "")
    if (length(lines) > 1L) {
      for (j in 2:length(lines)) cat(cont_prefix, lines[[j]], "\n", sep = "")
    }
    invisible(NULL)
  }

  print_year_detail <- function(df_year) {
    splits <- split(df_year, df_year$year_label)
    ord <- names(splits)[order(vapply(splits, function(x) x$tick[[1]], integer(1)))]
    first_block <- TRUE
    for (lbl in ord) {
      if (!first_block) cat("\n")
      cat(lbl, "\n", sep = "")
      d <- splits[[lbl]]
      for (i in seq_len(nrow(d))) {
        parts <- fmt_sentences(d$event_string[[i]])
        if (length(parts) == 0L) next
        cat_wrapped("- ", parts[[1L]], "  ")
        if (length(parts) > 1L) {
          for (k in 2:length(parts)) {
            cat_wrapped("  - ", parts[[k]], "    ")
          }
        }
      }
      first_block <- FALSE
    }
  }

  by_year <- split(ev, ev$year_int)
  years <- as.integer(names(by_year))
  years <- years[order(years)]

  render_timeline <- function() {
    if (!isTRUE(collapse_year_patterns) || length(years) <= 1L) {
      first_year <- TRUE
      for (yy in years) {
        if (!first_year) cat("\n")
        print_year_detail(by_year[[as.character(yy)]])
        first_year <- FALSE
      }
      return(invisible(NULL))
    }

    sig <- vapply(years, function(yy) {
      d <- by_year[[as.character(yy)]]
      paste(sprintf("%s|%s", d$within_year, d$template_string), collapse = " || ")
    }, character(1))

    run_start <- 1L
    first_year <- TRUE
    while (run_start <= length(years)) {
      run_end <- run_start
      while (run_end < length(years) && identical(sig[[run_end + 1L]], sig[[run_start]])) {
        run_end <- run_end + 1L
      }

      y0 <- years[[run_start]]
      if (!first_year) cat("\n")
      print_year_detail(by_year[[as.character(y0)]])
      if (run_end > run_start) {
        y1 <- years[[run_end]]
        cat(sprintf("Year %d-%d (same as Year %d)\n", y0, y1, y0))
      }
      run_start <- run_end + 1L
      first_year <- FALSE
    }
    invisible(NULL)
  }

  txt <- utils::capture.output(render_timeline(), type = "output")
  cat(paste(txt, collapse = "\n"), "\n", sep = "")
  if (!is.null(md_file)) {
    writeLines(txt, con = md_file, useBytes = TRUE)
  }

  invisible(state)
}

# Add a new cohort row and store its pop.
bp_add_cohort <- function(
  state,
  pop,
  stage,
  stream = "main",
  cycle_id = "cycle_1",
  source_cohort_id = NA_character_,
  selection_strategy = NA_character_,
  cross_strategy = NA_character_,
  duration_years = 0,
  active = TRUE
) {
  state$counters$cohort <- state$counters$cohort + 1L
  cohort_id <- sprintf("cohort_%07d", state$counters$cohort)

  dur_ticks <- years_to_ticks(state$time$dt, duration_years)
  done_tick <- as.integer(state$time$tick + dur_ticks)

  row <- data.frame(
    cohort_id = cohort_id,
    stage = as.character(stage),
    stream = as.character(stream),
    cycle_id = as.character(cycle_id),
    source_cohort_id = as.character(source_cohort_id),
    selection_strategy = as.character(selection_strategy),
    cross_strategy = as.character(cross_strategy),
    created_tick = as.integer(state$time$tick),
    done_tick = done_tick,
    available_tick = done_tick,
    closed_tick = NA_integer_,
    active = isTRUE(active),
    n_ind = pop_n_ind(pop),
    genotyped = FALSE,
    chips = "",
    stringsAsFactors = FALSE
  )

  state$cohorts <- rbind(state$cohorts, row)
  state$pops[[cohort_id]] <- pop
  state
}

# Propagate genotype-log availability from one source cohort to a copied subset cohort.
bp_inherit_genotypes_from_source <- function(state, new_cohort_id, source_ids) {
  ids <- unique(as.character(source_ids))
  ids <- ids[!is.na(ids) & nzchar(ids) & ids != "NA"]
  if (length(ids) != 1L) {
    return(state)
  }
  src <- ids[[1L]]
  src_rows <- state$genotype_log[state$genotype_log$cohort_id == src, , drop = FALSE]
  if (nrow(src_rows) == 0L) return(state)

  new_idx <- match(new_cohort_id, state$cohorts$cohort_id)
  if (is.na(new_idx)) return(state)
  new_n <- as.integer(state$cohorts$n_ind[new_idx])

  chips <- unique(as.character(src_rows$chip))
  for (ckey in chips) {
    if (any(state$genotype_log$cohort_id == new_cohort_id & state$genotype_log$chip == ckey)) next

    s2 <- src_rows[src_rows$chip == ckey, , drop = FALSE]
    j <- which.min(s2$available_tick)
    row <- s2[j, , drop = FALSE]
    row$cohort_id <- as.character(new_cohort_id)
    row$n_ind <- new_n
    state$genotype_log <- rbind(state$genotype_log, row)
  }

  bp_refresh_genotyped_flags(state)
}

# Return the most recently created cohort id.
bp_last_cohort_id <- function(state) {
  if (nrow(state$cohorts) == 0L) return(NA_character_)
  as.character(state$cohorts$cohort_id[nrow(state$cohorts)])
}

# Mark a cohort inactive.
bp_close_cohort <- function(state, cohort_id) {
  idx <- match(cohort_id, state$cohorts$cohort_id)
  if (is.na(idx)) return(state)
  state$cohorts$active[idx] <- FALSE
  state$cohorts$closed_tick[idx] <- as.integer(state$time$tick)
  state
}

#' Get Ready Cohorts
#'
#' Return cohorts available for use at `as_of_tick`, optionally filtered by
#' stage/stream and activity status.
#'
#' @param state Program state.
#' @param stage Optional stage filter.
#' @param stream Optional stream filter.
#' @param active_only Whether to return only active cohorts.
#' @param as_of_tick Tick used for availability filtering.
#'
#' @return Cohort metadata `data.frame`.
#' @export
bp_get_ready_cohorts <- function(state, stage = NULL, stream = NULL, active_only = TRUE, as_of_tick = state$time$tick) {
  df <- state$cohorts
  if (!is.null(stage)) {
    df <- df[df$stage %in% stage, , drop = FALSE]
  }
  if (!is.null(stream)) {
    df <- df[df$stream %in% stream, , drop = FALSE]
  }
  if (isTRUE(active_only)) {
    df <- df[df$active, , drop = FALSE]
  }
  df <- df[df$available_tick <= as.integer(as_of_tick), , drop = FALSE]
  df
}

# Select which source cohorts to process from a ready cohort table.
bp_select_source_rows <- function(state, ready, cfg) {
  if (nrow(ready) == 0L) return(ready)

  policy <- as.character(cfg$input_policy %||% "latest_one")
  if (policy == "all_ready") {
    return(ready)
  }

  ready_ord <- ready[order(ready$available_tick, ready$created_tick, decreasing = TRUE), , drop = FALSE]

  if (policy == "latest_one") {
    return(ready_ord[1, , drop = FALSE])
  }

  if (policy == "latest_n") {
    n <- as.integer(cfg$input_n %||% 1L)
    n <- max(1L, min(n, nrow(ready_ord)))
    return(ready_ord[seq_len(n), , drop = FALSE])
  }

  if (policy == "by_cycle") {
    cycles <- as.character(cfg$input_cycle_ids %||% cfg$cycle_id %||% NA_character_)
    cycles <- cycles[!is.na(cycles) & nzchar(cycles)]
    if (length(cycles) == 0L) {
      return(ready_ord[0, , drop = FALSE])
    }
    return(ready_ord[ready_ord$cycle_id %in% cycles, , drop = FALSE])
  }

  if (policy == "custom") {
    fn <- cfg$select_source_cohorts_fn
    if (!is.function(fn)) {
      stop("input_policy='custom' requires cfg$select_source_cohorts_fn", call. = FALSE)
    }
    out <- fn(state, ready_ord, cfg)
    if (is.null(out) || length(out) == 0L) {
      return(ready_ord[0, , drop = FALSE])
    }
    if (is.logical(out)) {
      if (length(out) != nrow(ready_ord)) {
        stop("custom selector logical output must match nrow(ready)", call. = FALSE)
      }
      return(ready_ord[out, , drop = FALSE])
    }
    if (is.numeric(out)) {
      idx <- as.integer(out)
      idx <- idx[idx >= 1L & idx <= nrow(ready_ord)]
      idx <- unique(idx)
      return(ready_ord[idx, , drop = FALSE])
    }
    ids <- as.character(out)
    return(ready_ord[ready_ord$cohort_id %in% ids, , drop = FALSE])
  }

  stop(sprintf("Unknown input_policy: %s", policy), call. = FALSE)
}

# Standardized behavior when a stage call has no eligible source cohorts.
bp_handle_no_ready <- function(cfg, fn_name, stage_label, context = "no available source cohorts") {
  msg <- sprintf("%s: %s for stage '%s'.", fn_name, context, stage_label)
  if (isTRUE(cfg$fail_if_no_ready %||% FALSE)) {
    stop(msg, call. = FALSE)
  }
  if (!isTRUE(cfg$silent %||% FALSE)) {
    cat(msg, "\n")
  }
}

#' Skip Event When Input Is Missing
#'
#' Standard helper for event verbs when an input cohort bundle is absent.
#'
#' @param state Program state.
#' @param input_obj Input bundle object (typically from [select_latest_available()]).
#' @param cfg Configuration list. Recognized flags:
#'   `fail_on_missing_input` (default `FALSE`) and
#'   `log_missing_input` (default `TRUE`).
#' @param event_name Optional event-function name override. If `NULL`, inferred
#'   from caller name.
#'
#' @return A list with `state` and logical `skip`.
#' @export
bp_skip_if_no_input <- function(state, input_obj, cfg, event_name = NULL) {
  if (!is.null(input_obj)) {
    return(list(state = state, skip = FALSE))
  }

  if (is.null(event_name)) {
    event_name <- tryCatch(as.character(sys.call(-1)[[1]]), error = function(e) "unknown_event")
  }
  event_name <- as.character(event_name %||% "unknown_event")
  msg <- sprintf("No input available for %s()", event_name)

  if (isTRUE(cfg$fail_on_missing_input %||% FALSE)) {
    stop(msg, call. = FALSE)
  }

  if (isTRUE(cfg$log_missing_input %||% TRUE)) {
    state <- bp_log_event(
      state = state,
      fn = event_name,
      event_type = "no_input_skip",
      stage = "",
      source_ids = character(0),
      output_id = NA_character_,
      event_string = msg,
      template_string = sprintf("No input for %s()", event_name),
      details = list()
    )
  }

  list(state = state, skip = TRUE)
}

#' Get Ready Pop Bundle
#'
#' Select source cohort(s), retrieve pop(s), and optionally merge them.
#'
#' @param state Program state.
#' @param stage Input stage(s).
#' @param stream Optional stream filter.
#' @param policy Source selection policy.
#' @param combine Whether to merge selected pops into one.
#' @param input_n Number of cohorts for `policy = "latest_n"`.
#' @param cycle_id Optional cycle id.
#' @param input_cycle_ids Optional cycle id vector.
#' @param select_source_cohorts_fn Custom source selector for `policy = "custom"`.
#' @param include_not_ready If `TRUE`, search active cohorts regardless of
#'   `available_tick`. Default `FALSE`.
#' @param silent Suppress no-ready messages.
#' @param fail_if_no_ready Error when no cohorts are eligible.
#'
#' @return `NULL` or a source bundle list with `pop`, `source_ids`, and metadata.
#' @export
get_ready_pop <- function(
  state,
  stage,
  stream = NULL,
  policy = "latest_one",
  combine = TRUE,
  input_n = NULL,
  cycle_id = NULL,
  input_cycle_ids = NULL,
  select_source_cohorts_fn = NULL,
  include_not_ready = FALSE,
  silent = FALSE,
  fail_if_no_ready = FALSE
) {
  cfg <- list(
    input_policy = policy,
    input_n = input_n,
    cycle_id = cycle_id,
    input_cycle_ids = input_cycle_ids,
    select_source_cohorts_fn = select_source_cohorts_fn,
    include_not_ready = isTRUE(include_not_ready),
    silent = isTRUE(silent),
    fail_if_no_ready = isTRUE(fail_if_no_ready)
  )
  ready <- if (isTRUE(include_not_ready)) {
    df <- state$cohorts
    if (!is.null(stage)) {
      df <- df[df$stage %in% stage, , drop = FALSE]
    }
    if (!is.null(stream)) {
      df <- df[df$stream %in% stream, , drop = FALSE]
    }
    df[df$active, , drop = FALSE]
  } else {
    bp_get_ready_cohorts(state, stage = stage, stream = stream)
  }
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "get_ready_pop", stage)
    return(NULL)
  }
  selected <- bp_select_source_rows(state, ready, cfg)
  if (nrow(selected) == 0L) {
    bp_handle_no_ready(cfg, "get_ready_pop", stage, context = "source selection policy returned no cohorts")
    return(NULL)
  }

  pops <- lapply(selected$cohort_id, function(cid) state$pops[[cid]])
  pop <- if (isTRUE(combine)) merge_pops(pops) else pops[[1L]]
  source_ids <- as.character(selected$cohort_id)
  cycle_values <- unique(as.character(selected$cycle_id))
  cycle_out <- if (length(cycle_values) == 1L) cycle_values else paste(cycle_values, collapse = ";")

  list(
    pop = pop,
    source_ids = source_ids,
    source_rows = selected,
    stage = as.character(stage),
    stream = if (is.null(stream)) as.character(selected$stream[[1L]]) else as.character(stream),
    cycle_id = cycle_out
  )
}

#' Select Latest Available Cohort Bundle
#'
#' Convenience wrapper around [get_ready_pop()] with `policy = "latest_one"`.
#'
#' @param state Program state.
#' @param stage Stage name(s).
#' @param stream Optional stream filter.
#' @param n Number of latest ready cohorts to select. Default `1L`.
#' @param combine Whether to merge selected pops.
#' @param include_not_ready If `TRUE`, search active cohorts regardless of
#'   `available_tick`. Default `FALSE`.
#' @param silent Suppress no-ready messages.
#' @param fail_if_no_ready Error when no cohorts are eligible.
#'
#' @return `NULL` or a source bundle list (`pop`, `source_ids`, ...).
#' @export
select_latest_available <- function(
  state,
  stage,
  stream = NULL,
  n = 1L,
  combine = TRUE,
  include_not_ready = FALSE,
  silent = FALSE,
  fail_if_no_ready = FALSE
) {
  n <- as.integer(n %||% 1L)
  if (!is.finite(n) || n < 1L) {
    stop("select_latest_available: n must be >= 1.", call. = FALSE)
  }
  if (!isTRUE(combine) && n > 1L) {
    stop("select_latest_available: combine=FALSE is ambiguous when n > 1; use combine=TRUE or n=1.", call. = FALSE)
  }

  get_ready_pop(
    state = state,
    stage = stage,
    stream = stream,
    policy = if (n == 1L) "latest_one" else "latest_n",
    combine = combine,
    input_n = if (n == 1L) NULL else n,
    include_not_ready = include_not_ready,
    silent = silent,
    fail_if_no_ready = fail_if_no_ready
  )
}

#' Select Current Stage Pop
#'
#' Convenience selector for quick inspection/debugging of the current cohort
#' at a stage.
#'
#' @param state Program state.
#' @param stage Stage name(s).
#' @param stream Optional stream filter.
#' @param policy Source selection policy (default `"latest_one"`).
#' @param combine Whether to merge selected pops.
#' @param bundle Return source bundle instead of pop.
#' @param include_not_ready If `TRUE`, search active cohorts regardless of
#'   `available_tick`. Default `FALSE`.
#' @param silent Suppress no-ready messages.
#' @param fail_if_no_ready Error when no cohorts are eligible.
#'
#' @return Pop object by default, or source bundle if `bundle = TRUE`.
#' @export
select_current <- function(
  state,
  stage,
  stream = NULL,
  policy = "latest_one",
  combine = TRUE,
  bundle = FALSE,
  include_not_ready = FALSE,
  silent = TRUE,
  fail_if_no_ready = FALSE
) {
  src <- get_ready_pop(
    state = state,
    stage = stage,
    stream = stream,
    policy = policy,
    combine = combine,
    include_not_ready = include_not_ready,
    silent = silent,
    fail_if_no_ready = fail_if_no_ready
  )
  if (is.null(src)) return(NULL)
  if (isTRUE(bundle)) return(src)
  src$pop
}

#' Add Output Stage Cohort
#'
#' Create a new cohort from an output pop, carrying defaults from source metadata.
#'
#' @param state Program state.
#' @param pop Output pop.
#' @param stage Output stage.
#' @param source Optional source bundle from [get_ready_pop()] or source
#'   cohort id(s) as character vector.
#' @param source_ids Optional explicit source cohort id vector. When provided,
#'   this overrides ids inferred from `source`.
#' @param ready_in_years Delay until cohort availability.
#' @param stream Optional output stream override.
#' @param cycle_id Optional output cycle id override.
#' @param active Whether the new cohort is active.
#' @param inherit_genotypes Whether to inherit genotype availability from source.
#' @param selection_strategy Optional human-readable selection mechanism text.
#' @param cross_strategy Optional human-readable crossing mechanism text.
#' @param cost_per_individual Optional cost logged once per individual in `pop`.
#' @param cost_event Event label used when `cost_per_individual` is provided.
#' @param cost_unit Unit label used when `cost_per_individual` is provided.
#'
#' @return Updated program state.
#' @export
put_stage_pop <- function(
  state,
  pop,
  stage,
  source = NULL,
  source_ids = NULL,
  ready_in_years = 0,
  stream = NULL,
  cycle_id = NULL,
  active = TRUE,
  inherit_genotypes = FALSE,
  selection_strategy = NA_character_,
  cross_strategy = NA_character_,
  cost_per_individual = NULL,
  cost_event = "cohort_creation",
  cost_unit = "individual"
) {
  source_is_bundle <- is.list(source) && !is.null(source$source_ids)
  source_id_vec <- if (!is.null(source_ids)) {
    as.character(source_ids)
  } else if (source_is_bundle) {
    as.character(source$source_ids)
  } else if (is.character(source)) {
    as.character(source)
  } else {
    character(0)
  }
  source_id_vec <- unique(source_id_vec[!is.na(source_id_vec) & nzchar(source_id_vec) & source_id_vec != "NA"])

  src_ids <- if (length(source_id_vec) == 0L) NA_character_ else paste(source_id_vec, collapse = ";")
  stream_val <- if (!is.null(stream)) stream else if (source_is_bundle && !is.null(source$stream)) source$stream else "main"
  cycle_val <- if (!is.null(cycle_id)) cycle_id else if (source_is_bundle && !is.null(source$cycle_id)) source$cycle_id else "cycle_1"

  state <- bp_add_cohort(
    state = state,
    pop = pop,
    stage = stage,
    stream = stream_val,
    cycle_id = cycle_val,
    source_cohort_id = src_ids,
    selection_strategy = selection_strategy,
    cross_strategy = cross_strategy,
    duration_years = ready_in_years,
    active = active
  )
  new_cohort_id <- bp_last_cohort_id(state)
  avail_tick <- state$cohorts$available_tick[match(new_cohort_id, state$cohorts$cohort_id)]
  yr_now <- bp_format_year(bp_tick_to_year(state, state$time$tick))
  yr_av <- bp_format_year(bp_tick_to_year(state, avail_tick))
  src_lbl <- bp_source_labels(state, source_id_vec, use = "created")
  ss_txt <- as.character(selection_strategy %||% NA_character_)
  cs_txt <- as.character(cross_strategy %||% NA_character_)
  extra <- c(
    if (!is.na(ss_txt) && nzchar(ss_txt)) paste0("Selection=", ss_txt) else NULL,
    if (!is.na(cs_txt) && nzchar(cs_txt)) paste0("Crossing=", cs_txt) else NULL
  )
  suffix <- if (length(extra) > 0L) paste0(" ", paste(extra, collapse = "; "), ".") else ""
  evt <- sprintf(
    "Year %s: Created %s from %s. Will be available Year %s.%s",
    yr_now, stage, src_lbl, yr_av, suffix
  )
  tpl <- sprintf(
    "Create %s from source_stage=%s",
    stage,
    paste(unique(vapply(source_id_vec, function(cid) {
      idx <- match(cid, state$cohorts$cohort_id)
      if (is.na(idx)) "unknown" else as.character(state$cohorts$stage[[idx]])
    }, character(1))), collapse = "+")
  )
  state <- bp_log_event(
    state = state,
    fn = "put_stage_pop",
    event_type = "stage_output",
    stage = stage,
    source_ids = source_id_vec,
    output_id = new_cohort_id,
    event_string = evt,
    template_string = tpl,
    details = list(
      ready_in_years = as.numeric(ready_in_years),
      selection_strategy = selection_strategy,
      cross_strategy = cross_strategy
    )
  )

  if (isTRUE(inherit_genotypes) && length(source_id_vec) > 0L) {
    state <- bp_inherit_genotypes_from_source(
      state = state,
      new_cohort_id = new_cohort_id,
      source_ids = source_id_vec
    )
  }
  if (!is.null(cost_per_individual)) {
    state <- bp_add_cost(
      state = state,
      stage = stage,
      cohort_id = new_cohort_id,
      event = as.character(cost_event %||% "cohort_creation"),
      unit = as.character(cost_unit %||% "individual"),
      n_units = pop_n_ind(pop),
      unit_cost = as.numeric(cost_per_individual)
    )
  }
  state
}

#' Add Stage Cost
#'
#' Append a cost event row to `state$cost_log`.
#'
#' @param state Program state.
#' @param event Event label.
#' @param n_units Number of units.
#' @param unit_cost Cost per unit.
#' @param stage Optional stage label.
#' @param unit Unit label.
#' @param cohort_id Optional cohort id.
#' @param n Alias for `n_units`.
#'
#' @return Updated program state.
#' @export
add_stage_cost <- function(
  state,
  event,
  n_units = NULL,
  unit_cost,
  stage = NULL,
  unit = "unit",
  cohort_id = NULL,
  n = NULL
) {
  if (is.null(n_units)) {
    n_units <- n
  }
  if (is.null(n_units)) {
    stop("add_stage_cost requires n_units (or n)", call. = FALSE)
  }
  cid <- cohort_id %||% bp_last_cohort_id(state)
  stg <- stage
  if (is.null(stg) || !nzchar(stg)) {
    idx <- match(cid, state$cohorts$cohort_id)
    stg <- if (is.na(idx)) "unknown" else state$cohorts$stage[[idx]]
  }
  bp_add_cost(
    state = state,
    stage = stg,
    cohort_id = cid,
    event = event,
    unit = unit,
    n_units = n_units,
    unit_cost = unit_cost
  )
}

#' Log Trial Phenotypes
#'
#' Log trial phenotypes from `pop_trial@pheno` into `state$phenotype_log`.
#'
#' @param state Program state.
#' @param pop_trial Pop with phenotype values in `@pheno`.
#' @param stage Stage label.
#' @param source Optional source bundle.
#' @param traits Trait index/indices.
#' @param n_loc Number of locations represented.
#' @param reps Number of replicates represented.
#' @param environment Environment code (`0` for aggregate).
#' @param p_value Optional environment offset.
#'
#' @return Updated program state.
#' @export
log_trial_pheno <- function(
  state,
  pop_trial,
  stage,
  source = NULL,
  traits = 1L,
  n_loc = 1L,
  reps = 1L,
  environment = 0L,
  p_value = NA_real_
) {
  cid <- bp_last_cohort_id(state)
  if (is.na(cid)) return(state)
  avail_tick <- state$cohorts$available_tick[match(cid, state$cohorts$cohort_id)]
  ph <- pop_trial@pheno
  tr <- as.integer(traits)
  if (is.null(dim(ph))) {
    ph <- matrix(ph, ncol = length(tr))
  }
  bp_record_pheno(
    state = state,
    cohort_id = cid,
    stage = stage,
    individual_id = pop_trial@id,
    traits = tr,
    pheno_matrix = ph,
    available_tick = avail_tick,
    n_loc = as.integer(n_loc),
    reps = as.integer(reps),
    environment = as.integer(environment),
    p_value = as.numeric(p_value)
  )
}

#' Close Source Cohorts
#'
#' Mark all source cohorts in a source bundle as inactive.
#'
#' @param state Program state.
#' @param source Source bundle from [get_ready_pop()].
#'
#' @return Updated program state.
#' @export
close_sources <- function(state, source) {
  if (is.null(source) || is.null(source$source_ids)) return(state)
  for (cid in as.character(source$source_ids)) {
    state <- bp_close_cohort(state, cid)
  }
  state
}

#' Advance Time
#'
#' Advance simulation time by integer ticks.
#'
#' @param state Program state.
#' @param n_ticks Number of ticks to advance.
#'
#' @return Updated program state.
#' @export
bp_advance_time <- function(state, n_ticks = 1L) {
  state$time$tick <- as.integer(state$time$tick + as.integer(n_ticks))
  state$time$t <- as.numeric(state$time$tick * state$time$dt)
  state <- bp_refresh_genotyped_flags(state)
  state
}

# Update cohort-level genotyped/chips fields based on available genotype records.
bp_refresh_genotyped_flags <- function(state) {
  if (nrow(state$cohorts) == 0L) return(state)
  chips_done <- state$genotype_log[state$genotype_log$available_tick <= state$time$tick, , drop = FALSE]

  state$cohorts$genotyped <- FALSE
  state$cohorts$chips <- ""

  if (nrow(chips_done) == 0L) return(state)

  by_cohort <- split(chips_done$chip, chips_done$cohort_id)
  for (cid in names(by_cohort)) {
    idx <- match(cid, state$cohorts$cohort_id)
    if (is.na(idx)) next
    uniq <- unique(as.character(by_cohort[[cid]]))
    state$cohorts$genotyped[idx] <- length(uniq) > 0L
    state$cohorts$chips[idx] <- paste(uniq, collapse = ";")
  }
  state
}

#' Record Phenotypes in Long Format
#'
#' Append phenotype observations to `state$phenotype_log`.
#'
#' This helper accepts a dense or sparse phenotype matrix with one row per
#' individual and one column per trait label in `traits`. Output is stored in
#' long format (`one row = one observed phenotype cell`).
#'
#' For sparse MET workflows, pass trait columns that encode trial/environment
#' identity (for example `yield_env1`, `yield_env2`, ...) and set
#' `drop_na = TRUE` so missing cells are not logged.
#'
#' @param state Program state.
#' @param cohort_id Output cohort id.
#' @param stage Output stage name.
#' @param individual_id Integer vector of individual IDs; length must match
#'   rows in `pheno_matrix`.
#' @param traits Trait labels, either integer trait codes or character names.
#'   If `NULL`, labels are inferred from `colnames(pheno_matrix)` when present,
#'   otherwise default to `trait1`, `trait2`, ...
#' @param pheno_matrix Numeric vector or matrix of phenotypes. If a vector is
#'   provided, it is treated as a one-column matrix.
#' @param available_tick Tick when these observations become available.
#' @param n_loc Number of locations represented by this record block.
#' @param reps Number of reps represented by this record block.
#' @param environment Environment id/label(s) for this record block. Can be a
#'   scalar or a vector. If a vector is provided, it is aligned to phenotype
#'   columns with the same recycling rules as `traits`. If `NULL`, defaults to
#'   `NA`.
#' @param p_value Optional environment mean/shift value for this record block.
#' @param drop_na Logical; if `TRUE`, rows with `NA` phenotype values are
#'   omitted from the appended log.
#'
#' @return Updated program state.
#' @export
bp_record_pheno <- function(
  state,
  cohort_id,
  stage,
  individual_id,
  traits = NULL,
  pheno_matrix,
  available_tick,
  n_loc,
  reps,
  environment = NULL,
  p_value = NA_real_,
  drop_na = TRUE
) {
  ph <- pheno_matrix
  if (is.null(dim(ph))) {
    ph <- matrix(ph, ncol = 1)
  }
  recycle_to_cols <- function(x, k, arg_name) {
    if (length(x) == 1L) return(rep(x, k))
    if (length(x) == k) return(x)
    if (k %% length(x) == 0L) return(rep(x, length.out = k))
    stop(
      sprintf("bp_record_pheno: length(%s)=%d is not compatible with ncol(pheno_matrix)=%d", arg_name, length(x), k),
      call. = FALSE
    )
  }

  tr <- traits
  if (is.null(tr)) {
    if (!is.null(colnames(ph)) && length(colnames(ph)) == ncol(ph)) {
      tr_labels <- as.character(colnames(ph))
    } else {
      tr_labels <- paste0("trait", seq_len(ncol(ph)))
    }
  } else if (is.numeric(tr)) {
    tr_labels <- paste0("trait", as.integer(recycle_to_cols(tr, ncol(ph), "traits")))
  } else {
    tr_labels <- as.character(recycle_to_cols(tr, ncol(ph), "traits"))
  }
  env_labels <- if (is.null(environment)) rep(NA, ncol(ph)) else recycle_to_cols(environment, ncol(ph), "environment")
  if (nrow(ph) != length(individual_id)) {
    stop("bp_record_pheno: nrow(pheno_matrix) must match length(individual_id)", call. = FALSE)
  }
  if (ncol(ph) != length(tr_labels)) {
    stop("bp_record_pheno: ncol(pheno_matrix) must match length(traits)", call. = FALSE)
  }
  p_vals <- recycle_to_cols(as.numeric(p_value), ncol(ph), "p_value")

  rows <- do.call(rbind, lapply(seq_along(tr_labels), function(k) {
    data.frame(
      cohort_id = as.character(cohort_id),
      stage = as.character(stage),
      individual_id = as.integer(individual_id),
      environment = as.character(env_labels[k]),
      trait = tr_labels[k],
      phenotype_value = as.numeric(ph[, k]),
      p_value = as.numeric(p_vals[k]),
      measured_tick = as.integer(state$time$tick),
      available_tick = as.integer(available_tick),
      n_loc = as.integer(n_loc),
      reps = as.integer(reps),
      stringsAsFactors = FALSE
    )
  }))
  if (isTRUE(drop_na) && nrow(rows) > 0L) {
    rows <- rows[!is.na(rows$phenotype_value), , drop = FALSE]
  }
  if (nrow(rows) == 0L) return(state)

  state$phenotype_log <- rbind(state$phenotype_log, rows)
  state
}

# Resolve explicit per-call trial environments with location and year effects.
bp_resolve_trial_env <- function(cfg, n_loc) {
  n_loc <- as.integer(n_loc)
  base_means <- if (!is.null(cfg$env_means)) {
    as.numeric(cfg$env_means)
  } else if (!is.null(cfg$env_mean_mu)) {
    mu <- as.numeric(cfg$env_mean_mu)
    if (length(mu) == 1L) rep(mu, n_loc) else mu
  } else {
    rep(0, n_loc)
  }
  if (length(base_means) != n_loc) {
    stop("env_means/env_mean_mu length must equal n_loc", call. = FALSE)
  }

  loc_sd <- as.numeric(cfg$env_mean_sd %||% 0)
  if (length(loc_sd) != 1L || !is.finite(loc_sd) || loc_sd < 0) {
    stop("env_mean_sd must be one non-negative numeric value", call. = FALSE)
  }
  year_sd <- as.numeric(cfg$env_year_sd %||% 0)
  if (length(year_sd) != 1L || !is.finite(year_sd) || year_sd < 0) {
    stop("env_year_sd must be one non-negative numeric value", call. = FALSE)
  }

  loc_dev <- stats::rnorm(n_loc, mean = 0, sd = loc_sd)
  year_eff <- stats::rnorm(1L, mean = 0, sd = year_sd)
  list(
    base_means = base_means,
    loc_dev = loc_dev,
    year_eff = as.numeric(year_eff),
    z_env = as.numeric(base_means + loc_dev + year_eff)
  )
}

# Convert latent environment deviation to AlphaSimR p-scale for each trait.
bp_env_p_from_latent <- function(simParam, traits, latent_env) {
  traits <- as.integer(traits)
  latent_env <- as.numeric(latent_env)
  out <- vapply(traits, function(trait_idx) {
    tr_obj <- simParam$traits[[trait_idx]]
    if ("envVar" %in% methods::slotNames(tr_obj)) {
      env_sd <- sqrt(as.numeric(tr_obj@envVar))
      if (!is.finite(env_sd) || env_sd <= 0) {
        return(0.5)
      }
      return(stats::pnorm(latent_env / env_sd))
    }
    0.5
  }, numeric(1))
  as.numeric(out)
}

# Pure helper: merge one or more AlphaSimR pops.
merge_pops <- function(pop_list) {
  if (length(pop_list) == 1L) return(pop_list[[1L]])
  if (all(vapply(pop_list, function(x) methods::is(x, "Pop"), logical(1)))) {
    return(AlphaSimR::mergePops(pop_list))
  }
  if (all(vapply(pop_list, is.data.frame, logical(1)))) {
    return(do.call(rbind, pop_list))
  }
  if (all(vapply(pop_list, is.matrix, logical(1)))) {
    return(do.call(rbind, pop_list))
  }
  stop("merge_pops requires homogeneous pop types (all Pop, all data.frame, or all matrix)", call. = FALSE)
}

#' Run Phenotyping Trial
#'
#' Generic field-trial runner with phenotype generation, phenotype logging,
#' output cohort creation, and cost logging.
#'
#' @param state Program state.
#' @param cfg Legacy configuration list (compatibility mode). Prefer explicit
#'   arguments with `pop`.
#' @param pop Pop object to phenotype (recommended mode).
#' @param output_stage Output stage name for created trial cohort.
#' @param input_cohorts Source cohort id(s), source bundle, or data frame with
#'   `cohort_id` column used for lineage logging.
#' @param selection_strategy Human-readable selection strategy text for event
#'   logging.
#' @param traits Trait index/vector. Default `1L`.
#' @param n_loc Number of locations represented. Default `1L`.
#' @param reps Replications per location. Default `1L`.
#' @param varE Residual variance passed to `AlphaSimR::setPheno`.
#' @param duration_years Delay until output cohort is available.
#' @param stream Output stream label (default `"main"`).
#' @param output_stream Optional output stream override.
#' @param cycle_id Optional output cycle id override.
#' @param cost_per_plot Cost per plot for logging.
#' @param trial_name Optional trial label for logging.
#' @param use_env_control Logical; when `TRUE`, use explicit environment
#'   generation (`onlyPheno=TRUE`) and aggregate line means.
#' @param env_means Optional fixed long-term latent environment means vector
#'   (length `n_loc`) on the normal scale.
#' @param env_mean_mu Long-term latent environment means when `env_means` is
#'   not provided. May be a scalar (recycled to all locations) or a vector of
#'   length `n_loc`.
#' @param env_mean_sd SD of location-specific deviation around the supplied
#'   long-term means for this trial call.
#' @param env_year_sd SD of one common year effect shared by all locations in
#'   this trial call.
#' @param log_per_environment Log per-environment rows in `phenotype_log`.
#' @param log_aggregate Log aggregate line-mean rows (`environment=0`).
#' @param inherit_genotypes Propagate genotype availability from source cohort(s)
#'   to output cohort when possible.
#' @param silent Suppress informational no-ready messages in legacy `cfg` mode.
#' @param fail_if_no_ready Error on no-ready cohorts in legacy `cfg` mode.
#'
#' @section Environment-control options:
#' To model explicit environment means, set one of `use_env_control = TRUE`,
#' `env_means`, `env_mean_mu`, `env_mean_sd`, or `env_year_sd`.
#' For each call, BPS constructs one latent environment value per location as
#' `env_means + location_deviation + year_effect`, where `location_deviation`
#' is independent across locations with SD `env_mean_sd` and `year_effect` is
#' one common draw shared across all locations with SD `env_year_sd`. Latent
#' environment values are converted to AlphaSimR `p` values with `pnorm()`, so
#' latent value `0` corresponds to the average environment (`p = 0.5`).
#'
#' A practical default is to keep most latent environment values in roughly
#' `[-2, 2]`. In practice this often means long-term means near `0`, with
#' `env_mean_sd` around `0.1` to `0.5` and `env_year_sd` around `0.1` to
#' `0.3`. Larger values push trials into more extreme parts of the GxE
#' reaction norm or create strongly correlated year effects.
#'
#' @section Logging behavior:
#' \describe{
#'   \item{`log_per_environment`}{Log each environment separately (default `TRUE` when environment control is used).}
#'   \item{`log_aggregate`}{Log aggregate line means as environment `0` (default `TRUE`).}
#' }
#'
#' @section Usage modes:
#' Recommended explicit mode:
#' \preformatted{
#' src <- select_latest_available(state, stage = "PYT", stream = "main")
#' sel <- AlphaSimR::selectInd(src$pop, nInd = 50, use = "pheno", simParam = state$sim$SP)
#' state <- run_phenotype_trial(
#'   state = state,
#'   pop = sel,
#'   output_stage = "AYT",
#'   input_cohorts = src$source_ids,
#'   selection_strategy = "Top phenotype from latest PYT",
#'   traits = 1L,
#'   n_loc = 4L,
#'   reps = 2L,
#'   varE = 1.0,
#'   duration_years = 1
#' )
#' }
#'
#' Legacy compatibility mode:
#' \preformatted{
#' state <- run_phenotype_trial(state, cfg = list(
#'   input_stage = "PYT",
#'   output_stage = "AYT",
#'   select_entries_fn = function(state, src, pop_in, cfg) seq_len(min(50, pop_n_ind(pop_in))),
#'   traits = 1L, n_loc = 4L, reps = 2L, varE = 1.0
#' ))
#' }
#'
#' @return Updated program state.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("AlphaSimR", quietly = TRUE)) {
#'   library(AlphaSimR)
#'   h <- quickHaplo(20, 2, 50)
#'   SP <- SimParam$new(h)
#'   SP$addTraitA(10)
#'
#'   state <- bp_init_state(SP = SP, dt = 0.25)
#'   parents <- newPop(h, simParam = SP)
#'   state <- put_stage_pop(state, parents, stage = "F5", ready_in_years = 0)
#'
#'   src <- select_latest_available(state, stage = "F5", stream = "main", combine = TRUE)
#'   state <- run_phenotype_trial(
#'     state = state,
#'     pop = src$pop,
#'     output_stage = "PYT",
#'     input_cohorts = src$source_ids,
#'     selection_strategy = "All entries from latest F5",
#'     traits = 1L,
#'     n_loc = 4L,
#'     reps = 2L,
#'     varE = 1.0,
#'     duration_years = 0.5,
#'     cost_per_plot = 20,
#'     env_means = c(-0.5, 0.0, 0.5, 1.0),
#'     env_mean_sd = 0.25,
#'     env_year_sd = 0.2
#'   )
#' }
#' }
#' @export
run_phenotype_trial <- function(
  state,
  cfg = NULL,
  pop = NULL,
  output_stage = NULL,
  input_cohorts = NULL,
  selection_strategy = NULL,
  traits = 1L,
  n_loc = 1L,
  reps = 1L,
  varE = NULL,
  duration_years = 1,
  stream = "main",
  output_stream = NULL,
  cycle_id = NULL,
  cost_per_plot = 10,
  trial_name = NULL,
  use_env_control = FALSE,
  env_means = NULL,
  env_mean_mu = NULL,
  env_mean_sd = NULL,
  env_year_sd = NULL,
  log_per_environment = TRUE,
  log_aggregate = TRUE,
  inherit_genotypes = TRUE,
  silent = FALSE,
  fail_if_no_ready = FALSE
) {
  # New explicit-pop mode.
  if (!is.null(pop)) {
    if (is.null(output_stage) || !nzchar(as.character(output_stage))) {
      stop("run_phenotype_trial: output_stage is required when pop is supplied.", call. = FALSE)
    }
    if (is.null(varE) || !is.finite(as.numeric(varE))) {
      stop("run_phenotype_trial: varE is required when pop is supplied.", call. = FALSE)
    }

    parse_source_ids <- function(x) {
      if (is.null(x) || length(x) == 0L) return(character(0))
      if (is.list(x) && !is.null(x$source_ids)) return(parse_source_ids(x$source_ids))
      if (is.data.frame(x) && "cohort_id" %in% names(x)) return(parse_source_ids(x$cohort_id))
      out <- as.character(x)
      out <- out[!is.na(out) & nzchar(out) & out != "NA"]
      unique(out)
    }
    source_ids <- parse_source_ids(input_cohorts)

    traits <- as.integer(traits %||% 1L)
    n_loc <- as.integer(n_loc %||% 1L)
    reps <- as.integer(reps %||% 1L)
    stage_name <- as.character(output_stage)
    pop_trial <- pop

    use_env <- isTRUE(use_env_control) ||
      !is.null(env_means) || !is.null(env_year_sd) ||
      !is.null(env_mean_mu) || !is.null(env_mean_sd)

    if (isTRUE(use_env)) {
      env_out <- bp_resolve_trial_env(
        cfg = list(
          env_means = env_means,
          env_mean_mu = env_mean_mu,
          env_mean_sd = env_mean_sd,
          env_year_sd = env_year_sd
        ),
        n_loc = n_loc
      )
      z_env <- env_out$z_env
      p_env <- matrix(NA_real_, nrow = n_loc, ncol = length(traits))

      env_pheno <- vector("list", n_loc)
      for (e in seq_len(n_loc)) {
        p_env[e, ] <- bp_env_p_from_latent(
          simParam = state$sim$SP,
          traits = traits,
          latent_env = z_env[e]
        )
        env_pheno[[e]] <- AlphaSimR::setPheno(
          pop_trial,
          varE = as.numeric(varE),
          reps = reps,
          traits = traits,
          p = p_env[e, ],
          onlyPheno = TRUE,
          simParam = state$sim$SP
        )
        if (is.null(dim(env_pheno[[e]]))) {
          env_pheno[[e]] <- matrix(env_pheno[[e]], ncol = length(traits))
        }
      }
      pheno_mean <- Reduce("+", env_pheno) / n_loc
      if (length(traits) == 1L) {
        pop_trial@pheno[] <- as.numeric(pheno_mean[, 1L])
      } else {
        pop_trial@pheno[,] <- pheno_mean
      }
    } else {
      env_pheno <- vector("list", n_loc)
      for (e in seq_len(n_loc)) {
        env_pheno[[e]] <- AlphaSimR::setPheno(
          pop_trial,
          varE = as.numeric(varE),
          reps = reps,
          traits = traits,
          onlyPheno = TRUE,
          simParam = state$sim$SP
        )
        if (is.null(dim(env_pheno[[e]]))) {
          env_pheno[[e]] <- matrix(env_pheno[[e]], ncol = length(traits))
        }
      }
      pheno_mean <- Reduce("+", env_pheno) / n_loc
      if (length(traits) == 1L) {
        pop_trial@pheno[] <- as.numeric(pheno_mean[, 1L])
      } else {
        pop_trial@pheno[,] <- pheno_mean
      }
      p_env <- rep(NA_real_, n_loc)
    }


    src_cycle <- {
      if (length(source_ids) == 0L) {
        as.character(cycle_id %||% "cycle_1")
      } else {
        rows <- state$cohorts[state$cohorts$cohort_id %in% source_ids, , drop = FALSE]
        cyc <- unique(as.character(rows$cycle_id))
        cyc <- cyc[!is.na(cyc) & nzchar(cyc) & cyc != "NA"]
        if (length(cyc) == 1L) cyc[[1L]] else as.character(cycle_id %||% "cycle_1")
      }
    }

    state <- bp_add_cohort(
      state = state,
      pop = pop_trial,
      stage = stage_name,
      stream = as.character(output_stream %||% stream %||% "main"),
      cycle_id = src_cycle,
      source_cohort_id = if (length(source_ids) == 0L) NA_character_ else paste(source_ids, collapse = ";"),
      selection_strategy = as.character(selection_strategy %||% "unspecified"),
      duration_years = as.numeric(duration_years %||% 1)
    )
    new_cohort_id <- bp_last_cohort_id(state)
    if (isTRUE(inherit_genotypes %||% TRUE)) {
      state <- bp_inherit_genotypes_from_source(
        state = state,
        new_cohort_id = new_cohort_id,
        source_ids = source_ids
      )
    }

    avail_tick <- state$cohorts$available_tick[match(new_cohort_id, state$cohorts$cohort_id)]
    src_label <- if (length(source_ids) == 0L) "unspecified source cohort(s)" else bp_source_labels(state, source_ids, use = "created")
    src_stage <- {
      if (length(source_ids) == 0L) {
        "unspecified"
      } else {
        rows <- state$cohorts[state$cohorts$cohort_id %in% source_ids, , drop = FALSE]
        paste(unique(as.character(rows$stage)), collapse = "+")
      }
    }
    yr_now <- bp_format_year(bp_tick_to_year(state, state$time$tick))
    yr_av <- bp_format_year(bp_tick_to_year(state, avail_tick))
    trait_txt <- paste(traits, collapse = ",")
    n_selected <- pop_n_ind(pop_trial)
    sel_desc <- as.character(selection_strategy %||% "unspecified")

    event_txt <- sprintf(
      "Year %s: Started a %s trial by selecting n=%d from %s by %s. The trial has %d locations with %d rep per location and takes %.2f years to complete and measures traits %s. Will be available Year %s.",
      yr_now, stage_name, n_selected, src_label, sel_desc, n_loc, reps,
      as.numeric(duration_years %||% 1), trait_txt, yr_av
    )
    tpl_txt <- sprintf(
      "Phenotype %s from %s by %s (%d loc x %d rep, traits %s, dur %.2f)",
      stage_name, src_stage, sel_desc, n_loc, reps, trait_txt, as.numeric(duration_years %||% 1)
    )
    state <- bp_log_event(
      state = state,
      fn = "run_phenotype_trial",
      event_type = "phenotyping",
      stage = stage_name,
      source_ids = source_ids,
      output_id = new_cohort_id,
      event_string = event_txt,
      template_string = tpl_txt,
      details = list(
        selection_strategy = sel_desc,
        traits = traits,
        n_loc = n_loc,
        reps = reps,
        varE = as.numeric(varE),
        duration_years = as.numeric(duration_years %||% 1),
        n_selected = n_selected
      )
    )

    if (isTRUE(use_env) && isTRUE(log_per_environment %||% TRUE)) {
      for (e in seq_len(n_loc)) {
        state <- bp_record_pheno(
          state = state,
          cohort_id = new_cohort_id,
          stage = stage_name,
          individual_id = pop_trial@id,
          traits = traits,
          pheno_matrix = env_pheno[[e]],
          available_tick = avail_tick,
          n_loc = n_loc,
          reps = reps,
          environment = e,
          p_value = p_env[e, ]
        )
      }
    }

    if (isTRUE(log_aggregate %||% TRUE)) {
      ph <- pop_trial@pheno
      if (is.null(dim(ph))) ph <- matrix(ph, ncol = length(traits))
      state <- bp_record_pheno(
        state = state,
        cohort_id = new_cohort_id,
        stage = stage_name,
        individual_id = pop_trial@id,
        traits = traits,
        pheno_matrix = ph,
        available_tick = avail_tick,
        n_loc = n_loc,
        reps = reps,
        environment = 0L,
        p_value = if (is.matrix(p_env)) colMeans(p_env) else mean(p_env)
      )
    }

    n_plots <- pop_n_ind(pop_trial) * n_loc * reps
    state <- bp_add_cost(
      state = state,
      stage = stage_name,
      cohort_id = new_cohort_id,
      event = "phenotype_trial",
      unit = "plot",
      n_units = n_plots,
      unit_cost = as.numeric(cost_per_plot %||% 10)
    )
    return(state)
  }

  if (is.null(cfg)) {
    stop("run_phenotype_trial: supply either cfg (legacy mode) or pop + explicit arguments.", call. = FALSE)
  }

  stage_label <- as.character(cfg$input_stage %||% "unknown")
  ready <- bp_get_ready_cohorts(state, stage = cfg$input_stage, stream = cfg$stream %||% NULL)
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "run_phenotype_trial", stage_label)
    return(state)
  }
  ready <- bp_select_source_rows(state, ready, cfg)
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "run_phenotype_trial", stage_label, context = "source selection policy returned no cohorts")
    return(state)
  }

  traits <- as.integer(cfg$traits %||% 1L)
  n_loc <- as.integer(cfg$n_loc %||% 1L)
  reps <- as.integer(cfg$reps %||% 1L)
  use_env_control <- isTRUE(cfg$use_env_control %||% FALSE) ||
    !is.null(cfg$env_means) || !is.null(cfg$env_year_sd) ||
    !is.null(cfg$env_mean_mu) || !is.null(cfg$env_mean_sd)

  for (i in seq_len(nrow(ready))) {
    src <- ready[i, , drop = FALSE]
    pop_in <- state$pops[[src$cohort_id]]

    idx <- if (is.function(cfg$select_entries_fn)) {
      as.integer(cfg$select_entries_fn(state, src, pop_in, cfg))
    } else {
      seq_len(pop_n_ind(pop_in))
    }
    if (length(idx) == 0L) next

    pop_trial <- pop_subset(pop_in, idx)
    if (isTRUE(use_env_control)) {
      env_out <- bp_resolve_trial_env(cfg, n_loc = n_loc)
      z_env <- env_out$z_env
      p_env <- matrix(NA_real_, nrow = n_loc, ncol = length(traits))

      env_pheno <- vector("list", n_loc)
      for (e in seq_len(n_loc)) {
        p_env[e, ] <- bp_env_p_from_latent(
          simParam = state$sim$SP,
          traits = traits,
          latent_env = z_env[e]
        )
        env_pheno[[e]] <- AlphaSimR::setPheno(
          pop_trial,
          varE = cfg$varE,
          reps = reps,
          traits = traits,
          p = p_env[e, ],
          onlyPheno = TRUE,
          simParam = state$sim$SP
        )
        if (is.null(dim(env_pheno[[e]]))) {
          env_pheno[[e]] <- matrix(env_pheno[[e]], ncol = length(traits))
        }
      }

      pheno_mean <- Reduce("+", env_pheno) / n_loc
      if (length(traits) == 1L) {
        pop_trial@pheno[] <- as.numeric(pheno_mean[, 1])
      } else {
        pop_trial@pheno[,] <- pheno_mean
      }
    } else {
      reps_eff <- as.integer(max(1L, reps * n_loc))
      pop_trial <- AlphaSimR::setPheno(
        pop_trial,
        varE = cfg$varE,
        reps = reps_eff,
        traits = traits,
        simParam = state$sim$SP
      )
      env_pheno <- NULL
      p_env <- rep(NA_real_, n_loc)
    }

    state <- bp_add_cohort(
      state = state,
      pop = pop_trial,
      stage = cfg$output_stage %||% cfg$trial_name,
      stream = cfg$output_stream %||% src$stream,
      cycle_id = src$cycle_id,
      source_cohort_id = src$cohort_id,
      duration_years = cfg$duration_years %||% 1
    )
    new_cohort_id <- bp_last_cohort_id(state)
    if (isTRUE(cfg$inherit_genotypes %||% TRUE)) {
      state <- bp_inherit_genotypes_from_source(
        state = state,
        new_cohort_id = new_cohort_id,
        source_ids = src$cohort_id
      )
    }

    avail_tick <- state$cohorts$available_tick[match(new_cohort_id, state$cohorts$cohort_id)]
    stage_name <- cfg$output_stage %||% cfg$trial_name
    n_selected <- pop_n_ind(pop_trial)
    sel_desc <- as.character(cfg$selection_strategy %||% if (is.function(cfg$select_entries_fn)) "select_entries_fn" else "all entries")
    src_label <- bp_source_labels(state, src$cohort_id, use = "created")
    yr_now <- bp_format_year(bp_tick_to_year(state, state$time$tick))
    yr_av <- bp_format_year(bp_tick_to_year(state, avail_tick))
    trait_txt <- paste(traits, collapse = ",")
    event_txt <- sprintf(
      "Year %s: Started a %s trial by selecting n=%d from %s by %s. The trial has %d locations with %d rep per location and takes %.2f years to complete and measures traits %s. Will be available Year %s.",
      yr_now,
      stage_name,
      n_selected,
      src_label,
      sel_desc,
      n_loc,
      reps,
      as.numeric(cfg$duration_years %||% 1),
      trait_txt,
      yr_av
    )
    tpl_txt <- sprintf(
      "Phenotype %s from %s by %s (%d loc x %d rep, traits %s, dur %.2f)",
      stage_name,
      as.character(src$stage),
      sel_desc,
      n_loc,
      reps,
      trait_txt,
      as.numeric(cfg$duration_years %||% 1)
    )
    state <- bp_log_event(
      state = state,
      fn = "run_phenotype_trial",
      event_type = "phenotyping",
      stage = stage_name,
      source_ids = src$cohort_id,
      output_id = new_cohort_id,
      event_string = event_txt,
      template_string = tpl_txt,
      details = list(
        selection_strategy = sel_desc,
        traits = traits,
        n_loc = n_loc,
        reps = reps,
        varE = as.numeric(cfg$varE),
        duration_years = as.numeric(cfg$duration_years %||% 1),
        n_selected = n_selected
      )
    )

    if (isTRUE(use_env_control) && isTRUE(cfg$log_per_environment %||% TRUE)) {
      for (e in seq_len(n_loc)) {
        state <- bp_record_pheno(
          state = state,
          cohort_id = new_cohort_id,
          stage = stage_name,
          individual_id = pop_trial@id,
          traits = traits,
          pheno_matrix = env_pheno[[e]],
          available_tick = avail_tick,
          n_loc = n_loc,
          reps = reps,
          environment = e,
          p_value = p_env[e, ]
        )
      }
    }

    if (isTRUE(cfg$log_aggregate %||% TRUE)) {
      ph <- pop_trial@pheno
      if (is.null(dim(ph))) ph <- matrix(ph, ncol = length(traits))
      state <- bp_record_pheno(
        state = state,
        cohort_id = new_cohort_id,
        stage = stage_name,
        individual_id = pop_trial@id,
        traits = traits,
        pheno_matrix = ph,
        available_tick = avail_tick,
        n_loc = n_loc,
        reps = reps,
        environment = 0L,
        p_value = if (is.matrix(p_env)) colMeans(p_env) else mean(p_env)
      )
    }

    n_plots <- pop_n_ind(pop_trial) * n_loc * reps
    state <- bp_add_cost(
      state,
      stage = stage_name,
      cohort_id = new_cohort_id,
      event = "phenotype_trial",
      unit = "plot",
      n_units = n_plots,
      unit_cost = cfg$cost_per_plot %||% 10
    )

    if (isTRUE(cfg$consume_input %||% TRUE)) {
      state <- bp_close_cohort(state, src$cohort_id)
    }
  }

  state
}

#' Run Genotyping
#'
#' Schedule/log genotyping events by cohort and chip, with genotyping costs.
#'
#' @param state Program state.
#' @param cfg Genotyping configuration list.
#'
#' @return Updated program state.
#' @export
run_genotyping <- function(state, cfg) {
  stage_label <- as.character(cfg$input_stage %||% "unknown")
  if (isTRUE(cfg$include_not_ready %||% FALSE)) {
    ready <- state$cohorts
    if (!is.null(cfg$input_stage)) {
      ready <- ready[ready$stage %in% cfg$input_stage, , drop = FALSE]
    }
    if (!is.null(cfg$stream)) {
      ready <- ready[ready$stream %in% cfg$stream, , drop = FALSE]
    }
    ready <- ready[ready$active, , drop = FALSE]
    ready <- ready[ready$created_tick <= as.integer(state$time$tick), , drop = FALSE]
  } else {
    ready <- bp_get_ready_cohorts(state, stage = cfg$input_stage, stream = cfg$stream %||% NULL)
  }
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "run_genotyping", stage_label)
    return(state)
  }
  ready <- bp_select_source_rows(state, ready, cfg)
  if (nrow(ready) == 0L) {
    bp_handle_no_ready(cfg, "run_genotyping", stage_label, context = "source selection policy returned no cohorts")
    return(state)
  }

  chip <- cfg$chip %||% state$sim$default_chip
  ckey <- chip_key(chip)
  dur_ticks <- years_to_ticks(state$time$dt, cfg$duration_years %||% 1)
  force <- isTRUE(cfg$force %||% FALSE)

  for (i in seq_len(nrow(ready))) {
    src <- ready[i, , drop = FALSE]
    n_ind <- src$n_ind
    already <- state$genotype_log$cohort_id == src$cohort_id & state$genotype_log$chip == ckey
    if (!force && any(already, na.rm = TRUE)) {
      next
    }

    row <- data.frame(
      cohort_id = src$cohort_id,
      chip = ckey,
      started_tick = as.integer(state$time$tick),
      done_tick = as.integer(state$time$tick + dur_ticks),
      available_tick = as.integer(state$time$tick + dur_ticks),
      n_ind = as.integer(n_ind),
      stringsAsFactors = FALSE
    )
    state$genotype_log <- rbind(state$genotype_log, row)

    yr_now <- bp_format_year(bp_tick_to_year(state, state$time$tick))
    yr_av <- bp_format_year(bp_tick_to_year(state, row$available_tick[[1]]))
    src_label <- bp_source_labels(state, src$cohort_id, use = "created")
    evt <- sprintf(
      "Year %s: Genotyped %s population using snpChip=%s (n=%d). Will be available Year %s.",
      yr_now, src_label, ckey, as.integer(n_ind), yr_av
    )
    tpl <- sprintf("Genotype stage=%s chip=%s", as.character(src$stage), ckey)
    state <- bp_log_event(
      state = state,
      fn = "run_genotyping",
      event_type = "genotyping",
      stage = as.character(src$stage),
      source_ids = src$cohort_id,
      output_id = src$cohort_id,
      event_string = evt,
      template_string = tpl,
      details = list(
        chip = ckey,
        n_ind = as.integer(n_ind),
        duration_years = as.numeric(cfg$duration_years %||% 1)
      )
    )

    state <- bp_add_cost(
      state,
      stage = src$stage,
      cohort_id = src$cohort_id,
      event = "genotyping",
      unit = "sample",
      n_units = n_ind,
      unit_cost = cfg$cost_per_sample %||% 15
    )
  }

  state <- bp_refresh_genotyped_flags(state)
  state
}

# Identify eligible training cohorts from a recent window with required genotype data.
bp_get_training_cohorts <- function(state, cfg) {
  ready <- bp_get_ready_cohorts(state, stage = cfg$from_stage, stream = cfg$stream %||% NULL)
  if (nrow(ready) == 0L) return(ready)

  lookback_ticks <- years_to_ticks(state$time$dt, cfg$lookback_years %||% 3)
  min_tick <- as.integer(state$time$tick - lookback_ticks)
  ready <- ready[ready$available_tick >= min_tick, , drop = FALSE]

  ckey <- chip_key(cfg$chip %||% state$sim$default_chip)
  if (nrow(ready) == 0L) return(ready)

  has_chip <- vapply(ready$cohort_id, function(cid) {
    any(state$genotype_log$cohort_id == cid & state$genotype_log$chip == ckey & state$genotype_log$available_tick <= state$time$tick)
  }, logical(1))
  ready <- ready[has_chip, , drop = FALSE]
  if (nrow(ready) == 0L) return(ready)

  training_policy <- as.character(cfg$training_policy %||% "all_ready")
  if (training_policy == "all_ready") {
    return(ready)
  }
  cfg2 <- cfg
  cfg2$input_policy <- training_policy
  bp_select_source_rows(state, ready, cfg2)
}

# Execute a user-provided hook with contextual error messages.
bp_call_user_fn <- function(fn, args, fn_label, stage_label) {
  tryCatch(
    do.call(fn, args),
    error = function(e) {
      stop(sprintf("%s failed for stage '%s': %s", fn_label, stage_label, conditionMessage(e)), call. = FALSE)
    }
  )
}

# Ensure one or more cohorts have available genotype records for the requested chip.
bp_assert_genotyped_cohorts <- function(state, cohort_ids, chip, stage_label = "unknown", context = "prediction") {
  ids <- unique(as.character(cohort_ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0L) {
    stop(sprintf("%s for stage '%s' requires cohort_ids for genotype validation", context, stage_label), call. = FALSE)
  }

  ckey <- chip_key(chip)
  bad <- vapply(ids, function(cid) {
    !any(
      state$genotype_log$cohort_id == cid &
        state$genotype_log$chip == ckey &
        state$genotype_log$available_tick <= state$time$tick
    )
  }, logical(1))

  if (any(bad)) {
    miss <- paste(ids[bad], collapse = ", ")
    stop(
      sprintf(
        "%s for stage '%s' requires genotyping for chip '%s'; missing cohorts: %s",
        context, stage_label, ckey, miss
      ),
      call. = FALSE
    )
  }
}

#' Predict EBV on a Pop
#'
#' Predict EBVs for a target pop using either a user function or AlphaSimR.
#'
#' @param pop Target pop.
#' @param model_entry Model entry from `state$gs_models`.
#' @param state Program state.
#' @param cfg Prediction configuration list.
#' @param stage_label Stage label used in errors.
#'
#' @section Key `cfg` fields:
#' \describe{
#'   \item{`cohort_ids`}{Cohort id(s) used for genotype-chip validation. Required when `require_genotyped = TRUE`.}
#'   \item{`chip`}{Optional chip override for validation/prediction compatibility.}
#'   \item{`require_genotyped`}{Whether to enforce genotype availability checks. Default `TRUE`.}
#'   \item{`predict_ebv_fn`}{Optional custom predictor function. Signature:
#'   `(target_pop, model_obj, state, cfg, model_entry)`. Must return numeric
#'   vector (`nInd`) or numeric matrix (`nInd x nTraits`).}
#' }
#'
#' @details
#' If `predict_ebv_fn` is not provided, this function calls
#' `AlphaSimR::setEBV(pop, solution = model_entry$model, simParam = state$sim$SP)`.
#' For custom predictors, returned values are written directly into `pop@ebv`
#' (single- or multi-trait).
#'
#' @return Pop with updated `@ebv`.
#' @export
predict_ebv_pop <- function(pop, model_entry, state, cfg, stage_label = "unknown") {
  require_genotyped <- isTRUE(cfg$require_genotyped %||% TRUE)
  if (require_genotyped) {
    chip_for_pred <- cfg$chip %||% model_entry$chip %||% state$sim$default_chip
    bp_assert_genotyped_cohorts(
      state = state,
      cohort_ids = cfg$cohort_ids %||% NULL,
      chip = chip_for_pred,
      stage_label = stage_label,
      context = "EBV prediction"
    )
  }

  predict_fn <- cfg$predict_ebv_fn %||% model_entry$predict_ebv_fn %||% NULL

  if (is.function(predict_fn)) {
    pred <- bp_call_user_fn(
      predict_fn,
      list(target_pop = pop, model_obj = model_entry$model, state = state, cfg = cfg, model_entry = model_entry),
      fn_label = "predict_ebv_fn",
      stage_label = stage_label
    )
    n <- pop_n_ind(pop)

    if (is.null(dim(pred))) {
      pred_vec <- as.numeric(pred)
      if (length(pred_vec) != n) {
        stop(sprintf("predict_ebv_fn returned vector length %d; expected %d", length(pred_vec), n), call. = FALSE)
      }
      if (anyNA(pred_vec)) {
        stop("predict_ebv_fn returned NA values", call. = FALSE)
      }
      pop@ebv <- matrix(pred_vec, ncol = 1)
    } else {
      pred_mat <- as.matrix(pred)
      if (nrow(pred_mat) != n) {
        stop(sprintf("predict_ebv_fn returned matrix with %d rows; expected %d", nrow(pred_mat), n), call. = FALSE)
      }
      if (ncol(pred_mat) < 1L) {
        stop("predict_ebv_fn returned matrix with zero columns", call. = FALSE)
      }
      if (anyNA(pred_mat)) {
        stop("predict_ebv_fn returned NA values", call. = FALSE)
      }
      storage.mode(pred_mat) <- "double"
      pop@ebv <- pred_mat
    }
    return(pop)
  }

  AlphaSimR::setEBV(pop, solution = model_entry$model, simParam = state$sim$SP)
}

#' Predict EBV for Cohorts in State
#'
#' Predict EBVs for one or more cohorts and persist them back into `state$pops`.
#'
#' @param state Program state.
#' @param cfg Prediction configuration list.
#'
#' @section Required `cfg` fields:
#' \describe{
#'   \item{`model_id`}{Model id from `state$gs_models`. If omitted, latest model is used.}
#' }
#'
#' @section Common optional `cfg` fields:
#' \describe{
#'   \item{`cohort_ids`}{Explicit cohort id(s) to score.}
#'   \item{`input_stage`}{Stage used when selecting cohorts via readiness/policy.}
#'   \item{`stream`}{Optional stream filter.}
#'   \item{`input_policy`}{Source selection policy when using `input_stage`.}
#'   \item{`chip`}{Optional chip override for validation/prediction compatibility.}
#'   \item{`require_genotyped`}{Whether to enforce genotype availability checks. Default `TRUE`.}
#'   \item{`predict_ebv_fn`}{Optional custom predictor function.}
#' }
#'
#' @return Updated program state with scored populations persisted.
#' @export
run_predict_ebv <- function(state, cfg) {
  model_id <- as.character(cfg$model_id %||% bp_latest_model_id(state))
  if (!nzchar(model_id) || is.null(state$gs_models[[model_id]])) {
    stop(sprintf("run_predict_ebv: model '%s' not found", model_id), call. = FALSE)
  }
  model_entry <- state$gs_models[[model_id]]

  if (!is.null(cfg$cohort_ids)) {
    ids <- unique(as.character(cfg$cohort_ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (length(ids) == 0L) return(state)
    keep <- state$cohorts$cohort_id %in% ids
    rows <- state$cohorts[keep, , drop = FALSE]
    if (nrow(rows) == 0L) return(state)
  } else {
    stage_label <- as.character(cfg$input_stage %||% "unknown")
    rows <- bp_get_ready_cohorts(state, stage = cfg$input_stage, stream = cfg$stream %||% NULL)
    if (nrow(rows) == 0L) {
      bp_handle_no_ready(cfg, "run_predict_ebv", stage_label)
      return(state)
    }
    rows <- bp_select_source_rows(state, rows, cfg)
    if (nrow(rows) == 0L) {
      bp_handle_no_ready(cfg, "run_predict_ebv", stage_label, context = "source selection policy returned no cohorts")
      return(state)
    }
  }

  for (i in seq_len(nrow(rows))) {
    cid <- as.character(rows$cohort_id[i])
    pop <- state$pops[[cid]]
    cfg_pred <- utils::modifyList(as.list(cfg), list(cohort_ids = cid))
    pop2 <- predict_ebv_pop(
      pop = pop,
      model_entry = model_entry,
      state = state,
      cfg = cfg_pred,
      stage_label = as.character(rows$stage[i] %||% "unknown")
    )
    state$pops[[cid]] <- pop2
  }

  state
}

#' Train GP Model
#'
#' Train and store a genomic prediction model from eligible cohorts.
#'
#' @param state Program state.
#' @param cfg Training configuration list.
#'
#' @section Required `cfg` fields:
#' \describe{
#'   \item{`from_stage`}{Stage used as training source cohorts.}
#' }
#'
#' @section Common optional `cfg` fields:
#' \describe{
#'   \item{`chip`}{Genotyping chip key/index. Default `state$sim$default_chip`.}
#'   \item{`trait`}{Trait index for default RRBLUP model.}
#'   \item{`lookback_years`}{Training cohort lookback window.}
#'   \item{`training_policy`}{Subset policy over eligible training cohorts.}
#'   \item{`model_id`}{Stored model id. Auto-generated if omitted.}
#'   \item{`train_model_fn`}{Custom trainer `(train_pop, state, cfg) -> model_object`.}
#'   \item{`predict_ebv_fn`}{Optional custom predictor stored alongside model.}
#' }
#'
#' @details
#' By default, this function trains `AlphaSimR::RRBLUP(...)` on merged training
#' cohorts filtered to the requested chip. The stored model entry includes model
#' object, trained tick, chip key, trait index, and source cohorts.
#'
#' @return Updated program state.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("AlphaSimR", quietly = TRUE)) {
#'   library(AlphaSimR)
#'   h <- quickHaplo(20, 2, 50)
#'   SP <- SimParam$new(h)
#'   SP$addTraitA(10)
#'
#'   state <- bp_init_state(SP = SP, dt = 1)
#'   pop <- newPop(h, simParam = SP)
#'   pop <- setPheno(pop, varE = 1, reps = 1, traits = 1, simParam = SP)
#'   state <- put_stage_pop(state, pop, stage = "PYT", ready_in_years = 0)
#'
#'   state <- run_genotyping(state, list(input_stage = "PYT", chip = 1L, duration_years = 0))
#'
#'   state <- run_train_gp_model(
#'     state,
#'     list(
#'       from_stage = "PYT",
#'       chip = 1L,
#'       trait = 1L,
#'       lookback_years = 3,
#'       model_id = "rrblup_pyt"
#'     )
#'   )
#'   names(state$gs_models)
#' }
#' }
#' @export
run_train_gp_model <- function(state, cfg) {
  train_cohorts <- bp_get_training_cohorts(state, cfg)
  if (nrow(train_cohorts) == 0L) {
    bp_handle_no_ready(cfg, "run_train_gp_model", as.character(cfg$from_stage %||% "unknown"), context = "no eligible training cohorts")
    return(state)
  }

  pops <- lapply(train_cohorts$cohort_id, function(cid) state$pops[[cid]])
  train_pop <- merge_pops(pops)

  chip_raw <- cfg$chip %||% state$sim$default_chip
  ckey <- chip_key(chip_raw)
  cidx <- chip_index(state, chip_raw)
  trait <- as.integer(cfg$trait %||% 1L)
  stage_label <- as.character(cfg$from_stage %||% "unknown")

  if (is.function(cfg$train_model_fn)) {
    model <- bp_call_user_fn(
      cfg$train_model_fn,
      list(train_pop = train_pop, state = state, cfg = cfg),
      fn_label = "train_model_fn",
      stage_label = stage_label
    )
  } else {
    model <- AlphaSimR::RRBLUP(train_pop, traits = trait, use = "pheno", snpChip = cidx, simParam = state$sim$SP)
  }
  if (is.null(model)) {
    stop("run_train_gp_model: model object is NULL", call. = FALSE)
  }

  state$counters$model <- state$counters$model + 1L
  model_id <- as.character(cfg$model_id %||% sprintf("gp_%03d", state$counters$model))

  state$gs_models[[model_id]] <- list(
    model = model,
    model_id = model_id,
    model_name = as.character(cfg$model_name %||% if (is.function(cfg$train_model_fn)) "custom_train_model_fn" else "AlphaSimR::RRBLUP"),
    predict_ebv_fn = cfg$predict_ebv_fn %||% NULL,
    trained_tick = as.integer(state$time$tick),
    chip = ckey,
    trait = trait,
    source_cohorts = train_cohorts$cohort_id
  )

  yr_now <- bp_format_year(bp_tick_to_year(state, state$time$tick))
  src_label <- bp_source_labels(state, train_cohorts$cohort_id, use = "created")
  n_total <- pop_n_ind(train_pop)
  evt <- sprintf(
    "Year %s: Trained GS model %s on %s (n_total=%d, trait=%d, chip=%s).",
    yr_now,
    state$gs_models[[model_id]]$model_name,
    src_label,
    as.integer(n_total),
    as.integer(trait),
    ckey
  )
  tpl <- sprintf(
    "Train model=%s from stage=%s (trait=%d, chip=%s)",
    state$gs_models[[model_id]]$model_name,
    as.character(cfg$from_stage %||% "unknown"),
    as.integer(trait),
    ckey
  )
  state <- bp_log_event(
    state = state,
    fn = "run_train_gp_model",
    event_type = "train_model",
    stage = as.character(cfg$from_stage %||% "unknown"),
    source_ids = train_cohorts$cohort_id,
    output_id = model_id,
    event_string = evt,
    template_string = tpl,
    details = list(
      model_id = model_id,
      model_name = state$gs_models[[model_id]]$model_name,
      n_total = as.integer(n_total),
      trait = as.integer(trait),
      chip = ckey
    )
  )

  state
}

#' Get Latest Model ID
#'
#' Return the id of the most recently trained model in `state$gs_models`.
#'
#' @param state Program state.
#'
#' @return Character scalar model id, or `NULL`.
#' @export
bp_latest_model_id <- function(state) {
  if (length(state$gs_models) == 0L) return(NULL)
  ord <- order(vapply(state$gs_models, function(x) x$trained_tick, numeric(1)), decreasing = TRUE)
  names(state$gs_models)[ord][1]
}

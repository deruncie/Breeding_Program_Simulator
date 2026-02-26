# S3 print method and helpers for bp_state.

#' State Overview
#'
#' Build a compact overview of a `bp_state` object.
#'
#' @param state A `bp_state` object.
#'
#' @return Named list with summary statistics.
#' @export
bp_state_overview <- function(state) {
  stopifnot(inherits(state, "bp_state"))

  cohorts <- state$cohorts
  total_cohorts <- if (is.null(cohorts)) 0L else nrow(cohorts)
  active_cohorts <- if (total_cohorts == 0L) 0L else sum(cohorts$active, na.rm = TRUE)
  ready_cohorts <- if (total_cohorts == 0L) {
    0L
  } else {
    sum(cohorts$available_tick <= as.integer(state$time$tick) & cohorts$active, na.rm = TRUE)
  }

  list(
    tick = as.integer(state$time$tick),
    dt = as.numeric(state$time$dt),
    year = as.numeric(bp_tick_to_year(state, state$time$tick)),
    total_cohorts = as.integer(total_cohorts),
    active_cohorts = as.integer(active_cohorts),
    ready_cohorts = as.integer(ready_cohorts),
    n_models = as.integer(length(state$gs_models)),
    n_varieties = as.integer(nrow(state$outputs$varieties %||% data.frame())),
    n_cost_events = as.integer(nrow(state$cost_log %||% data.frame())),
    n_pheno_rows = as.integer(nrow(state$phenotype_log %||% data.frame())),
    n_geno_rows = as.integer(nrow(state$genotype_log %||% data.frame())),
    n_event_rows = as.integer(nrow(state$event_log %||% data.frame()))
  )
}

#' Print `bp_state`
#'
#' Pretty-print a compact summary for a breeding-program simulation state.
#'
#' @param x A `bp_state` object.
#' @param ... Optional args:
#'   `max_stages` (default `8`) and `recent_events` (default `5`).
#'
#' @return Invisibly returns `x`.
#' @export
print.bp_state <- function(x, ...) {
  stopifnot(inherits(x, "bp_state"))
  dots <- list(...)
  max_stages <- as.integer(dots$max_stages %||% 8L)
  recent_events <- as.integer(dots$recent_events %||% 5L)

  ov <- bp_state_overview(x)
  cat("Breeding Program State\n")
  cat(
    sprintf(
      "Time: tick=%d, dt=%.3f years, year=%.2f\n",
      ov$tick, ov$dt, ov$year
    )
  )
  cat(
    sprintf(
      "Cohorts: total=%d, active=%d, ready=%d\n",
      ov$total_cohorts, ov$active_cohorts, ov$ready_cohorts
    )
  )
  cat(
    sprintf(
      "Models: %d | Varieties: %d\n",
      ov$n_models, ov$n_varieties
    )
  )
  cat(
    sprintf(
      "Logs: events=%d, pheno_rows=%d, geno_rows=%d, cost_events=%d\n",
      ov$n_event_rows, ov$n_pheno_rows, ov$n_geno_rows, ov$n_cost_events
    )
  )

  cohorts <- x$cohorts
  if (!is.null(cohorts) && nrow(cohorts) > 0L) {
    act <- cohorts[cohorts$active, c("stage"), drop = FALSE]
    if (nrow(act) > 0L) {
      tab <- sort(table(act$stage), decreasing = TRUE)
      tab <- tab[seq_len(min(length(tab), max_stages))]
      cat("Active by stage:\n")
      for (nm in names(tab)) {
        cat(sprintf("  - %s: %d\n", nm, as.integer(tab[[nm]])))
      }
    }
  }

  ev <- x$event_log
  if (!is.null(ev) && nrow(ev) > 0L && recent_events > 0L) {
    cat("Recent events:\n")
    k <- min(nrow(ev), recent_events)
    d <- ev[(nrow(ev) - k + 1L):nrow(ev), , drop = FALSE]
    for (i in seq_len(nrow(d))) {
      cat(sprintf("  - %s\n", as.character(d$event_string[[i]])))
    }
  }

  invisible(x)
}

# Simple plotting helpers for an existing simulation `state`.
#
# Usage in an interactive R session:
#   source("R/operators.R")
#   source("R/readable_api.R")
#   source("R/monitoring.R")
#   source("examples/plot_state_metrics_simple.R")
#   # after you create `state`:
  m <- bp_collect_metrics(state)
  bp_plot_mean_gv_origin(m)
  bp_plot_mean_gv_available(m)
  bp_plot_metric_available(m, metric = "mean_gv")
  bp_plot_metric_origin(m, metric = "mean_gv")

library(ggplot2)

# Extract standard cohort metrics from a state object.
bp_collect_metrics <- function(
  state,
  stages = c("PD_DH_INPUT", "PYT", "AYT", "EYT", "Variety"),
  trait = 1L,
  origin_stage = "PD_DH_INPUT",
  include_inactive = TRUE
) {
  ticks_per_year <- as.integer(round(1 / state$time$dt))
  bp_extract_cohort_metrics(
    state = state,
    stages = stages,
    trait = trait,
    origin_stage = origin_stage,
    include_inactive = include_inactive,
    ticks_per_year = ticks_per_year
  )
}

# Plot any metric summarized by origin year.
bp_plot_metric_origin <- function(metrics_df, metric = "mean_gv") {
  df <- bp_summarize_metric_by_year(
    metrics_df,
    metric = metric,
    year_col = "origin_year"
  )
  ggplot(df, aes(x = year, y = value, color = stage, group = stage)) +
    geom_line() +
    geom_point() +
    labs(
      title = paste(metric, "by Origin Year"),
      x = "Origin Year",
      y = metric
    )
}

# Plot any metric summarized by available year.
bp_plot_metric_available <- function(metrics_df, metric = "mean_gv") {
  df <- bp_summarize_metric_by_year(
    metrics_df,
    metric = metric,
    year_col = "available_year"
  )
  ggplot(df, aes(x = year, y = value, color = stage, group = stage)) +
    geom_line() +
    geom_point() +
    labs(
      title = paste(metric, "by Available Year"),
      x = "Available Year",
      y = metric
    )
}

# Convenience wrappers for the most common debugging plots.
bp_plot_mean_gv_origin <- function(metrics_df) {
  bp_plot_metric_origin(metrics_df, metric = "mean_gv")
}

bp_plot_mean_gv_available <- function(metrics_df) {
  bp_plot_metric_available(metrics_df, metric = "mean_gv")
}

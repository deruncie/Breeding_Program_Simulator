library(shiny)
library(ggplot2)

find_repo_root <- function(start = getwd(), max_up = 5L) {
  p <- normalizePath(start, mustWork = TRUE)
  for (i in 0:max_up) {
    cand <- if (i == 0L) p else normalizePath(file.path(p, paste(rep("..", i), collapse = "/")), mustWork = TRUE)
    if (file.exists(file.path(cand, "R", "readable_api.R")) &&
        file.exists(file.path(cand, "examples", "AdvanceYear_GSTP_readable_structured.R"))) {
      return(cand)
    }
  }
  stop("Could not locate repository root from current working directory", call. = FALSE)
}

repo_root <- find_repo_root()
source(file.path(repo_root, "R", "operators.R"))
source(file.path(repo_root, "R", "readable_api.R"))
source(file.path(repo_root, "R", "readable_wrappers.R"))
source(file.path(repo_root, "R", "monitoring.R"))

# Load the GSTP example functions without executing its top-level run block.
ex_env <- new.env(parent = globalenv())
old_wd <- getwd()
setwd(repo_root)
tryCatch(
  source(file.path(repo_root, "examples", "AdvanceYear_GSTP_readable_structured.R"), local = ex_env),
  finally = setwd(old_wd)
)
run_gstp_loop_demo <- ex_env$run_gstp_loop_demo

metric_choices <- c("mean_gv", "var_gv", "max_gv", "cor_ebv_gv", "h2", "H2")

ui <- fluidPage(
  titlePanel("Breeding Program Monitoring Dashboard"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n_years", "Years", value = 14, min = 4, max = 40, step = 1),
      numericInput("trait", "Trait", value = 1, min = 1, step = 1),
      checkboxInput("include_inactive", "Include inactive cohorts", value = TRUE),
      checkboxInput("make_plots", "Run scheme plotting in backend", value = FALSE),
      actionButton("run_btn", "Run Simulation", class = "btn-primary"),
      hr(),
      checkboxGroupInput(
        "stages",
        "Stages",
        choices = c("DH_PIPE", "PYT", "AYT", "EYT", "Variety"),
        selected = c("DH_PIPE", "PYT", "AYT", "EYT", "Variety")
      ),
      selectInput("metric", "Metric", choices = metric_choices, selected = "mean_gv")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Selected Metric",
          fluidRow(
            column(6, plotOutput("plot_origin", height = "360px")),
            column(6, plotOutput("plot_available", height = "360px"))
          ),
          fluidRow(
            column(6, h4("Origin-Year Summary"), tableOutput("tbl_origin")),
            column(6, h4("Available-Year Summary"), tableOutput("tbl_available"))
          )
        ),
        tabPanel(
          "All Metrics",
          uiOutput("all_metric_plots")
        ),
        tabPanel(
          "Cohort Metrics",
          h4("Per-cohort metrics"),
          tableOutput("metrics_table")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  sim_state <- eventReactive(input$run_btn, {
    run_gstp_loop_demo(n_years = as.integer(input$n_years), make_plots = isTRUE(input$make_plots))
  }, ignoreNULL = FALSE)

  metrics_df <- reactive({
    st <- sim_state()
    ticks_per_year <- as.integer(round(1 / st$time$dt))

    bp_extract_cohort_metrics(
      state = st,
      stages = input$stages,
      trait = as.integer(input$trait),
      origin_stage = "DH_PIPE",
      include_inactive = isTRUE(input$include_inactive),
      ticks_per_year = ticks_per_year
    )
  })

  origin_summary <- reactive({
    bp_summarize_metric_by_year(
      metrics_df(),
      metric = input$metric,
      year_col = "origin_year"
    )
  })

  available_summary <- reactive({
    bp_summarize_metric_by_year(
      metrics_df(),
      metric = input$metric,
      year_col = "available_year"
    )
  })

  output$plot_origin <- renderPlot({
    df <- origin_summary()
    ggplot(df, aes(x = year, y = value, color = stage, group = stage)) +
      geom_line() +
      geom_point() +
      labs(title = paste(input$metric, "by Origin Year"), x = "Origin Year", y = input$metric)
  })

  output$plot_available <- renderPlot({
    df <- available_summary()
    ggplot(df, aes(x = year, y = value, color = stage, group = stage)) +
      geom_line() +
      geom_point() +
      labs(title = paste(input$metric, "by Available Year"), x = "Available Year", y = input$metric)
  })

  output$tbl_origin <- renderTable({
    origin_summary()
  }, digits = 4)

  output$tbl_available <- renderTable({
    available_summary()
  }, digits = 4)

  output$metrics_table <- renderTable({
    df <- metrics_df()
    keep <- c(
      "cohort_id", "stage", "origin_cohort_id", "origin_year", "available_year",
      "mean_gv", "var_gv", "max_gv", "cor_ebv_gv", "h2", "H2"
    )
    df[order(df$available_year, df$stage), keep, drop = FALSE]
  }, digits = 4)

  output$all_metric_plots <- renderUI({
    tagList(lapply(metric_choices, function(m) {
      fluidRow(
        column(6, plotOutput(outputId = paste0("plot_origin_", m), height = "260px")),
        column(6, plotOutput(outputId = paste0("plot_available_", m), height = "260px"))
      )
    }))
  })

  for (m in metric_choices) {
    local({
      mm <- m
      output[[paste0("plot_origin_", mm)]] <- renderPlot({
        df <- bp_summarize_metric_by_year(metrics_df(), metric = mm, year_col = "origin_year")
        ggplot(df, aes(x = year, y = value, color = stage, group = stage)) +
          geom_line() +
          geom_point() +
          labs(title = paste(mm, "by Origin Year"), x = "Origin Year", y = mm)
      })

      output[[paste0("plot_available_", mm)]] <- renderPlot({
        df <- bp_summarize_metric_by_year(metrics_df(), metric = mm, year_col = "available_year")
        ggplot(df, aes(x = year, y = value, color = stage, group = stage)) +
          geom_line() +
          geom_point() +
          labs(title = paste(mm, "by Available Year"), x = "Available Year", y = mm)
      })
    })
  }
}

shinyApp(ui, server)

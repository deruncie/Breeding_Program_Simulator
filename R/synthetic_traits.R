# Synthetic-trait definitions and population-local score storage.

#' Define a Synthetic Trait
#'
#' @param name Unique synthetic-trait name.
#' @param traits Biological trait indices or names consumed by `fun`.
#' @param fun Function receiving the selected trait matrix as its first
#'   argument and returning one numeric value per individual.
#' @param args Optional named list of additional arguments passed to `fun`.
#' @param linear Whether `fun` is linear. Linear synthetic GVs are evaluated
#'   exactly; nonlinear traits use Monte Carlo integration.
#' @param score_level Phenotype scoring level. Currently `"aggregate"` is the
#'   standard supported level.
#' @param missing_component How to handle a requested component trait with no
#'   values. `"propagate"` leaves missing values as `NA`. `"zero_if_all_missing"`
#'   replaces any component column that is entirely `NA` with zero before
#'   evaluating `fun`; this is useful for index traits where unobserved or
#'   unpredicted components should contribute zero.
#'
#' @return A `bp_synthetic_trait` object.
#' @export
bp_synthetic_trait <- function(
  name,
  traits,
  fun,
  args = list(),
  linear = FALSE,
  score_level = c("aggregate", "environment"),
  missing_component = c("propagate", "zero_if_all_missing")
) {
  name <- as.character(name)[1L]
  if (is.na(name) || !nzchar(name)) stop("name must be non-empty.", call. = FALSE)
  if (length(traits) == 0L) stop("traits must contain at least one trait.", call. = FALSE)
  if (!is.function(fun)) stop("fun must be a function.", call. = FALSE)
  if (!is.list(args)) stop("args must be a list.", call. = FALSE)
  out <- list(
    name = name,
    traits = traits,
    fun = fun,
    args = args,
    linear = isTRUE(linear),
    score_level = match.arg(score_level),
    missing_component = match.arg(missing_component)
  )
  class(out) <- "bp_synthetic_trait"
  out
}

#' Register Synthetic Traits
#'
#' @param state Program state.
#' @param synthetic_traits One definition or a list of definitions.
#'
#' @return Updated program state.
#' @export
bp_register_synthetic_traits <- function(state, synthetic_traits) {
  defs <- bp_normalize_synthetic_traits(synthetic_traits)
  if (is.null(state$sim$synthetic_traits)) state$sim$synthetic_traits <- list()
  for (def in defs) state$sim$synthetic_traits[[def$name]] <- def
  state
}

bp_normalize_synthetic_traits <- function(x) {
  if (is.null(x)) return(list())
  if (inherits(x, "bp_synthetic_trait")) return(list(x))
  if (!is.list(x)) stop("synthetic_traits must be a definition or list of definitions.", call. = FALSE)
  if (length(x) == 0L) return(list())
  if (all(vapply(x, function(z) inherits(z, "bp_synthetic_trait"), logical(1)))) return(unname(x))
  stop("Every synthetic trait must be created by bp_synthetic_trait().", call. = FALSE)
}

bp_resolve_synthetic_traits <- function(state, x) {
  if (is.null(x)) return(list())
  if (inherits(x, "bp_synthetic_trait") || is.list(x) && all(vapply(x, function(z) inherits(z, "bp_synthetic_trait"), logical(1)))) {
    return(bp_normalize_synthetic_traits(x))
  }
  names_requested <- as.character(x)
  registry <- state$sim$synthetic_traits %||% list()
  missing <- setdiff(names_requested, names(registry))
  if (length(missing) > 0L) {
    stop(sprintf("Unknown synthetic trait(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  unname(registry[names_requested])
}

#' Extract a Biological Trait Matrix
#'
#' @param pop AlphaSimR `Pop` object.
#' @param use Matrix source: `gv`, `pheno`, or `ebv`.
#' @param traits Optional biological trait indices or names.
#'
#' @return Numeric matrix.
#' @export
bp_trait_matrix <- function(pop, use = c("gv", "pheno", "ebv"), traits = NULL) {
  use <- match.arg(use)
  bp_report_matrix(pop, traits = traits, slot = use)
}

bp_synthetic_missing_component <- function(def) {
  def$missing_component %||% "propagate"
}

bp_apply_synthetic_missing_rule <- function(mat, def) {
  mat <- as.matrix(mat)
  if (identical(bp_synthetic_missing_component(def), "zero_if_all_missing") && ncol(mat) > 0L) {
    all_missing <- vapply(seq_len(ncol(mat)), function(j) all(is.na(mat[, j])), logical(1))
    if (any(all_missing)) mat[, all_missing] <- 0
  }
  mat
}

bp_synthetic_component_matrix <- function(pop, use = c("gv", "pheno", "ebv"), def) {
  use <- match.arg(use)
  full <- bp_report_matrix(pop, traits = NULL, slot = use)
  n <- pop_n_ind(pop)
  wanted <- def$traits
  out <- matrix(NA_real_, nrow = n, ncol = length(wanted))
  out_names <- if (is.numeric(wanted)) {
    colnames(full)[as.integer(wanted)] %||% paste0("trait", as.integer(wanted))
  } else {
    as.character(wanted)
  }
  colnames(out) <- out_names
  if (ncol(full) == 0L) return(bp_apply_synthetic_missing_rule(out, def))

  idx <- if (is.numeric(wanted)) {
    as.integer(wanted)
  } else {
    match(as.character(wanted), colnames(full))
  }
  present <- !is.na(idx) & idx >= 1L & idx <= ncol(full)
  if (any(present)) {
    out[, present] <- full[, idx[present], drop = FALSE]
  }
  bp_apply_synthetic_missing_rule(out, def)
}

bp_eval_synthetic_trait <- function(values, synthetic_trait) {
  def <- synthetic_trait
  mat <- bp_apply_synthetic_missing_rule(values, def)
  out <- do.call(def$fun, c(list(mat), def$args))
  out <- as.numeric(out)
  if (length(out) != nrow(mat)) {
    stop(sprintf("Synthetic trait '%s' returned %d values; expected %d.", def$name, length(out), nrow(mat)), call. = FALSE)
  }
  out
}

bp_eval_synthetic_trial_matrix <- function(values, synthetic_trait, measured_traits, pop) {
  def <- synthetic_trait
  mat <- as.matrix(values)
  measured_traits <- as.integer(measured_traits)
  wanted_global <- if (is.numeric(def$traits)) {
    as.integer(def$traits)
  } else {
    all_names <- colnames(pop@gv)
    match(as.character(def$traits), all_names)
  }
  wanted <- match(wanted_global, measured_traits)
  if (anyNA(wanted) && !identical(bp_synthetic_missing_component(def), "zero_if_all_missing")) {
    stop(
      sprintf("Synthetic trait '%s' requires biological traits not measured in this trial.", def$name),
      call. = FALSE
    )
  }
  component_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = length(def$traits))
  present <- !is.na(wanted)
  if (any(present)) component_mat[, present] <- mat[, wanted[present], drop = FALSE]
  bp_eval_synthetic_trait(component_mat, def)
}

bp_score_synthetic_trial <- function(pop, definitions, measured_traits, env_pheno = NULL) {
  defs <- bp_normalize_synthetic_traits(definitions)
  if (length(defs) == 0L) return(list(pop = pop, aggregate = NULL, environments = NULL))
  aggregate <- matrix(NA_real_, nrow = pop@nInd, ncol = length(defs))
  colnames(aggregate) <- vapply(defs, `[[`, character(1), "name")
  env_scores <- if (!is.null(env_pheno)) vector("list", length(env_pheno)) else NULL
  if (!is.null(env_scores)) {
    env_scores <- lapply(env_scores, function(x) {
      out <- matrix(NA_real_, nrow = pop@nInd, ncol = length(defs))
      colnames(out) <- colnames(aggregate)
      out
    })
  }
  aggregate_pheno <- pop@pheno[, as.integer(measured_traits), drop = FALSE]
  for (j in seq_along(defs)) {
    def <- defs[[j]]
    if (!is.null(env_pheno)) {
      per_env <- vapply(
        env_pheno,
        function(x) bp_eval_synthetic_trial_matrix(x, def, measured_traits, pop),
        numeric(pop@nInd)
      )
      for (e in seq_along(env_scores)) env_scores[[e]][, j] <- per_env[, e]
      aggregate[, j] <- if (identical(def$score_level, "environment")) {
        rowMeans(per_env)
      } else {
        bp_eval_synthetic_trial_matrix(aggregate_pheno, def, measured_traits, pop)
      }
    } else {
      aggregate[, j] <- bp_eval_synthetic_trial_matrix(aggregate_pheno, def, measured_traits, pop)
    }
    pop <- bp_set_synthetic_values(pop, def$name, aggregate[, j], type = "pheno")
  }
  list(pop = pop, aggregate = aggregate, environments = env_scores)
}

#' Evaluate a Synthetic Trait
#'
#' @param pop AlphaSimR `Pop` object.
#' @param synthetic_trait Synthetic-trait definition.
#' @param use Biological source matrix: `gv`, `pheno`, or `ebv`.
#'
#' @return Numeric vector with one value per individual.
#' @export
bp_synthetic_values <- function(pop, synthetic_trait, use = c("gv", "pheno", "ebv")) {
  use <- match.arg(use)
  def <- bp_normalize_synthetic_traits(synthetic_trait)[[1L]]
  mat <- bp_synthetic_component_matrix(pop, use = use, def = def)
  bp_eval_synthetic_trait(mat, def)
}

bp_synthetic_store_key <- function(type) {
  type <- match.arg(type, c("pheno", "ebv", "gv"))
  paste0("bps_synthetic_", type)
}

#' Store Synthetic Trait Values
#'
#' @param pop AlphaSimR `Pop` object.
#' @param name Synthetic-trait name.
#' @param values Numeric vector with one value per individual.
#' @param type Storage type: `pheno`, `ebv`, or `gv`.
#'
#' @return Updated population.
#' @export
bp_set_synthetic_values <- function(pop, name, values, type = c("pheno", "ebv", "gv")) {
  if (!methods::is(pop, "Pop")) stop("pop must be an AlphaSimR Pop object.", call. = FALSE)
  key <- bp_synthetic_store_key(match.arg(type))
  name <- as.character(name)[1L]
  values <- as.numeric(values)
  if (length(values) != pop@nInd) stop("values must have one value per individual.", call. = FALSE)
  current <- pop@misc[[key]]
  if (is.null(current) || !is.matrix(current) || nrow(current) != pop@nInd) {
    current <- matrix(numeric(0), nrow = pop@nInd, ncol = 0L)
  }
  if (name %in% colnames(current)) {
    current[, name] <- values
  } else {
    current <- cbind(current, stats::setNames(data.frame(values, check.names = FALSE), name))
    current <- as.matrix(current)
  }
  pop@misc[[key]] <- current
  pop
}

#' Get Stored Synthetic Trait Values
#'
#' @param pop AlphaSimR `Pop` object.
#' @param name Synthetic-trait name.
#' @param type Storage type: `pheno`, `ebv`, or `gv`.
#' @param missing_ok Return `NULL` instead of error when unavailable.
#'
#' @return Numeric vector or `NULL`.
#' @export
bp_get_stored_synthetic_values <- function(pop, name, type = c("pheno", "ebv", "gv"), missing_ok = FALSE) {
  key <- bp_synthetic_store_key(match.arg(type))
  mat <- pop@misc[[key]]
  if (is.null(mat) || !is.matrix(mat) || !(name %in% colnames(mat))) {
    if (isTRUE(missing_ok)) return(NULL)
    stop(sprintf("Stored synthetic %s '%s' is unavailable.", match.arg(type), name), call. = FALSE)
  }
  as.numeric(mat[, name])
}

bp_draw_correlated_noise <- function(n, covariance) {
  covariance <- as.matrix(covariance)
  eig <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  tol <- max(1, max(abs(eig$values))) * 1e-10
  if (any(eig$values < -tol)) stop("varE must be positive semidefinite.", call. = FALSE)
  root <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values))
  matrix(stats::rnorm(n * nrow(covariance)), nrow = n) %*% t(root)
}

bp_subset_covariance <- function(varE, traits, n_traits_total) {
  if (is.null(varE)) return(matrix(0, nrow = length(traits), ncol = length(traits)))
  if (is.null(dim(varE))) {
    x <- as.numeric(varE)
    if (length(x) == 1L) return(diag(x, nrow = length(traits)))
    if (length(x) == length(traits)) return(diag(x, nrow = length(traits)))
    if (is.numeric(traits) && length(x) >= max(as.integer(traits))) return(diag(x[as.integer(traits)], nrow = length(traits)))
    stop("varE vector is incompatible with synthetic trait components.", call. = FALSE)
  }
  mat <- as.matrix(varE)
  if (all(dim(mat) == length(traits))) return(mat)
  if (is.numeric(traits) && nrow(mat) >= max(as.integer(traits))) {
    return(mat[as.integer(traits), as.integer(traits), drop = FALSE])
  }
  if (is.character(traits) && !is.null(rownames(mat)) && all(traits %in% rownames(mat))) {
    return(mat[traits, traits, drop = FALSE])
  }
  if (all(dim(mat) == n_traits_total)) {
    idx <- if (is.numeric(traits)) as.integer(traits) else match(traits, colnames(mat))
    return(mat[idx, idx, drop = FALSE])
  }
  stop("varE matrix is incompatible with synthetic trait components.", call. = FALSE)
}

#' Estimate Genetic Values for a Synthetic Trait
#'
#' Linear traits are evaluated directly. Nonlinear traits estimate expected
#' performance across random TPE trials and plant residuals using common random
#' numbers across individuals.
#'
#' @param pop AlphaSimR `Pop` object.
#' @param synthetic_trait Synthetic-trait definition.
#' @param state Program state containing `state$sim$SP`.
#' @param varE Residual variance vector or covariance matrix.
#' @param n_trials Number of random TPE trials. Default `100`.
#' @param n_plants_per_trial Number of residual plant draws per trial. Default
#'   `10`.
#' @param seed Optional random seed.
#'
#' @return Numeric vector with one estimated synthetic GV per individual.
#' @export
bp_get_synthetic_gv <- function(
  pop,
  synthetic_trait,
  state,
  varE = NULL,
  n_trials = 100L,
  n_plants_per_trial = 10L,
  seed = NULL
) {
  def <- bp_normalize_synthetic_traits(synthetic_trait)[[1L]]
  if (isTRUE(def$linear)) return(bp_synthetic_values(pop, def, use = "gv"))
  if (!is.null(seed)) set.seed(seed)
  n_trials <- as.integer(n_trials)
  n_plants_per_trial <- as.integer(n_plants_per_trial)
  if (n_trials < 1L || n_plants_per_trial < 1L) stop("n_trials and n_plants_per_trial must be positive.", call. = FALSE)

  base_gv <- bp_trait_matrix(pop, use = "gv", traits = def$traits)
  k <- ncol(base_gv)
  trait_idx <- if (is.numeric(def$traits)) {
    as.integer(def$traits)
  } else {
    match(as.character(def$traits), colnames(pop@gv))
  }
  if (anyNA(trait_idx)) stop("Synthetic trait components were not found in pop@gv.", call. = FALSE)
  Sigma <- bp_subset_covariance(varE, def$traits, ncol(pop@gv))
  SP <- state$sim$SP

  env_sd <- numeric(k)
  gxe <- matrix(0, nrow = pop@nInd, ncol = k)
  for (j in seq_len(k)) {
    tr <- SP$traits[[trait_idx[[j]]]]
    if (!is.null(tr) && "envVar" %in% methods::slotNames(tr) && !is.null(tr@envVar)) {
      env_sd[[j]] <- sqrt(as.numeric(tr@envVar))
    }
    if (length(pop@gxe) >= trait_idx[[j]] && !is.null(pop@gxe[[trait_idx[[j]]]])) {
      gxe[, j] <- as.numeric(pop@gxe[[trait_idx[[j]]]])
    }
  }

  total <- numeric(pop@nInd)
  n_draws <- n_trials * n_plants_per_trial
  for (trial in seq_len(n_trials)) {
    p <- stats::runif(1L)
    z <- stats::qnorm(p, sd = env_sd)
    trial_gv <- base_gv + sweep(gxe, 2L, z, `*`)
    noise <- bp_draw_correlated_noise(n_plants_per_trial, Sigma)
    for (plant in seq_len(n_plants_per_trial)) {
      x <- sweep(trial_gv, 2L, noise[plant, ], `+`)
      total <- total + bp_eval_synthetic_trait(x, def)
    }
  }
  total / n_draws
}

#' Materialize Synthetic Genetic Values
#'
#' Compute synthetic genetic values and cache them in
#' `pop@misc$bps_synthetic_gv`. Existing cached values are reused unless
#' `overwrite = TRUE`.
#'
#' @param pop AlphaSimR `Pop` object.
#' @param synthetic_traits Synthetic-trait definitions or registered names.
#' @param state Program state containing registered traits and `state$sim$SP`.
#' @param varE Residual variance vector or covariance matrix.
#' @param n_trials Number of random TPE trials.
#' @param n_plants_per_trial Number of residual plant draws per trial.
#' @param seed Optional random seed.
#' @param overwrite Whether to replace existing cached values.
#'
#' @return Population with cached synthetic genetic values.
#' @export
bp_materialize_synthetic_gv <- function(
  pop,
  synthetic_traits,
  state,
  varE = NULL,
  n_trials = 100L,
  n_plants_per_trial = 10L,
  seed = NULL,
  overwrite = FALSE
) {
  defs <- bp_resolve_synthetic_traits(state, synthetic_traits)
  for (def in defs) {
    cached <- bp_get_stored_synthetic_values(pop, def$name, type = "gv", missing_ok = TRUE)
    if (is.null(cached) || isTRUE(overwrite)) {
      values <- bp_get_synthetic_gv(
        pop = pop,
        synthetic_trait = def,
        state = state,
        varE = varE,
        n_trials = n_trials,
        n_plants_per_trial = n_plants_per_trial,
        seed = seed
      )
      pop <- bp_set_synthetic_values(pop, def$name, values, type = "gv")
    }
  }
  pop
}

#' Select Individuals by a Synthetic Trait
#'
#' @param pop AlphaSimR `Pop` object.
#' @param n_select Number of individuals to select.
#' @param synthetic_trait Synthetic-trait definition.
#' @param use Synthetic source: `ebv`, `pheno`, or `gv`.
#' @param state Optional state, required for nonlinear synthetic GV evaluation.
#' @param varE Residual covariance used for nonlinear synthetic GV evaluation.
#' @param prediction For EBV selection, prefer direct stored predictions or
#'   derive scores from component EBVs.
#' @param simParam AlphaSimR simulation parameters.
#' @param ... Additional arguments passed to `AlphaSimR::selectInd()`.
#'
#' @return Selected population.
#' @export
bp_select_synthetic <- function(
  pop,
  n_select,
  synthetic_trait,
  use = c("ebv", "pheno", "gv"),
  state = NULL,
  varE = NULL,
  prediction = c("auto", "direct", "derived"),
  simParam = NULL,
  ...
) {
  use <- match.arg(use)
  prediction <- match.arg(prediction)
  def <- bp_normalize_synthetic_traits(synthetic_trait)[[1L]]
  values <- if (use == "ebv") {
    direct <- bp_get_stored_synthetic_values(pop, def$name, "ebv", missing_ok = TRUE)
    if (prediction == "direct" && is.null(direct)) stop("Direct synthetic EBV is unavailable.", call. = FALSE)
    if (!is.null(direct) && prediction != "derived") direct else bp_synthetic_values(pop, def, "ebv")
  } else if (use == "pheno") {
    stored <- bp_get_stored_synthetic_values(pop, def$name, "pheno", missing_ok = TRUE)
    if (!is.null(stored)) stored else bp_synthetic_values(pop, def, "pheno")
  } else {
    stored <- bp_get_stored_synthetic_values(pop, def$name, "gv", missing_ok = TRUE)
    if (!is.null(stored)) stored else bp_get_synthetic_gv(pop, def, state = state, varE = varE)
  }
  if (anyNA(values)) stop("Synthetic selection values contain NA.", call. = FALSE)
  source_fn <- function(candidate_pop, ...) {
    idx <- match(candidate_pop@id, pop@id)
    matrix(values[idx], ncol = 1L)
  }
  AlphaSimR::selectInd(
    pop,
    nInd = as.integer(n_select),
    trait = function(x, ...) x,
    use = source_fn,
    simParam = simParam %||% if (!is.null(state)) state$sim$SP else NULL,
    ...
  )
}

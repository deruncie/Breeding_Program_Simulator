#' Install the BPS Codex Skills
#'
#' Install the Codex skills bundled with BPS for the current user or for one
#' project. The installed skills remain tied to the version shipped with the
#' installed BPS package.
#'
#' @details
#' After installation, restart Codex and open a new thread before invoking the
#' skills. A thread that was already open may retain the skill list with which
#' it started.
#'
#' @param scope Install for the current `"user"` or one `"project"`.
#' @param project Project root used when `scope = "project"`. By default, use
#'   the current working directory.
#' @param overwrite Replace existing BPS skill directories.
#'
#' @return Invisibly, the installed skill directories.
#' @export
bp_install_codex_skills <- function(
    scope = c("user", "project"),
    project = getwd(),
    overwrite = FALSE) {
  scope <- match.arg(scope)
  source_root <- system.file(
    "codex-skills",
    package = "BreedingProgramSimulator"
  )

  if (!nzchar(source_root)) {
    stop("The installed BPS package does not contain its Codex skills.")
  }

  destination_root <- if (scope == "user") {
    path.expand("~/.agents/skills")
  } else {
    file.path(path.expand(project), ".agents", "skills")
  }
  dir.create(destination_root, recursive = TRUE, showWarnings = FALSE)

  source_dirs <- list.dirs(source_root, recursive = FALSE, full.names = TRUE)
  installed <- character()
  skipped <- character()

  for (source_dir in source_dirs) {
    skill_name <- basename(source_dir)
    destination <- file.path(destination_root, skill_name)

    if (dir.exists(destination)) {
      if (!isTRUE(overwrite)) {
        skipped <- c(skipped, skill_name)
        next
      }
      unlink(destination, recursive = TRUE, force = TRUE)
    }

    copied <- file.copy(
      source_dir,
      destination_root,
      recursive = TRUE,
      copy.mode = TRUE,
      copy.date = TRUE
    )
    if (!isTRUE(copied)) {
      stop(sprintf("Could not install Codex skill '%s'.", skill_name))
    }
    installed <- c(installed, destination)
  }

  if (length(skipped) > 0L) {
    warning(
      sprintf(
        "Skills already exist and were not replaced: %s. Use overwrite = TRUE to replace them.",
        paste(skipped, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(installed) > 0L) {
    package_version <- as.character(
      utils::packageVersion("BreedingProgramSimulator")
    )
    skill_prompts <- paste0("  $", basename(source_dirs), collapse = "\n")
    message(
      sprintf(
        paste0(
          "Installed BPS %s Codex skills in %s\n\n",
          "Next steps:\n",
          "1. Restart Codex (or reload the IDE extension).\n",
          "2. Open a new thread; an already-open thread may keep its original skill list.\n",
          "3. Invoke either skill by name:\n%s"
        ),
        package_version,
        destination_root,
        skill_prompts
      )
    )
  }

  invisible(installed)
}

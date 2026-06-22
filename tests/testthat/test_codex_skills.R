test_that("bundled Codex skills are complete", {
  skill_root <- system.file(
    "codex-skills",
    package = "BreedingProgramSimulator"
  )
  skill_names <- c(
    "breeding-scheme-drafter",
    "breeding-experiment-designer"
  )

  expect_true(nzchar(skill_root))
  expect_true(all(file.exists(file.path(skill_root, skill_names, "SKILL.md"))))
  expect_true(all(file.exists(
    file.path(skill_root, skill_names, "agents", "openai.yaml")
  )))
})

test_that("Codex skills install at project scope", {
  project <- tempfile("bps-skill-project-")
  dir.create(project)
  on.exit(unlink(project, recursive = TRUE, force = TRUE), add = TRUE)

  expect_message(
    installed <- bp_install_codex_skills(
      scope = "project",
      project = project
    ),
    "Restart Codex"
  )

  expect_length(installed, 2L)
  expect_true(all(file.exists(file.path(installed, "SKILL.md"))))

  expect_warning(
    second_install <- bp_install_codex_skills(
      scope = "project",
      project = project
    ),
    "already exist"
  )
  expect_length(second_install, 0L)

  expect_message(
    replaced <- bp_install_codex_skills(
      scope = "project",
      project = project,
      overwrite = TRUE
    ),
    "\\$breeding-scheme-drafter"
  )
  expect_length(replaced, 2L)
})

test_that("bundled scheme style guide matches the package guide", {
  bundled <- system.file(
    "codex-skills",
    "breeding-scheme-drafter",
    "references",
    "script-style-guide.md",
    package = "BreedingProgramSimulator"
  )
  package_guide <- system.file(
    "style-guide",
    "BREEDING_SCHEME_SCRIPT_STYLE_GUIDE.md",
    package = "BreedingProgramSimulator"
  )

  expect_identical(
    readLines(bundled, warn = FALSE),
    readLines(package_guide, warn = FALSE)
  )
})

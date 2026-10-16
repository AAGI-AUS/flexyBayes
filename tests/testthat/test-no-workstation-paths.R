# =============================================================================
# No workstation path in the repository (review of 2026-09-11).
#
# `data-raw/asreml_reference.R` set .libPaths() to an absolute path under the
# maintainer's home directory. `data-raw/` is excluded from the tarball, so it
# was not an installed-package defect, but it made the ASReml reference
# recipe unrunnable by anyone else, and it put a workstation path in a public
# repository.
#
# The estate-level style rule that forbids this (I13) is enforced by a
# pre-commit grep over `R/`, `tests/`, `vignettes/` and `inst/`. It did not
# fire, because `data-raw/` is not one of the directories it reads. A guard
# protects the altitude it is written at and no other, so this package now
# carries its own, over every source directory instead of four of them.
# =============================================================================

.wp_root <- function() {
  root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = FALSE
  )
  if (file.exists(file.path(root, "DESCRIPTION"))) root else NA_character_
}

test_that("no source file carries an absolute home-directory path", {
  root <- .wp_root()
  if (is.na(root)) {
    testthat::skip("source tree not reachable (installed-package test run)")
  }

  dirs <- c("R", "tests", "vignettes", "inst", "data-raw", "tools", "man")
  dirs <- dirs[dir.exists(file.path(root, dirs))]
  files <- unlist(lapply(
    dirs,
    function(d) {
      list.files(
        file.path(root, d),
        pattern = "\\.(R|r|Rmd|rmd|Rd|md|yaml|yml|sh|txt|csv)$",
        recursive = TRUE,
        full.names = TRUE
      )
    }
  ))
  # This file names the pattern it searches for, so it cannot search itself.
  files <- files[!grepl("test-no-workstation-paths\\.R$", files)]

  # Assembled instead of written whole, for the same reason.
  pattern <- paste0("/", "Users", "/[A-Za-z0-9._-]+/")

  hits <- unlist(lapply(files, function(f) {
    lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) {
      character(0)
    })
    idx <- grep(pattern, lines)
    if (!length(idx)) {
      return(NULL)
    }
    paste0(sub(paste0("^", root, "/"), "", f), ":", idx)
  }))

  expect_equal(hits, NULL)
})

test_that("the ASReml reference recipe reads its library from the environment", {
  root <- .wp_root()
  if (is.na(root)) {
    testthat::skip("source tree not reachable (installed-package test run)")
  }
  path <- file.path(root, "data-raw", "asreml_reference.R")
  if (!file.exists(path)) {
    testthat::skip("data-raw/ not present in this tree")
  }
  src <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(src, "FLEXYBAYES_ASREML_LIB", fixed = TRUE)
})

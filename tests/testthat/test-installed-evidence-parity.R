# =============================================================================
# The installed capability CSV (inst/validation/execution_grid/
# capability_matrix.csv) is a copy of .fb_capability_matrix(), kept in step
# by tools/generate_capability_matrix.R using the row builder in
# tools/capability_grid_matrix.R -- the same one tools/execution_grid.R
# drives its own cells from. Before 2026-09-16 the copy was frozen at commit
# 4c81984 (0.9.3): it still said INLA refused the sectioned residual and the
# combined MET model, four model classes short of the current registry,
# while README.md, inst/KNOWN_ISSUES.md and .fb_capability_matrix() itself
# had all moved on and nothing read the CSV to notice.
#
# This file is the installed-evidence oracle: it reads the CSV the way a
# user does, through system.file(), and fails when it disagrees with the
# live registry on cell, model_class, backend or verdict.
# =============================================================================

.iep_pkg_file <- function(...) {
  candidate <- testthat::test_path("..", "..", ...)
  if (file.exists(candidate)) {
    return(normalizePath(candidate))
  }
  NA_character_
}

# The row builder lives in tools/, outside the package namespace, so it is
# sourced into its own environment rather than exported.
.iep_grid_builder <- function() {
  path <- .iep_pkg_file("tools", "capability_grid_matrix.R")
  testthat::skip_if(
    is.na(path), "tools/capability_grid_matrix.R is not in the installed tree"
  )
  env <- new.env()
  source(path, local = env)
  env
}

# Semantic-column parity: every row of `csv` names a real registry cell at
# the verdict the registry itself records. Returns a character vector of
# mismatch descriptions; empty means the two agree completely.
.iep_diff <- function(csv, cap, m_calls, builder) {
  expected <- builder$.fb_capability_grid_cells(cap, m_calls)
  cols <- c("cell", "model_class", "backend", "verdict")
  left <- expected[order(expected$cell), cols, drop = FALSE]
  right <- csv[order(csv$cell), cols, drop = FALSE]
  rownames(left) <- NULL
  rownames(right) <- NULL

  if (!identical(nrow(left), nrow(right))) {
    return(sprintf(
      "row count differs: registry-derived %d, installed CSV %d",
      nrow(left), nrow(right)
    ))
  }
  if (!identical(left$cell, right$cell)) {
    return("cell sets differ between the registry and the installed CSV")
  }
  bad <- !mapply(identical, left$verdict, right$verdict)
  if (!any(bad)) {
    return(character(0))
  }
  sprintf(
    "%s: registry says '%s', installed CSV says '%s'",
    left$cell[bad], left$verdict[bad], right$verdict[bad]
  )
}

test_that("the installed capability CSV agrees with the live registry", {
  csv_path <- .iep_pkg_file(
    "inst", "validation", "execution_grid", "capability_matrix.csv"
  )
  skip_if(is.na(csv_path), "the installed capability CSV is not in the tree")
  builder <- .iep_grid_builder()

  csv <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  cap <- flexyBayes:::.fb_capability_matrix()
  m_calls <- builder$.fb_capability_grid_m_calls()

  expect_identical(
    .iep_diff(csv, cap, m_calls, builder),
    character(0),
    label = paste0(
      "inst/validation/execution_grid/capability_matrix.csv disagrees with ",
      ".fb_capability_matrix() -- re-run tools/generate_capability_matrix.R"
    )
  )
})

test_that("the parity check rejects a CSV with one changed verdict", {
  # The falsifier. Without the row-by-row verdict comparison above, a
  # check that only counted rows would pass a CSV carrying the right
  # cells at the wrong verdicts.
  csv_path <- .iep_pkg_file(
    "inst", "validation", "execution_grid", "capability_matrix.csv"
  )
  skip_if(is.na(csv_path), "the installed capability CSV is not in the tree")
  builder <- .iep_grid_builder()

  csv <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  cap <- flexyBayes:::.fb_capability_matrix()
  m_calls <- builder$.fb_capability_grid_m_calls()

  mutated <- csv
  flip <- if (identical(mutated$verdict[[1L]], "fits")) "refuses" else "fits"
  mutated$verdict[[1L]] <- flip

  diffs <- .iep_diff(mutated, cap, m_calls, builder)
  expect_gt(length(diffs), 0L)
})

# =============================================================================
# The installed capability CSV (inst/validation/execution_grid/
# capability_matrix.csv) is a copy of .fb_capability_matrix(), kept in step
# by tools/generate_capability_matrix.R using the row builder in
# tools/capability_grid_matrix.R, the same one tools/execution_grid.R
# drives its own cells from. Before 2026-09-16 the copy was frozen at commit
# 4c81984 (0.9.3): it still said INLA refused the sectioned residual and the
# combined MET model, four model classes short of the current registry,
# while README.md, inst/KNOWN_ISSUES.md and .fb_capability_matrix() itself
# had all moved on and nothing read the CSV to notice.
#
# Two checks guard it. The first reads the CSV in the source tree and builds
# the expected rows with the grid's own row builder in tools/, so it also
# catches a missing row. It skips where tools/ is absent, as in an installed
# package. The second reads the CSV through system.file(), as an installed
# user does, and checks every shipped row against the live registry on cell,
# model_class, backend and verdict, so it runs under R CMD check.
# =============================================================================

.iep_pkg_file <- function(...) {
  candidate <- testthat::test_path("..", "..", ...)
  if (file.exists(candidate)) {
    return(normalizePath(candidate))
  }
  NA_character_
}

# The row builder lives in tools/, outside the package namespace, so it is
# sourced into its own environment instead of exported.
.iep_grid_builder <- function() {
  path <- .iep_pkg_file("tools", "capability_grid_matrix.R")
  testthat::skip_if(
    is.na(path),
    "tools/capability_grid_matrix.R is not in the installed tree"
  )
  env <- new.env()
  source(path, local = env)
  env
}

# Semantic-column parity: every row of `csv` names a real registry cell at
# the verdict the registry itself records. Returns a character vector of
# mismatch descriptions. Empty means the two agree completely.
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
      nrow(left),
      nrow(right)
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
    left$cell[bad],
    left$verdict[bad],
    right$verdict[bad]
  )
}

test_that("the installed capability CSV agrees with the live registry", {
  csv_path <- .iep_pkg_file(
    "inst",
    "validation",
    "execution_grid",
    "capability_matrix.csv"
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
      ".fb_capability_matrix(), re-run tools/generate_capability_matrix.R"
    )
  )
})

test_that("the parity check rejects a CSV with one changed verdict", {
  # The falsifier. Without the row-by-row verdict comparison above, a
  # check that only counted rows would pass a CSV carrying the right
  # cells at the wrong verdicts.
  csv_path <- .iep_pkg_file(
    "inst",
    "validation",
    "execution_grid",
    "capability_matrix.csv"
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

# Each shipped row checked against the registry alone: the model class must be
# a registry row, the cell id must be that class's slug and the backend, and
# the verdict must be the registry's for that backend. Returns a character
# vector of mismatch descriptions. Empty means every shipped row agrees.
.iep_installed_diff <- function(csv, cap) {
  slug <- function(x) {
    gsub("(^-|-$)", "", gsub("[^a-z0-9]+", "-", tolower(x)))
  }
  row_of <- match(csv$model_class, cap$model_class)
  expected_cell <- paste0("cap-", slug(csv$model_class), "-", csv$backend)
  expected_verdict <- vapply(
    seq_len(nrow(csv)),
    function(i) {
      if (is.na(row_of[[i]]) || !csv$backend[[i]] %in% names(cap)) {
        return(NA_character_)
      }
      as.character(cap[[csv$backend[[i]]]][[row_of[[i]]]])
    },
    character(1L)
  )
  c(
    sprintf("%s: model class not in the registry", csv$cell[is.na(row_of)]),
    sprintf(
      "%s: cell id should be %s",
      csv$cell[csv$cell != expected_cell],
      expected_cell[csv$cell != expected_cell]
    ),
    sprintf(
      "%s: registry says '%s', installed CSV says '%s'",
      csv$cell[!is.na(expected_verdict) & csv$verdict != expected_verdict],
      expected_verdict[
        !is.na(expected_verdict) &
          csv$verdict != expected_verdict
      ],
      csv$verdict[!is.na(expected_verdict) & csv$verdict != expected_verdict]
    )
  )
}

.iep_installed_csv <- function() {
  path <- system.file(
    "validation",
    "execution_grid",
    "capability_matrix.csv",
    package = "flexyBayes"
  )
  expect_true(nzchar(path))
  utils::read.csv(path, stringsAsFactors = FALSE)
}

test_that("every row of the installed capability CSV is the registry's", {
  csv <- .iep_installed_csv()

  expect_gt(nrow(csv), 0L)
  expect_identical(
    .iep_installed_diff(csv, flexyBayes:::.fb_capability_matrix()),
    character(0)
  )
})

test_that("the installed-row check rejects one changed verdict", {
  csv <- .iep_installed_csv()
  csv$verdict[[1L]] <- if (identical(csv$verdict[[1L]], "fits")) {
    "refuses"
  } else {
    "fits"
  }

  expect_gt(
    length(.iep_installed_diff(csv, flexyBayes:::.fb_capability_matrix())),
    0L
  )
})

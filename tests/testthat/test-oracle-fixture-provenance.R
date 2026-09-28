# =============================================================================
# Banked ASReml evidence names the commit it came from (review finding A21).
#
# The first tranche built its tarball with `R CMD build .` from the working
# directory while stamping the record with `git rev-parse HEAD`. Those agree
# only when the tree is clean, and it was not: five fixtures carry
# `source_sha = 5e4e0ef` for a run whose code was what later became `9b1d93c`.
# The content was probably identical, and the stored evidence gives a reader
# no way to establish that, which is the entire purpose of storing it.
#
# `run_tranche.sh` now refuses a dirty tree and builds from
# `git archive <sha>`. This file is the assertion that checks it, and it
# is two-sided by construction:
#
#   * every banked fixture carries the provenance fields, and its
#     `source_sha` is declared in ANCHORS.md;
#   * and the checking function itself is shown to reject a record that is
#     missing a field, carries an undeclared SHA, or was run dirty -- so a
#     checker that accepted everything could not pass this file.
#
# The git-ancestry half runs only where git and the repository are available,
# which is the maintainer's machine and not a CRAN or CI check of a tarball.
# =============================================================================

.ofp_dir <- function() {
  d <- system.file("validation", "asreml-4.2", package = "flexyBayes")
  if (!nzchar(d) || !dir.exists(d)) {
    d <- testthat::test_path("..", "..", "inst", "validation", "asreml-4.2")
  }
  d
}

.ofp_anchors <- function(dir) {
  path <- file.path(dir, "ANCHORS.md")
  if (!file.exists(path)) {
    return(character(0))
  }
  txt <- readLines(path, warn = FALSE)
  # One anchor per line, as a 40-character hex SHA in backticks at the start
  # of a table row. Read rather than parsed loosely: an anchor that is only
  # mentioned in prose is not a declaration.
  m <- regmatches(txt, regexpr("`[0-9a-f]{40}`", txt))
  gsub("`", "", m, fixed = TRUE)
}

# The verdict function. Returns character(0) for a record that is sound, and
# one message per problem otherwise.
.ofp_verdict <- function(record, anchors) {
  msgs <- character(0)
  need <- c("scenario", "source_sha", "tarball_sha256", "versions", "provenance")
  missing <- setdiff(need, names(record))
  if (length(missing)) {
    msgs <- c(msgs, paste0("missing field(s): ", paste(missing, collapse = ", ")))
    return(msgs)
  }
  p <- record$provenance
  for (f in c("scenarios_sha256", "compare_sha256", "runner_sha256", "dirty",
              "run_command", "blas", "rng_kind", "platform")) {
    if (is.null(p[[f]])) {
      msgs <- c(msgs, paste0("provenance is missing `", f, "`"))
    }
  }
  if (isTRUE(p$dirty)) {
    msgs <- c(msgs, "the tree was dirty when this evidence was produced")
  }
  if (!grepl("^[0-9a-f]{40}$", record$source_sha %||% "")) {
    msgs <- c(msgs, paste0("source_sha is not a full SHA: ", record$source_sha))
  } else if (!record$source_sha %in% anchors) {
    msgs <- c(
      msgs,
      paste0("source_sha ", substr(record$source_sha, 1L, 12L),
             " is not declared in ANCHORS.md")
    )
  }
  if (!nzchar(record$tarball_sha256 %||% "")) {
    msgs <- c(msgs, "tarball_sha256 is empty")
  }
  msgs
}

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a


test_that("every banked ASReml fixture names a declared commit and a clean tree", {
  dir <- .ofp_dir()
  skip_if_not(dir.exists(dir))
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  skip_if(length(files) == 0L)

  anchors <- .ofp_anchors(dir)
  expect_gt(length(anchors), 0L)

  for (f in files) {
    record <- readRDS(f)
    expect_identical(
      .ofp_verdict(record, anchors),
      character(0),
      label = basename(f)
    )
  }
})


test_that("the provenance check rejects the records it exists to reject", {
  # The falsifier. Without this block a verdict function that returned
  # character(0) unconditionally would pass the test above.
  good <- list(
    scenario = "01_iid_intercept",
    source_sha = strrep("a", 40L),
    tarball_sha256 = strrep("b", 64L),
    versions = list(r = "R version 4.5.2"),
    provenance = list(
      scenarios_sha256 = strrep("c", 64L),
      compare_sha256 = strrep("d", 64L),
      runner_sha256 = strrep("e", 64L),
      dirty = FALSE,
      run_command = "run_tranche.sh <work-dir> 01_iid_intercept",
      blas = "libRblas.dylib",
      lapack = "libRlapack.dylib",
      rng_kind = "Mersenne-Twister / Inversion / Rejection",
      platform = "aarch64-apple-darwin20"
    )
  )
  anchors <- strrep("a", 40L)

  expect_identical(.ofp_verdict(good, anchors), character(0))

  # Undeclared commit.
  bad <- good
  bad$source_sha <- strrep("f", 40L)
  expect_match(.ofp_verdict(bad, anchors), "not declared", all = FALSE)

  # A short or absent SHA, which is what a stamped-but-unbuilt record looks
  # like when someone abbreviates it.
  bad <- good
  bad$source_sha <- "5e4e0ef"
  expect_match(.ofp_verdict(bad, anchors), "not a full SHA", all = FALSE)

  # Run on a dirty tree, which the runner now refuses but can be forced.
  bad <- good
  bad$provenance$dirty <- TRUE
  expect_match(.ofp_verdict(bad, anchors), "dirty", all = FALSE)

  # A record from the old runner, which recorded no provenance block at all.
  bad <- good
  bad$provenance <- NULL
  expect_match(.ofp_verdict(bad, anchors), "missing field", all = FALSE)

  # A provenance block that forgot the comparator's own hash: the manifest
  # would pin the package and not the comparison.
  bad <- good
  bad$provenance$scenarios_sha256 <- NULL
  expect_match(.ofp_verdict(bad, anchors), "scenarios_sha256", all = FALSE)
})


test_that("each declared anchor is an ancestor of HEAD", {
  # Maintainer-side only: a CRAN check of a tarball has no repository.
  skip_on_cran()
  skip_on_ci()
  skip_if_not(nzchar(Sys.which("git")))
  root <- tryCatch(
    normalizePath(testthat::test_path("..", "..")),
    error = function(e) NA_character_
  )
  skip_if(is.na(root) || !dir.exists(file.path(root, ".git")))

  anchors <- .ofp_anchors(file.path(root, "inst", "validation", "asreml-4.2"))
  skip_if(length(anchors) == 0L)

  for (sha in anchors) {
    status <- suppressWarnings(system2(
      "git",
      c("-C", shQuote(root), "merge-base", "--is-ancestor", sha, "HEAD"),
      stdout = FALSE, stderr = FALSE
    ))
    expect_identical(
      status, 0L,
      label = paste0("anchor ", substr(sha, 1L, 12L), " reachable from HEAD")
    )
  }
})

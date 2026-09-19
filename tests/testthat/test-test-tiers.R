# test-test-tiers.R -- the long test tier and the release criterion that
# proves it ran.
#
# skip_if_not_long_tier() (helper-test-tiers.R) holds a block back unless
# FLEXYBAYES_RUN_LONG is "true". Criterion D11 of tools/definition_of_done.R
# reads a suite record and fails when a block was held back or a tagged block
# is missing. The D11 blocks source the tool from the source tree and skip
# where tools/ is absent, as in a built tarball.

.tiers_root <- function() {
  root <- testthat::test_path("..", "..")
  if (!file.exists(file.path(root, "tools", "definition_of_done.R"))) {
    testthat::skip("tools/definition_of_done.R not present in test env")
  }
  normalizePath(root)
}

.tiers_dod <- function() {
  env <- new.env(parent = globalenv())
  sys.source(file.path(.tiers_root(), "tools", "definition_of_done.R"),
             envir = env)
  env
}

# A record shaped as tools/run_test_suite.R writes it: every tagged block of
# this tree plus one untagged block, none skipped.
.tiers_record <- function(dod, path) {
  tagged <- withr::with_dir(.tiers_root(), dod$.dod_long_tier_blocks())
  tab <- data.frame(
    file = c(tagged$file, "test-other.R"),
    test = c(tagged$test, "an untagged block"),
    nb = 1L, failed = 0L, skipped = FALSE, error = FALSE, warning = 0L,
    passed = 1L, real = 0.1, skip_reason = ""
  )
  utils::write.table(tab, path, sep = "\t", row.names = FALSE, quote = FALSE)
  tab
}

test_that("the long tier holds a block back unless FLEXYBAYES_RUN_LONG is true", {
  for (value in c("", "TRUE", "1", "yes")) {
    withr::local_envvar(FLEXYBAYES_RUN_LONG = value)
    reason <- tryCatch(
      {
        skip_if_not_long_tier()
        "ran"
      },
      skip = function(s) conditionMessage(s)
    )
    expect_match(reason, "FLEXYBAYES_RUN_LONG", fixed = TRUE)
  }
  withr::local_envvar(FLEXYBAYES_RUN_LONG = "true")
  expect_true(skip_if_not_long_tier())
})

test_that("D11 passes a record in which every block ran", {
  dod <- .tiers_dod()
  path <- withr::local_tempfile(fileext = ".tsv")
  .tiers_record(dod, path)
  res <- withr::with_dir(.tiers_root(), dod$.dod_tier_check(path))
  expect_true(res$pass)
})

test_that("D11 fails a record in which the long tier was held back", {
  dod <- .tiers_dod()
  path <- withr::local_tempfile(fileext = ".tsv")
  tab <- .tiers_record(dod, path)
  tab$skipped[[1L]] <- TRUE
  tab$skip_reason[[1L]] <- "Reason: FLEXYBAYES_RUN_LONG != \"true\"; long tier"
  utils::write.table(tab, path, sep = "\t", row.names = FALSE, quote = FALSE)
  res <- withr::with_dir(.tiers_root(), dod$.dod_tier_check(path))
  expect_false(res$pass)
  expect_match(res$detail, "^1 of ")
})

test_that("D11 fails a record that lacks a tagged block, and a missing record", {
  dod <- .tiers_dod()
  path <- withr::local_tempfile(fileext = ".tsv")
  tab <- .tiers_record(dod, path)
  if (nrow(tab) > 1L) {
    utils::write.table(tab[-1L, ], path, sep = "\t", row.names = FALSE,
                       quote = FALSE)
    res <- withr::with_dir(.tiers_root(), dod$.dod_tier_check(path))
    expect_false(res$pass)
    expect_match(res$detail, "; 1 of ")
  }
  res <- withr::with_dir(.tiers_root(), dod$.dod_tier_check(""))
  expect_false(res$pass)
})

test_that("the skip ledger files the long tier under the time budget", {
  p <- file.path(.tiers_root(), "tools", "skip_ledger.R")
  env <- new.env(parent = baseenv())
  sys.source(p, envir = env)
  expect_identical(
    env$.skip_classify_site("skip_if_not_long_tier", ""),
    "cran_ci_time_budget"
  )
})

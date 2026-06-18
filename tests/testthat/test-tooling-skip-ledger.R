# test-tooling-skip-ledger.R -- smoke tests for tools/skip_ledger.R.
#
# The ledger generator lives at flexyBayes/tools/skip_ledger.R (excluded
# from the package tarball by .Rbuildignore line 36, `^tools$`).  These
# tests source it from the workspace tree and exercise the classifier
# + Markdown writer against a small fixture corpus.  When the tarball
# layout removes tools/, every test in this file skips cleanly.

.skip_ledger_path <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "tools", "skip_ledger.R"),
    "tools/skip_ledger.R",
    "flexyBayes/tools/skip_ledger.R"
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) NA_character_ else hit[1L]
}

.source_skip_ledger <- function() {
  p <- .skip_ledger_path()
  if (is.na(p)) {
    testthat::skip("tools/skip_ledger.R not present in test env")
  }
  env <- new.env(parent = baseenv())
  source(p, local = env)
  env
}

test_that("skip_ledger.R sources and exposes build_skip_ledger + helpers", {
  env <- .source_skip_ledger()
  expect_true(is.function(env$build_skip_ledger))
  expect_true(is.function(env$.skip_classify_site))
  expect_true(is.function(env$.parse_one_test_file))
})

test_that(".skip_classify_site routes flake-prefix tokens with top priority", {
  env <- .source_skip_ledger()
  expect_identical(
    env$.skip_classify_site(
      "skip",
      "[flake-test-infrastructure] suite-order race on flag"
    ),
    "test_infrastructure_flake"
  )
  expect_identical(
    env$.skip_classify_site(
      "skip",
      "[flake-stochastic-rng] TF RNG leak"
    ),
    "stochastic_rng_dependent"
  )
})

test_that(".skip_classify_site routes audit's five reason classes", {
  env <- .source_skip_ledger()
  expect_identical(
    env$.skip_classify_site("skip_on_cran", ""),
    "cran_ci_time_budget"
  )
  expect_identical(
    env$.skip_classify_site("skip_if_no_greta", ""),
    "unavailable_optional_backend"
  )
  expect_identical(
    env$.skip_classify_site("skip_if_not_installed", "greta"),
    "unavailable_optional_backend"
  )
  expect_identical(
    env$.skip_classify_site(
      "skip",
      "FLEXYBAYES_RUN_STRESS != true; skip Tier-3 stress"
    ),
    "stress_only"
  )
  expect_identical(
    env$.skip_classify_site(
      "skip",
      "defer to v0.3.9 once isolation pass lands"
    ),
    "deliberately_deferred_feature"
  )
})

test_that(".parse_one_test_file returns the expected columns + types", {
  env <- .source_skip_ledger()
  fixture <- tempfile(fileext = ".R")
  writeLines(
    c(
      "test_that('fixture skip site', {",
      "  testthat::skip_if_not_installed('greta')",
      "  testthat::skip('[flake-stochastic-rng] fixture flake')",
      "  expect_true(TRUE)",
      "})"
    ),
    fixture
  )
  on.exit(unlink(fixture))

  dt <- env$.parse_one_test_file(fixture)
  expect_s3_class(dt, "data.table")
  expect_setequal(
    names(dt),
    c("file", "line", "call", "reason", "reason_class")
  )
  expect_equal(nrow(dt), 2L)
  expect_true("stochastic_rng_dependent" %in% dt$reason_class)
  expect_true("unavailable_optional_backend" %in% dt$reason_class)
})

test_that("build_skip_ledger writes a valid Markdown ledger", {
  env <- .source_skip_ledger()
  tdir <- tempfile()
  dir.create(file.path(tdir, "tests", "testthat"), recursive = TRUE)
  on.exit(unlink(tdir, recursive = TRUE))

  writeLines(
    c(
      "test_that('a', { testthat::skip_on_cran(); expect_true(TRUE) })",
      "test_that('b', { testthat::skip_if_not_installed('INLA'); expect_true(TRUE) })"
    ),
    file.path(tdir, "tests", "testthat", "test-fixture.R")
  )

  out_path <- file.path(tdir, "inst", "skip-ledger.md")
  dt <- env$build_skip_ledger(
    test_dir = file.path(tdir, "tests", "testthat"),
    output_path = out_path
  )
  expect_true(file.exists(out_path))
  body <- readLines(out_path)
  expect_match(body[1L], "^# Skip ledger$")
  expect_true(any(grepl("CRAN / CI time budget", body)))
  expect_true(any(grepl("Unavailable optional backend", body)))
  expect_equal(nrow(dt), 2L)
})

test_that(".skip_ledger_reason_classes lists all seven classes", {
  env <- .source_skip_ledger()
  expect_setequal(
    names(env$.skip_ledger_reason_classes),
    c(
      "unavailable_optional_backend",
      "cran_ci_time_budget",
      "stress_only",
      "external_toolchain_unavailable",
      "deliberately_deferred_feature",
      "test_infrastructure_flake",
      "stochastic_rng_dependent"
    )
  )
})

# ---------------------------------------------------------------- #
# v0.4.0 section-11 gate: taxonomy lock over the real suite.        #
# ---------------------------------------------------------------- #

test_that("every real skip site classifies within the seven-class taxonomy", {
  env <- .source_skip_ledger()
  expect_true(is.function(env$.assert_taxonomy_lock))

  # Resolve the live testthat directory (the dir holding this file).
  candidates <- c(
    testthat::test_path(),
    "tests/testthat",
    "flexyBayes/tests/testthat"
  )
  test_dir <- candidates[dir.exists(candidates)][1L]
  skip_if(is.na(test_dir), "testthat dir not resolvable")

  dt <- env$.assert_taxonomy_lock(test_dir) # stops on any out-of-class site
  expect_true(nrow(dt) >= 1L)
  expect_true(all(
    dt$reason_class %in%
      names(env$.skip_ledger_reason_classes)
  ))
})

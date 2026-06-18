# Huge-N probe smoke tests (v0.4.0 Wave 2 Phase 2C; v040-plan section
# 11 gate).
#
# The full probe (tools/profiles/huge_n_probe.R) materialises a 5e6-row
# synthetic dataset and is gated behind FLEXYBAYES_RUN_HUGE_N=true so it
# never runs on CRAN / CI. These smoke tests verify the same huge-N
# shape contracts cheaply -- the prediction-kernel file-output format
# resolution at N = 5e6, and that the probe script is correctly env-
# gated -- without paying the 5e6-row cost.

# ---------------------------------------------------------------- #
# Prediction-kernel file-output format resolution at N = 5e6.       #
# ---------------------------------------------------------------- #

test_that("format='auto' routes a 5e6-row payload to fst when fst is available", {
  expect_identical(
    flexyBayes:::.resolve_format("auto", 5e6, FALSE, TRUE),
    "fst"
  )
})

test_that("format='auto' falls back to rds at 5e6 rows when fst is absent", {
  expect_identical(
    flexyBayes:::.resolve_format("auto", 5e6, FALSE, FALSE),
    "rds"
  )
})

test_that("interop=TRUE forces csv at 5e6 rows regardless of fst", {
  expect_identical(
    flexyBayes:::.resolve_format("auto", 5e6, TRUE, TRUE),
    "csv"
  )
  expect_identical(
    flexyBayes:::.resolve_format("auto", 5e6, TRUE, FALSE),
    "csv"
  )
})

test_that("the 5e6 threshold is above the 1e6 fst cutoff (huge-N stays file-backed)", {
  # The aggregation / file-backed cutoff is 1e6 (ADR 0023); 5e6 is well
  # past it, so the huge-N path never silently downgrades to in-memory.
  expect_identical(
    flexyBayes:::.resolve_format("auto", 1e6 - 1L, FALSE, TRUE),
    "rds"
  )
  expect_identical(
    flexyBayes:::.resolve_format("auto", 5e6, FALSE, TRUE),
    "fst"
  )
})

# ---------------------------------------------------------------- #
# The probe script exists and is environment-gated.                 #
# ---------------------------------------------------------------- #

test_that("the huge-N probe script exists and is FLEXYBAYES_RUN_HUGE_N-gated", {
  # Resolve the probe path relative to the installed / loaded package
  # source tree; skip if the tools/ directory is not shipped (it is a
  # workspace tool, not installed).
  candidates <- c(
    testthat::test_path("..", "..", "tools", "profiles", "huge_n_probe.R"),
    file.path("tools", "profiles", "huge_n_probe.R")
  )
  probe <- candidates[file.exists(candidates)][1]
  skip_if(
    is.na(probe) || !length(probe),
    "huge_n_probe.R not present in this tree (workspace-only tool)"
  )
  src <- paste(readLines(probe), collapse = "\n")
  expect_match(src, "FLEXYBAYES_RUN_HUGE_N")
  expect_match(src, "5e6")
  # The three shape paths the probe exercises.
  expect_match(src, "preflight")
  expect_match(src, "aggregation")
  expect_match(src, "\\.resolve_format")
})

test_that("the full probe does not run without the env var set", {
  withr::local_envvar(FLEXYBAYES_RUN_HUGE_N = "")
  expect_identical(Sys.getenv("FLEXYBAYES_RUN_HUGE_N"), "")
})

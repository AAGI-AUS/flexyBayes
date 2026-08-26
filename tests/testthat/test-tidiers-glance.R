# Tests for glance() on all three fit shapes (D15 / WP-D, found by
# WP-G2's stress run; see reports/WP-G2.md Handoffs).
#
# Before this fix, glance.flexybayes_inla() hard-refused, and the PARENT
# glance.flexybayes() (reachable by calling it directly, bypassing S3
# dispatch, as a diagnostic probe would) crashed with "arguments imply
# differing number of rows: 1, 0" -- INLA's model_info has no n_params
# key and call_info's chains/n_samples are NULL, both feeding
# data.frame() unevenly against the 1-row numeric/character columns.
# glance() now returns the same 10-column one-row shape on INLA, brms
# and an aggregated (INLA) fit alike, with sampler-specific columns
# (chains, samples, max_rhat, min_ess) NA on the two INLA-backed shapes
# rather than guessed at.

.glance_cols <- c(
  "nobs", "npar", "logLik", "family", "link",
  "chains", "samples", "max_rhat", "min_ess", "run_time"
)

test_that("glance() returns the one-row shape on a plain INLA fit", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(20260826L)
  d <- data.frame(y = stats::rnorm(60), x = stats::rnorm(60),
                   g = factor(rep(1:6, 10)))
  fit <- suppressMessages(flexybayes(
    y ~ x + (1 | g), data = d, backend = "inla", verbose = FALSE
  ))
  g <- suppressMessages(glance(fit))
  expect_s3_class(g, "data.frame")
  expect_identical(nrow(g), 1L)
  expect_named(g, .glance_cols)
  expect_identical(g$nobs, 60L)
  expect_identical(g$family, "gaussian")
  expect_true(is.na(g$logLik))
  expect_true(is.na(g$chains))
  expect_true(is.na(g$samples))
  expect_true(is.na(g$max_rhat))
  expect_true(is.na(g$min_ess))
  expect_true(is.numeric(g$run_time) && !is.na(g$run_time))
})

test_that("glance() returns the one-row shape on a brms fit (regression guard)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(20260826L)
  d <- data.frame(y = stats::rnorm(60), x = stats::rnorm(60),
                   g = factor(rep(1:6, 10)))
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ x + (1 | g),
    data = d,
    backend = "brms",
    n_samples = 200,
    warmup = 200,
    chains = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  g <- glance(fit)
  expect_s3_class(g, "data.frame")
  expect_identical(nrow(g), 1L)
  expect_named(g, .glance_cols)
  expect_identical(g$nobs, 60L)
  expect_identical(g$chains, 1)
  expect_identical(g$samples, 200)
  expect_true(is.numeric(g$logLik) && !is.na(g$logLik))
})

test_that("glance() returns the one-row shape on an aggregated (INLA) fit", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(20260826L)
  d <- data.frame(
    y = stats::rnorm(2000L),
    env = factor(rep(1:5, 400L)),
    geno = factor(rep(1:20, 100L))
  )
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~geno,
    data = d,
    backend = "inla",
    aggregate = TRUE,
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")
  g <- suppressMessages(glance(fit))
  expect_s3_class(g, "data.frame")
  expect_identical(nrow(g), 1L)
  expect_named(g, .glance_cols)
  expect_identical(g$nobs, 2000L)
  expect_true(is.na(g$chains))
  expect_true(is.na(g$max_rhat))
})

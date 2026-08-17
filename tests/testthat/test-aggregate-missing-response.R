# =============================================================================
# The aggregated Gaussian path and missing responses.
#
# The aggregated path compresses rows into per-cell means and counts. A missing
# response has no sufficient statistic to contribute, and the compression did
# not account for that: the cell mean came back NA and a downstream test failed
# with "missing value where TRUE/FALSE needed".
#
# The failure was reachable from an ordinary call -- a randomised complete
# block design with a lost plot, fitted with the default aggregate = "auto" --
# and it surfaced only when the augmentation layer began retaining
# missing-response rows instead of letting them be dropped upstream.
#
# The repair is to fall through to the per-row path, NOT to teach the
# compression to skip NA rows. Skipping them would make the aggregated path
# mean complete-case deletion while the per-row path means augmentation, so
# the same call would answer two different questions depending on a
# compression ratio the user never chose.
# =============================================================================

.rcbd_with_hole <- function(n_row = 8L, n_col = 6L, n_missing = 7L,
                            seed = 1L) {
  set.seed(seed)
  d <- expand.grid(col_num = seq_len(n_col), row_num = seq_len(n_row))
  d$blk <- factor(rep(seq_len(6L), length.out = nrow(d)))
  b <- stats::rnorm(6L, 0, sqrt(1.5))
  d$y <- 20 + b[as.integer(d$blk)] + stats::rnorm(nrow(d), 0, sqrt(0.8))
  d$y[sample.int(nrow(d), n_missing)] <- NA
  d
}

test_that("a missing response does not break the default aggregate = auto", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .rcbd_with_hole()
  expect_no_error(suppressMessages(flexybayes(
    fixed = y ~ 1, random = ~ blk, residual = ~ units,
    data = d, backend = "inla", na_action = "augment", verbose = FALSE
  )))
})

test_that("auto and explicit aggregate = FALSE agree when a response is missing", {
  # Falling through must give the per-row answer, not a third thing.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .rcbd_with_hole()
  f_auto <- suppressMessages(flexybayes(
    fixed = y ~ 1, random = ~ blk, residual = ~ units, data = d,
    backend = "inla", na_action = "augment", verbose = FALSE
  ))
  f_row <- suppressMessages(flexybayes(
    fixed = y ~ 1, random = ~ blk, residual = ~ units, data = d,
    backend = "inla", na_action = "augment", aggregate = FALSE,
    verbose = FALSE
  ))
  h1 <- f_auto$inla$summary.hyperpar[, "mean"]
  h2 <- f_row$inla$summary.hyperpar[, "mean"]
  # Both calls take the same route, so this is testing that the fallthrough
  # lands on the per-row path rather than a third behaviour. The tolerance is
  # not tight because INLA is not bitwise reproducible across runs under
  # threading -- the same call here differs in the fourth decimal -- so a
  # 1e-6 expectation would fail on the engine's own noise rather than on
  # anything this test is about.
  expect_equal(unname(h1), unname(h2), tolerance = 1e-2)
  expect_identical(
    f_auto$extras$backend_decision$path,
    f_row$extras$backend_decision$path
  )
})

test_that("aggregate = TRUE refuses rather than compressing a missing response", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .rcbd_with_hole()
  expect_error(
    suppressMessages(flexybayes(
      fixed = y ~ 1, random = ~ blk, residual = ~ units, data = d,
      backend = "inla", na_action = "augment", aggregate = TRUE,
      verbose = FALSE
    )),
    "missing values"
  )
})

test_that("the aggregated path is untouched when nothing is missing", {
  # The guard must not close a route that was working: a complete response
  # still aggregates, and the fit still reports the aggregated path.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .rcbd_with_hole(n_missing = 0L)
  fit <- suppressMessages(flexybayes(
    fixed = y ~ 1, random = ~ blk, residual = ~ units,
    data = d, backend = "inla", verbose = FALSE
  ))
  path <- fit$extras$backend_decision$path %||% ""
  expect_match(paste(path, collapse = " "), "aggregat")
})

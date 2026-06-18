# Regression tests for the per-row INLA post-fit method surface.
#
# Guards two bugs found by the 2026-06-02 corner-to-corner stress run:
#   * fitted()/residuals() on an INLA fit dispatched to *.default and
#     silently returned NULL (the object carries no $glm slot). They now
#     read INLA's posterior-mean fitted values.
#   * logLik() had no flexybayes_inla method, so it (and glance()) errored
#     with "no applicable method". It now returns an honest NA logLik
#     rather than fabricating INLA's marginal likelihood as a conditional one.

skip_if_no_inla <- function() skip_if_not_installed("INLA")

test_that("fitted() returns posterior-mean fitted values for an INLA fit", {
  skip_if_no_inla()
  set.seed(404L)
  d <- data.frame(x = rnorm(50))
  d$y <- 1 + 0.7 * d$x + rnorm(50, 0, 0.6)
  fit <- fb_inla(y ~ x, data = d)

  fv <- fitted(fit)
  expect_type(fv, "double")
  expect_length(fv, nrow(d))
  expect_true(all(is.finite(fv)))
  # not the silent-NULL of the pre-fix default dispatch
  expect_false(is.null(fv))
})

test_that("residuals() are response residuals that reconstruct y", {
  skip_if_no_inla()
  set.seed(405L)
  d <- data.frame(x = rnorm(50))
  d$y <- 2 - 0.5 * d$x + rnorm(50, 0, 0.5)
  fit <- fb_inla(y ~ x, data = d)

  rv <- residuals(fit)
  expect_length(rv, nrow(d))
  expect_true(all(is.finite(rv)))
  # observed = fitted + residual, exactly
  expect_equal(fitted(fit) + rv, d$y, tolerance = 1e-8)
})

test_that("fitted() on a non-identity (Poisson) INLA fit is on the response scale", {
  skip_if_no_inla()
  set.seed(406L)
  d <- data.frame(x = rnorm(60))
  d$y <- rpois(60, exp(0.5 + 0.3 * d$x))
  fit <- fb_inla(y ~ x, data = d, family = "poisson")

  fv <- fitted(fit)
  expect_length(fv, nrow(d))
  expect_true(all(fv > 0))   # response-scale rates, strictly positive
})

test_that("logLik() on an INLA fit returns an honest NA, not a fabricated value", {
  skip_if_no_inla()
  set.seed(407L)
  d <- data.frame(x = rnorm(40))
  d$y <- 0.3 + d$x + rnorm(40, 0, 0.7)
  fit <- fb_inla(y ~ x, data = d)

  expect_message(ll <- logLik(fit), "marginal")
  expect_s3_class(ll, "logLik")
  expect_true(is.na(as.numeric(ll)))
  expect_equal(attr(ll, "df"), 2L)      # intercept + slope
  expect_equal(attr(ll, "nobs"), nrow(d))
})

# =============================================================================
# Bayesian labels on the standard-R compatibility view.
#
# summary.flexybayes_glm() computed pd = Phi(|beta/sd|) correctly, then
# reported 2 * (1 - pd) in a column headed `Pr(>|z|)`, and a footnote called
# THAT the "approximate posterior probability of direction". It is a doubled
# posterior sign-tail area formatted to look like a frequentist p-value:
# for a coefficient with pd = 0.6, the old column read 0.8.
#
# confint.flexybayes_glm() was documented as quantile-based and returned
# beta +/- z * sd, a normal approximation. The two coincide only for a
# symmetric posterior.
#
# These are value-and-label tests. The previous tests asserted on the shape
# of the returned matrix, which was correct throughout.
# =============================================================================

# A flexybayes_glm carrying known draws, so every reported quantity has an
# exactly computable target. Column 2 is deliberately centred near zero so
# pd is far from 1 and the old and new quantities diverge visibly.
.mk_glm <- function(n_draws = 4000L, seed = 5L, with_draws = TRUE) {
  set.seed(seed)
  draws <- cbind(
    `(Intercept)` = stats::rnorm(n_draws, mean = 2, sd = 0.5),
    x = stats::rnorm(n_draws, mean = 0.15, sd = 0.5)
  )
  structure(
    list(
      coefficients = colMeans(draws),
      family = stats::gaussian(),
      df.residual = 50L
    ),
    class = c("flexybayes_glm", "list"),
    posterior_vcov = stats::cov(draws),
    posterior_draws = if (with_draws) draws else NULL
  )
}

test_that("summary() no longer reports a doubled sign-tail area as Pr(>|z|)", {
  g <- .mk_glm()
  s <- suppressMessages(summary(g))
  expect_false("Pr(>|z|)" %in% colnames(s))
  expect_false("z value" %in% colnames(s))
  expect_true("pd" %in% colnames(s))
})

test_that("pd is the probability of direction, counted from the draws", {
  g <- .mk_glm()
  draws <- attr(g, "posterior_draws")
  s <- suppressMessages(summary(g))
  pd_reported <- s[, "pd"]
  pd_target <- pmax(colMeans(draws > 0), colMeans(draws < 0))
  expect_equal(unname(pd_reported), unname(pd_target))

  # pd is a probability of direction, so it cannot fall below one half.
  expect_true(all(pd_reported >= 0.5))
  expect_true(all(pd_reported <= 1))

  # And it is materially different from what the old column reported for
  # a coefficient whose posterior straddles zero.
  se <- sqrt(diag(attr(g, "posterior_vcov")))
  old_column <- 2 * (1 - stats::pnorm(abs(g$coefficients / se)))
  expect_gt(abs(pd_reported[["x"]] - old_column[["x"]]), 0.1)
})

test_that("confint() returns posterior quantiles when draws are present", {
  g <- .mk_glm()
  draws <- attr(g, "posterior_draws")
  ci <- confint(g)
  target <- t(apply(draws, 2L, stats::quantile, probs = c(0.025, 0.975)))
  # Compare the numeric content: `ci` carries an interval_basis attribute
  # that unname() does not strip and expect_equal() would otherwise flag.
  expect_equal(as.numeric(ci), as.numeric(target))
  expect_identical(rownames(ci), colnames(draws))
  expect_identical(attr(ci, "interval_basis"), "posterior_quantile")
})

test_that("a skewed posterior separates the quantile interval from the normal one", {
  # The guard that matters: on a symmetric posterior the old normal
  # approximation and the quantile interval agree, so a test built on one
  # would pass either way. A skewed posterior tells them apart.
  set.seed(9L)
  draws <- cbind(b = stats::rlnorm(6000L, meanlog = 0, sdlog = 0.9))
  g <- structure(
    list(coefficients = c(b = mean(draws[, 1L])),
         family = stats::gaussian(), df.residual = 30L),
    class = c("flexybayes_glm", "list"),
    posterior_vcov = stats::cov(draws),
    posterior_draws = draws
  )
  ci <- confint(g)
  se <- sqrt(diag(attr(g, "posterior_vcov")))
  normal_lower <- g$coefficients - stats::qnorm(0.975) * se

  # The lognormal posterior is positive everywhere, so its lower quantile
  # must be positive; the normal approximation runs below zero.
  expect_gt(ci[1L, 1L], 0)
  expect_lt(normal_lower, 0)
})

test_that("without draws both methods fall back and say so", {
  g <- .mk_glm(with_draws = FALSE)
  s <- suppressMessages(summary(g))
  expect_true("pd_normal_approx" %in% colnames(s))
  expect_false("pd" %in% colnames(s))

  ci <- confint(g)
  expect_identical(attr(ci, "interval_basis"), "normal_approximation")
  se <- sqrt(diag(attr(g, "posterior_vcov")))
  z <- stats::qnorm(0.975)
  expect_equal(
    as.numeric(ci),
    as.numeric(cbind(g$coefficients - z * se, g$coefficients + z * se))
  )
})

test_that("confint honours the level argument", {
  g <- .mk_glm()
  draws <- attr(g, "posterior_draws")
  ci <- confint(g, level = 0.5)
  target <- t(apply(draws, 2L, stats::quantile, probs = c(0.25, 0.75)))
  expect_equal(as.numeric(ci), as.numeric(target))
  expect_identical(colnames(ci), c("25%", "75%"))
})

test_that("a real brms fit carries its fixed-effect draws through", {
  # Integration check: the live path that actually produces a
  # flexybayes_glm must attach the draws, or every test above is
  # exercising a fixture the package never builds.
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(42L)
  n <- 40L
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1.5 + 0.8 * d$x + stats::rnorm(n)
  fit <- suppressMessages(flexybayes(
    y ~ x, data = d, backend = "brms",
    n_samples = 400, warmup = 400, chains = 2,
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  draws <- attr(fit$glm, "posterior_draws")
  expect_true(is.matrix(as.matrix(draws)))
  expect_identical(colnames(draws), names(fit$glm$coefficients))
  expect_identical(attr(confint(fit$glm), "interval_basis"),
                   "posterior_quantile")
})

# vcov() on an INLA fit: repeatable, and equal to the covariance of INLA's
# own posterior draws of the fixed effects.

.vcov_inla_data <- function() {
  set.seed(1)
  d <- data.frame(
    site = factor(rep(1:6, each = 12)),
    trt = factor(rep(c("a", "b", "c"), 24)),
    x = stats::rnorm(72)
  )
  d$y <- 1 +
    c(0, 0.5, -0.3)[d$trt] +
    0.4 * d$x +
    stats::rnorm(6)[d$site] +
    stats::rnorm(72)
  d$count <- stats::rpois(
    72,
    exp(0.5 + 0.3 * d$x + 0.5 * stats::rnorm(6)[d$site])
  )
  d
}

test_that("DEF-005: two vcov() calls on one INLA fit return the same matrix", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )

  fit <- flexybayes(
    y ~ trt + x,
    random = ~site,
    data = .vcov_inla_data(),
    backend = "inla",
    verbose = FALSE
  )

  set.seed(42)
  rng_before <- .Random.seed
  v1 <- vcov(fit)
  expect_identical(.Random.seed, rng_before)
  v2 <- vcov(fit)
  expect_identical(v1, v2)
  expect_identical(vcov(fit, n_samples = 50L), v1)
  expect_identical(dimnames(v1), list(names(coef(fit)), names(coef(fit))))
  expect_true(isSymmetric(v1))
})

test_that("DEF-005: vcov() on an INLA fit matches INLA's posterior draws", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )

  d <- .vcov_inla_data()
  fit <- flexybayes(
    y ~ trt + x,
    random = ~site,
    data = d,
    backend = "inla",
    verbose = FALSE
  )
  v <- vcov(fit)
  sf <- fit$inla$summary.fixed

  # Slope and contrast SDs from INLA's marginals agree to 1e-4 on this fit.
  slopes <- rownames(sf) != "(Intercept)"
  expect_equal(
    sqrt(diag(v))[slopes],
    sf$sd[slopes],
    tolerance = 0.01,
    ignore_attr = TRUE
  )

  # Reference: 2,000 draws straight from INLA's sampler. Over 20 repeats of
  # 4,000 draws the largest gap was 0.049 in a correlation and 7.7% in a
  # variance, or 0.07 and 11% scaled to 2,000 draws.
  s <- INLA::inla.posterior.sample(2000L, fit$inla, num.threads = "1:1")
  rows <- paste0(rownames(sf), ":1")
  draws <- t(vapply(s, function(z) z$latent[rows, 1L], numeric(length(rows))))
  scale <- fit$extras$response_standardisation$scale
  expect_lt(max(abs(stats::cor(draws) - stats::cov2cor(v))), 0.1)
  expect_lt(
    max(abs(diag(stats::cov(draws)) * scale^2 / diag(v) - 1)),
    0.15
  )

  # A Poisson fit standardises the covariate x. Its SD on the original
  # scale agrees with INLA's marginal to 3e-4.
  fit_pois <- flexybayes(
    count ~ trt + x,
    random = ~site,
    data = d,
    family = stats::poisson(),
    backend = "inla",
    verbose = FALSE
  )
  expect_false(is.null(fit_pois$extras$covariate_standardisation))
  expect_equal(
    sqrt(vcov(fit_pois)["x", "x"]),
    fit_pois$inla$summary.fixed["x", "sd"],
    tolerance = 0.01
  )
})

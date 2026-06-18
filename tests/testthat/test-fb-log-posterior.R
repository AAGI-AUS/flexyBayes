# Tests for the C4 producer: fb_log_posterior().
#
# The greta path is the real producer and needs a working greta /
# TensorFlow stack, so its tests are gated by skip_if_no_greta(). The
# abstain paths (brms / INLA / default) and the input-contract guards need
# no backend and run everywhere.

# ---- abstain paths (no backend required) ------------------------------------

test_that("the brms backend abstains with a classed, informative condition", {
  brms_fit <- structure(
    list(),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
  cnd <- tryCatch(fb_log_posterior(brms_fit), error = function(e) e)
  expect_s3_class(cnd, "fb_c4_unavailable")
  expect_match(conditionMessage(cnd), "brms backend")
})

test_that("the INLA backend abstains with a classed, informative condition", {
  inla_fit <- structure(list(), class = c("flexybayes_inla", "list"))
  cnd <- tryCatch(fb_log_posterior(inla_fit), error = function(e) e)
  expect_s3_class(cnd, "fb_c4_unavailable")
  expect_match(conditionMessage(cnd), "INLA")
})

test_that("a non-flexyBayes object abstains via the default method", {
  cnd <- tryCatch(fb_log_posterior(list(a = 1)), error = function(e) e)
  expect_s3_class(cnd, "fb_c4_unavailable")
})

test_that("a greta-classed fit without a retained model abstains", {
  ## A greta-classed fit that lost its model graph must abstain, not crash.
  no_model <- structure(
    list(greta = list(model = NULL)),
    class = c("flexybayes", "list")
  )
  cnd <- tryCatch(fb_log_posterior(no_model), error = function(e) e)
  expect_s3_class(cnd, "fb_c4_unavailable")
})

# ---- transform helper (pure, no backend) ------------------------------------

test_that("the natural -> free transform mirrors greta's default bijectors", {
  ## Unbounded: identity.
  expect_equal(
    .fb_natural_to_free(c(-2, 0, 3), -Inf, Inf),
    c(-2, 0, 3)
  )
  ## Lower-bounded at 0: log(x). Out-of-support -> NaN.
  expect_equal(
    suppressWarnings(.fb_natural_to_free(c(1, exp(1), 2), 0, Inf)),
    c(0, 1, log(2))
  )
  expect_true(is.nan(suppressWarnings(.fb_natural_to_free(-1, 0, Inf))))
  expect_true(is.nan(suppressWarnings(.fb_natural_to_free(0, 0, Inf))))
  ## Upper-bounded: log(upper - x).
  expect_equal(
    suppressWarnings(.fb_natural_to_free(c(0, 0.5), -Inf, 1)),
    c(log(1), log(0.5))
  )
  ## Both finite: logit on (lower, upper).
  expect_equal(
    suppressWarnings(.fb_natural_to_free(0.5, 0, 1)),
    stats::qlogis(0.5)
  )
  expect_true(is.nan(suppressWarnings(.fb_natural_to_free(1.5, 0, 1))))
})

# ---- greta path: the real producer ------------------------------------------

test_that("the greta producer matches the analytic conjugate log-posterior", {
  skip_if_no_greta()
  suppressMessages(library(greta))
  withr::local_seed(11L)

  ## Conjugate Gaussian mean with known variance: the joint
  ## log p(y | mu) + log p(mu) is analytic, an independent oracle the
  ## producer's log-density must match up to an additive constant.
  n <- 40L
  s0 <- 2
  m0 <- 0
  t0 <- 5
  y_obs <- stats::rnorm(n, 1.3, s0)
  mu <- greta::normal(m0, t0)
  yd <- greta::as_data(y_obs)
  greta::distribution(yd) <- greta::normal(mu, s0)
  m <- greta::model(mu)
  fit <- suppressMessages(suppressWarnings(fb_greta(
    fb_from_greta(m, canonical_names = c(mu = "mu")),
    n_samples = 300L, warmup = 300L, chains = 1L,
    verbose = FALSE, mcmc_verbose = FALSE
  )))

  producer <- fb_log_posterior(fit)
  expect_s3_class(producer, "fb_log_posterior_producer")
  expect_true(is.function(producer))
  expect_equal(attr(producer, "parameter_names"), "mu")
  ## Marginal likelihood unknown for a posterior -> NA, honestly.
  expect_true(is.na(attr(producer, "log_normalizer")))
  ## mu is unbounded -> NA support both ends.
  expect_true(is.na(attr(producer, "support_lower")))
  expect_true(is.na(attr(producer, "support_upper")))
  ## Draws supplied to seed the proposal.
  expect_true(is.matrix(attr(producer, "draws")))
  expect_equal(ncol(attr(producer, "draws")), 1L)

  joint_R <- function(v) {
    sum(stats::dnorm(y_obs, v, s0, log = TRUE)) +
      stats::dnorm(v, m0, t0, log = TRUE)
  }
  grid <- seq(0, 3, by = 0.25)
  pv <- producer(matrix(grid, ncol = 1L))
  av <- vapply(grid, joint_R, numeric(1L))

  ## Up to an additive constant: correlation 1, constant offset.
  expect_equal(stats::cor(pv, av), 1, tolerance = 1e-8)
  expect_lt(stats::sd(pv - av), 1e-8)
})

test_that("the greta producer is vectorised and domain-safe on constraints", {
  skip_if_no_greta()
  suppressMessages(library(greta))
  withr::local_seed(21L)

  ## mu unbounded + sigma > 0 (a truncated-normal prior -> log transform):
  ## exercises the constrained path and the natural-scale truncation
  ## normaliser, again against an analytic oracle.
  n <- 50L
  y_obs <- stats::rnorm(n, 1.0, 1.5)
  mu <- greta::normal(0, 5)
  sigma <- greta::normal(0, 5, truncation = c(0, Inf))
  yd <- greta::as_data(y_obs)
  greta::distribution(yd) <- greta::normal(mu, sigma)
  m <- greta::model(mu, sigma)
  fit <- suppressMessages(suppressWarnings(fb_greta(
    fb_from_greta(m, canonical_names = c(mu = "mu", sigma = "sigma")),
    n_samples = 300L, warmup = 300L, chains = 1L,
    verbose = FALSE, mcmc_verbose = FALSE
  )))

  producer <- fb_log_posterior(fit)
  expect_equal(attr(producer, "parameter_names"), c("mu", "sigma"))
  ## sigma is lower-bounded at 0; mu is unbounded.
  expect_equal(attr(producer, "support_lower"), c(NA_real_, 0))
  expect_equal(attr(producer, "support_upper"), c(NA_real_, NA_real_))

  joint_R <- function(muv, sigv) {
    if (sigv <= 0) {
      return(-Inf)
    }
    sum(stats::dnorm(y_obs, muv, sigv, log = TRUE)) +
      stats::dnorm(muv, 0, 5, log = TRUE) +
      (stats::dnorm(sigv, 0, 5, log = TRUE) -
         log(1 - stats::pnorm(0, 0, 5)))
  }
  gr <- expand.grid(
    mu = seq(0.4, 1.6, 0.4),
    sigma = seq(0.8, 2.2, 0.4)
  )
  pv <- producer(as.matrix(gr))
  av <- mapply(joint_R, gr$mu, gr$sigma)
  expect_equal(stats::cor(pv, av), 1, tolerance = 1e-8)
  expect_lt(stats::sd(pv - av), 1e-8)

  ## Domain safety: sigma <= 0 returns -Inf, never an error, and a mixed
  ## batch keeps the valid row finite.
  expect_identical(producer(matrix(c(1.0, -1.0), nrow = 1L)), -Inf)
  expect_identical(producer(matrix(c(1.0, 0.0), nrow = 1L)), -Inf)
  mixed <- producer(matrix(c(1.0, 1.5, 1.0, -1.0), nrow = 2L, byrow = TRUE))
  expect_true(is.finite(mixed[1L]))
  expect_identical(mixed[2L], -Inf)

  ## Vectorisation: a batch equals row-by-row evaluation.
  rows <- matrix(c(1.0, 1.5, 0.5, 2.0, 1.2, 1.1), nrow = 3L, byrow = TRUE)
  batch <- producer(rows)
  single <- vapply(
    seq_len(nrow(rows)),
    function(i) producer(rows[i, , drop = FALSE]),
    numeric(1L)
  )
  expect_equal(batch, single)

  ## The 2-row zero probe proxymix runs at construction must not error.
  expect_silent(producer(matrix(0, nrow = 2L, ncol = 2L)))
})

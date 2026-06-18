# Tests for fb_brms() -- ADR 0014.
#
# Contract:
#   - Accepts the five-shape brms corpus: fixed-only Gaussian;
#     fixed + one RI; fixed + crossed RIs; binomial RI (single-
#     column Bernoulli); Poisson RI.
#   - Refuses everything outside the corpus at ingest with a
#     structured message (random slopes, smoothers, GP, autocor,
#     LHS-pipe addition forms).
#   - Inherits ADR 0006 (backend = c("greta", "inla", "auto") +
#     decision trace), ADR 0011 (review_code), ADR 0005 (canonical-
#     name auto-resolution via the existing greta + INLA mappers).
#   - Returns class "flexybayes" (greta or auto-fell-back) or
#     "flexybayes_inla" (INLA accepted); per ADR 0014 sec.7 no
#     subclass distinguishes brms-entry from asreml-entry fits.
#
# Greta is required for the greta-path fits; INLA for the auto-
# accept + cross-engine subtests (skip_if_not_installed guards).

# ---------------------------------------------------------------- #
# Test fixtures                                                    #
# ---------------------------------------------------------------- #

mk_brms_gaussian_data <- function() {
  set.seed(20260523L)
  n <- 50L
  g1 <- factor(rep(letters[1:5], each = 10L))
  g2 <- factor(rep(LETTERS[1:5], times = 10L))
  x <- rnorm(n)
  b_g1 <- rnorm(5L, sd = 1.5)[as.integer(g1)]
  b_g2 <- rnorm(5L, sd = 0.8)[as.integer(g2)]
  y <- 1 + 0.5 * x + b_g1 + b_g2 + rnorm(n, sd = 0.3)
  data.frame(y = y, x = x, g1 = g1, g2 = g2)
}

mk_brms_binomial_data <- function() {
  set.seed(20260523L)
  n <- 60L
  g <- factor(rep(letters[1:6], each = 10L))
  x <- rnorm(n)
  eta <- 0.3 + 0.5 * x + rnorm(6L, sd = 0.5)[as.integer(g)]
  y <- rbinom(n, size = 1L, prob = plogis(eta))
  data.frame(y = y, x = x, g = g)
}

mk_brms_poisson_data <- function() {
  set.seed(20260523L)
  n <- 60L
  g <- factor(rep(letters[1:6], each = 10L))
  x <- rnorm(n)
  eta <- 0.2 + 0.4 * x + rnorm(6L, sd = 0.3)[as.integer(g)]
  y <- rpois(n, lambda = exp(eta))
  data.frame(y = y, x = x, g = g)
}


# ---------------------------------------------------------------- #
# (1) Corpus shape 1: fixed-only Gaussian                          #
# ---------------------------------------------------------------- #

test_that("fb() corpus 1: fixed-only Gaussian (brms grammar forced)", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_brms_gaussian_data()
  # A bar-free formula is ASReml grammar by default; `syntax = "brms"`
  # forces the brms ingest so this still exercises the brms corpus path.
  fit <- suppressMessages(fb(
    y ~ x,
    data = d,
    backend = "greta",
    syntax = "brms",
    n_samples = 60L,
    warmup = 60L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
  expect_identical(fit$extras$fb_terms$source, "brms")
  expect_true("(Intercept)" %in% names(coef(fit)))
  expect_true("x" %in% names(coef(fit)))
})


# ---------------------------------------------------------------- #
# (2) Corpus shape 2: Gaussian + one random intercept              #
# ---------------------------------------------------------------- #

test_that("fb_brms() corpus 2: Gaussian + one random intercept", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_brms_gaussian_data()
  fit <- suppressMessages(fb(
    y ~ x + (1 | g1),
    data = d,
    backend = "greta",
    n_samples = 60L,
    warmup = 60L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
  expect_identical(fit$extras$fb_terms$source, "brms")
  expect_length(fit$extras$fb_terms$random_terms, 1L)
  expect_identical(fit$extras$fb_terms$random_terms[[1L]]$var, "g1")
})


# ---------------------------------------------------------------- #
# (3) Corpus shape 3: Gaussian + crossed random intercepts         #
# ---------------------------------------------------------------- #

test_that("fb_brms() corpus 3: crossed random intercepts", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_brms_gaussian_data()
  fit <- suppressMessages(fb(
    y ~ x + (1 | g1) + (1 | g2),
    data = d,
    backend = "greta",
    n_samples = 60L,
    warmup = 60L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
  expect_length(fit$extras$fb_terms$random_terms, 2L)
  rg <- vapply(
    fit$extras$fb_terms$random_terms,
    function(t) t$var,
    character(1)
  )
  expect_setequal(rg, c("g1", "g2"))
})


# ---------------------------------------------------------------- #
# (4) Corpus shape 4: Binomial RI (single-column Bernoulli)        #
# ---------------------------------------------------------------- #

test_that("fb_brms() corpus 4: binomial RI as single-column Bernoulli", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_brms_binomial_data()
  fit <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = d,
    family = "binomial",
    backend = "greta",
    n_samples = 60L,
    warmup = 60L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
  expect_identical(fit$extras$fb_terms$family, "binomial")
  expect_identical(fit$extras$fb_terms$link, "logit")
})


# ---------------------------------------------------------------- #
# (5) Corpus shape 5: Poisson RI                                    #
# ---------------------------------------------------------------- #

test_that("fb_brms() corpus 5: Poisson random-intercept GLMM", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_brms_poisson_data()
  fit <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = d,
    family = "poisson",
    backend = "greta",
    n_samples = 60L,
    warmup = 60L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
  expect_identical(fit$extras$fb_terms$family, "poisson")
  expect_identical(fit$extras$fb_terms$link, "log")
})


# ---------------------------------------------------------------- #
# (6) ADR 0020 refusal: correlated random slopes (x | g)            #
# ---------------------------------------------------------------- #
#
# Uncorrelated random slopes (x || g) are now supported per ADR 0020
# (covered by test-random-slopes-uncor.R). The correlated form
# (x | g) continues to refuse, but with the new structured
# `flexybayes_correlated_slope_unsupported` condition carrying
# precise slots (deferral_target = "v0.3"; workaround = "(x || g)";
# grouping_factor; slope_variable).

test_that("fb_brms() refuses correlated random slopes (x | g) with ADR 0020 structured condition", {
  d <- mk_brms_gaussian_data()
  err <- tryCatch(
    fb(y ~ x + (x | g1), data = d, verbose = FALSE),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_correlated_slope_unsupported")
  expect_identical(err$deferral_target, "a future release")
  expect_identical(err$workaround, "(x || g)")
  expect_identical(err$grouping_factor, "g1")
  expect_identical(err$slope_variable, "x")
  msg <- conditionMessage(err)
  expect_true(grepl("Correlated random slopes", msg))
  expect_true(grepl("future release", msg))
  expect_true(grepl("\\(x \\|\\| g\\)", msg))
})


# ---------------------------------------------------------------- #
# (7) Non-corpus refusal: smoother                                  #
# ---------------------------------------------------------------- #

test_that("fb_brms() refuses brms smoothers with structured message", {
  d <- mk_brms_gaussian_data()
  err <- tryCatch(
    fb(y ~ s(x) + (1 | g1), data = d, verbose = FALSE),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("does not yet support", err))
  expect_true(grepl("smoother|s\\(", err))
})


# ---------------------------------------------------------------- #
# (8) Non-corpus refusal: LHS-pipe addition form                   #
# ---------------------------------------------------------------- #

test_that("fb_brms() refuses LHS-pipe addition forms cleanly", {
  d <- mk_brms_gaussian_data()
  err <- tryCatch(
    fb(y | weights(x) ~ x + (1 | g1), data = d, verbose = FALSE),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("LHS addition forms", err))
  expect_true(grepl("future expansion|brms ingest", err))
})


# ---------------------------------------------------------------- #
# (9) backend = "auto" routes to INLA on accept                    #
# ---------------------------------------------------------------- #

test_that("fb(..., backend = 'auto') routes to INLA when LGM-feasible", {
  testthat::skip_if_not_installed("INLA")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_brms_gaussian_data()
  fit <- suppressMessages(fb(
    y ~ x + (1 | g1),
    data = d,
    backend = "auto",
    verbose = FALSE
  ))
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "inla")
  expect_identical(bd$path, "auto_accept")
  expect_s3_class(fit, "flexybayes_inla")
})


# ---------------------------------------------------------------- #
# (10) review_code = TRUE returns the deferred token; proceed()    #
# ---------------------------------------------------------------- #

test_that("fb(..., review_code = TRUE) returns <flexybayes_review>; proceed() fits", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_brms_gaussian_data()
  rev <- suppressMessages(fb(
    y ~ x + (1 | g1),
    data = d,
    n_samples = 50L,
    warmup = 50L,
    chains = 1L,
    review_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(rev, "flexybayes_review")
  expect_true(is.character(rev$code) && nchar(rev$code) > 0L)
  expect_identical(rev$backend, "greta")
  expect_identical(rev$ir$source, "brms")
  fit <- suppressMessages(proceed(rev))
  expect_s3_class(fit, "flexybayes")
  expect_true("(Intercept)" %in% names(coef(fit)))
  # Second proceed() returns the cached fit.
  fit2 <- proceed(rev)
  expect_identical(coef(fit), coef(fit2))
})


# ---------------------------------------------------------------- #
# (11) backend_decision(fit) uniform shape across paths            #
# ---------------------------------------------------------------- #

test_that("backend_decision() returns uniform-shape trace on fb_brms() fits", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_brms_gaussian_data()
  fit <- suppressMessages(fb(
    y ~ x + (1 | g1),
    data = d,
    backend = "greta",
    n_samples = 50L,
    warmup = 50L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  bd <- backend_decision(fit)
  expect_setequal(
    names(bd),
    c(
      "backend",
      "path",
      "gate_checks",
      "reason",
      "preflight_summary",
      "representation_plan",
      "rejected_routes",
      "routing_policy_version"
    )
  )
  expect_identical(bd$backend, "greta")
  expect_identical(bd$path, "explicit_greta")
  expect_null(bd$gate_checks)
  # ADR 0024 v0.3.6+ four new fields: NULL on the small-data fast
  # path (no preflight, no representation plan); empty list of
  # rejected routes on an explicit user request (policy-bypass
  # semantics); routing_policy_version is the live constant.
  expect_null(bd$preflight_summary)
  expect_null(bd$representation_plan)
  expect_identical(bd$rejected_routes, list())
  expect_type(bd$routing_policy_version, "character")
  expect_length(bd$routing_policy_version, 1L)
})


# ---------------------------------------------------------------- #
# (12) review_code under backend != "greta" raises clean refusal   #
# ---------------------------------------------------------------- #

test_that("fb(review_code = TRUE, backend != 'greta') raises structured refusal", {
  d <- mk_brms_gaussian_data()
  err <- tryCatch(
    fb(
      y ~ x + (1 | g1),
      data = d,
      backend = "inla",
      review_code = TRUE,
      verbose = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("review_code", err, fixed = TRUE))
  expect_true(grepl("future ADR|greta", err))
})


# ---------------------------------------------------------------- #
# (13) triangulate(fit_greta, fit_inla) auto-resolves canonical    #
#      names via the ADR 0005 registry (no name_map supplied)      #
# ---------------------------------------------------------------- #

test_that("triangulate() on fb_brms() greta + INLA fits resolves canonical names automatically", {
  skip_if_greta_backend_unusable()
  testthat::skip_if_not_installed("INLA")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_brms_gaussian_data()
  fit_g <- suppressMessages(fb(
    y ~ x + (1 | g1),
    data = d,
    backend = "greta",
    n_samples = 80L,
    warmup = 80L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fit_i <- suppressMessages(fb(
    y ~ x + (1 | g1),
    data = d,
    backend = "inla",
    verbose = FALSE
  ))
  tri <- triangulate(fit_g, fit_i)
  expect_true(all(
    c("(Intercept)", "x", "sd_g1", "sigma") %in%
      tri$common
  ))
})


# ---------------------------------------------------------------- #
# (14) Invalid backend value raises match.arg error                #
# ---------------------------------------------------------------- #

test_that("fb_brms() rejects invalid backend value", {
  d <- mk_brms_gaussian_data()
  err <- tryCatch(
    fb(y ~ x + (1 | g1), data = d, backend = "stan", verbose = FALSE),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("'arg'.*should be one of|match.arg", err, perl = TRUE))
})

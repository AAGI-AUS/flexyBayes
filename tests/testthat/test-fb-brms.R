# Tests for fb_brms().
#
# Contract:
#   - Accepts the five-shape brms corpus: fixed-only Gaussian;
#     fixed + one RI; fixed + crossed RIs; binomial RI (single-
#     column Bernoulli); Poisson RI.
#   - Refuses everything outside the corpus at ingest with a
#     structured message (random slopes, smoothers, GP, autocor,
#     LHS-pipe addition forms).
#   - Takes the `backend` choice (c("inla", "brms", "auto")) with its
#     decision trace, `review_code`, and canonical-name auto-resolution
#     via the existing brms + INLA mappers.
#   - Returns class "flexybayes_brms" (brms accepted) or
#     "flexybayes_inla" (INLA accepted); no subclass
#     distinguishes brms-entry from asreml-entry fits.
#
# brms is required for the brms-path fits; INLA for the auto-
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
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_not_long_tier()
  d <- mk_brms_gaussian_data()
  # A bar-free formula is ASReml grammar by default; `syntax = "brms"`
  # forces the brms ingest so this still exercises the brms corpus path.
  fit <- suppressMessages(flexybayes(
    y ~ x,
    data = d,
    backend = "brms",
    syntax = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
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
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_not_long_tier()
  d <- mk_brms_gaussian_data()
  # The assertions are structural, so the low-ESS warning of the one-chain
  # 500-draw budget is muffled (see helper-ess-warnings.R).
  fit <- .muffle_ess_warnings(suppressMessages(flexybayes(
    y ~ x + (1 | g1),
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  expect_s3_class(fit, "flexybayes")
  expect_identical(fit$extras$fb_terms$source, "brms")
  expect_length(fit$extras$fb_terms$random_terms, 1L)
  expect_identical(fit$extras$fb_terms$random_terms[[1L]]$var, "g1")
})


# ---------------------------------------------------------------- #
# (3) Corpus shape 3: Gaussian + crossed random intercepts         #
# ---------------------------------------------------------------- #

test_that("fb_brms() corpus 3: crossed random intercepts", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_not_long_tier()
  d <- mk_brms_gaussian_data()
  fit <- suppressMessages(flexybayes(
    y ~ x + (1 | g1) + (1 | g2),
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
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
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_not_long_tier()
  d <- mk_brms_binomial_data()
  fit <- suppressMessages(flexybayes(
    y ~ x + (1 | g),
    data = d,
    family = "binomial",
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
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
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_not_long_tier()
  d <- mk_brms_poisson_data()
  # The assertions are structural, so the low-ESS warning of the one-chain
  # 500-draw budget is muffled (see helper-ess-warnings.R).
  fit <- .muffle_ess_warnings(suppressMessages(flexybayes(
    y ~ x + (1 | g),
    data = d,
    family = "poisson",
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  expect_s3_class(fit, "flexybayes")
  expect_identical(fit$extras$fb_terms$family, "poisson")
  expect_identical(fit$extras$fb_terms$link, "log")
})


# ---------------------------------------------------------------- #
# (6) Correlated random slopes (x | g) through the fb() surface      #
# ---------------------------------------------------------------- #
#
# Both bar forms reach the engines: (x || g) is two independent
# variances, (x | g) the same two with the correlation between them
# freed. The single pipe used to refuse with a structured deferral.

test_that("fb() carries correlated random slopes (x | g) to the brms formula", {
  d <- mk_brms_gaussian_data()
  fb_ir <- fb_from_brms(y ~ x + (x | g1), data = d)
  expect_identical(fb_ir$random_terms[[1]]$type, "simple_slope_cor")
  expect_identical(fb_ir$random_terms[[1]]$var, "g1")
  expect_identical(fb_ir$random_terms[[1]]$slope_var, "x")
  txt <- gsub(
    "\\s+",
    " ",
    paste(deparse(.fb_to_brms_formula(fb_ir)), collapse = " ")
  )
  expect_match(txt, "(x | g1)", fixed = TRUE)
})


# ---------------------------------------------------------------- #
# (7) Non-corpus refusal: smoother                                  #
# ---------------------------------------------------------------- #

test_that("fb_brms() refuses brms smoothers with structured message", {
  d <- mk_brms_gaussian_data()
  # Typed since 0.9.2: the smoother used to reach the
  # engine and surface as `could not find function "s"`, on brms only
  # after a completed sampling run.
  err <- expect_error(
    flexybayes(y ~ s(x) + (1 | g1), data = d, verbose = FALSE),
    class = "flexybayes_refusal_fixed_smoother_not_supported"
  )
  expect_match(conditionMessage(err), "s\\(x\\)")
  expect_match(conditionMessage(err), "spl\\(x\\)")
  expect_match(conditionMessage(err), "backend = \"inla\"")
})


# ---------------------------------------------------------------- #
# (8) Non-corpus refusal: LHS-pipe addition form                   #
# ---------------------------------------------------------------- #

test_that("fb_brms() refuses LHS-pipe addition forms cleanly", {
  d <- mk_brms_gaussian_data()
  err <- tryCatch(
    flexybayes(y | weights(x) ~ x + (1 | g1), data = d, verbose = FALSE),
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
  fit <- suppressMessages(flexybayes(
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
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_not_long_tier()
  d <- mk_brms_gaussian_data()
  rev <- suppressMessages(flexybayes(
    y ~ x + (1 | g1),
    data = d,
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    review_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(rev, "flexybayes_review")
  expect_true(is.character(rev$code) && nchar(rev$code) > 0L)
  expect_identical(rev$backend, "stan_via_brms")
  expect_identical(rev$ir$source, "brms")
  # Structural test (does proceed() fit and cache?), the 500-draw
  # budget above still leaves residual ESS-only warnings (no R-hat /
  # divergence) on this particular model, so muffle those specifically
  # instead of push the budget even higher for a shape assertion.
  fit <- .muffle_ess_warnings(suppressMessages(
    suppressWarnings(proceed(rev), classes = "flexybayes_deprecated")
  ))
  expect_s3_class(fit, "flexybayes")
  expect_true("(Intercept)" %in% names(coef(fit)))
  # Second proceed() returns the cached fit (and prints a note to that
  # effect, suppress it, same as the first proceed() call above).
  fit2 <- suppressMessages(
    suppressWarnings(proceed(rev), classes = "flexybayes_deprecated")
  )
  expect_identical(coef(fit), coef(fit2))
})


# ---------------------------------------------------------------- #
# (11) backend_decision(fit) uniform shape across paths            #
# ---------------------------------------------------------------- #

test_that("backend_decision() returns uniform-shape trace on fb_brms() fits", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_brms_gaussian_data()
  # The assertions are structural, so the low-ESS warning of the one-chain
  # 500-draw budget is muffled (see helper-ess-warnings.R).
  fit <- .muffle_ess_warnings(suppressMessages(flexybayes(
    y ~ x + (1 | g1),
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
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
  expect_identical(bd$backend, "brms")
  expect_identical(bd$path, "explicit_brms")
  expect_null(bd$gate_checks)
  # Four fields added in v0.3.6: NULL on the small-data fast
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
# (12) review_code under backend = "inla" raises clean refusal      #
# ---------------------------------------------------------------- #

test_that("fb(review_code = TRUE, backend = 'inla') raises structured refusal", {
  d <- mk_brms_gaussian_data()
  err <- tryCatch(
    flexybayes(
      y ~ x + (1 | g1),
      data = d,
      backend = "inla",
      review_code = TRUE,
      verbose = FALSE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_review_code_backend_unsupported")
  expect_true(grepl("review_code", conditionMessage(err), fixed = TRUE))
  expect_true(grepl("brms", conditionMessage(err), fixed = TRUE))
})


# ---------------------------------------------------------------- #
# (13) [removed] canonical-name alignment across brms and INLA,     #
#      covered by test-as-draws-generics.R ("both engines answer    #
#      with the same parameter tokens")                             #
# ---------------------------------------------------------------- #

# ---------------------------------------------------------------- #
# (14) Invalid backend value raises match.arg error                #
# ---------------------------------------------------------------- #

test_that("fb_brms() rejects invalid backend value", {
  # 0.9.3: an unrecognised backend name is intercepted by
  # .check_known_backend_name() before match.arg() runs, so the error is
  # a typed flexyBayes refusal naming the active engines.
  d <- mk_brms_gaussian_data()
  err <- tryCatch(
    flexybayes(y ~ x + (1 | g1), data = d, backend = "stan", verbose = FALSE),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_unknown_backend_refusal")
  expect_identical(err$reason_code, "unknown_backend")
})

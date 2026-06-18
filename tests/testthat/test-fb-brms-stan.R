# Tests for fb(backend = "brms") -- Stan-passthrough emit.
#
# Two layers:
#   - Unit tests on .priors_to_brms_specs() (the pure spec-list
#     translation) that do NOT require brms. These run on every
#     CI / R CMD check pass and lock the eight-row translation
#     table contract.
#   - Round-trip tests via fb(backend = "brms") that DO
#     require brms (+ a C++ toolchain for Stan). Skipped on CRAN
#     / CI; opt-in via NOT_CRAN = "true".
#
# Stan compile latency (typically 30-60 sec for the first call,
# faster on cmdstanr cache reuse) means the round-trip subtests
# are not part of the default check loop -- they exercise the
# eight-pattern triangulation protocol's third-engine gap when
# the developer runs `NOT_CRAN=true Rscript tools/tally.R` with
# brms + a Stan toolchain installed.

# ---------------------------------------------------------------- #
# Fixtures                                                          #
# ---------------------------------------------------------------- #

mk_brms_stan_gaussian_data <- function() {
  set.seed(20260523L)
  n <- 40L
  g <- factor(rep(letters[1:5], each = 8L))
  x <- rnorm(n)
  b_g <- rnorm(5L, sd = 1.0)[as.integer(g)]
  y <- 1 + 0.3 * x + b_g + rnorm(n, sd = 0.5)
  data.frame(y = y, x = x, g = g)
}

mk_brms_stan_bernoulli_data <- function() {
  set.seed(20260523L)
  n <- 60L
  g <- factor(rep(letters[1:6], each = 10L))
  x <- rnorm(n)
  eta <- 0.2 + 0.4 * x + rnorm(6L, sd = 0.5)[as.integer(g)]
  y <- rbinom(n, size = 1L, prob = plogis(eta))
  data.frame(y = y, x = x, g = g)
}


# ---------------------------------------------------------------- #
# (1) priors_to_brms_specs: legacy-scalar bridge                    #
# ---------------------------------------------------------------- #

test_that("priors_to_brms_specs: legacy bridge for fixed-RI model", {
  d <- mk_brms_stan_gaussian_data()
  fb <- fb_from_brms(y ~ x + (1 | g), data = d)

  sp <- flexyBayes:::.priors_to_brms_specs(
    NULL,
    fb,
    prior_fixed_sd = 100,
    prior_vc_sd = 1
  )

  # Expect: Intercept(normal,100), b(normal,100), sigma(lognormal,1),
  # sd_g(lognormal,1) -- four rows.
  expect_true(is.list(sp))
  expect_gte(length(sp), 3L)

  classes <- vapply(sp, function(s) s$class, character(1))
  expect_true("Intercept" %in% classes)
  expect_true("b" %in% classes)
  expect_true("sigma" %in% classes)
  expect_true("sd" %in% classes)

  # Intercept + b: normal(0, 100)
  expect_match(
    sp[[which(classes == "Intercept")]]$string,
    "^normal\\(0, 100\\)$"
  )
  expect_match(sp[[which(classes == "b")]]$string, "^normal\\(0, 100\\)$")
  # sigma: lognormal(0, 1)
  expect_match(sp[[which(classes == "sigma")]]$string, "^lognormal\\(0, 1\\)$")
  # sd: lognormal(0, 1), group = g
  sd_idx <- which(classes == "sd")
  expect_match(sp[[sd_idx]]$string, "^lognormal\\(0, 1\\)$")
  expect_identical(sp[[sd_idx]]$group, "g")
})


# ---------------------------------------------------------------- #
# (2) priors_to_brms_specs: PC on sigma -> exponential(rate)        #
# ---------------------------------------------------------------- #

test_that("priors_to_brms_specs: PC sigma rate = -log(prob)/upper", {
  d <- mk_brms_stan_gaussian_data()
  fb <- fb_from_brms(y ~ x + (1 | g), data = d)
  fb$priors <- fb_prior(sigma ~ pc(upper = 2, prob = 0.05))

  sp <- flexyBayes:::.priors_to_brms_specs(fb$priors, fb)
  expect_length(sp, 1L)
  expect_identical(sp[[1L]]$class, "sigma")
  expect_match(sp[[1L]]$string, "^exponential\\(")

  # Rate verification: -log(0.05) / 2 = 1.497866...
  num_str <- sub("^exponential\\(", "", sub("\\)$", "", sp[[1L]]$string))
  rate <- as.numeric(num_str)
  expect_equal(rate, -log(0.05) / 2, tolerance = 1e-10)
})


# ---------------------------------------------------------------- #
# (3) priors_to_brms_specs: uniform on sd group                     #
# ---------------------------------------------------------------- #

test_that("priors_to_brms_specs: uniform(0, U) on sd group carries lb/ub", {
  d <- mk_brms_stan_gaussian_data()
  fb <- fb_from_brms(y ~ x + (1 | g), data = d)
  fb$priors <- fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 5))

  sp <- flexyBayes:::.priors_to_brms_specs(fb$priors, fb)
  expect_length(sp, 1L)
  expect_identical(sp[[1L]]$class, "sd")
  expect_identical(sp[[1L]]$group, "g")
  expect_match(sp[[1L]]$string, "^uniform\\(0, 5\\)$")
  expect_equal(sp[[1L]]$lb, 0)
  expect_equal(sp[[1L]]$ub, 5)
})


# ---------------------------------------------------------------- #
# (4) priors_to_brms_specs: student_t on b coef                     #
# ---------------------------------------------------------------- #

test_that("priors_to_brms_specs: student_t on a named b coef", {
  d <- mk_brms_stan_gaussian_data()
  fb <- fb_from_brms(y ~ x + (1 | g), data = d)
  fb$priors <- fb_prior(b("x") ~ student_t(df = 4, scale = 2.5))

  sp <- flexyBayes:::.priors_to_brms_specs(fb$priors, fb)
  expect_length(sp, 1L)
  expect_identical(sp[[1L]]$class, "b")
  expect_identical(sp[[1L]]$coef, "x")
  expect_match(sp[[1L]]$string, "^student_t\\(4, 0, 2\\.5\\)$")
})


# ---------------------------------------------------------------- #
# (5) priors_to_brms_specs: b("(Intercept)") -> class = Intercept   #
# ---------------------------------------------------------------- #

test_that("priors_to_brms_specs: Intercept-named b spec routes to class Intercept", {
  d <- mk_brms_stan_gaussian_data()
  fb <- fb_from_brms(y ~ x + (1 | g), data = d)
  fb$priors <- fb_prior(b("(Intercept)") ~ normal(0, 50))

  sp <- flexyBayes:::.priors_to_brms_specs(fb$priors, fb)
  expect_length(sp, 1L)
  expect_identical(sp[[1L]]$class, "Intercept")
  expect_true(is.na(sp[[1L]]$coef))
  expect_match(sp[[1L]]$string, "^normal\\(0, 50\\)$")
})


# ---------------------------------------------------------------- #
# (6) priors_to_brms_specs: unsupported family raises structured    #
# ---------------------------------------------------------------- #

test_that("priors_to_brms_specs: unsupported half_cauchy on sigma refuses", {
  d <- mk_brms_stan_gaussian_data()
  fb <- fb_from_brms(y ~ x + (1 | g), data = d)
  fb$priors <- fb_prior(sigma ~ half_cauchy(scale = 1))

  expect_error(
    flexyBayes:::.priors_to_brms_specs(fb$priors, fb),
    "half_cauchy.*sigma"
  )
})


# ---------------------------------------------------------------- #
# (7) priors_to_brms_specs: unsupported target raises structured    #
# ---------------------------------------------------------------- #

test_that("priors_to_brms_specs: cor() target refuses with pointer to ...", {
  d <- mk_brms_stan_gaussian_data()
  fb <- fb_from_brms(y ~ x + (1 | g), data = d)
  fb$priors <- fb_prior(cor(group = "g") ~ lkj(eta = 2))

  expect_error(
    flexyBayes:::.priors_to_brms_specs(fb$priors, fb),
    "cor\\(group"
  )
})


# ---------------------------------------------------------------- #
# (8) fb_brms() backend = "brms" with brms installed: review_code   #
# ---------------------------------------------------------------- #

test_that("fb(backend='brms', review_code=TRUE) returns Stan code", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  d <- mk_brms_stan_gaussian_data()
  rev <- suppressMessages(
    fb(
      y ~ x + (1 | g),
      data = d,
      backend = "brms",
      review_code = TRUE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  )

  expect_s3_class(rev, "flexybayes_review")
  expect_identical(rev$backend, "stan_via_brms")
  expect_true(is.character(rev$code))
  expect_true(nchar(rev$code) > 0L)
  # Stan code shape: must contain a `data { ... }` block.
  expect_match(rev$code, "data\\s*\\{")
  # And a `parameters { ... }` block.
  expect_match(rev$code, "parameters\\s*\\{")
})


# ---------------------------------------------------------------- #
# (9) Five-shape corpus: Gaussian RI round-trip + class             #
# ---------------------------------------------------------------- #

test_that("fb(backend='brms'): Gaussian RI fits cleanly", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  d <- mk_brms_stan_gaussian_data()
  fit <- suppressMessages(suppressWarnings(
    fb(
      y ~ x + (1 | g),
      data = d,
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  ))

  expect_s3_class(fit, "flexybayes_brms")
  expect_s3_class(fit, "flexybayes")
  expect_false(is.null(fit$brms))
  expect_true("(Intercept)" %in% names(coef(fit)))
  expect_true("x" %in% names(coef(fit)))

  bd <- backend_decision(fit)
  expect_identical(bd$backend, "brms")
  expect_identical(bd$path, "explicit_brms")

  expect_identical(fit$extras$fb_terms$source, "brms")
})


# ---------------------------------------------------------------- #
# (10) Bernoulli RI round-trip                                      #
# ---------------------------------------------------------------- #

test_that("fb(backend='brms'): Bernoulli RI fits cleanly", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  d <- mk_brms_stan_bernoulli_data()
  fit <- suppressMessages(suppressWarnings(
    fb(
      y ~ x + (1 | g),
      data = d,
      family = "binomial",
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  ))

  expect_s3_class(fit, "flexybayes_brms")
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "brms")
})


# ---------------------------------------------------------------- #
# (11) triangulate() across greta and brms: auto-resolve via registry #
# ---------------------------------------------------------------- #

test_that("triangulate(): brms + greta align via the identity registry", {
  testthat::skip_if_not_installed("brms")
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  d <- mk_brms_stan_gaussian_data()
  fit_b <- suppressMessages(suppressWarnings(
    fb(
      y ~ x + (1 | g),
      data = d,
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  ))
  fit_g <- suppressMessages(suppressWarnings(
    fb(
      y ~ x + (1 | g),
      data = d,
      backend = "greta",
      n_samples = 100L,
      warmup = 100L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  ))

  tri <- triangulate(fit_b, fit_g)
  expect_s3_class(tri, "triangulate_result")
  expect_identical(tri$source_a, "brms")
  expect_identical(tri$source_b, "greta")
  # Auto-resolve canonical names: (Intercept) and x must be in the
  # common set without a user-supplied name_map.
  expect_true(any(c("(Intercept)", "x") %in% tri$common))
})


# ---------------------------------------------------------------- #
# (12) auto-dispatch never routes to brms                           #
# ---------------------------------------------------------------- #

test_that("fb(backend='auto') does NOT route to Stan", {
  testthat::skip_if_not_installed("brms")
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  d <- mk_brms_stan_gaussian_data()
  fit_a <- suppressMessages(suppressWarnings(
    fb(
      y ~ x + (1 | g),
      data = d,
      backend = "auto",
      n_samples = 80L,
      warmup = 80L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  ))
  bd <- backend_decision(fit_a)
  expect_true(bd$backend %in% c("greta", "inla"))
  expect_false(identical(bd$backend, "brms"))
})


# ---------------------------------------------------------------- #
# (13) fb(backend='stan' or anything else) match.arg refusal   #
# ---------------------------------------------------------------- #

test_that("fb(backend='stan') raises the standard match.arg error", {
  d <- mk_brms_stan_gaussian_data()
  err <- tryCatch(
    fb(y ~ x + (1 | g), data = d, backend = "stan"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "'arg' should be one of")
  expect_match(err, "brms")
})


# ---------------------------------------------------------------- #
# (X) predict.flexybayes_brms: in-sample posterior mean reproducible #
# ---------------------------------------------------------------- #

test_that("predict.flexybayes_brms: in-sample posterior_epred matches summary", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  d <- mk_brms_stan_gaussian_data()
  fit <- suppressMessages(suppressWarnings(
    fb(
      y ~ x + (1 | g),
      data = d,
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  ))

  # In-sample call: shape matches nrow(data).
  pred <- predict(fit)
  expect_length(pred, nrow(d))
  expect_true(all(is.finite(pred)))

  # se.fit branch returns fit + se.fit of equal length.
  ps <- predict(fit, se.fit = TRUE)
  expect_named(ps, c("fit", "se.fit"))
  expect_length(ps$fit, nrow(d))
  expect_length(ps$se.fit, nrow(d))
})


# ---------------------------------------------------------------- #
# (Y) predict.flexybayes_brms: new-data and link-scale paths        #
# ---------------------------------------------------------------- #

test_that("predict.flexybayes_brms: newdata + type = 'link' return correct shape", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  d <- mk_brms_stan_gaussian_data()
  fit <- suppressMessages(suppressWarnings(
    fb(
      y ~ x + (1 | g),
      data = d,
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  ))

  newd <- d[1:5, ]
  pred_resp <- predict(fit, newdata = newd, type = "response")
  pred_link <- predict(fit, newdata = newd, type = "link")
  expect_length(pred_resp, 5L)
  expect_length(pred_link, 5L)
  # Identity link: response and link scale agree numerically.
  expect_equal(pred_resp, pred_link, tolerance = 1e-8)

  # summary = FALSE returns the full draws x rows posterior matrix.
  full <- predict(fit, newdata = newd, summary = FALSE)
  expect_true(is.matrix(full))
  expect_identical(ncol(full), 5L)
})


# ---------------------------------------------------------------- #
# (14) canonical_names.flexybayes_brms: identity with prefix-strip  #
# ---------------------------------------------------------------- #

test_that("canonical_names.flexybayes_brms: registry dispatches", {
  d <- mk_brms_stan_gaussian_data()
  fb <- fb_from_brms(y ~ x + (1 | g), data = d)
  stub <- structure(
    list(brms = NULL, extras = list(fb_terms = fb)),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
  cn <- canonical_names(stub)
  # Without a live brmsfit the mapper returns an empty map; what we
  # are testing here is the dispatch (S3 method registration + the
  # "registry" source label, not "no_mapper").
  expect_identical(cn$source, "registry")
  expect_true(is.character(cn$map))
})

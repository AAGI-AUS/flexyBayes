# Tests for prior_summary() -- user-facing accessor for the
# resolved prior on a flexyBayes fit.
#
# Contract:
#   - S3 generic + methods for flexybayes / flexybayes_inla /
#     flexybayes_direct_greta, plus a default that refuses with a
#     structured message.
#   - Returns a prior_summary_flexybayes object carrying `kind`
#     (one of "fb_prior" / "legacy_scalar" / "no_prior_recorded" /
#     "unknown_shape"), `backend`, and (when applicable) the
#     fb_prior + auto-default origin metadata.
#   - The print method labels auto-default vs user-supplied priors
#     and flags fb_greta() fits as declaration-only.



# ---------------------------------------------------------------- #
# (a) Default-fired uniform prior -- auto origin                    #
# ---------------------------------------------------------------- #

test_that("prior_summary() on a default-prior fit reports auto-uniform", {
  skip_if_no_greta()
  set.seed(20260524L)
  d <- data.frame(
    y = rnorm(30, 50, 5),
    x = rnorm(30),
    g = factor(rep(1:5, 6))
  )
  fit <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = d,
    backend = "greta",
    n_samples = 50L,
    warmup = 50L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  ps <- prior_summary(fit)
  expect_s3_class(ps, "prior_summary_flexybayes")
  expect_identical(ps$kind, "fb_prior")
  expect_identical(ps$backend, "greta")
  expect_identical(ps$default_origin, "auto")
  expect_false(isTRUE(ps$declaration_only))
  # The auto-default carries the family-aware scale attributes.
  expect_true(!is.null(ps$default_basis))
  expect_true(inherits(ps$fb_prior, "fb_prior"))
  # Print returns invisibly.
  out <- capture.output(print(ps))
  expect_true(
    any(grepl("auto-default", out, ignore.case = TRUE)) ||
      any(grepl("uniform", out, ignore.case = TRUE))
  )
})


# ---------------------------------------------------------------- #
# (b) User-supplied fb_prior() -- user origin                       #
# ---------------------------------------------------------------- #

test_that("prior_summary() on a user-supplied fb_prior fit reports user origin", {
  skip_if_no_greta()
  set.seed(20260524L)
  d <- data.frame(
    y = rnorm(30, 50, 5),
    x = rnorm(30),
    g = factor(rep(1:5, 6))
  )
  priors <- fb_prior(
    sigma ~ half_normal(scale = 5),
    sd(group = "g") ~ half_normal(scale = 5),
    b("x") ~ normal(mean = 0, sd = 10)
  )
  fit <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = d,
    prior = priors,
    backend = "greta",
    n_samples = 50L,
    warmup = 50L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  ps <- prior_summary(fit)
  expect_identical(ps$kind, "fb_prior")
  expect_identical(ps$default_origin, "user")
  expect_null(ps$default_basis)
  out <- capture.output(print(ps))
  expect_true(any(grepl("user-supplied", out, ignore.case = TRUE)))
})


# ---------------------------------------------------------------- #
# (c) Legacy scalar bridge -- prior_vc_sd passed explicitly         #
# ---------------------------------------------------------------- #

test_that("prior_summary() on a legacy-scalar fit reports the bridge values", {
  skip_if_no_greta()
  set.seed(20260524L)
  d <- data.frame(
    y = rnorm(30, 50, 5),
    x = rnorm(30),
    g = factor(rep(1:5, 6))
  )
  fit <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = d,
    prior_vc_sd = 1, # opt back into legacy default
    prior_fixed_sd = 100,
    backend = "greta",
    n_samples = 50L,
    warmup = 50L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  ps <- prior_summary(fit)
  expect_identical(ps$kind, "legacy_scalar")
  expect_equal(ps$fixed_sd, 100)
  expect_equal(ps$vc_sd, 1)
  out <- capture.output(print(ps))
  expect_true(any(grepl("legacy scalar", out, ignore.case = TRUE)))
})


# ---------------------------------------------------------------- #
# (d) Default-method refusal on a non-flexybayes object             #
# ---------------------------------------------------------------- #

test_that("prior_summary() default method refuses unknown classes", {
  expect_error(
    prior_summary(lm(mpg ~ wt, data = mtcars)),
    "prior_summary"
  )
})


# ---------------------------------------------------------------- #
# (e) fb_greta() declaration-only flag                              #
# ---------------------------------------------------------------- #

test_that("prior_summary() flags fb_greta() fits as declaration-only", {
  skip_if_no_greta()
  # Build a minimal direct-greta fit.
  d <- data.frame(y = rnorm(30, 50, 5), x = rnorm(30))
  m <- local({
    y <- greta::as_data(d$y)
    x <- greta::as_data(d$x)
    b0 <- greta::normal(0, 100)
    b1 <- greta::normal(0, 100)
    sigma <- greta::uniform(0, 5 * sd(d$y))
    mu <- b0 + b1 * x
    greta::distribution(y) <- greta::normal(mu, sigma)
    greta::model(b0, b1, sigma)
  })
  # v0.5.0: canonical_names attaches at IR build via fb_from_greta();
  # the greta pin fits the IR.
  ir_dg <- suppressMessages(fb_from_greta(
    m,
    data = d,
    canonical_names = c(b0 = "(Intercept)", b1 = "x", sigma = "sigma")
  ))
  fit_dg <- suppressMessages(fb_greta(
    ir_dg,
    n_samples = 50L,
    warmup = 50L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  ps <- prior_summary(fit_dg)
  expect_s3_class(ps, "prior_summary_flexybayes")
  expect_identical(ps$backend, "greta-direct")
  expect_true(isTRUE(ps$declaration_only))
  out <- capture.output(print(ps))
  expect_true(any(grepl("declaration", out, ignore.case = TRUE)))
})

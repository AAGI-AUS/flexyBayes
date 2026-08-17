# The combined multi-environment-trial model: interaction random effects
# and a sectioned residual in one fit.
#
# Until this file existed the two halves were only ever tested apart. The
# README and inst/KNOWN_ISSUES.md said so plainly -- the capability row
# read `emits`, meaning the emitted Stan program carried both blocks and
# nobody had run it -- and the August 2026 external audit recorded the
# combination as unverified for exactly that reason.
#
# Two assertions close it, in the order the fix specification requires:
#
#   (a) EMIT. The reconstructed brms formula carries the interaction
#       group-level term AND the distributional predictor on sigma, and
#       the Stan program brms writes from it carries both the group-level
#       standard-deviation vector and the sigma design matrix. Half a
#       model that samples beautifully is still half a model, so this is
#       asserted before anything is fitted.
#
#   (b) FIT. One live fit on simulated multi-environment data samples
#       with acceptable diagnostics, and the per-environment residual
#       variances come back in the order they were simulated in.
#
# The recovery check is against the REALISED draw, not the nominal
# generating parameters. On a single dataset of this size an environment
# can easily draw a residual sample two standard errors from its own
# variance, and a test that asserted the nominal value would be asserting
# the luck of the seed.

.met_combined_data <- function(
  n_gen = 10L,
  n_env = 4L,
  n_rep = 3L,
  seed = 20260815L
) {
  set.seed(seed)
  env_mean <- c(5.0, 6.0, 4.5, 7.0)[seq_len(n_env)]
  resid_sd <- c(0.4, 0.8, 1.2, 1.6)[seq_len(n_env)]
  d <- expand.grid(
    rep = seq_len(n_rep),
    gen = factor(sprintf("G%02d", seq_len(n_gen))),
    env = factor(sprintf("E%d", seq_len(n_env)))
  )
  g <- stats::rnorm(n_gen, 0, 1.0)
  ge <- stats::rnorm(n_gen * n_env, 0, 0.7)
  cell <- (as.integer(d$env) - 1L) * n_gen + as.integer(d$gen)
  e <- stats::rnorm(nrow(d), 0, resid_sd[as.integer(d$env)])
  d$yield <- env_mean[as.integer(d$env)] + g[as.integer(d$gen)] +
    ge[cell] + e
  # What the fit can actually be judged against: the residual variance
  # this draw realised, environment by environment.
  attr(d, "realised_resid_var") <- vapply(
    split(e, d$env),
    stats::var,
    numeric(1L)
  )
  d
}

test_that("the combined model emits both blocks, in the formula and in Stan", {
  skip_if_not_installed("brms")
  d <- .met_combined_data()

  fb <- fb_from_asreml(
    fixed = yield ~ env,
    random = ~ gen + gen:env,
    residual = ~ dsum(~ units | env),
    data = d
  )
  bf <- flexyBayes:::.fb_to_brms_formula(fb)
  main <- paste(deparse(bf$formula), collapse = " ")
  sigma <- paste(
    vapply(bf$pforms, function(p) paste(deparse(p), collapse = " "),
           character(1L)),
    collapse = " "
  )
  expect_match(main, "(1 | gen)", fixed = TRUE)
  expect_match(main, "(1 | gen:env)", fixed = TRUE)
  expect_match(sigma, "sigma ~ 0 + env", fixed = TRUE)

  # The Stan program is the thing that runs, so it gets its own
  # assertion: `sd_2` is the interaction term's standard deviation and
  # `X_sigma` / `b_sigma` are the residual predictor's design matrix and
  # coefficients. A program with one and not the other is a different
  # model from the one the user wrote.
  code <- paste(utils::capture.output(print(flexybayes(
    yield ~ env, random = ~ gen + gen:env,
    residual = ~ dsum(~ units | env), data = d,
    backend = "brms", return_code = TRUE, verbose = FALSE
  ))), collapse = "\n")
  expect_true(grepl("sd_1", code, fixed = TRUE))
  expect_true(grepl("sd_2", code, fixed = TRUE))
  expect_true(grepl("X_sigma", code, fixed = TRUE))
  expect_true(grepl("b_sigma", code, fixed = TRUE))
  # `sigma = exp(...)` rather than a scalar residual scale: the
  # sectioning replaces the residual, it does not sit beside it.
  expect_false(grepl("real<lower=0> sigma;", code, fixed = TRUE))
})

test_that("the combined model samples with acceptable diagnostics", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  skip_if_not_installed("posterior")
  d <- .met_combined_data()

  fit <- suppressWarnings(suppressMessages(flexybayes(
    fixed = yield ~ env,
    random = ~ gen + gen:env,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "brms",
    chains = 2L,
    warmup = 1000L,
    n_samples = 1000L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  expect_s3_class(fit, "flexybayes_brms")

  conv <- fit$extras$convergence
  rhat <- conv$gelman$psrf[, "Point est."]
  rhat <- rhat[is.finite(rhat)]
  ess <- conv$n_eff[is.finite(conv$n_eff)]
  expect_true(length(rhat) > 0L)
  expect_lt(max(rhat), 1.05)
  expect_gt(min(ess), 100)
  expect_identical(as.integer(conv$n_divergent %||% 0L), 0L)

  # Both variance components are present and separately estimated.
  vc <- fit$extras$variance_comps$component
  expect_true("sd_gen" %in% vc)
  expect_true("sd_gen:env" %in% vc)
  expect_false("sigma" %in% vc)

  # And the residual is reported per environment, in the order the
  # environments were simulated to be noisy in.
  tab <- flexyBayes:::.brms_residual_by_level_table(fit)
  expect_identical(tab$level, levels(d$env))
  expect_identical(
    order(tab$variance),
    order(attr(d, "realised_resid_var"))
  )
})

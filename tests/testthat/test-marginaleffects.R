# Integration tests for the marginaleffects support
# (get_coef / set_coef / get_vcov / get_predict / insight::get_data).
#
# predictions(), avg_slopes() and avg_comparisons() are exercised
# end-to-end on real greta- and INLA-backed fits, with numeric sanity
# against the known data-generating process. Closes the E-emmeans /
# marginaleffects coverage gap: before these, the fit class was not in
# marginaleffects' support list and no get_predict method existed.

suppressPackageStartupMessages(library(testthat))

.mfx_data <- function(seed = 5L, N = 150L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(sample(letters[1:3], N, replace = TRUE)),
    x = rnorm(N)
  )
  truef <- c(a = 0, b = 0.8, c = -0.5)
  d$y <- 1 + truef[as.character(d$f)] + 0.4 * d$x + rnorm(N, 0, 0.6)
  d
}

.mfx_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}


# ---------------------------------------------------------------- #
# INLA backend (deterministic, cheap)                               #
# ---------------------------------------------------------------- #

test_that("marginaleffects accessors work on an INLA fit", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("marginaleffects")
  skip_on_cran() # INLA posterior sampling: heavy + core-limited under --as-cran
  skip_on_ci()
  .mfx_silence()
  d <- .mfx_data()
  fit <- suppressMessages(fb(
    y ~ f + x,
    data = d,
    backend = "inla",
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  cf <- marginaleffects::get_coef(fit)
  expect_setequal(names(cf), c("(Intercept)", "fb", "fc", "x"))
  v <- marginaleffects::get_vcov(fit)
  expect_true(is.matrix(v) && all(dim(v) == length(cf)))
  # set_coef round-trips through get_predict (delta-method contract).
  fit2 <- marginaleffects::set_coef(
    fit,
    stats::setNames(rep(0, length(cf)), names(cf))
  )
  gp <- marginaleffects::get_predict(fit2, newdata = d)
  expect_true(all(gp$estimate == 0))
})

test_that("predictions() / avg_slopes() on an INLA fit recover the DGP", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("marginaleffects")
  skip_on_cran() # INLA posterior sampling: heavy + core-limited under --as-cran
  skip_on_ci()
  .mfx_silence()
  d <- .mfx_data()
  fit <- suppressMessages(fb(
    y ~ f + x,
    data = d,
    backend = "inla",
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  pr <- as.data.frame(marginaleffects::predictions(fit))
  expect_identical(nrow(pr), nrow(d))
  expect_true(all(is.finite(pr$estimate)))

  sl <- as.data.frame(marginaleffects::avg_slopes(fit))
  xs <- sl[sl$term == "x", "estimate"]
  expect_length(xs, 1L)
  expect_true(abs(xs - 0.4) < 0.2) # true x slope 0.4

  cmp <- as.data.frame(marginaleffects::avg_comparisons(fit, variables = "f"))
  est <- stats::setNames(cmp$estimate, cmp$contrast)
  expect_true(abs(est[["b - a"]] - 0.8) < 0.4)
  expect_true(abs(est[["c - a"]] - (-0.5)) < 0.4)
})


# ---------------------------------------------------------------- #
# greta backend (MCMC -- gated)                                     #
# ---------------------------------------------------------------- #

test_that("predictions() / avg_slopes() work on a greta fit", {
  skip_if_greta_backend_unusable()
  skip_if_not_installed("marginaleffects")
  skip_on_cran()
  skip_on_ci()
  .mfx_silence()
  d <- .mfx_data()
  fit <- suppressMessages(flexybayes(
    fixed = y ~ f + x,
    data = d,
    backend = "greta",
    n_samples = 400L,
    warmup = 400L,
    chains = 2L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  pr <- as.data.frame(marginaleffects::predictions(fit))
  expect_identical(nrow(pr), nrow(d))
  expect_true(all(is.finite(pr$estimate)))

  sl <- as.data.frame(marginaleffects::avg_slopes(fit))
  xs <- sl[sl$term == "x", "estimate"]
  expect_length(xs, 1L)
  expect_true(abs(xs - 0.4) < 0.25)
})


# ---------------------------------------------------------------- #
# Cross-engine consistency of the integration                       #
# ---------------------------------------------------------------- #

test_that("INLA and greta avg_comparisons agree on the factor effect", {
  skip_if_not_installed("INLA")
  skip_if_greta_backend_unusable()
  skip_if_not_installed("marginaleffects")
  skip_on_cran()
  skip_on_ci()
  .mfx_silence()
  d <- .mfx_data()
  fi <- suppressMessages(fb(
    y ~ f + x,
    data = d,
    backend = "inla",
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fg <- suppressMessages(flexybayes(
    fixed = y ~ f + x,
    data = d,
    backend = "greta",
    n_samples = 500L,
    warmup = 500L,
    chains = 2L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  ci <- as.data.frame(marginaleffects::avg_comparisons(fi, variables = "f"))
  cg <- as.data.frame(marginaleffects::avg_comparisons(fg, variables = "f"))
  ei <- stats::setNames(ci$estimate, ci$contrast)
  eg <- stats::setNames(cg$estimate, cg$contrast)
  common <- intersect(names(ei), names(eg))
  expect_true(length(common) >= 1L)
  # Two independent engines on the same data + matched fixed-effect
  # priors should agree to within MC tolerance.
  expect_true(max(abs(ei[common] - eg[common])) < 0.15)
})

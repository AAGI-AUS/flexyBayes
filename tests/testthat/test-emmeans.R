# Integration tests for the emmeans support (recover_data / emm_basis).
#
# emmeans(fit, ~ f) and contrast() are exercised end-to-end on real
# greta- and INLA-backed fits, with numeric sanity against the known
# data-generating process. Closes the E-emmeans coverage gap: before
# these, the exported methods had no standalone test and were in fact
# non-functional (no nbasis on the basis, no INLA model interface).

suppressPackageStartupMessages(library(testthat))

# Fixed-seed design with a 3-level factor and a continuous covariate.
.emm_data <- function(seed = 5L, N = 150L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(sample(letters[1:3], N, replace = TRUE)),
    x = rnorm(N)
  )
  truef <- c(a = 0, b = 0.8, c = -0.5)
  d$y <- 1 + truef[as.character(d$f)] + 0.4 * d$x + rnorm(N, 0, 0.6)
  d
}

.emm_silence <- function() {
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

test_that("emmeans(fit, ~ f) works on an INLA fit and recovers cell means", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran() # INLA posterior sampling: heavy + core-limited under --as-cran
  skip_on_ci()
  .emm_silence()
  d <- .emm_data()
  fit <- suppressMessages(fb(
    y ~ f + x,
    data = d,
    backend = "inla",
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  emm <- emmeans::emmeans(fit, ~f)
  s <- as.data.frame(summary(emm))
  expect_identical(nrow(s), 3L)
  expect_setequal(as.character(s$f), c("a", "b", "c"))
  # Cell means at mean(x): a ~ 1.0, b ~ 1.8, c ~ 0.5 (DGP), loose band.
  m <- stats::setNames(s$emmean, s$f)
  expect_true(m["b"] > m["a"]) # b above a (true gap +0.8)
  expect_true(m["a"] > m["c"]) # a above c (true gap +0.5)
  expect_true(all(abs(m - c(a = 1, b = 1.8, c = 0.5)) < 0.4))
  expect_true(all(is.finite(s$SE)) && all(s$SE > 0))
})

test_that("emmeans pairwise contrasts on an INLA fit recover the DGP gaps", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran() # INLA posterior sampling: heavy + core-limited under --as-cran
  skip_on_ci()
  .emm_silence()
  d <- .emm_data()
  fit <- suppressMessages(fb(
    y ~ f + x,
    data = d,
    backend = "inla",
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  ct <- as.data.frame(summary(emmeans::contrast(
    emmeans::emmeans(fit, ~f),
    "pairwise"
  )))
  expect_identical(nrow(ct), 3L) # a-b, a-c, b-c
  est <- stats::setNames(ct$estimate, ct$contrast)
  # a - b ~ -0.8 (true), a - c ~ +0.5 (true); generous band.
  expect_true(abs(est[["a - b"]] - (-0.8)) < 0.4)
  expect_true(abs(est[["a - c"]] - (0.5)) < 0.4)
  expect_true(all(is.finite(ct$SE)))
})


# ---------------------------------------------------------------- #
# greta backend (MCMC -- gated)                                     #
# ---------------------------------------------------------------- #

test_that("emmeans(fit, ~ f) works on a greta fit (over-parameterised basis)", {
  skip_if_greta_backend_unusable()
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .emm_silence()
  d <- .emm_data()
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

  emm <- emmeans::emmeans(fit, ~f)
  s <- as.data.frame(summary(emm))
  expect_identical(nrow(s), 3L)
  m <- stats::setNames(s$emmean, s$f)
  # The over-parameterised (intercept + all levels) greta basis is
  # rank-deficient; emmeans must still return finite, estimable EMMs
  # via the non-estimability basis.
  expect_true(all(is.finite(m)))
  expect_true(m["b"] > m["a"] && m["a"] > m["c"])
  expect_true(all(abs(m - c(a = 1, b = 1.8, c = 0.5)) < 0.5))
})


# ---------------------------------------------------------------- #
# Honest scope boundary                                             #
# ---------------------------------------------------------------- #

test_that("emmeans refuses informatively when the design cannot be reconciled", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  .emm_silence()
  # A model the integration cannot reconcile is refused, not silently
  # mis-mapped. We construct the refusal directly via the internal
  # design-matrix builder with a deliberately mismatched coef basis,
  # so the test does not depend on which formulae happen to reconcile.
  d <- .emm_data()
  trms <- stats::delete.response(stats::terms(y ~ f + x))
  expect_error(
    flexyBayes:::.fb_fixef_model_matrix(trms, d, c("(Intercept)", "WRONG"), d),
    "could not reconcile|posterior draws"
  )
})

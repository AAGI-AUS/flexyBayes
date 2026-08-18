# =============================================================================
# The posterior generics on a flexyBayes fit.
#
# A Bayesian reaches for `posterior::as_draws_df(fit)`, not for a
# package-specific extractor. These tests hold three things:
#
#   * the same canonical parameter tokens come back from both engines --
#     (Intercept), the fixed-effect terms, sigma, sd_<group> -- so a
#     workflow written against one engine reads the other unchanged;
#   * the variance components are on the standard-deviation scale, and the
#     values are the same ones `.fb_canonical_draws()` returns, not a
#     second extraction that could drift from it;
#   * the Suggests guard refuses with a written message naming the package
#     and the install command, never a bare requireNamespace() failure.
#
# The INLA draws are SAMPLED from the fitted approximation, so any test
# comparing two extractions has to seed the sampling step. The fitting is
# deterministic; the sampling is not.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


# ---------------------------------------------------------------- #
# Fixtures                                                          #
# ---------------------------------------------------------------- #

.adg_data <- function(seed = 20260817L, n = 60L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(rep(c("a", "b"), length.out = n)),
    g = factor(rep(letters[1:5], length.out = n)),
    x = stats::rnorm(n)
  )
  b_g <- stats::rnorm(5L, sd = 0.8)[as.integer(d$g)]
  d$y <- 1 + 0.6 * (d$f == "b") + 0.4 * d$x + b_g + stats::rnorm(n, sd = 0.5)
  d
}

.adg_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

# Both engines are fitted at most once per file; testthat gives each file
# its own process, so a file-local cache is safe and keeps the Stan
# compile off the clock for every block after the first.
.adg_cache <- new.env(parent = emptyenv())

.adg_inla_fit <- function() {
  if (is.null(.adg_cache$inla)) {
    .adg_cache$inla <- suppressMessages(fb(
      y ~ f + x,
      random = ~g,
      data = .adg_data(),
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .adg_cache$inla
}

.adg_brms_fit <- function() {
  if (is.null(.adg_cache$brms)) {
    .adg_cache$brms <- suppressMessages(suppressWarnings(fb(
      y ~ f + x + (1 | g),
      data = .adg_data(),
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .adg_cache$brms
}

# The canonical tokens this model earns on either engine: an intercept,
# one fixed-effect term per non-reference level plus the covariate, the
# residual SD, and the group SD.
.ADG_TOKENS <- c("(Intercept)", "fb", "x", "sigma", "sd_g")


# ---------------------------------------------------------------- #
# 1. The draws object, per engine                                   #
# ---------------------------------------------------------------- #

test_that("as_draws_df() returns canonical names on an INLA fit", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("posterior")
  skip_on_cran()
  skip_on_ci()
  .adg_silence()
  fit <- .adg_inla_fit()

  set.seed(4242L)
  d <- posterior::as_draws_df(fit, n_draws = 50L)

  expect_s3_class(d, "draws_df")
  expect_setequal(posterior::variables(d), .ADG_TOKENS)
  expect_equal(posterior::ndraws(d), 50L)
  # A standard deviation is positive. A precision draw handed over
  # untransformed would still be positive, so the scale is checked
  # against the seam itself in section 2 -- this is the cheap guard.
  expect_true(all(as.numeric(d[["sd_g"]]) > 0))
  expect_true(all(as.numeric(d[["sigma"]]) > 0))
})

test_that("as_draws_df() returns canonical names on a brms fit", {
  skip_if_not_installed("brms")
  skip_if_not_installed("posterior")
  skip_on_cran()
  skip_on_ci()
  .adg_silence()
  fit <- .adg_brms_fit()

  d <- posterior::as_draws_df(fit)

  expect_s3_class(d, "draws_df")
  expect_setequal(posterior::variables(d), .ADG_TOKENS)
  expect_equal(posterior::ndraws(d), 200L)
  expect_true(all(as.numeric(d[["sd_g"]]) > 0))
})

test_that("both engines answer with the same parameter tokens", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  skip_if_not_installed("posterior")
  skip_on_cran()
  skip_on_ci()
  .adg_silence()

  set.seed(4242L)
  di <- posterior::as_draws_df(.adg_inla_fit(), n_draws = 50L)
  db <- posterior::as_draws_df(.adg_brms_fit())

  # Cross-engine name equality is the whole point of the canonical view:
  # a workflow written against one engine reads the other unchanged.
  expect_identical(
    sort(posterior::variables(di)),
    sort(posterior::variables(db))
  )
})


# ---------------------------------------------------------------- #
# 2. The values are the seam's own                                  #
# ---------------------------------------------------------------- #

test_that("the draws_df columns are the seam's values, unaltered", {
  # The exact identity, held where it can be held exactly: the seam is
  # mocked, so what comes back is what went in, column for column. On a
  # live INLA fit two extractions cannot be compared this way -- see the
  # next block for why.
  skip_if_not_installed("posterior")
  ref <- list(
    `(Intercept)` = c(1.5, 2.5, 3.5),
    x = c(-1, 0, 1),
    sigma = c(0.4, 0.5, 0.6),
    sd_g = c(0.7, 0.8, 0.9)
  )
  local_mocked_bindings(
    .fb_canonical_draws = function(fit, ...) ref,
    .package = "flexyBayes"
  )
  d <- flexyBayes:::.fb_door_draws_df(
    structure(list(), class = "flexybayes"),
    n_draws = 3L
  )

  expect_setequal(posterior::variables(d), names(ref))
  for (nm in names(ref)) {
    expect_identical(as.numeric(d[[nm]]), ref[[nm]], label = nm)
  }
})

test_that("brms draws match the seam exactly, INLA hyperparameters do", {
  # What is reproducible on each engine, and nothing more.
  #
  # A brms fit carries its draws; two extractions are byte-equal. An
  # INLA fit is SAMPLED from the fitted approximation, and only half of
  # that sampling follows R's random stream: the hyperparameter draws
  # do, so they repeat under the same set.seed(); the latent field --
  # the intercept and the fixed-effect terms -- is drawn by INLA's own
  # generator, seeded at random unless inla.posterior.sample()'s `seed`
  # argument is given (?INLA::inla.posterior.sample). So the fixed
  # effects are compared where they can be, on their posterior mean,
  # against Monte-Carlo error rather than for equality.
  skip_if_not_installed("posterior")
  skip_on_cran()
  skip_on_ci()
  .adg_silence()

  if (requireNamespace("brms", quietly = TRUE)) {
    fitb <- .adg_brms_fit()
    db <- posterior::as_draws_df(fitb)
    rawb <- flexyBayes:::.fb_canonical_draws(fitb)
    for (nm in .ADG_TOKENS) {
      expect_identical(as.numeric(db[[nm]]), as.numeric(rawb[[nm]]), label = nm)
    }
  }

  skip_if_not_installed("INLA")
  fit <- .adg_inla_fit()
  set.seed(20260817L)
  d <- posterior::as_draws_df(fit, n_draws = 300L)
  set.seed(20260817L)
  raw <- flexyBayes:::.fb_canonical_draws(fit, n_samples = 300L)

  for (nm in c("sigma", "sd_g")) {
    expect_identical(as.numeric(d[[nm]]), as.numeric(raw[[nm]]), label = nm)
  }
  for (nm in c("(Intercept)", "fb", "x")) {
    expect_equal(
      mean(as.numeric(d[[nm]])),
      mean(as.numeric(raw[[nm]])),
      tolerance = 0.05,
      label = nm
    )
  }

  # And the SD scale is a transform of the precision the engine stored,
  # not the precision itself: the two live on reciprocal-square-root
  # scales, so one cannot be mistaken for the other at these magnitudes.
  prec_mean <- fit$inla$summary.hyperpar[
    "Precision for the Gaussian observations", "mean"
  ]
  expect_lt(abs(median(as.numeric(d[["sigma"]])) - 1 / sqrt(prec_mean)), 0.1)
})

test_that("as_draws() and as_draws_matrix() derive from the same frame", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("posterior")
  skip_on_cran()
  skip_on_ci()
  .adg_silence()
  fit <- .adg_inla_fit()

  set.seed(31337L)
  a <- posterior::as_draws(fit, n_draws = 25L)
  set.seed(31337L)
  m <- posterior::as_draws_matrix(fit, n_draws = 25L)

  expect_s3_class(a, "draws_df")
  expect_s3_class(m, "draws_matrix")
  expect_setequal(posterior::variables(a), .ADG_TOKENS)
  expect_setequal(colnames(m), .ADG_TOKENS)
  expect_equal(
    as.numeric(m[, "sd_g"]),
    as.numeric(a[["sd_g"]]),
    tolerance = 1e-12
  )
})

test_that("n_draws is honoured on the INLA path and validated", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("posterior")
  skip_on_cran()
  skip_on_ci()
  .adg_silence()
  fit <- .adg_inla_fit()

  expect_equal(
    posterior::ndraws(posterior::as_draws_df(fit, n_draws = 15L)),
    15L
  )
  expect_error(
    posterior::as_draws_df(fit, n_draws = c(10L, 20L)),
    "single positive number"
  )
  expect_error(posterior::as_draws_df(fit, n_draws = 0L), "single positive")
})


# ---------------------------------------------------------------- #
# 3. The Suggests guard refuses by name, not by namespace           #
# ---------------------------------------------------------------- #

test_that("posterior absent refuses with a written message", {
  # Exercised without uninstalling posterior, by mocking the named
  # availability predicate the guard consults. The message has to name
  # the package and the install command; R's own "there is no package
  # called" is what the guard exists to prevent.
  local_mocked_bindings(
    .fb_posterior_available = function() FALSE,
    .package = "flexyBayes"
  )
  err <- tryCatch(
    flexyBayes:::.fb_door_draws_df(
      structure(list(), class = "flexybayes"),
      n_draws = 10L
    ),
    condition = function(e) e
  )
  msg <- conditionMessage(err)
  expect_match(msg, "posterior", fixed = TRUE)
  expect_match(msg, "install.packages(\"posterior\")", fixed = TRUE)
  expect_match(msg, "fb_as_draws_simple(fit)", fixed = TRUE)
  expect_false(grepl("there is no package called", msg, fixed = TRUE))
})

test_that("the guard fires before any draw is extracted", {
  # A guard that only fired after the extraction would let a fit with no
  # posterior mask the missing dependency.
  local_mocked_bindings(
    .fb_posterior_available = function() FALSE,
    .package = "flexyBayes"
  )
  expect_error(
    flexyBayes:::.fb_door_draws_df(
      structure(list(), class = "flexybayes"),
      n_draws = "not a count"
    ),
    "install.packages"
  )
})

test_that("a fit carrying no canonical parameter refuses by name", {
  # An empty draws list is what a fit with nothing to read returns. The
  # seam is mocked rather than an engine coerced into producing one:
  # the branch under test is the guard, not the extraction.
  skip_if_not_installed("posterior")
  local_mocked_bindings(
    .fb_canonical_draws = function(fit, ...) list(),
    .package = "flexyBayes"
  )
  err <- tryCatch(
    flexyBayes:::.fb_door_draws_df(
      structure(list(), class = "flexybayes"),
      n_draws = 10L
    ),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_fit_lacks_posterior_draws")
  expect_match(conditionMessage(err), "as_draws()", fixed = TRUE)
})

test_that("ragged draw vectors are refused rather than reshaped", {
  # posterior::as_draws_df() would recycle or error obscurely. A draws
  # object needs one value per parameter per draw, and the refusal names
  # the parameters whose lengths disagree.
  skip_if_not_installed("posterior")
  local_mocked_bindings(
    .fb_canonical_draws = function(fit, ...) {
      list(`(Intercept)` = c(1, 2, 3), sigma = c(1, 2))
    },
    .package = "flexyBayes"
  )
  expect_error(
    flexyBayes:::.fb_door_draws_df(
      structure(list(), class = "flexybayes"),
      n_draws = 10L
    ),
    "unequal length"
  )
})


# ---------------------------------------------------------------- #
# 4. The convergence surface is engine-native (H2 smoke)            #
# ---------------------------------------------------------------- #

test_that("summary(fit)$converge is present and engine-native on INLA", {
  # H2: summary(fit)$converge is the Bayesian diagnostic surface. Built
  # in Phase 1; asserted here so the Bayesian door's own tests carry a
  # tripwire on it. A Laplace approximation has no R-hat and none is
  # invented for it.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .adg_silence()
  cv <- summary(.adg_inla_fit())$converge

  expect_type(cv, "list")
  expect_true("engine" %in% names(cv))
  expect_match(cv$engine, "Laplace")
  expect_true(all(c("mode_status", "kld_max") %in% names(cv)))
  expect_false(any(grepl("rhat", names(cv), ignore.case = TRUE)))
})

test_that("summary(fit)$converge is present and engine-native on brms", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .adg_silence()
  cv <- summary(.adg_brms_fit())$converge

  expect_type(cv, "list")
  expect_match(cv$engine, "Stan")
  expect_true(all(
    c("max_rhat", "min_ess_bulk", "min_ess_tail", "n_divergent") %in% names(cv)
  ))
  expect_true(is.finite(cv$max_rhat))
})

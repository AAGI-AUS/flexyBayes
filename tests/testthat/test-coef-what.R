# =============================================================================
# coef(what = ) and ranef() -- the ASReml-shaped coefficient accessors.
#
# `coef(asreml_fit)$random` is how an ASReml user extracts BLUPs. The
# equivalent here is `coef(fit, what = "random")`, and `ranef(fit)` is the
# name the mixed-model ecosystem uses for the same thing.
#
# Two properties are load-bearing:
#
#   * the DEFAULT is unchanged. `coef(fit)` is still the fixed-effect named
#     numeric on every engine, because emmeans, marginaleffects, predict(),
#     vcov() dimnames and the tidiers all read it as one. The
#     characterisation file test-characterisation-accessors-090.R freezes
#     that contract from the other side; what is asserted here is that the
#     new argument did not disturb it.
#
#   * the random table has the SAME six columns on both engines, although
#     the engines record their random effects in unrelated shapes -- INLA a
#     data frame per term keyed by `ID`, brms a three-dimensional array per
#     grouping factor.
#
# summary(fit)$random and coef(fit, what = "random") are one construction,
# so the two accessors are asserted identical rather than separately
# correct.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


# ---------------------------------------------------------------- #
# Fixtures                                                          #
# ---------------------------------------------------------------- #

.cw_data <- function(seed = 20260817L, n = 60L) {
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

.cw_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

.cw_cache <- new.env(parent = emptyenv())

.cw_inla_fit <- function() {
  if (is.null(.cw_cache$inla)) {
    .cw_cache$inla <- suppressMessages(fb(
      y ~ f + x,
      random = ~g,
      data = .cw_data(),
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .cw_cache$inla
}

.cw_brms_fit <- function() {
  if (is.null(.cw_cache$brms)) {
    .cw_cache$brms <- suppressMessages(suppressWarnings(fb(
      y ~ f + x,
      random = ~g,
      data = .cw_data(),
      backend = "brms",
      aggregate = FALSE,
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .cw_cache$brms
}

.CW_RANDOM_COLS <- c(
  "group", "level", "estimate", "std.error", "conf.low", "conf.high"
)

.CW_MISSING_COLS <- c(
  "row", "estimate", "std.error", "conf.low", "conf.high"
)


# ---------------------------------------------------------------- #
# 1. The default is the old default                                 #
# ---------------------------------------------------------------- #

test_that("coef() with no `what` is the fixed-effect vector on INLA", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fit <- .cw_inla_fit()

  beta <- stats::coef(fit)
  expect_true(is.numeric(beta))
  expect_null(dim(beta))
  expect_identical(
    beta,
    stats::setNames(
      fit$inla$summary.fixed$mean,
      rownames(fit$inla$summary.fixed)
    )
  )
  # The explicit spelling of the default returns the same object.
  expect_identical(stats::coef(fit, what = "fixed"), beta)
  # And the emmeans / marginaleffects alignment survives.
  expect_identical(dimnames(stats::vcov(fit))[[1L]], names(beta))
})

test_that("coef() with no `what` is the fixed-effect vector on brms", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fit <- .cw_brms_fit()

  beta <- stats::coef(fit)
  expect_identical(beta, fit$glm$coefficients)
  expect_identical(stats::coef(fit, what = "fixed"), beta)
})

test_that("an unrecognised `what` is refused rather than defaulted", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fit <- .cw_inla_fit()
  expect_error(stats::coef(fit, what = "blups"), "should be one of")
})


# ---------------------------------------------------------------- #
# 2. what = "random" -- one frame per grouping factor               #
# ---------------------------------------------------------------- #

test_that("coef(what = 'random') has one row per level on INLA", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fit <- .cw_inla_fit()

  re <- stats::coef(fit, what = "random")
  expect_true(is.list(re))
  expect_identical(names(re), "g")
  expect_s3_class(re$g, "data.frame")
  expect_identical(names(re$g), .CW_RANDOM_COLS)
  expect_identical(nrow(re$g), 5L)
  expect_identical(sort(re$g$level), sort(levels(.cw_data()$g)))
  expect_true(all(re$g$group == "g"))
  expect_true(all(is.finite(re$g$estimate)))
  expect_true(all(re$g$std.error > 0))
  expect_true(all(re$g$conf.low <= re$g$conf.high))
  # The values are INLA's own, not a re-derivation.
  expect_equal(
    re$g$estimate,
    as.numeric(fit$inla$summary.random$g$mean),
    tolerance = 1e-12
  )
})

test_that("coef(what = 'random') has one row per level on brms", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fit <- .cw_brms_fit()

  re <- stats::coef(fit, what = "random")
  expect_true(is.list(re))
  expect_identical(names(re), "g")
  expect_identical(names(re$g), .CW_RANDOM_COLS)
  expect_identical(nrow(re$g), 5L)
  expect_true(all(re$g$std.error > 0))
  expect_true(all(re$g$conf.low <= re$g$conf.high))
  # brms::ranef() is the source, unmodified.
  arr <- brms::ranef(fit$brms)$g
  expect_equal(
    re$g$estimate,
    as.numeric(arr[, "Estimate", "Intercept"]),
    tolerance = 1e-12
  )
})

test_that("the summary slot and the accessor are one construction", {
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fits <- list()
  if (requireNamespace("INLA", quietly = TRUE)) {
    fits$inla <- .cw_inla_fit()
  }
  if (requireNamespace("brms", quietly = TRUE)) {
    fits$brms <- .cw_brms_fit()
  }
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    s <- suppressMessages(utils::capture.output(
      out <- summary(fits[[nm]])
    ))
    expect_identical(
      out$random,
      stats::coef(fits[[nm]], what = "random"),
      label = nm
    )
  }
})


# ---------------------------------------------------------------- #
# 3. what = "all"                                                   #
# ---------------------------------------------------------------- #

test_that("what = 'all' is the named list of the other three", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fit <- .cw_inla_fit()

  all <- stats::coef(fit, what = "all")
  expect_identical(names(all), c("fixed", "random", "missing"))
  expect_identical(all$fixed, stats::coef(fit))
  expect_identical(all$random, stats::coef(fit, what = "random"))
  expect_identical(all$missing, stats::coef(fit, what = "missing"))
})

test_that("what = 'missing' is a typed frame on complete data", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fit <- .cw_inla_fit()

  mv <- stats::coef(fit, what = "missing")
  expect_s3_class(mv, "data.frame")
  expect_identical(names(mv)[seq_along(.CW_MISSING_COLS)], .CW_MISSING_COLS)
  expect_identical(nrow(mv), 0L)
})


# ---------------------------------------------------------------- #
# 4. ranef()                                                        #
# ---------------------------------------------------------------- #

test_that("ranef() is coef(what = 'random')", {
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fits <- list()
  if (requireNamespace("INLA", quietly = TRUE)) {
    fits$inla <- .cw_inla_fit()
  }
  if (requireNamespace("brms", quietly = TRUE)) {
    fits$brms <- .cw_brms_fit()
  }
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    expect_identical(
      flexyBayes::ranef(fits[[nm]]),
      stats::coef(fits[[nm]], what = "random"),
      label = nm
    )
  }
})

test_that("the method is registered on nlme's generic too", {
  skip_if_not_installed("nlme")
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fit <- .cw_inla_fit()
  # nlme's is the generic lme4 and brms re-export; a session with any of
  # them attached masks the local one, and this registration is what keeps
  # a bare ranef(fit) working there.
  expect_identical(nlme::ranef(fit), stats::coef(fit, what = "random"))
})


# ---------------------------------------------------------------- #
# 5. The callers the default protects                               #
# ---------------------------------------------------------------- #

test_that("emmeans still runs through the unchanged coef() default", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fit <- .cw_inla_fit()

  emm <- suppressMessages(emmeans::emmeans(fit, ~f))
  s <- as.data.frame(summary(emm))
  expect_identical(nrow(s), 2L)
  expect_true(all(is.finite(s$SE)))
  expect_true(all(s$SE > 0))
})

test_that("predict() on newdata still reads the fixed vector", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .cw_silence()
  fit <- .cw_inla_fit()

  nd <- .cw_data()[1:4, , drop = FALSE]
  pred <- stats::predict(fit, newdata = nd)
  expect_true(is.numeric(pred))
  expect_identical(length(pred), 4L)
  expect_true(all(is.finite(pred)))
})

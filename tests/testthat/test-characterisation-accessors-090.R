# test-characterisation-accessors-090.R -- the 0.9.0 accessor surface,
# frozen before the ASReml-hands UX work touches it.
#
# These are characterisation tests, not contract tests: they record what
# coef(), summary(), predict() and the emmeans seam actually do today, so
# that a later change to any of them is a visible diff rather than a
# silent one. Where a test states a shape the coming work is expected to
# replace, the block says so in a comment naming what will supersede it.
#
# Two behaviours recorded here are not what the surface documentation
# suggests, and both matter to whoever edits next:
#
#   * coef() reaches its value by a different route per engine. On a brms
#     fit it is object$glm$coefficients. On an INLA fit there is no $glm
#     slot at all -- coef.flexybayes_inla() reads $inla$summary.fixed.
#     Anything that reasons about "the $glm coefficients" holds on one
#     engine only.
#
#   * summary() returned the engine's own summary object rather than a
#     shared one -- a four-slot list on INLA, brms's brmssummary on brms,
#     neither carrying a $varcomp. Section 3 has been updated to record
#     the 0.9.1 return, the shared summary.flexybayes object; the
#     engine's own summary stays reachable at $extras$summary.
#
# The asreml na.method() fixture the missingness work needs lives in
# helper-asreml-shapes.R and is exercised here without asreml installed.

suppressPackageStartupMessages(library(testthat))


# ---------------------------------------------------------------- #
# Fixtures                                                          #
# ---------------------------------------------------------------- #

# One dataset for every fit in this file: a two-level treatment factor
# (the smallest design emmeans can be asked about), a five-level
# grouping factor for the brms random intercept, and one covariate.
.char090_data <- function(seed = 20260817L, n = 60L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(rep(c("a", "b"), length.out = n)),
    g = factor(rep(letters[1:5], length.out = n)),
    x = stats::rnorm(n)
  )
  b_g <- stats::rnorm(5L, sd = 0.8)[as.integer(d$g)]
  d$y <- 1 + 0.6 * (d$f == "b") + 0.4 * d$x + b_g +
    stats::rnorm(n, sd = 0.5)
  d
}

.char090_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

# Both engines are fitted at most once per file. testthat runs each file
# in its own process, so a file-local cache is safe, and it keeps the
# brms Stan compile off the clock for every block after the first.
.char090_cache <- new.env(parent = emptyenv())

.char090_inla_fit <- function() {
  if (is.null(.char090_cache$inla)) {
    .char090_cache$inla <- suppressMessages(fb(
      y ~ f + x,
      data = .char090_data(),
      backend = "inla",
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .char090_cache$inla
}

.char090_brms_fit <- function() {
  if (is.null(.char090_cache$brms)) {
    .char090_cache$brms <- suppressMessages(suppressWarnings(fb(
      y ~ f + x + (1 | g),
      data = .char090_data(),
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .char090_cache$brms
}


# ---------------------------------------------------------------- #
# 1. The recorded asreml na.method() shape                          #
# ---------------------------------------------------------------- #

test_that("the recorded na.method() fixture is a bare list keyed x/y", {
  # asreml is deliberately absent from DESCRIPTION. The fixture in
  # helper-asreml-shapes.R is the oracle, recorded once from the live
  # licensed 4.2.0.392 install, and this block runs with asreml
  # uninstalled.
  m <- .asreml_na_method_recorded("explicit")

  expect_type(m, "list")
  # Detection is by shape, never by class: the value carries no class
  # attribute of its own, so an is() test would have nothing to match.
  expect_identical(class(m), "list")
  expect_identical(names(m), c("x", "y"))
  expect_identical(m$y, "include")
  expect_identical(m$x, "fail")
  expect_identical(setdiff(names(attributes(m)), "names"), character(0))
})

test_that("an unsupplied na.method() argument arrives as its whole vector", {
  # The trap the normaliser has to survive: na.method() does not reduce
  # an unsupplied argument to a scalar, so a partially specified call
  # carries a length-3 character vector whose FIRST element is the
  # effective policy.
  partial <- .asreml_na_method_recorded("y_include")
  expect_identical(partial$y, "include")
  expect_length(partial$x, 3L)
  expect_identical(partial$x[[1L]], "fail")

  full <- .asreml_na_method_recorded("default")
  expect_length(full$y, 3L)
  expect_length(full$x, 3L)
  expect_identical(full$y[[1L]], "include")
  expect_identical(full$x[[1L]], "fail")

  # The recorded defaults agree with the recorded no-argument call.
  expect_identical(.asreml_na_method_defaults()$y, full$y)
  expect_identical(.asreml_na_method_defaults()$x, full$x)
})

test_that("every recorded na.method() case keeps the same two-slot shape", {
  cases <- c(
    "explicit", "default", "y_include", "y_omit", "y_fail",
    "x_omit", "x_include"
  )
  for (case in cases) {
    m <- .asreml_na_method_recorded(case)
    expect_identical(names(m), c("x", "y"), info = case)
    expect_true(is.character(m$x) && is.character(m$y), info = case)
  }
  expect_identical(.asreml_na_method_version(), "4.2.0.392")
})


# ---------------------------------------------------------------- #
# 2. coef() -- the default is the fixed-effect named numeric        #
# ---------------------------------------------------------------- #

test_that("coef() on an INLA fit is the fixed-effect posterior means", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .char090_silence()
  fit <- .char090_inla_fit()
  beta <- stats::coef(fit)

  expect_true(is.numeric(beta))
  expect_null(dim(beta))
  expect_identical(names(beta), c("(Intercept)", "fb", "x"))

  # The route, not just the value: an INLA fit carries no $glm slot, so
  # coef.flexybayes_inla() reads INLA's own fixed-effect summary. Any
  # rewrite that assumes $glm$coefficients is the single source will
  # return NULL here.
  expect_null(fit$glm)
  expect_identical(
    beta,
    stats::setNames(
      fit$inla$summary.fixed$mean,
      rownames(fit$inla$summary.fixed)
    )
  )

  # vcov() dimnames follow coef() names -- the emmeans / marginaleffects
  # seam depends on this alignment.
  expect_identical(dimnames(stats::vcov(fit)), list(names(beta), names(beta)))
})

test_that("coef() on a brms fit is identical to fit$glm$coefficients", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .char090_silence()
  fit <- .char090_brms_fit()
  beta <- stats::coef(fit)

  expect_identical(beta, fit$glm$coefficients)
  expect_true(is.numeric(beta))
  expect_null(dim(beta))
  expect_identical(names(beta), c("(Intercept)", "fb", "x"))
})


# ---------------------------------------------------------------- #
# 3. summary() -- invisible return, one object on both engines      #
# ---------------------------------------------------------------- #
#
# Updated at 0.9.1. The 0.9.0 baseline this block recorded was two
# incomparable objects -- INLA's own four-slot list (fixed, random,
# hyperpar, fitted) and brms's `brmssummary`, each identical to
# `fit$extras$summary`, neither carrying a `$varcomp`. Both engines now
# return the same eleven-slot `summary.flexybayes` object, and the
# engine's own summary stays where it was, on `$extras$summary`, because
# seven internal readers reach for it there.
#
# The invisibility assertions are unchanged: they were true before and
# are true now. The full shape contract lives in
# test-summary-asreml-shape.R; what is recorded here is the diff from
# the frozen baseline.

.CHAR090_SUMMARY_SLOTS <- c(
  "fixed", "varcomp", "random", "missing", "converge",
  "n_design", "n_observed", "na_action", "model", "engine", "call"
)

test_that("summary() on an INLA fit returns the shared summary object", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .char090_silence()
  fit <- .char090_inla_fit()
  captured <- utils::capture.output(res <- withVisible(summary(fit)))

  expect_false(res$visible)
  expect_true(length(captured) > 0L)

  expect_identical(class(res$value), c("summary.flexybayes", "list"))
  expect_identical(names(res$value), .CHAR090_SUMMARY_SLOTS)
  expect_identical(res$value$engine, "inla")
  # The gap the baseline recorded is closed.
  expect_true(is.data.frame(res$value$varcomp))

  # And INLA's own four-slot summary is still where it was.
  expect_identical(
    names(fit$extras$summary),
    c("fixed", "random", "hyperpar", "fitted")
  )
})

test_that("summary() on a brms fit returns the shared summary object", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .char090_silence()
  fit <- .char090_brms_fit()
  captured <- utils::capture.output(res <- withVisible(summary(fit)))

  expect_false(res$visible)
  expect_true(length(captured) > 0L)

  expect_identical(class(res$value), c("summary.flexybayes", "list"))
  expect_identical(names(res$value), .CHAR090_SUMMARY_SLOTS)
  expect_identical(res$value$engine, "brms")
  expect_true(is.data.frame(res$value$varcomp))

  # brms's own summary is still where it was, and still brms's.
  expect_identical(class(fit$extras$summary), "brmssummary")
  expect_true(all(
    c("fixed", "random", "nobs") %in% names(fit$extras$summary)
  ))
})


# ---------------------------------------------------------------- #
# 4. predict() -- the no-classify path                              #
# ---------------------------------------------------------------- #

test_that("predict() with no newdata returns the fixed-effect predictor", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .char090_silence()
  d <- .char090_data()
  fit <- .char090_inla_fit()
  pred <- stats::predict(fit)

  expect_true(is.numeric(pred))
  expect_length(pred, nrow(d))
  # Unnamed, and no dim: several downstream consumers index it
  # positionally.
  expect_null(names(pred))
  expect_null(dim(pred))

  # On the identity link the two `type` values coincide, because
  # predict.flexybayes_inla() returns the linear predictor either way.
  expect_identical(
    stats::predict(fit, type = "link"),
    stats::predict(fit, type = "response")
  )

  # The value is X %*% bhat over the fit's own data -- fixed effects
  # only, which is why it is not fitted() on a fit carrying random
  # terms. WP-C adds `classify =` alongside this path and must leave it
  # untouched.
  trms <- stats::delete.response(stats::terms(stats::formula(fit)))
  x_mat <- stats::model.matrix(trms, fit$data)
  expect_equal(pred, as.numeric(x_mat %*% stats::coef(fit)))

  # se.fit = TRUE returns a two-element list, fit unnamed and se.fit
  # named by row.
  se <- stats::predict(fit, se.fit = TRUE)
  expect_identical(names(se), c("fit", "se.fit"))
  expect_length(se$fit, nrow(d))
  expect_length(se$se.fit, nrow(d))
  expect_true(all(se$se.fit > 0))
})


# ---------------------------------------------------------------- #
# 5. emmeans on a two-level factor                                  #
# ---------------------------------------------------------------- #

test_that("emmeans(fit, ~ f) runs on a two-level factor, INLA", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .char090_silence()
  fit <- .char090_inla_fit()

  emm <- emmeans::emmeans(fit, ~f)
  s <- as.data.frame(summary(emm))
  expect_identical(nrow(s), 2L)
  expect_setequal(as.character(s$f), c("a", "b"))
  expect_true(all(is.finite(s$emmean)))
  expect_true(all(is.finite(s$SE)) && all(s$SE > 0))
  # b is the raised level in the data-generating process.
  m <- stats::setNames(s$emmean, s$f)
  expect_true(m[["b"]] > m[["a"]])
})

test_that("emmeans(fit, ~ f) runs on a two-level factor, brms", {
  skip_if_not_installed("brms")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .char090_silence()
  fit <- .char090_brms_fit()

  emm <- emmeans::emmeans(fit, ~f)
  s <- as.data.frame(summary(emm))
  expect_identical(nrow(s), 2L)
  expect_setequal(as.character(s$f), c("a", "b"))
  expect_true(all(is.finite(s$emmean)))
  expect_true(all(is.finite(s$SE)) && all(s$SE > 0))
  m <- stats::setNames(s$emmean, s$f)
  expect_true(m[["b"]] > m[["a"]])
})

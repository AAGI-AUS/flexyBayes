# =============================================================================
# The engine and the representation a re-fit comes back on.
#
# update() re-issues the recorded call, which carries the requested engine
# (`backend`) and representation (`aggregate`). Without them an identity
# update() of a Stan fit would come back as an aggregated INLA fit.
#
# Five identity re-fits, one per route, plus the override. The record holds
# the request ("auto"), not the engine it resolved to, so a changed model
# is routed afresh. `verbose` is recorded but not carried.
#
# An INLA fit records no sampler settings. The last section checks that an
# auto re-fit whose INLA attempt fails still reaches brms with valid ones.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


.uec_data <- function(seed = 90210L, n = 60L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(rep(c("a", "b"), length.out = n)),
    g = factor(rep(letters[1:5], length.out = n)),
    x = stats::rnorm(n)
  )
  b_g <- stats::rnorm(5L, sd = 0.8)[as.integer(d$g)]
  d$y <- 1 + 0.6 * (d$f == "b") + 0.4 * d$x + b_g + stats::rnorm(n, sd = 0.5)
  # A count response for the aggregated fit, drawn after `y` so the Gaussian
  # fits see the same data. A Gaussian model is fitted per row.
  d$count <- stats::rpois(n, 3)
  d
}

.uec_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

.uec_cache <- new.env(parent = emptyenv())

# Named backend, per-row: the route whose identity re-fit used to leave
# Stan entirely.
.uec_brms_fit <- function() {
  if (is.null(.uec_cache$brms)) {
    .uec_cache$brms <- suppressMessages(suppressWarnings(flexybayes(
      y ~ f + x,
      random = ~g,
      data = .uec_data(),
      backend = "brms",
      aggregate = FALSE,
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .uec_cache$brms
}

# The continuous covariate keeps this fit per row.
.uec_inla_fit <- function() {
  if (is.null(.uec_cache$inla)) {
    .uec_cache$inla <- suppressMessages(flexybayes(
      y ~ f + x,
      random = ~g,
      data = .uec_data(),
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .uec_cache$inla
}

# The sharp per-row case: nothing about the model holds it per-row, only
# `aggregate = FALSE` does. Drop the recorded `aggregate` and the re-fit
# compresses, which is a different representation and a different object
# shape under the same name.
.uec_inla_forced_fit <- function() {
  if (is.null(.uec_cache$inla_forced)) {
    .uec_cache$inla_forced <- suppressMessages(flexybayes(
      y ~ f,
      random = ~g,
      data = .uec_data(),
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .uec_cache$inla_forced
}

# No continuous covariate and no named backend, which is what auto-routes
# to the aggregated Gaussian emit.
.uec_agg_fit <- function() {
  if (is.null(.uec_cache$agg)) {
    .uec_cache$agg <- suppressMessages(flexybayes(
      count ~ f,
      random = ~g,
      data = .uec_data(),
      family = "poisson",
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .uec_cache$agg
}


# ---------------------------------------------------------------- #
# 1. The record carries the request, not the resolution             #
# ---------------------------------------------------------------- #

test_that("the record carries the requested backend, not the resolved one", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .uec_silence()

  agg <- .uec_agg_fit()
  ci <- agg$extras$call_info
  # The fit ran on INLA, through the aggregated emit. What the call asked
  # for was neither of those things, and the record says what the call
  # asked for.
  expect_s3_class(agg, "flexybayes_aggregated")
  expect_identical(ci$backend, "auto")
  expect_identical(ci$aggregate, "auto")
  expect_false(ci$verbose)

  per_row <- .uec_inla_fit()
  ci2 <- per_row$extras$call_info
  expect_identical(ci2$backend, "inla")
  expect_false(ci2$aggregate)
  expect_false(ci2$verbose)
})


# ---------------------------------------------------------------- #
# 2. Five identity re-fits, one per route, plus the override         #
# ---------------------------------------------------------------- #

test_that("an identity update() of a brms fit comes back on brms", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  skip_if_not_long_tier()
  .uec_silence()
  fit <- .uec_brms_fit()
  expect_s3_class(fit, "flexybayes_brms")

  refit <- suppressMessages(suppressWarnings(stats::update(fit)))
  # The three assertions the defect broke: the engine, the object shape,
  # and the representation that follows from both.
  expect_s3_class(refit, "flexybayes_brms")
  expect_false(inherits(refit, "flexybayes_inla"))
  expect_false(inherits(refit, "flexybayes_aggregated"))
  expect_identical(class(refit), class(fit))
})

test_that("an identity update() of an auto aggregated fit stays aggregated", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .uec_silence()
  fit <- .uec_agg_fit()
  expect_s3_class(fit, "flexybayes_aggregated")

  refit <- suppressMessages(stats::update(fit))
  expect_s3_class(refit, "flexybayes_aggregated")
  expect_identical(class(refit), class(fit))
})

test_that("an identity update() of a per-row INLA fit stays per-row INLA", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .uec_silence()
  fit <- .uec_inla_fit()
  expect_s3_class(fit, "flexybayes_inla")
  expect_false(inherits(fit, "flexybayes_aggregated"))

  refit <- suppressMessages(stats::update(fit))
  expect_s3_class(refit, "flexybayes_inla")
  expect_false(inherits(refit, "flexybayes_aggregated"))
  expect_identical(class(refit), class(fit))
  # The unified variance-component table has to survive the re-fit on
  # whichever representation the fit comes back on.
  expect_true(is.data.frame(suppressMessages(summary(refit))$varcomp))
})

test_that("a per-row fit held there by aggregate = FALSE is not compressed", {
  # The same route as above, on a model the aggregation gate would
  # otherwise compress. Aggregated and per-row fits return the same
  # summary object, so only the class shows the difference.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .uec_silence()
  fit <- .uec_inla_forced_fit()
  expect_s3_class(fit, "flexybayes_inla")
  expect_false(inherits(fit, "flexybayes_aggregated"))

  refit <- suppressMessages(stats::update(fit))
  expect_false(inherits(refit, "flexybayes_aggregated"))
  expect_s3_class(refit, "flexybayes_inla")
  expect_identical(class(refit), class(fit))
  expect_true(is.data.frame(suppressMessages(summary(refit))$varcomp))
})

test_that("a backend override in ... still wins over the record", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  skip_if_not_long_tier()
  .uec_silence()
  fit <- .uec_inla_fit()

  # The sampler settings travel with the override because the INLA record
  # carries none.
  refit <- suppressMessages(suppressWarnings(stats::update(
    fit,
    backend = "brms",
    n_samples = 200L,
    warmup = 200L,
    chains = 1L
  )))
  expect_s3_class(refit, "flexybayes_brms")
})


# ---------------------------------------------------------------- #
# 3. `verbose` is recorded and deliberately not carried               #
# ---------------------------------------------------------------- #

test_that("verbose is recorded but does not travel with the re-fit", {
  # update() reproduces the model and its engine, representation,
  # missingness and prior, not the display settings of the first fit.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .uec_silence()
  fit <- .uec_inla_fit()

  # Recorded, so a reader of the fit can see what the call asked for.
  expect_false(fit$extras$call_info$verbose)

  # Not carried: the re-fit follows the current call's default, and what
  # is held invariant is the answer, not the console.
  refit <- suppressMessages(stats::update(fit))
  expect_identical(class(refit), class(fit))
  expect_true(refit$extras$call_info$verbose)
  expect_identical(
    prior_summary(refit)$kind,
    prior_summary(fit)$kind
  )
  expect_identical(
    suppressMessages(summary(refit))$varcomp$component,
    suppressMessages(summary(fit))$varcomp$component
  )

  # And a caller who wants quiet says so, on the spot.
  quiet <- utils::capture.output(
    suppressMessages(qfit <- stats::update(fit, verbose = FALSE))
  )
  expect_false(
    any(grepl("-- flexyBayes: INLA fit", quiet, fixed = TRUE))
  )
  expect_s3_class(qfit, "flexybayes_inla")
})


# ---------------------------------------------------------------- #
# 4. The fallback no longer forwards sampler settings that are not  #
#    there                                                          #
# ---------------------------------------------------------------- #

test_that("an auto fallback to brms omits absent sampler settings", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  skip_if_not_installed("testthat", "3.2.0")
  skip_on_cran()
  skip_on_ci()
  .uec_silence()

  # The INLA failure is mocked. The real trigger, a bounded uniform whose
  # upper bound cuts into the posterior, crashes the inla binary
  # (inst/KNOWN_ISSUES.md).
  testthat::local_mocked_bindings(
    inla = function(...) stop("forced INLA failure"),
    .package = "INLA"
  )
  captured <- NULL
  testthat::local_mocked_bindings(
    brm = function(...) {
      captured <<- list(...)
      # Signalling out of the mock keeps the test to the argument check.
      stop("captured")
    },
    .package = "brms"
  )

  # NULL sampler settings are what an INLA record holds, and what
  # update() therefore re-issues.
  suppressMessages(try(
    flexybayes(
      y ~ f + x,
      random = ~g,
      data = .uec_data(),
      backend = "auto",
      aggregate = FALSE,
      n_samples = NULL,
      warmup = NULL,
      chains = NULL,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ),
    silent = TRUE
  ))

  expect_false(is.null(captured))
  # brms builds `iter` by adding two settings together and refuses the
  # result when either is absent ("Cannot coerce 'iter' to a single
  # numeric value"). What it must receive is one usable number.
  expect_identical(length(captured$iter), 1L)
  expect_true(is.finite(captured$iter))
  expect_gt(captured$iter, captured$warmup)
  expect_identical(length(captured$chains), 1L)
  expect_true(is.finite(captured$chains))
})

test_that("an auto fallback hands the emit its arguments unevaluated", {
  # The fall-through emit is reached through do.call(), which evaluates
  # its arguments. `the_call` would re-run flexybayes() from inside its
  # own dispatch, and the formulas would take the dispatch frame as
  # their environment.
  #
  # The family here has no INLA likelihood, so the gate refuses and the
  # fit completes on the brms fallback.
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  skip_if_not_long_tier()
  .uec_silence()
  d <- .uec_data()
  d$y <- pmax(d$y, 0)

  fit <- suppressMessages(suppressWarnings(flexybayes(
    y ~ f + x,
    random = ~g,
    data = d,
    family = "hurdle_gamma",
    backend = "auto",
    aggregate = FALSE,
    n_samples = 100L,
    warmup = 100L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))

  expect_s3_class(fit, "flexybayes_brms")
  expect_true(nzchar(fit$extras$backend_decision$path %||% ""))
  expect_true(is.call(fit$extras$the_call))
  expect_s3_class(fit$extras$call_info$fixed, "formula")
  expect_s3_class(fit$extras$call_info$random, "formula")
})

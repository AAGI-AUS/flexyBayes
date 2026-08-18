# =============================================================================
# The engine and the representation a re-fit comes back on.
#
# update() re-issues the recorded call, and until 0.9.1 the record held
# neither the engine the call asked for nor the representation it ran
# under. Both therefore came back as the formal defaults -- backend =
# "auto", aggregate = "auto" -- so an identity update() of a Stan fit
# returned an aggregated INLA fit: a different inference engine and a
# different representation, under the same name and with nothing said.
#
# Five identity re-fits are held here, one per route, plus the override.
# What is recorded is the REQUEST ("auto"), not the engine the request
# resolved to: the recorded value is a policy, and a re-fit whose model
# has changed has to be free to route the changed model.
#
# `verbose` is recorded beside those two and deliberately not carried.
# The cut is what the argument can get wrong: a silently switched engine
# is a different answer under the same name, while a banner printed where
# the first fit printed none is a banner. update() reproduces the model
# and the policy behind it, not the display settings of the session that
# first ran it.
#
# The last section is the fallback these NULL-carrying INLA records used
# to break. An INLA fit records no sampler settings -- the nested Laplace
# approximation runs none -- so an auto re-fit whose INLA attempt fails
# reached the brms emit with all three NULL, and brms refused the `iter`
# it built from them. The user saw a brms argument error where an INLA
# failure had happened.
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
    .uec_cache$brms <- suppressMessages(suppressWarnings(fb(
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

# The continuous covariate keeps this one per-row rather than letting the
# aggregation gate compress it.
.uec_inla_fit <- function() {
  if (is.null(.uec_cache$inla)) {
    .uec_cache$inla <- suppressMessages(fb(
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
    .uec_cache$inla_forced <- suppressMessages(fb(
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
    .uec_cache$agg <- suppressMessages(fb(
      y ~ f,
      random = ~g,
      data = .uec_data(),
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
  # otherwise compress. This is where the missing `aggregate` bit: the
  # re-fit came back <flexybayes_aggregated>, a different representation
  # from the one the user pinned. (An aggregated fit now returns the same
  # unified summary object, so the class assertions below are what
  # catches the drift; before that they were the only thing that could.)
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
  .uec_silence()
  fit <- .uec_inla_fit()

  # The sampler settings travel with the override because the INLA record
  # carries none -- this is the caller doing deliberately what the
  # fallback below has to do on its own.
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
  # The cut this file is drawn along. update() reproduces the model and
  # the policy behind it -- engine, representation, missingness, prior --
  # and not the display settings of the session that first ran the fit. A
  # silently switched engine is a different answer under the same name. A
  # banner printed where the first fit printed none is a banner.
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

  # The INLA failure is forced rather than provoked. The failure that
  # found this -- a bounded uniform whose upper bound cuts into the
  # posterior -- segfaults the inla binary, which is recorded in
  # inst/KNOWN_ISSUES.md and has no business in a test suite.
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
    fb(
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
  # The fall-through emit is now reached through do.call(), which
  # re-evaluates everything it is handed unless told not to. Four of the
  # arguments are language objects: `the_call` is a call, so evaluating
  # it would re-run flexybayes() from inside its own dispatch, and
  # `fixed` / `random` / `residual` are formulas, which would come back
  # re-created against the dispatch frame rather than the caller's.
  #
  # This fit takes the gate-refusal fallback -- a multi-stratum
  # interaction random term routes to brms -- and completes, which is
  # what forces the arguments the mocked test above never reaches.
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .uec_silence()
  d <- .uec_data()
  d$h <- factor(rep(c("p", "q"), length.out = nrow(d)))

  fit <- suppressMessages(suppressWarnings(fb(
    y ~ f + x,
    random = ~ g:h,
    data = d,
    backend = "auto",
    aggregate = FALSE,
    n_samples = 100L,
    warmup = 100L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))

  expect_s3_class(fit, "flexybayes_brms")
  expect_identical(
    fit$extras$backend_decision$path,
    "auto_multistratum_to_brms"
  )
  expect_true(is.call(fit$extras$the_call))
  expect_s3_class(fit$extras$call_info$fixed, "formula")
  expect_s3_class(fit$extras$call_info$random, "formula")
})

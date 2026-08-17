# test-sampler-passthrough.R -- `seed` and `control` reach Stan.
#
# Until 0.9.0 flexybayes() had no `...` and no sampler arguments beyond
# n_samples / warmup / chains, so `set.seed()` before a call did not pin
# Stan's stream and `adapt_delta` was unreachable without dropping to
# brms::brm() directly. The gap was load-bearing: the package's own
# heterogeneous-variance oracle measured 0.27 percent run-to-run drift on
# the numbers it quotes.
#
# Two levels of assertion. The call-level test intercepts brms::brm() and
# reads the arguments it was handed, which runs without a compiler. The
# live test fits the same small model twice under one seed and compares
# the draws, which is the property a reader of a quoted posterior cares
# about.

mk_passthrough_data <- function(n_group = 30L, n_rep = 6L) {
  set.seed(20260815L)
  g <- factor(rep(seq_len(n_group), each = n_rep))
  x <- stats::rnorm(n_group * n_rep)
  u <- stats::rnorm(n_group, sd = 0.7)
  y <- 1 + 0.5 * x + u[as.integer(g)] + stats::rnorm(n_group * n_rep, sd = 0.4)
  data.frame(y = y, x = x, g = g)
}

# ---------------------------------------------------------------- #
# Call level -- the arguments arrive at brms::brm()                 #
# ---------------------------------------------------------------- #

test_that("seed and control reach the brms::brm() call", {
  skip_if_not_installed("brms")
  skip_if_not_installed("testthat", "3.2.0")

  captured <- NULL
  testthat::local_mocked_bindings(
    brm = function(...) {
      captured <<- list(...)
      # Signalling out of the mock keeps the test to the argument check:
      # nothing downstream of brm() is under test here.
      stop("captured")
    },
    .package = "brms"
  )

  d <- mk_passthrough_data()
  suppressMessages(try(
    flexybayes(
      fixed = y ~ x,
      random = ~g,
      data = d,
      backend = "brms",
      n_samples = 50L,
      warmup = 50L,
      chains = 1L,
      seed = 4321L,
      control = list(adapt_delta = 0.97, max_treedepth = 12L),
      verbose = FALSE,
      mcmc_verbose = FALSE
    ),
    silent = TRUE
  ))

  expect_false(is.null(captured))
  expect_identical(captured$seed, 4321L)
  expect_identical(captured$control$adapt_delta, 0.97)
  expect_identical(captured$control$max_treedepth, 12L)
})

test_that("an unset seed maps onto brms's own default rather than NULL", {
  skip_if_not_installed("brms")
  skip_if_not_installed("testthat", "3.2.0")

  captured <- NULL
  testthat::local_mocked_bindings(
    brm = function(...) {
      captured <<- list(...)
      stop("captured")
    },
    .package = "brms"
  )

  d <- mk_passthrough_data()
  suppressMessages(try(
    flexybayes(
      fixed = y ~ x,
      random = ~g,
      data = d,
      backend = "brms",
      n_samples = 50L,
      warmup = 50L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ),
    silent = TRUE
  ))

  # brms::brm(seed = NULL) is not the same as omitting the argument;
  # brms's own default is NA. The mapping matters because a NULL would
  # error inside rstan rather than drawing a seed.
  expect_true(is.na(captured$seed))
  expect_null(captured$control)
})

# ---------------------------------------------------------------- #
# The engines that do not consume them say so                      #
# ---------------------------------------------------------------- #

test_that("the INLA path reports seed and control as no-ops", {
  skip_if_not_installed("INLA")
  d <- mk_passthrough_data()

  msgs <- testthat::capture_messages(
    fit <- flexybayes(
      fixed = y ~ x,
      random = ~g,
      data = d,
      backend = "inla",
      seed = 99L,
      verbose = FALSE
    )
  )
  expect_s3_class(fit, "flexybayes_inla")
  expect_true(any(grepl("INLA backend does not use", msgs, fixed = TRUE)))
  expect_true(any(grepl("draws no random numbers", msgs, fixed = TRUE)))
})

test_that("the sampler-argument note names the engine and is silenceable", {
  expect_message(
    flexyBayes:::.note_sampler_args_ignored("INLA", seed = 7L),
    "INLA backend does not use"
  )
  expect_message(
    flexyBayes:::.note_sampler_args_ignored(
      "INLA",
      seed = 7L,
      control = list(adapt_delta = 0.9)
    ),
    "`seed` and `control` are"
  )
  expect_message(
    flexyBayes:::.note_sampler_args_ignored("INLA"),
    NA
  )
  withr::local_options(flexyBayes.silence_sampler_arg_note = TRUE)
  expect_message(
    flexyBayes:::.note_sampler_args_ignored("INLA", seed = 7L),
    NA
  )
})

# ---------------------------------------------------------------- #
# Live -- the same seed returns the same posterior                  #
# ---------------------------------------------------------------- #

test_that("two brms fits under one seed return identical draws", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- mk_passthrough_data()

  # Live brms fits are wrapped as the rest of the suite wraps them: brms's
  # ESS advisory fires on the group-SD of a modest design and would
  # otherwise leave the suite non-clean. The diagnostics are asserted
  # below instead of being taken on the advisory's word.
  fit_one <- suppressWarnings(suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 2L,
    seed = 20260815L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  fit_two <- suppressWarnings(suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 2L,
    seed = 20260815L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))

  draws_one <- posterior::as_draws_matrix(fit_one$brms$fit)
  draws_two <- posterior::as_draws_matrix(fit_two$brms$fit)

  expect_identical(dim(draws_one), dim(draws_two))
  expect_identical(colnames(draws_one), colnames(draws_two))
  # Bit-for-bit, not "close": a pinned seed that only narrows the spread
  # is not a pinned seed.
  expect_equal(
    as.numeric(draws_one),
    as.numeric(draws_two),
    tolerance = 0
  )

  # The seed is recorded on the fit, so update() repeats it.
  expect_identical(fit_one$extras$call_info$seed, 20260815L)

  # The fit is usable, not merely repeatable: identical draws from two
  # broken chains would satisfy the comparison above on its own.
  rhat <- brms::rhat(fit_one$brms)
  expect_lt(max(rhat[is.finite(rhat)]), 1.05)
})

# Tests for fb_engine() -- the inference-engine constructor noun
# (v0.4.0 Wave 1 Phase 1D).

suppressPackageStartupMessages(library(testthat))


# ---------------------------------------------------------------- #
# Construction + shape                                              #
# ---------------------------------------------------------------- #

test_that("fb_engine() returns the documented classed-object shape", {
  e <- fb_engine("brms", chains = 4L)
  expect_s3_class(e, "fb_engine")
  expect_true(is_fb_engine(e))
  expect_identical(class(e), c("fb_engine", "list"))
  expect_setequal(names(e), c("name", "paradigm", "toolchain_status", "opts"))
  expect_identical(e$name, "brms")
  expect_identical(e$opts, list(chains = 4L))
})

test_that("fb_engine() maps each engine to its inference paradigm", {
  expect_identical(fb_engine("inla")$paradigm, "laplace")
  expect_identical(fb_engine("brms")$paradigm, "mcmc")
})

test_that("fb_engine() reports toolchain status from package availability", {
  e <- fb_engine("brms")
  expect_true(
    e$toolchain_status %in%
      c("ready", "requires_install", "unavailable")
  )
  # brms is a real dependency; status reflects whether it is installed here.
  expected <- if (requireNamespace("brms", quietly = TRUE)) {
    "ready"
  } else {
    "requires_install"
  }
  expect_identical(e$toolchain_status, expected)
})

test_that("fb_engine() merges opts = list(...) and ... options", {
  e <- fb_engine("brms", opts = list(n_samples = 100L), chains = 2L)
  expect_identical(e$opts, list(n_samples = 100L, chains = 2L))
})


# ---------------------------------------------------------------- #
# Validation                                                        #
# ---------------------------------------------------------------- #

test_that("fb_engine() rejects unknown engines and 'auto'", {
  expect_error(fb_engine("stan"), "unknown engine")
  # 'auto' is a routing directive, not an engine.
  expect_error(fb_engine("auto"), "routing directive")
})

test_that("fb_engine() rejects unrecognised and unnamed options", {
  expect_error(fb_engine("brms", foo = 1), "unrecognised option")
  expect_error(fb_engine("brms", opts = list(3)), "must be named")
  expect_error(fb_engine("brms", opts = "x"), "must be a named list")
})

test_that("fb_engine() rejects a non-string name", {
  expect_error(fb_engine(1L), "non-empty single string")
  expect_error(fb_engine(c("inla", "brms")), "non-empty single string")
})


# ---------------------------------------------------------------- #
# print                                                             #
# ---------------------------------------------------------------- #

test_that("print.fb_engine() shows name, paradigm, status, and opts", {
  out <- utils::capture.output(print(fb_engine("brms", chains = 4L)))
  expect_true(any(grepl("<fb_engine> brms \\(mcmc,", out)))
  expect_true(any(grepl("chains = 4", out)))
})


# ---------------------------------------------------------------- #
# Internal resolvers                                                #
# ---------------------------------------------------------------- #

test_that(".resolve_engine_string() resolves fb_engine and passes strings", {
  expect_identical(
    flexyBayes:::.resolve_engine_string(fb_engine("inla")),
    "inla"
  )
  expect_identical(flexyBayes:::.resolve_engine_string("brms"), "brms")
  expect_identical(
    flexyBayes:::.resolve_engine_string(c("brms", "inla", "auto")),
    c("brms", "inla", "auto")
  )
})

test_that(".fb_engine_opts() returns opts for fb_engine, NULL otherwise", {
  expect_identical(
    flexyBayes:::.fb_engine_opts(fb_engine("brms", chains = 2L)),
    list(chains = 2L)
  )
  expect_null(flexyBayes:::.fb_engine_opts("brms"))
  expect_null(flexyBayes:::.fb_engine_opts(fb_engine("brms")))
})


# ---------------------------------------------------------------- #
# Consumption by the fitting verbs                                  #
# ---------------------------------------------------------------- #

test_that("flexybayes(backend = fb_engine('inla')) routes to INLA", {
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE
  )
  set.seed(1L)
  d <- data.frame(
    y = rnorm(40),
    x = rnorm(40),
    g = factor(rep(1:5, length.out = 40))
  )
  fit <- suppressMessages(flexybayes(
    fixed = y ~ x + (1 | g),
    data = d,
    backend = fb_engine("inla"),
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_identical(backend_decision(fit)$backend, "inla")
})

test_that("fb_engine() opts override the sampler controls on brms", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  )
  set.seed(1L)
  d <- data.frame(y = rnorm(40), x = rnorm(40))
  # Structural test (does the draw count match the overridden sampler
  # controls?), not a numeric-agreement one. A 60-draw / chains=1L /
  # unset-seed budget produced a genuine R-hat warning (not just ESS,
  # which .muffle_ess_warnings() deliberately leaves unmuffled) --
  # chains must stay at 1L for the exact draw-count assertion below,
  # so the fix is a larger, seeded budget that actually converges
  # rather than widening the muffle pattern.
  fit <- suppressMessages(flexybayes(
    fixed = y ~ x,
    data = d,
    backend = fb_engine("brms", chains = 1L),
    n_samples = 300L,
    warmup = 300L,
    seed = 20260523L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  # 300 post-warmup samples x 1 chain -> 300 posterior draws.
  expect_identical(nrow(posterior::as_draws_matrix(fit$brms)), 300L)
})

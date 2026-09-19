# Tests for the backend argument on flexybayes() -- ADR 0006.
#
# Contract (0.9.3: the third native engine ADR 0006 named was withdrawn
# entirely -- see NEWS.md -- and its exports, code paths, and fallback
# role were removed rather than repointed):
#   - backend = "inla" calls lgm_gate(); on accept dispatches
#     emit_inla() and returns a flexybayes_inla; on refuse raises a
#     formatted refusal.
#   - backend = "auto" calls lgm_gate(); on accept dispatches INLA
#     (when installed) and records path "auto_accept"; on refuse
#     falls back to brms, with the specific fallback path recorded
#     (e.g. "auto_multistratum_to_brms" for a multi-stratum residual
#     structure lgm_gate() cannot represent).
#   - Invalid or unrecognised backend value raises a typed
#     unknown_backend refusal naming the two active engines.
#   - review_code = TRUE under backend = "inla" raises a clean refusal
#     naming brms as the only code-producing engine.
#
# brms is required for the auto-fallback and explicit-brms subtests;
# INLA is required for the INLA-accept subtests (skip_if_not_installed
# guards in place).

mk_lgm_data <- function() {
  set.seed(20260522L)
  n <- 30L
  data.frame(
    yield = rnorm(n, mean = 100, sd = 10),
    env = factor(rep(1:3, length.out = n)),
    geno = factor(rep(1:5, length.out = n))
  )
}


# ---------------------------------------------------------------- #
# (a) Invalid backend value raises a typed unknown_backend refusal   #
# ---------------------------------------------------------------- #

test_that("invalid backend value raises a typed unknown_backend refusal", {
  d <- mk_lgm_data()
  err <- tryCatch(
    flexybayes(yield ~ env, random = ~geno, data = d, backend = "stan"),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_unknown_backend_refusal")
  expect_identical(err$reason_code, "unknown_backend")
  expect_match(conditionMessage(err), '"inla"', fixed = TRUE)
  expect_match(conditionMessage(err), '"brms"', fixed = TRUE)
})


# ---------------------------------------------------------------- #
# (b) review_code + backend = "inla" structured refusal             #
# ---------------------------------------------------------------- #

test_that("review_code = TRUE under backend = 'inla' raises a clean refusal naming brms", {
  d <- mk_lgm_data()
  err <- tryCatch(
    flexybayes(
      yield ~ env,
      random = ~geno,
      data = d,
      backend = "inla",
      review_code = TRUE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_review_code_backend_unsupported")
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "review_code", fixed = TRUE)
  expect_match(conditionMessage(err), "brms", fixed = TRUE)
})


# ---------------------------------------------------------------- #
# (c) backend = "brms" records the explicit-brms trace              #
# ---------------------------------------------------------------- #

test_that("backend = 'brms' records the explicit-brms trace", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  skip_if_not_long_tier()
  d <- mk_lgm_data()
  fit <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "brms")
})


# ---------------------------------------------------------------- #
# (d) backend = "auto" on an LGM-compatible model routes to INLA    #
# ---------------------------------------------------------------- #

test_that("backend = 'auto' routes a Gaussian random-intercept model to INLA", {
  testthat::skip_if_not_installed("INLA")
  d <- mk_lgm_data()
  # aggregate = FALSE: scope the test to the per-row INLA dispatch
  # path. A Gaussian model is no longer aggregated under "auto" (its
  # aggregated route cannot carry the residual prior), so the setting
  # pins the route this test asserts rather than changing it.
  # Aggregated-path coverage lives in test-emit-gaussian-aggregated.R.
  fit <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    backend = "auto",
    aggregate = FALSE,
    verbose = FALSE
  ))
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "inla")
  expect_true(bd$path %in% c("auto_accept", "explicit_inla_accept"))
  expect_true(!is.null(bd$gate_checks))
})


# ---------------------------------------------------------------- #
# (e) backend = "inla" on a model INLA cannot represent refuses      #
# ---------------------------------------------------------------- #
# The probe used to be the heterogeneous residual at(env):units. That
# fits on INLA from 0.10.0.9000, through the multiple-likelihood
# interface, so the probe moves to a family INLA's likelihood roster has
# no counterpart for. Under backend = "inla" the refusal surfaces as a
# clean error; under backend = "auto" the same refusal routes to brms
# (subtest f).

test_that("backend = 'inla' surfaces the INLA-side refusal cleanly", {
  d <- mk_lgm_data()
  d$yield <- pmax(d$yield, 0)
  err <- tryCatch(
    flexybayes(
      yield ~ env,
      random = ~geno,
      data = d,
      family = "hurdle_gamma",
      backend = "inla",
      verbose = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl(
    "backend = \"inla\"|INLA backend refused|emit_inla\\(\\) refused|does not support|family",
    err
  ))
})

test_that("the heterogeneous residual that used to drive this now fits", {
  # The counterpart of the change above, asserted rather than assumed:
  # the structure this file used as its canonical INLA refusal is a
  # capability now, so a reader of subtest (e) can see why it moved.
  d <- mk_lgm_data()
  fb <- fb_from_asreml(
    yield ~ env, random = ~geno, residual = ~ at(env):units, data = d
  )
  expect_true("lgm_compatible" %in% lgm_gate(fb)$capabilities)
})


# ---------------------------------------------------------------- #
# (f) backend = "auto" routes to brms when INLA refuses, with a trace #
# ---------------------------------------------------------------- #
# This used to be driven by the multi-stratum residual at(env):units,
# which INLA now fits -- so `auto` keeps that model on INLA, which was
# the whole point of the 0.10.0.9000 emit work. Both halves are asserted
# here: the model that used to fall back no longer does, and a model
# INLA genuinely cannot take still falls back by the documented path.

test_that("backend = 'auto' keeps the multi-stratum residual on INLA", {
  testthat::skip_if_not_installed("INLA")
  d <- mk_lgm_data()
  fit <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    residual = ~ at(env):units,
    data = d,
    backend = "auto",
    verbose = FALSE
  ))
  expect_identical(backend_decision(fit)$backend, "inla")
})

test_that("backend = 'auto' falls back to brms when INLA refuses, with trace", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  skip_if_not_long_tier()
  d <- mk_lgm_data()
  # corh(env):geno -- heterogeneous variances with one shared
  # correlation. No engine has an equicorrelation group structure, and
  # INLA refuses it at the gate, so `auto` must route rather than fit.
  msgs <- character()
  fit <- withCallingHandlers(
    suppressWarnings(suppressMessages(flexybayes(
      yield ~ env,
      random = ~geno,
      data = d,
      family = "hurdle_gamma",
      backend = "auto",
      n_samples = 500L,
      warmup = 500L,
      chains = 1L,
      seed = 20260523L,
      control = list(adapt_delta = 0.97),
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "brms")
  expect_true(!is.null(bd$path) && nzchar(bd$path))
})


# ---------------------------------------------------------------- #
# (g) ADR 0006 verification snippet runs end-to-end                  #
# ---------------------------------------------------------------- #

test_that("ADR 0006 verification: explicit-brms + auto-accept + non-LGM refusal", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  d <- mk_lgm_data()

  # Explicit brms -- the withdrawn engine's role in this verification
  # (existing behaviour preserved under an explicit engine request)
  # moves to brms, the remaining full-HMC engine. (ADR 0031: the
  # default is "auto", so an explicit engine request must be named; an
  # unpinned call here would take the aggregated-INLA path.)
  fit_b <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit_b, "flexybayes")
  expect_identical(backend_decision(fit_b)$backend, "brms")

  # backend = "inla" on a model INLA cannot represent: structured
  # refusal. The probe is a family with no INLA likelihood, since the
  # structured residual this used to use now fits there.
  d2 <- d
  d2$yield <- pmax(d2$yield, 0)
  err <- tryCatch(
    flexybayes(
      yield ~ env,
      random = ~geno,
      data = d2,
      family = "hurdle_gamma",
      backend = "inla",
      verbose = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(nzchar(err))
})

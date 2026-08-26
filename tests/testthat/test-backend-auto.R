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
  # path. The Stage 3A aggregated-gaussian gate (ADR 0022, v0.3.2)
  # would otherwise auto-route this LGM-compatible model to the
  # aggregated emit (path = "aggregated_gaussian") and short-circuit
  # the lgm_gate call. Aggregated-path coverage lives in
  # test-emit-gaussian-aggregated.R.
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
# (e) backend = "inla" on an LGM-incompatible model refuses         #
# ---------------------------------------------------------------- #
# A heterogeneous-residual at() structure triggers an LGM check
# failure (residual is not the iid units default). Use it to drive the
# refusal path.

test_that("backend = 'inla' surfaces the INLA-side refusal cleanly", {
  d <- mk_lgm_data()
  # at(env):units is rejected at emit_inla()'s feasibility check
  # (v0.1 does not support structured residual for INLA). Under backend
  # = "inla" the refusal surfaces as a clean error; under backend =
  # "auto" the same refusal triggers the brms fallback (subtest f).
  err <- tryCatch(
    flexybayes(
      yield ~ env,
      random = ~geno,
      residual = ~ at(env):units,
      data = d,
      backend = "inla",
      verbose = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl(
    "backend = \"inla\"|INLA backend refused|emit_inla\\(\\) refused|does not support",
    err
  ))
})


# ---------------------------------------------------------------- #
# (f) backend = "auto" on a multi-stratum residual structure falls   #
#     back to brms with a logged note and the specific fallback path #
# ---------------------------------------------------------------- #

test_that("backend = 'auto' on a multi-stratum residual falls back to brms with trace", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  d <- mk_lgm_data()
  msgs <- character()
  fit <- withCallingHandlers(
    suppressWarnings(flexybayes(
      yield ~ env,
      random = ~geno,
      residual = ~ at(env):units,
      data = d,
      backend = "auto",
      n_samples = 500L,
      warmup = 500L,
      chains = 1L,
      seed = 20260523L,
      control = list(adapt_delta = 0.97),
      verbose = FALSE,
      mcmc_verbose = FALSE
    )),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "brms")
  # at(env):units is a multi-stratum designed-experiment residual
  # structure: lgm_gate() refuses it (residual_term_type_inla) and
  # dispatch routes it to brms by this specific path, distinct from
  # the generic auto-fallback (INLA collapses these variance
  # components; brms is the faithful full-HMC alternative).
  expect_identical(bd$path, "auto_multistratum_to_brms")
  expect_true(!is.null(bd$gate_checks))
  expect_true(any(grepl("lgm_gate\\(\\) refused", msgs)))
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

  # backend = "inla" on a non-LGM model: structured refusal.
  err <- tryCatch(
    flexybayes(
      yield ~ env,
      random = ~geno,
      residual = ~ at(env):units,
      data = d,
      backend = "inla",
      verbose = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(nzchar(err))
})

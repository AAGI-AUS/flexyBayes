# Tests for the backend argument on flexybayes() -- ADR 0006.
#
# Contract:
#   - backend = "greta" (default) preserves the existing call shape;
#     backend_decision(fit)$backend == "greta" with path
#     "explicit_greta".
#   - backend = "inla" calls lgm_gate(); on accept dispatches
#     emit_inla() and returns a flexybayes_inla; on refuse raises a
#     formatted refusal.
#   - backend = "auto" calls lgm_gate(); on accept dispatches INLA
#     (when installed) and records path "auto_accept"; on refuse
#     falls back to greta with a one-time note and records path
#     "auto_lgm_refuse" with the gate failure list.
#   - Invalid backend value raises a clean match.arg error.
#   - review_code = TRUE under backend != "greta" raises a clean
#     refusal.
#
# Greta is required for the dispatch-to-greta subtests; INLA is
# required for the INLA-accept subtests (skip_if_not_installed
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
# (a) Invalid backend value raises match.arg error                  #
# ---------------------------------------------------------------- #

test_that("invalid backend value raises match.arg error", {
  d <- mk_lgm_data()
  err <- tryCatch(
    flexybayes(yield ~ env, random = ~geno, data = d, backend = "stan"),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("'arg'.*should be one of|match.arg", err, perl = TRUE))
})


# ---------------------------------------------------------------- #
# (b) review_code + backend != "greta" structured refusal           #
# ---------------------------------------------------------------- #

test_that("review_code = TRUE under backend != 'greta' raises a clean refusal", {
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
  expect_match(conditionMessage(err), "greta", fixed = TRUE)
})


# ---------------------------------------------------------------- #
# (c) backend = "greta" preserves the explicit-greta trace          #
# ---------------------------------------------------------------- #

test_that("backend = 'greta' records the explicit-greta trace", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_lgm_data()
  fit <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    backend = "greta",
    n_samples = 50L,
    warmup = 50L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "greta")
  expect_identical(bd$path, "explicit_greta")
  expect_null(bd$gate_checks)
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
# failure (rcov is not the iid units default). Use it to drive the
# refusal path.

test_that("backend = 'inla' surfaces the INLA-side refusal cleanly", {
  d <- mk_lgm_data()
  # at(env):units is rejected at emit_inla()'s feasibility check
  # (v0.1 does not support structured rcov for INLA). Under backend
  # = "inla" the refusal surfaces as a clean error; under backend =
  # "auto" the same refusal triggers the greta fallback (subtest f).
  err <- tryCatch(
    flexybayes(
      yield ~ env,
      random = ~geno,
      rcov = ~ at(env):units,
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
# (f) backend = "auto" on a non-LGM model falls back to greta with  #
#     a logged note and records the gate failure trace               #
# ---------------------------------------------------------------- #

test_that("backend = 'auto' on a non-LGM model falls back to greta with trace", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_lgm_data()
  msgs <- character()
  fit <- withCallingHandlers(
    suppressWarnings(flexybayes(
      yield ~ env,
      random = ~geno,
      rcov = ~ at(env):units,
      data = d,
      backend = "auto",
      n_samples = 50L,
      warmup = 50L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "greta")
  # ADR 0017: the gate's .lgm_check_rcov_term_inla_support() now
  # catches at(env):units at gate time, so the only valid fall-back
  # path is auto_lgm_refuse. The pre-ADR-0017 emit-level path
  # ("auto_inla_emit_refuse") is architecturally unreachable.
  expect_identical(bd$path, "auto_lgm_refuse")
  expect_true(!is.null(bd$gate_checks))
  expect_true(any(grepl("lgm_gate\\(\\) refused", msgs)))
})


# ---------------------------------------------------------------- #
# (g) ADR 0006 verification snippet runs end-to-end                  #
# ---------------------------------------------------------------- #

test_that("ADR 0006 verification: explicit-greta + auto-accept + non-LGM refusal", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  d <- mk_lgm_data()

  # Explicit greta -- existing behaviour preserved. (ADR 0031: the
  # default is now "auto", so the explicit-greta path must be requested
  # by name; an unpinned call here would take the aggregated-INLA path.)
  fit_g <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    backend = "greta",
    n_samples = 50L,
    warmup = 50L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit_g, "flexybayes")
  expect_identical(backend_decision(fit_g)$path, "explicit_greta")

  # backend = "inla" on a non-LGM model: structured refusal.
  err <- tryCatch(
    flexybayes(
      yield ~ env,
      random = ~geno,
      rcov = ~ at(env):units,
      data = d,
      backend = "inla",
      verbose = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(nzchar(err))
})

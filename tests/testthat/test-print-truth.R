# The "Representation:" and "Engine:" lines (ADR 0030 Decision 4) on
# print.flexybayes(), summary.flexybayes() and print.flexybayes_aggregated():
#
#   (a) snapshot of print() on an INLA Laplace fit -> Engine: INLA Laplace
#   (b) snapshot of print() on a brms / Stan HMC fit -> Engine: brms / Stan HMC
#   (c) snapshot of print() on an aggregated-exact fit ->
#       Representation: aggregated_exact (compression 5:1)
#   (d) .engine_label_for_fit() returns the raw backend string for any
#       engine name other than "inla" or "brms"
#
# The helpers are tested on synthesised fit-shaped objects, so the tests
# need no engine and do not depend on engine versions.

# Helper: build a minimal flexybayes-shaped object carrying just the
# slots the truth-display helpers consume ($exactness +
# $extras$backend_decision + $extras$aggregation_meta).
.test_truth_fit <- function(
  backend,
  path,
  exactness = "exact",
  N = NULL,
  K = NULL
) {
  am <- if (!is.null(N) && !is.null(K)) {
    list(N = N, K = K)
  } else {
    NULL
  }
  structure(
    list(
      exactness = exactness,
      extras = list(
        backend_decision = list(backend = backend, path = path),
        aggregation_meta = am
      )
    ),
    class = "flexybayes"
  )
}

# ---------------------------------------------------------------- #
# (a) unrecognised engine name falls through verbatim                #
# ---------------------------------------------------------------- #

test_that("Engine: falls through to the raw backend string for an unrecognised engine", {
  fit <- .test_truth_fit(backend = "future_engine", path = "explicit_future_engine")
  repr <- flexyBayes:::.repr_label_for_fit(fit, fit$extras$backend_decision)
  engine <- flexyBayes:::.engine_label_for_fit(fit, fit$extras$backend_decision)
  expect_equal(repr, "exact (model, not inference)")
  expect_equal(engine, "future_engine")
})

# ---------------------------------------------------------------- #
# (b) INLA Laplace fit                                               #
# ---------------------------------------------------------------- #

test_that("Engine: INLA Laplace for an INLA-backed fit", {
  fit <- .test_truth_fit(backend = "inla", path = "explicit_inla_accept")
  repr <- flexyBayes:::.repr_label_for_fit(fit, fit$extras$backend_decision)
  engine <- flexyBayes:::.engine_label_for_fit(fit, fit$extras$backend_decision)
  expect_equal(repr, "exact (model, not inference)")
  expect_equal(engine, "INLA Laplace")
})

# ---------------------------------------------------------------- #
# (c) aggregated-exact fit                                           #
# ---------------------------------------------------------------- #

test_that("Representation: aggregated_exact (compression 5:1) for aggregated fit", {
  fit <- .test_truth_fit(
    backend = "inla",
    path = "aggregated_inla",
    exactness = "aggregated_exact",
    N = 1000L,
    K = 200L
  )
  repr <- flexyBayes:::.repr_label_for_fit(fit, fit$extras$backend_decision)
  engine <- flexyBayes:::.engine_label_for_fit(fit, fit$extras$backend_decision)
  expect_equal(repr, "aggregated_exact (compression 5:1)")
  expect_equal(engine, "INLA Laplace (aggregated)")
})

# ---------------------------------------------------------------- #
# Bonus: brms / Stan HMC label                                       #
# ---------------------------------------------------------------- #

test_that("Engine: brms / Stan HMC for a brms-backed fit", {
  fit <- .test_truth_fit(backend = "brms", path = "explicit_brms")
  engine <- flexyBayes:::.engine_label_for_fit(fit, fit$extras$backend_decision)
  expect_equal(engine, "brms / Stan HMC")
})

# ---------------------------------------------------------------- #
# Bonus: adjacent lines + no Exact.: line on print.flexybayes()      #
# ---------------------------------------------------------------- #

test_that("print.flexybayes emits Representation:/Engine: adjacent (no Exact.:)", {
  # Synthesise a fuller fit so print.flexybayes() does not stumble on
  # the convergence + Params lines. Only the truth-display block is
  # checked.
  fit <- structure(
    list(
      exactness = "exact",
      extras = list(
        backend_decision = list(backend = "brms", path = "explicit_brms"),
        call_info = list(
          fixed = y ~ x,
          random = NULL,
          residual = NULL,
          chains = 2,
          n_samples = 100,
          warmup = 50
        ),
        model_info = list(
          family = "gaussian",
          link = "identity",
          n_params = 3L,
          n_fixed = 2L,
          n_random = 1L
        ),
        run_time = 1.0,
        convergence = list(gelman = NULL, n_eff = NULL),
        aggregation_meta = NULL
      )
    ),
    class = "flexybayes"
  )
  out <- utils::capture.output(print(fit))
  has_repr <- grep("Representation:", out, fixed = TRUE)
  has_engine <- grep("Engine:", out, fixed = TRUE)
  expect_length(has_repr, 1L)
  expect_length(has_engine, 1L)
  expect_equal(has_engine[1L] - has_repr[1L], 1L)
  # And no legacy single-line Exact.: rendering.
  expect_length(grep("Exact\\. :", out), 0L)
})

# test-plan-live-parity.R -- the plan surface must never name a route
# dispatch would not take.
#
# The 2026-08-16 adversarial review found `plan = TRUE` announcing
# `Backend chosen: inla / Path: aggregated_inla` on a
# `residual = ~ dsum(~ units | env)` model in the same breath as
# `Gate outcome: refuse_structural (residual_term_type_inla)`. The live
# call routed to brms, as it should. Two independent defects produced
# that: `.fb_aggregation_plan()` treated a sectioned residual as
# aggregatable, and `fb_plan()`'s aggregation override did not require
# the INLA gate to have accepted the model.
#
# These tests pin the contract from both ends -- the refused shape must
# plan and fit the same way, and a genuinely aggregation-eligible model
# must still plan and fit the aggregated INLA route.

# ---------------------------------------------------------------- #
# Fixtures                                                           #
# ---------------------------------------------------------------- #

# A sectioned-residual (dsum) shape: three environments with genuinely
# different residual spread, replicated. This is the reviewer's probe.
.parity_dsum_data <- function(n_rep = 20L, seed = 20260816L) {
  set.seed(seed)
  d <- data.frame(env = factor(rep(c("e1", "e2", "e3"), each = n_rep)))
  d$y <- stats::rnorm(
    nrow(d),
    mean = c(1, 3, 6)[as.integer(d$env)],
    sd = c(0.3, 1.1, 2.2)[as.integer(d$env)]
  )
  d
}

# A model the aggregated INLA route genuinely fits: factor-only fixed
# and random terms, 1000 rows compressing to 100 cells.
.parity_aggregatable_data <- function(n = 1000L, seed = 20260816L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(sample(paste0("f", seq_len(5L)), n, replace = TRUE)),
    g = factor(sample(paste0("g", seq_len(20L)), n, replace = TRUE))
  )
  d$y <- stats::rnorm(n)
  d
}


# ---------------------------------------------------------------- #
# (1) The plan names the backend the live fit uses                   #
# ---------------------------------------------------------------- #

test_that("plan for a dsum residual names brms, never aggregated_inla", {
  d <- .parity_dsum_data()
  p <- flexybayes(
    y ~ env,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "auto",
    plan = TRUE
  )

  expect_s3_class(p, "fb_plan")
  expect_identical(p$backend_chosen, "brms")
  expect_false(identical(p$path, "aggregated_inla"))
  # The gate verdict on the plan is the one the live INLA route reports.
  expect_identical(p$gate_outcome, "refuse_structural")
  expect_identical(p$gate_primary_rule, "residual_term_type_inla")
})

test_that("a sectioned residual is not aggregatable at plan level", {
  d <- .parity_dsum_data()
  p <- flexybayes(
    y ~ env,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "auto",
    plan = TRUE
  )
  expect_false(isTRUE(p$aggregation$eligible))
  expect_true(
    "structured_residual_not_aggregatable" %in% p$aggregation$reason_codes
  )
})

test_that("explicit backend = 'inla' refuses the dsum residual typed", {
  skip_if_not_installed("INLA")
  d <- .parity_dsum_data()
  err <- tryCatch(
    flexybayes(
      y ~ env,
      residual = ~ dsum(~ units | env),
      data = d,
      backend = "inla",
      verbose = FALSE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal")
  expect_s3_class(err, "flexybayes_lgm_residual_term_type_inla")
  expect_identical(err$reason_code, "inla_gate_refused")
  # The refusal names the engine that does represent the model.
  expect_match(conditionMessage(err), "backend = \"brms\"", fixed = TRUE)
})

test_that("aggregate = TRUE on a dsum residual refuses typed", {
  skip_if_not_installed("INLA")
  d <- .parity_dsum_data()
  err <- tryCatch(
    flexybayes(
      y ~ env,
      residual = ~ dsum(~ units | env),
      data = d,
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_aggregate_refusal")
})


# ---------------------------------------------------------------- #
# (2) The live fit under `auto` takes the route the plan named       #
# ---------------------------------------------------------------- #

test_that("live auto fit of the dsum shape lands on brms as planned", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")
  d <- .parity_dsum_data()
  p <- flexybayes(
    y ~ env,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "auto",
    plan = TRUE
  )
  fit <- suppressMessages(flexybayes(
    y ~ env,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "auto",
    n_samples = 400L,
    warmup = 400L,
    chains = 2L,
    verbose = FALSE
  ))
  expect_identical(fit$extras$backend_decision$backend, p$backend_chosen)
  expect_identical(fit$extras$backend_decision$backend, "brms")
})


# ---------------------------------------------------------------- #
# (3) Positive control -- the aggregated INLA route still runs        #
# ---------------------------------------------------------------- #

test_that("an aggregation-eligible model still plans aggregated_inla", {
  skip_if_not_installed("INLA")
  d <- .parity_aggregatable_data()
  p <- flexybayes(
    y ~ f,
    random = ~g,
    data = d,
    backend = "auto",
    plan = TRUE
  )
  expect_identical(p$gate_outcome, "accept")
  expect_true(isTRUE(p$aggregation$eligible))
  expect_identical(p$backend_chosen, "inla")
  expect_identical(p$path, "aggregated_inla")
})

test_that("an aggregation-eligible model still fits the aggregated route", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .parity_aggregatable_data()
  fit <- suppressMessages(flexybayes(
    y ~ f,
    random = ~g,
    data = d,
    backend = "auto",
    aggregate = TRUE,
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")
  expect_identical(fit$extras$backend_decision$backend, "inla")
})

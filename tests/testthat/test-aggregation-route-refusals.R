# =============================================================================
# The aggregated-route entry refusals in .maybe_aggregate_gaussian() raise
# typed conditions: missing response, backend without an aggregated emit, and
# no resolvable route. Each also carries the family class
# `flexybayes_aggregate_refusal`, so a handler for the family catches all.
#
# A missing response (no sufficient statistic for the cell) and an
# unavailable route (no aggregated emit reachable) have separate reason codes.
# =============================================================================

.agg_route_trial <- function(seed = 20260819L) {
  set.seed(seed)
  d <- expand.grid(
    rep = seq_len(20L),
    trt = factor(seq_len(5L)),
    blk = factor(seq_len(4L))
  )
  d$y <- 10 + as.numeric(d$trt) * 0.5 + stats::rnorm(nrow(d), 0, 0.3)
  # These route tests need a model the GATE refuses on INLA so that no
  # aggregated route resolves. Any us() is that model: estimated correlation
  # matrices are fitted on brms only, and the gate refuses them
  # (.lgm_check_correlated_covariance()). Six levels, which also exceed
  # INLA's iidKd range, date from when us() fitted on INLA up to five. The
  # mechanism under test is unchanged: gate refusal above, plan-eligibility
  # gate below.
  d$g2 <- factor(rep(seq_len(6L), length.out = nrow(d)))
  d
}

# ---------------------------------------------------------------- #
# 1. Site 1: a missing response has no sufficient statistic        #
# ---------------------------------------------------------------- #

test_that("aggregate = TRUE on a missing response raises a typed refusal", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()
  d$y[c(3L, 17L)] <- NA_real_

  expect_error(
    suppressMessages(flexybayes(
      fixed = y ~ trt,
      random = ~blk,
      residual = ~units,
      data = d,
      backend = "inla",
      na_action = "augment",
      aggregate = TRUE,
      verbose = FALSE
    )),
    class = "flexybayes_refusal_aggregation_response_incomplete"
  )
})

test_that("the missing-response refusal keeps the family class and wording", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()
  d$y[c(3L, 17L)] <- NA_real_

  err <- tryCatch(
    suppressMessages(flexybayes(
      fixed = y ~ trt,
      random = ~blk,
      residual = ~units,
      data = d,
      backend = "inla",
      na_action = "augment",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  # A caller pattern-matching the aggregated-route family still catches it.
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_s3_class(err, "flexybayes_refusal")
  # Users match on this text, so it is pinned.
  expect_match(
    conditionMessage(err),
    "the response has missing values",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "per-cell sufficient", fixed = TRUE)
  expect_identical(err$reason_code, "aggregation_response_incomplete")
})

# ---------------------------------------------------------------- #
# 2. Site 2: a backend with no aggregated emit at all              #
# ---------------------------------------------------------------- #
#
# This is the capability matrix's one `n/a` cell: aggregation on brms. The
# execution grid asserts that the request does not silently produce the
# structure, which means it must refuse, and refuse typed.

test_that("aggregate = TRUE on brms raises a typed refusal", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()

  expect_error(
    suppressMessages(flexybayes(
      fixed = y ~ trt,
      random = ~blk,
      data = d,
      backend = "brms",
      aggregate = TRUE,
      verbose = FALSE
    )),
    class = "flexybayes_refusal_aggregation_route_unavailable"
  )
})

test_that("the brms aggregation refusal names the backend and its wording", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()

  err <- tryCatch(
    suppressMessages(flexybayes(
      fixed = y ~ trt,
      random = ~blk,
      data = d,
      backend = "brms",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_identical(err$backend, "brms")
  expect_match(
    conditionMessage(err),
    "is not supported on backend = \"brms\"",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "switch backend", fixed = TRUE)
})

# ---------------------------------------------------------------- #
# 3. Site 3: no aggregated backend resolves for this model         #
# ---------------------------------------------------------------- #
#
# Reached without mocking: on backend = "auto" the resolver consults
# lgm_gate(), and an unstructured-covariance random effect is refused there,
# so no aggregated route resolves. This guard sits ABOVE the plan-eligibility
# gate, so the plan refusal cannot pre-empt it.

test_that("aggregate = TRUE with no resolvable route raises a typed refusal", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()

  expect_error(
    suppressMessages(flexybayes(
      fixed = y ~ trt,
      random = ~ us(g2):blk,
      data = d,
      backend = "auto",
      aggregate = TRUE,
      verbose = FALSE
    )),
    class = "flexybayes_refusal_aggregation_route_unavailable"
  )
})

test_that("the no-route refusal keeps the family class and its wording", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()

  err <- tryCatch(
    suppressMessages(flexybayes(
      fixed = y ~ trt,
      random = ~ us(g2):blk,
      data = d,
      backend = "auto",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_match(
    conditionMessage(err),
    "no active aggregated backend",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(err),
    "INLA refused or unavailable",
    fixed = TRUE
  )
  expect_identical(err$reason_code, "aggregation_route_unavailable")
})

# ---------------------------------------------------------------- #
# 4. The message text did not move when the class arrived           #
# ---------------------------------------------------------------- #
#
# Callers match on the message text, so the whole string is pinned. A
# deliberate rewording must edit this test.

test_that("the two route refusals carry their 0.9.3 message verbatim", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()

  err_brms <- tryCatch(
    suppressMessages(flexybayes(
      fixed = y ~ trt,
      random = ~blk,
      data = d,
      backend = "brms",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_identical(
    conditionMessage(err_brms),
    paste0(
      "`aggregate = TRUE` is not supported on backend = \"brms\" (the ",
      "aggregated path is wired for inla only). Pass ",
      "aggregate = FALSE or switch backend."
    )
  )

  err_route <- tryCatch(
    suppressMessages(flexybayes(
      fixed = y ~ trt,
      random = ~ us(g2):blk,
      data = d,
      backend = "auto",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_identical(
    conditionMessage(err_route),
    paste0(
      "`aggregate = TRUE` refused: no active aggregated backend (INLA ",
      "refused or unavailable). Pass aggregate = FALSE for the ",
      "per-row path."
    )
  )
})

test_that("the missing-response refusal carries its message verbatim", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()
  d$y[c(3L, 17L)] <- NA_real_

  err <- tryCatch(
    suppressMessages(flexybayes(
      fixed = y ~ trt,
      random = ~blk,
      residual = ~units,
      data = d,
      backend = "inla",
      na_action = "augment",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_identical(
    conditionMessage(err),
    paste0(
      "`aggregate = TRUE` refused: the response has missing values, and ",
      "the aggregated path compresses rows into per-cell sufficient ",
      "statistics that are not defined for a missing response. Pass ",
      "aggregate = FALSE (or leave it at \"auto\") for the per-row path, ",
      "which carries a missing response as a latent quantity."
    )
  )
})

# ---------------------------------------------------------------- #
# 5. Registry and documentation surfaces                            #
# ---------------------------------------------------------------- #

test_that("both codes are registered and reachable from fb_refusals()", {
  for (code in c(
    "aggregation_response_incomplete",
    "aggregation_route_unavailable"
  )) {
    entry <- flexyBayes:::.lookup_refusal(code)
    expect_type(entry, "list")
    expect_identical(entry$since_version, "0.9.2")

    tbl <- fb_refusals(reason_code = code)
    expect_s3_class(tbl, "fb_refusals_table")
    expect_identical(nrow(tbl), 1L)
    expect_match(tbl$description, "aggregate = TRUE", fixed = TRUE)
  }
})

test_that("the two codes describe two different boundaries", {
  tbl <- fb_refusals(
    reason_code = c(
      "aggregation_response_incomplete",
      "aggregation_route_unavailable"
    )
  )
  expect_identical(nrow(tbl), 2L)
  expect_match(
    tbl$description[tbl$reason_code == "aggregation_response_incomplete"],
    "sufficient",
    fixed = TRUE
  )
  expect_match(
    tbl$description[tbl$reason_code == "aggregation_route_unavailable"],
    "no aggregated emit route",
    fixed = TRUE
  )
})

# ---------------------------------------------------------------- #
# 6. Site 4: a family or link the aggregated route cannot take     #
# ---------------------------------------------------------------- #
#
# .fb_aggregation_plan() checks family and link before term shapes, so a
# family or link reason is always plan$reason_codes[1].
#
# A non-canonical link (binomial(link = "probit"), poisson(link = "sqrt"))
# raises inla_gate_refused, as on the per-row route. The aggregated links
# (identity, logit, log) are INLA's defaults for the three families.
# A family outside {gaussian, binomial, poisson} raises
# aggregation_route_unavailable. Gamma, for example, fits on INLA per row
# but has no aggregated emit.

.agg_link_trial <- function(seed = 1L) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(6L), each = 8L)))
  eta <- 0.3 + stats::rnorm(6L, 0, 0.3)[as.integer(d$g)]
  d$y_bin <- stats::rbinom(48L, 1L, stats::plogis(eta))
  d$y_poi <- stats::rpois(48L, exp(eta))
  d$y_gam <- stats::rgamma(48L, shape = 2, rate = 2 * exp(-eta))
  d
}

test_that("aggregate = TRUE with a non-canonical link raises inla_gate_refused", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_link_trial()

  d$y <- d$y_bin
  expect_error(
    suppressMessages(flexybayes(
      y ~ 1 + (1 | g),
      data = d,
      family = stats::binomial(link = "probit"),
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    )),
    class = "flexybayes_refusal_inla_gate_refused"
  )

  d$y <- d$y_poi
  expect_error(
    suppressMessages(flexybayes(
      y ~ 1 + (1 | g),
      data = d,
      family = stats::poisson(link = "sqrt"),
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    )),
    class = "flexybayes_refusal_inla_gate_refused"
  )
})

test_that("aggregate = TRUE with a family the aggregated route cannot take raises aggregation_route_unavailable", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_link_trial()
  d$y <- d$y_gam

  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ 1 + (1 | g),
      data = d,
      family = "gamma",
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_aggregation_route_unavailable")
  expect_s3_class(err, "flexybayes_refusal")
  expect_identical(err$reason_code, "aggregation_route_unavailable")
  expect_match(conditionMessage(err), "gamma", fixed = TRUE)

  # The family scope check does not depend on which backend was asked
  # for, so "auto" reaches the same reason.
  err_auto <- tryCatch(
    suppressMessages(flexybayes(
      y ~ 1 + (1 | g),
      data = d,
      family = "gamma",
      backend = "auto",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(
    err_auto,
    "flexybayes_refusal_aggregation_route_unavailable"
  )
})

test_that("flexybayes_stream() with a non-canonical link or family refuses typed", {
  d <- .agg_link_trial()

  d_bin <- d
  d_bin$y <- d_bin$y_bin
  expect_error(
    flexybayes_stream(
      fixed = y ~ 1,
      random = ~g,
      source = d_bin,
      family = stats::binomial(link = "probit"),
      fit = FALSE,
      verbose = FALSE
    ),
    class = "flexybayes_refusal_aggregation_route_unavailable"
  )

  d_poi <- d
  d_poi$y <- d_poi$y_poi
  expect_error(
    flexybayes_stream(
      fixed = y ~ 1,
      random = ~g,
      source = d_poi,
      family = stats::poisson(link = "sqrt"),
      fit = FALSE,
      verbose = FALSE
    ),
    class = "flexybayes_refusal_aggregation_route_unavailable"
  )

  d_gam <- d
  d_gam$y <- d_gam$y_gam
  expect_error(
    flexybayes_stream(
      fixed = y ~ 1,
      random = ~g,
      source = d_gam,
      family = "gamma",
      fit = FALSE,
      verbose = FALSE
    ),
    class = "flexybayes_refusal_aggregation_route_unavailable"
  )
})

test_that("the default link still aggregates on both entries", {
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE
  )
  d <- .agg_link_trial()
  d$y <- d$y_bin

  fit <- suppressMessages(flexybayes(
    y ~ 1 + (1 | g),
    data = d,
    family = "binomial",
    backend = "inla",
    aggregate = TRUE,
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")

  agg <- flexybayes_stream(
    fixed = y ~ 1,
    random = ~g,
    source = d,
    family = "binomial",
    fit = FALSE,
    verbose = FALSE
  )
  expect_s3_class(agg, "fb_aggregated")
})

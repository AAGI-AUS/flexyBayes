# =============================================================================
# The three aggregated-route entry refusals raise typed conditions.
#
# The execution grid's capability-matrix section (tools/execution_grid.R,
# section M) reached three refusals in .maybe_aggregate_gaussian() that raised
# a bare stop(): the missing-response guard, the backend-capability guard, and
# the no-resolvable-route guard. Each is a documented refusal -- the function's
# own contract says it "raises a typed condition" -- but a caller could not
# match any of them by class, only by message text.
#
# The two refusals lower down the same function already carried
# `flexybayes_aggregate_refusal`, so the family class is retained beneath the
# per-code class here: a handler written against the family keeps catching
# every aggregated-route refusal rather than only some of them.
#
# The three sites split across two boundaries, not one. A missing response is
# a property of the DATA (no sufficient statistic exists for the cell); an
# unavailable route is a property of the ENGINE ROSTER (no aggregated emit is
# reachable). They carry separate reason codes accordingly.
# =============================================================================

.agg_route_trial <- function(seed = 20260819L) {
  set.seed(seed)
  d <- expand.grid(
    rep = seq_len(20L),
    trt = factor(seq_len(5L)),
    blk = factor(seq_len(4L))
  )
  d$y <- 10 + as.numeric(d$trt) * 0.5 + stats::rnorm(nrow(d), 0, 0.3)
  d$g2 <- factor(rep(seq_len(2L), length.out = nrow(d)))
  d
}

# ---------------------------------------------------------------- #
# 1. Site 1 -- a missing response has no sufficient statistic       #
# ---------------------------------------------------------------- #

test_that("aggregate = TRUE on a missing response raises a typed refusal", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()
  d$y[c(3L, 17L)] <- NA_real_

  expect_error(
    suppressMessages(flexybayes(
      fixed = y ~ trt, random = ~ blk, residual = ~ units, data = d,
      backend = "inla", na_action = "augment", aggregate = TRUE,
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
      fixed = y ~ trt, random = ~ blk, residual = ~ units, data = d,
      backend = "inla", na_action = "augment", aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  # A caller pattern-matching the aggregated-route family still catches it.
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_s3_class(err, "flexybayes_refusal")
  # The message text is the pre-0.9.2 text verbatim -- users match on it.
  expect_match(conditionMessage(err), "the response has missing values",
               fixed = TRUE)
  expect_match(conditionMessage(err), "per-cell sufficient", fixed = TRUE)
  expect_identical(err$reason_code, "aggregation_response_incomplete")
})

# ---------------------------------------------------------------- #
# 2. Site 2 -- a backend with no aggregated emit at all             #
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
      fixed = y ~ trt, random = ~ blk, data = d,
      backend = "brms", aggregate = TRUE, verbose = FALSE
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
      fixed = y ~ trt, random = ~ blk, data = d,
      backend = "brms", aggregate = TRUE, verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_identical(err$backend, "brms")
  expect_match(conditionMessage(err),
               "is not supported on backend = \"brms\"", fixed = TRUE)
  expect_match(conditionMessage(err), "switch backend", fixed = TRUE)
})

# ---------------------------------------------------------------- #
# 3. Site 3 -- no aggregated backend resolves for this model        #
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
      fixed = y ~ trt, random = ~ us(g2):blk, data = d,
      backend = "auto", aggregate = TRUE, verbose = FALSE
    )),
    class = "flexybayes_refusal_aggregation_route_unavailable"
  )
})

test_that("the no-route refusal keeps the family class and its wording", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()

  err <- tryCatch(
    suppressMessages(flexybayes(
      fixed = y ~ trt, random = ~ us(g2):blk, data = d,
      backend = "auto", aggregate = TRUE, verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_match(conditionMessage(err), "no active aggregated backend",
               fixed = TRUE)
  expect_match(conditionMessage(err), "greta aggregated path is quarantined",
               fixed = TRUE)
  expect_identical(err$reason_code, "aggregation_route_unavailable")
})

# ---------------------------------------------------------------- #
# 4. The message text did not move when the class arrived           #
# ---------------------------------------------------------------- #
#
# Typing these refusals is an additive change: the condition class is new,
# the message body is the pre-0.9.2 text verbatim. Callers already in the
# field match on the text, so the whole string is pinned here rather than a
# fragment of it. A deliberate re-wording must edit this test.

test_that("the two route refusals carry their pre-0.9.2 message verbatim", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .agg_route_trial()

  err_brms <- tryCatch(
    suppressMessages(flexybayes(
      fixed = y ~ trt, random = ~ blk, data = d,
      backend = "brms", aggregate = TRUE, verbose = FALSE
    )),
    error = function(e) e
  )
  expect_identical(
    conditionMessage(err_brms),
    paste0(
      "`aggregate = TRUE` is not supported on backend = \"brms\" (the ",
      "aggregated path is wired for inla only since the greta ",
      "quarantine). Pass aggregate = FALSE or switch backend."
    )
  )

  err_route <- tryCatch(
    suppressMessages(flexybayes(
      fixed = y ~ trt, random = ~ us(g2):blk, data = d,
      backend = "auto", aggregate = TRUE, verbose = FALSE
    )),
    error = function(e) e
  )
  expect_identical(
    conditionMessage(err_route),
    paste0(
      "`aggregate = TRUE` refused: no active aggregated backend (INLA ",
      "refused or unavailable; the greta aggregated path is quarantined). ",
      "Pass aggregate = FALSE for the per-row path."
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
      fixed = y ~ trt, random = ~ blk, residual = ~ units, data = d,
      backend = "inla", na_action = "augment", aggregate = TRUE,
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
  for (code in c("aggregation_response_incomplete",
                 "aggregation_route_unavailable")) {
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
  tbl <- fb_refusals(reason_code = c("aggregation_response_incomplete",
                                     "aggregation_route_unavailable"))
  expect_identical(nrow(tbl), 2L)
  expect_match(tbl$description[tbl$reason_code ==
                                 "aggregation_response_incomplete"],
               "sufficient", fixed = TRUE)
  expect_match(tbl$description[tbl$reason_code ==
                                 "aggregation_route_unavailable"],
               "no aggregated emit route", fixed = TRUE)
})

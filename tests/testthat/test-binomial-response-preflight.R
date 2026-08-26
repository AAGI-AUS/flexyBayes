# Tests for the non-binary binomial-response preflight refusal
# (D14 / WP-D, found by WP-G's test sweep).
#
# Before this fix, family = "binomial" with a response outside {0, 1}
# reached the engine unrefused on every route and failed there raw:
# INLA with a subprocess death wrapped into an unrelated
# inla_program_failed diagnosis, brms with an untyped condition, and
# the aggregated route with its own untyped stop(). One typed
# preflight refusal (binomial_response_not_binary,
# flexybayes_binomial_response_not_binary_refusal) now fires before
# any of the three, at the single dispatch choke point every route
# passes through.

test_that(".refuse_non_binary_binomial_response() is a no-op for a valid 0/1 response", {
  fb <- list(family = "binomial", response = "y")
  d <- data.frame(y = c(0, 1, 1, 0, NA))
  expect_null(flexyBayes:::.refuse_non_binary_binomial_response(fb, d))
})

test_that(".refuse_non_binary_binomial_response() is a no-op for a non-binomial family", {
  fb <- list(family = "gaussian", response = "y")
  d <- data.frame(y = c(2, 3, 4, 5))
  expect_null(flexyBayes:::.refuse_non_binary_binomial_response(fb, d))
})

test_that(".refuse_non_binary_binomial_response() refuses a count response, naming the column and values", {
  fb <- list(family = "binomial", response = "y")
  d <- data.frame(y = c(0L, 1L, 2L, 3L, 5L))
  err <- tryCatch(
    flexyBayes:::.refuse_non_binary_binomial_response(fb, d),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_binomial_response_not_binary")
  expect_s3_class(err, "flexybayes_binomial_response_not_binary_refusal")
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "`y`", fixed = TRUE)
  expect_match(conditionMessage(err), "2, 3, 5", fixed = TRUE)
  expect_match(conditionMessage(err), "trials = ", fixed = TRUE)
  expect_identical(err$column, "y")
  expect_identical(err$offending_values, c(2L, 3L, 5L))
})

test_that(".refuse_non_binary_binomial_response() refuses a non-numeric response by class", {
  fb <- list(family = "binomial", response = "y")
  d <- data.frame(y = factor(c("yes", "no", "yes")))
  err <- tryCatch(
    flexyBayes:::.refuse_non_binary_binomial_response(fb, d),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_binomial_response_not_binary")
  expect_match(conditionMessage(err), "class `factor`", fixed = TRUE)
})

# ---------------------------------------------------------------- #
# Live: both engines, and the aggregated route                     #
# ---------------------------------------------------------------- #

.mk_count_binomial_data <- function(n = 60L, seed = 20260826L) {
  set.seed(seed)
  data.frame(
    y = sample(0:5, n, replace = TRUE),
    x = stats::rnorm(n),
    g = factor(rep(seq_len(6L), length.out = n))
  )
}

test_that("a live INLA call refuses before the engine runs", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .mk_count_binomial_data()
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ x + (1 | g),
      data = d,
      family = "binomial",
      backend = "inla",
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_binomial_response_not_binary")
})

test_that("a live brms call refuses before the engine runs", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- .mk_count_binomial_data()
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ x + (1 | g),
      data = d,
      family = "binomial",
      backend = "brms",
      n_samples = 50,
      warmup = 50,
      chains = 1,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_binomial_response_not_binary")
})

test_that("the aggregated route (aggregate = TRUE) refuses with the typed class, not the old untyped stop()", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(20260826L)
  d <- data.frame(
    y = sample(0:5, 2000L, replace = TRUE),
    env = factor(rep(1:5, 400L)),
    geno = factor(rep(1:20, 100L))
  )
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ env,
      random = ~geno,
      data = d,
      family = "binomial",
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_binomial_response_not_binary")
  # Not the pre-fix ad hoc condition, which carried no flexybayes_refusal
  # class at all.
  expect_false(inherits(err, "flexybayes_aggregate_refusal"))
})

test_that("flexybayes_stream()'s own trials = route is unaffected", {
  # The new guard runs only in .dispatch_backend(), the non-streaming
  # entry's choke point; flexybayes_stream() calls the aggregated emit
  # functions directly and has its own trials = contract, which this
  # item must not break.
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(20260826L)
  n_trials <- sample(3:10, 400L, replace = TRUE)
  d <- data.frame(
    y = stats::rbinom(400L, n_trials, 0.4),
    ntrial = n_trials,
    env = factor(rep(1:5, 80L)),
    geno = factor(rep(1:20, 20L))
  )
  fit <- tryCatch(
    suppressMessages(flexybayes_stream(
      y ~ env,
      random = ~geno,
      source = d,
      family = "binomial",
      trials = "ntrial",
      backend = "inla",
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_false(inherits(fit, "error"),
    label = if (inherits(fit, "error")) conditionMessage(fit) else "stream fit"
  )
})

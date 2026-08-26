# test-stream-overflow-and-hooks.R -- two small structural guards.
#
# (a) flexybayes_stream() offers only inla as a backend choice, and a
#     caller who names a withdrawn or otherwise unrecognised backend gets
#     the unknown-backend refusal by name rather than match.arg's
#     "'arg' should be one of".
# (b) The single-file and in-memory streaming sources refuse a row count
#     past R's integer limit instead of recording it as NA. The row count
#     is mocked; nothing here allocates billions of rows.
#
# (0.9.3: a third section here tested the S3 shim that registered
# flexyBayes methods on a later load of the withdrawn native engine's
# package -- that shim registrar and its two .onLoad() hook
# registrations are deleted entirely (see NEWS.md), so that section is
# removed rather than rewritten; there is no successor mechanism, since
# no active engine needs a deferred-load S3 shim.)

# ---------------------------------------------------------------- #
# (a) The streaming entry point takes only inla                     #
# ---------------------------------------------------------------- #

test_that("flexybayes_stream() defaults to inla and refuses an unrecognised backend by name", {
  expect_identical(
    eval(formals(flexybayes_stream)$backend),
    "inla"
  )

  err <- tryCatch(
    flexybayes_stream(
      y ~ env,
      source = data.frame(y = 1, env = factor("a")),
      backend = "stan"
    ),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_unknown_backend")
  expect_s3_class(err, "flexybayes_unknown_backend_refusal")
  expect_identical(err$backend, "stan")
  expect_match(conditionMessage(err), "not a recognised flexyBayes engine")
  expect_match(conditionMessage(err), '"inla"', fixed = TRUE)
  # The failure must not be argument matching -- that is the untyped
  # outcome this change exists to avoid.
  expect_false(grepl("should be one of", conditionMessage(err)))
})

# ---------------------------------------------------------------- #
# (b) Row counts past the integer limit refuse                      #
# ---------------------------------------------------------------- #

test_that("a row count past 2^31 - 1 refuses instead of recording NA", {
  big <- 5e9
  err <- tryCatch(
    flexyBayes:::.fb_checked_row_count(big, "This .fst source"),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_row_count_exceeds_integer")
  expect_identical(err$n_rows, big)
  expect_match(conditionMessage(err), "5000000000")
  expect_match(conditionMessage(err), "2147483647")
  expect_match(conditionMessage(err), "Partition the source")

  # What the old code did with the same input, pinned so the regression
  # is visible: a silent NA.
  expect_true(is.na(suppressWarnings(as.integer(big))))
})

test_that("row counts inside the integer range pass through unchanged", {
  expect_identical(
    flexyBayes:::.fb_checked_row_count(1e6, "This .fst source"),
    1000000L
  )
  expect_identical(
    flexyBayes:::.fb_checked_row_count(.Machine$integer.max, "src"),
    .Machine$integer.max
  )
  # A generator source has no knowable row count and must stay NA rather
  # than refuse.
  expect_identical(
    flexyBayes:::.fb_checked_row_count(NA_integer_, "src"),
    NA_integer_
  )
})

test_that("an in-memory source with a mocked row count refuses", {
  # The row count is mocked at the source constructor rather than
  # materialised. The point of the guard is that nothing downstream ever
  # sees an NA where a count belongs.
  df <- data.frame(y = stats::rnorm(10L), g = factor(rep(1:2, 5L)))
  testthat::local_mocked_bindings(
    .fb_checked_row_count = function(n_rows, context) {
      stop(flexyBayes:::.fb_refusal_condition(
        reason_code = "row_count_exceeds_integer",
        message = paste0(context, " reports too many rows."),
        n_rows = n_rows
      ))
    }
  )
  err <- tryCatch(
    flexyBayes:::.fb_stream_source(df, chunk_rows = 5L),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_row_count_exceeds_integer")
})

test_that("the aggregation plan routes its row count through the guard", {
  src <- deparse(body(flexyBayes:::.fb_aggregation_plan))
  expect_true(any(grepl(".fb_checked_row_count", src, fixed = TRUE)))
  expect_false(any(grepl("as.integer(fb_dataset$n_rows)", src, fixed = TRUE)))
})


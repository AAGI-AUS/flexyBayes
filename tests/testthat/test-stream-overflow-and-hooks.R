# test-stream-overflow-and-hooks.R -- three small structural guards.
#
# (a) flexybayes_stream() no longer offers greta as a backend choice, and
#     a caller who passes it gets the quarantine refusal by name rather
#     than match.arg's "'arg' should be one of".
# (b) The single-file and in-memory streaming sources refuse a row count
#     past R's integer limit instead of recording it as NA. The row count
#     is mocked; nothing here allocates billions of rows.
# (c) The greta S3 shim registers on a later greta load as well as on an
#     earlier one, so the registration is not order-dependent. greta is
#     never loaded by this file.

# ---------------------------------------------------------------- #
# (a) The streaming entry point drops greta                         #
# ---------------------------------------------------------------- #

test_that("flexybayes_stream() defaults to inla and refuses greta by name", {
  expect_identical(
    eval(formals(flexybayes_stream)$backend),
    "inla"
  )

  err <- tryCatch(
    flexybayes_stream(
      y ~ env,
      source = data.frame(y = 1, env = factor("a")),
      backend = "greta"
    ),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_backend_quarantined")
  expect_s3_class(err, "flexybayes_backend_quarantined_refusal")
  expect_identical(err$backend, "greta")
  expect_match(conditionMessage(err), "quarantined")
  expect_match(conditionMessage(err), "aggregated emit")
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

# ---------------------------------------------------------------- #
# (c) The greta shim registers whichever order the loads happen in  #
# ---------------------------------------------------------------- #

test_that("a greta onLoad hook is registered so a later load gets the shim", {
  # isNamespaceLoaded() answers only for the instant .onLoad() ran, so a
  # user loading greta after flexyBayes would never have got the shim.
  # The hook covers that direction. greta is not loaded here: the test
  # asserts the registration, not the shim's effect.
  hooks <- getHook(packageEvent("greta", "onLoad"))
  expect_true(length(hooks) >= 1L)
  expect_true(any(vapply(hooks, is.function, logical(1))))

  # The hook body calls the shared registrar, so both routes register the
  # same methods.
  bodies <- vapply(
    hooks,
    function(h) paste(deparse(body(h)), collapse = " "),
    character(1)
  )
  expect_true(any(grepl(".register_greta_shim", bodies, fixed = TRUE)))
})

test_that("the shim registrar refuses to load greta itself", {
  # asNamespace() LOADS a namespace that is not already loaded, so
  # without its own guard the registrar would pull in greta, reticulate
  # and TensorFlow -- the exact cost the isNamespaceLoaded() policy
  # exists to avoid. The guard is asserted at the source, because a
  # runtime assertion depends on whether some earlier test file in the
  # suite has loaded greta already.
  src <- paste(
    deparse(body(flexyBayes:::.register_greta_shim)),
    collapse = " "
  )
  expect_true(grepl("isNamespaceLoaded(\"greta\")", src, fixed = TRUE))
  expect_false(grepl("requireNamespace(\"greta\")", src, fixed = TRUE))
  # The guard precedes the asNamespace() call it protects.
  expect_lt(
    regexpr("isNamespaceLoaded", src, fixed = TRUE),
    regexpr("asNamespace(\"greta\")", src, fixed = TRUE)
  )
  expect_silent(flexyBayes:::.register_greta_shim("flexyBayes"))
})

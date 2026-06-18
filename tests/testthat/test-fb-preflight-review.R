# Tests for the Stage 2 MVP preflight integration into the review-code
# path of flexybayes() / fb_brms() (ADR 0021 / v0.3.0).
#
# Two flows are exercised:
#   - small-n (< 1e5): review token's $preflight slot stays NULL;
#     print() output unchanged from v0.2
#   - large-n (>= 1e5): preflight runs, attaches to the review object,
#     print() shows the "Preflight summary" header above the cat_code()
#     prompt, and a forced-tight ceiling raises the typed condition
#     before any backend code runs.

test_that("review small-n: preflight slot is NULL; print unchanged", {
  df <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    g = factor(rep(letters[1:5], 10))
  )
  rev <- flexybayes(y ~ x, random = ~g, data = df, review_code = TRUE)

  expect_s3_class(rev, "flexybayes_review")
  expect_null(rev$preflight)

  out <- capture.output(print(rev))
  expect_false(any(grepl("Preflight summary", out, fixed = TRUE)))
  expect_true(any(grepl("<flexybayes_review>", out, fixed = TRUE)))
})

test_that("review large-n: preflight attaches; print shows the summary", {
  N <- 2e5L
  df <- data.frame(
    y = rnorm(N),
    g = factor(sample.int(50L, N, replace = TRUE))
  )
  rev <- flexybayes(y ~ 1, random = ~g, data = df, review_code = TRUE)

  expect_s3_class(rev, "flexybayes_review")
  expect_s3_class(rev$preflight, "fb_preflight")
  expect_null(rev$preflight$refusal)

  out <- capture.output(print(rev))
  expect_true(any(grepl("Preflight summary", out, fixed = TRUE)))
  expect_true(any(grepl("(1 | g)", out, fixed = TRUE)))
  expect_true(any(grepl("ceiling", out, fixed = TRUE)))
})

test_that("review large-n: tight ceiling raises preflight refusal before code emit", {
  N <- 2e5L
  df <- data.frame(
    y = rnorm(N),
    g = factor(sample.int(50L, N, replace = TRUE))
  )
  withr::local_options(flexyBayes.preflight_ceiling_gb = 0.0001)

  err <- tryCatch(
    flexybayes(y ~ 1, random = ~g, data = df, review_code = TRUE),
    flexybayes_preflight_refusal = function(c) c,
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_preflight_refusal")
  expect_identical(err$reason_code, "design_memory_exceeds_ceiling")
})

test_that("review small-n via fb_brms: preflight slot stays NULL", {
  df <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    g = factor(rep(letters[1:5], 10))
  )
  rev <- fb(y ~ x + (1 | g), data = df, review_code = TRUE)

  expect_s3_class(rev, "flexybayes_review")
  expect_null(rev$preflight)
})

# test-aggregated-labels.R -- the aggregated header names the fit's own
# family.
#
# `print()` and `summary()` on a `<flexybayes_aggregated>` object
# hard-coded "aggregated-gaussian", so a binomial or Poisson fit from the
# count aggregator printed the wrong family on the surface the package
# uses to signal exactness (2026-08-16 adversarial review, P2-1).

test_that("the label is read off the fit's family", {
  expect_identical(
    flexyBayes:::.agg_fit_label(list(family = "gaussian")),
    "aggregated-gaussian"
  )
  expect_identical(
    flexyBayes:::.agg_fit_label(list(family = "binomial")),
    "aggregated-binomial"
  )
  expect_identical(
    flexyBayes:::.agg_fit_label(list(family = "poisson")),
    "aggregated-poisson"
  )
})

test_that("a fit with no recorded family degrades to the bare label", {
  expect_identical(flexyBayes:::.agg_fit_label(list()), "aggregated")
  expect_identical(
    flexyBayes:::.agg_fit_label(list(family = NA_character_)),
    "aggregated"
  )
  expect_identical(
    flexyBayes:::.agg_fit_label(list(family = c("a", "b"))),
    "aggregated"
  )
})

test_that("print and summary on a live aggregated fit carry the family", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(20260816L)
  n <- 1000L
  d <- data.frame(
    f = factor(sample(paste0("f", seq_len(5L)), n, replace = TRUE)),
    g = factor(sample(paste0("g", seq_len(20L)), n, replace = TRUE))
  )
  d$y <- stats::rnorm(n)
  fit <- suppressMessages(flexybayes(
    y ~ f,
    random = ~g,
    data = d,
    backend = "inla",
    aggregate = TRUE,
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")

  print_out <- utils::capture.output(print(fit))
  expect_true(any(grepl(
    "[flexyBayes / aggregated-gaussian]", print_out, fixed = TRUE
  )))

  summary_out <- utils::capture.output(summary(fit))
  expect_true(any(grepl(
    "[flexyBayes / aggregated-gaussian]", summary_out, fixed = TRUE
  )))
})

test_that("a Poisson aggregated fit prints aggregated-poisson", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(20260816L)
  n <- 1000L
  d <- data.frame(
    f = factor(sample(paste0("f", seq_len(5L)), n, replace = TRUE)),
    g = factor(sample(paste0("g", seq_len(20L)), n, replace = TRUE))
  )
  d$y <- stats::rpois(n, lambda = 3)
  fit <- tryCatch(
    suppressMessages(flexybayes(
      y ~ f,
      random = ~g,
      data = d,
      family = "poisson",
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  # The engine is installed, so a fit error here is a regression, not a
  # reason to skip.
  expect_false(inherits(fit, "error"),
               label = if (inherits(fit, "error")) {
                 conditionMessage(fit)
               } else {
                 "aggregated Poisson fit"
               })
  skip_if(inherits(fit, "error"))
  print_out <- utils::capture.output(print(fit))
  expect_true(any(grepl(
    "[flexyBayes / aggregated-poisson]", print_out, fixed = TRUE
  )))
})

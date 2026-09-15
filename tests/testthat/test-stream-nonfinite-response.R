# flexybayes_stream() reads its chunks directly and does not pass through the
# missing-response layer that refuses a NaN or infinite response on every
# flexybayes() route. On a Poisson stream a NaN count fitted without a message
# and an infinite one ended in an untyped INLA program failure. The chunk
# reader now refuses both by name before any statistic is accumulated.

.stream_nonfinite_data <- function(value) {
  set.seed(3L)
  d <- data.frame(
    env = factor(sample(1:4, 400L, replace = TRUE)),
    geno = factor(sample(1:10, 400L, replace = TRUE))
  )
  d$y <- stats::rpois(400L, 3)
  d$y[7L] <- value
  d
}

test_that("a streamed count model refuses a NaN or infinite response by name", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  for (value in list(Inf, NaN)) {
    for (fit in c(TRUE, FALSE)) {
      err <- tryCatch(
        suppressMessages(flexybayes_stream(
          y ~ env,
          random = ~geno,
          source = .stream_nonfinite_data(value),
          family = "poisson",
          fit = fit,
          verbose = FALSE
        )),
        error = function(e) e
      )
      expect_s3_class(err, "flexybayes_refusal_response_not_finite")
      expect_match(conditionMessage(err), as.character(value), fixed = TRUE)
    }
  }
})

test_that("a streamed count model refuses a missing response as aggregate = TRUE does", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  # A missing count used to fit without a message: the cell holding it
  # summed to NA and dropped out of the fit.
  for (fit in c(TRUE, FALSE)) {
    err <- tryCatch(
      suppressMessages(flexybayes_stream(
        y ~ env,
        random = ~geno,
        source = .stream_nonfinite_data(NA_real_),
        family = "poisson",
        fit = fit,
        verbose = FALSE
      )),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_refusal_aggregation_response_incomplete")
    expect_s3_class(err, "flexybayes_aggregate_refusal")
  }
})

test_that("a finite streamed count model is unaffected", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  agg <- flexybayes_stream(
    y ~ env,
    random = ~geno,
    source = .stream_nonfinite_data(2),
    family = "poisson",
    fit = FALSE,
    verbose = FALSE
  )
  expect_identical(as.numeric(agg$N), 400)
})

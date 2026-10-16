# The aggregated header names the fit's own family, so a binomial or
# Poisson fit from the count aggregator does not print as
# "aggregated-gaussian".

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
  # A Gaussian model is fitted per row (its aggregated route cannot carry
  # the residual prior), so this label is exercised on a Bernoulli model.
  d$y <- stats::rbinom(n, 1L, 0.4)
  fit <- suppressMessages(flexybayes(
    y ~ f,
    random = ~g,
    data = d,
    family = "binomial",
    backend = "inla",
    aggregate = TRUE,
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")

  print_out <- utils::capture.output(print(fit))
  expect_true(any(grepl(
    "[flexyBayes / aggregated-binomial]",
    print_out,
    fixed = TRUE
  )))

  summary_out <- utils::capture.output(summary(fit))
  expect_true(any(grepl(
    "[flexyBayes / aggregated-binomial]",
    summary_out,
    fixed = TRUE
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
  expect_false(
    inherits(fit, "error"),
    label = if (inherits(fit, "error")) {
      conditionMessage(fit)
    } else {
      "aggregated Poisson fit"
    }
  )
  skip_if(inherits(fit, "error"))
  print_out <- utils::capture.output(print(fit))
  expect_true(any(grepl(
    "[flexyBayes / aggregated-poisson]",
    print_out,
    fixed = TRUE
  )))
})

test_that("print on a flexybayes_stream() fit names the backend and the aggregated path", {
  # flexybayes_stream() bypasses .dispatch_backend(), so the emit itself
  # must set $exactness and $extras$backend_decision. Without them the
  # banner reads "(path = )" and "(unknown engine)".
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(20260826L)
  n <- 1000L
  d <- data.frame(
    f = factor(sample(paste0("f", seq_len(5L)), n, replace = TRUE)),
    g = factor(sample(paste0("g", seq_len(20L)), n, replace = TRUE))
  )
  # A streamed Gaussian fit refuses, so the banner is checked on a Poisson
  # stream.
  d$y <- stats::rpois(n, 3)
  fit <- suppressMessages(flexybayes_stream(
    y ~ f,
    random = ~g,
    source = d,
    family = "poisson",
    backend = "inla",
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")
  expect_identical(fit$exactness, "aggregated_exact")
  expect_identical(fit$extras$backend_decision$backend, "inla")
  expect_identical(fit$extras$backend_decision$path, "aggregated_count")

  print_out <- utils::capture.output(print(fit))
  expect_true(any(grepl(
    "backend:    inla (path = aggregated_count )",
    print_out,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Engine:         INLA Laplace",
    print_out,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Representation: aggregated_exact",
    print_out,
    fixed = TRUE
  )))
  expect_false(any(grepl("(unknown engine)", print_out, fixed = TRUE)))
})

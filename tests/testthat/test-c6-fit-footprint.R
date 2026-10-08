# What an INLA fit costs to keep, save or send to another process, and the
# factorisations vcov() and the posterior sampler run per configuration.

.c6_fit <- function(n, seed = 1L) {
  set.seed(seed)
  d <- data.frame(
    site = factor(sample(sprintf("s%02d", 1:40), n, replace = TRUE)),
    x = stats::rnorm(n)
  )
  d$y <- stats::rnorm(40)[d$site] + 0.5 * d$x + stats::rnorm(n)
  suppressWarnings(suppressMessages(flexybayes(
    y ~ x,
    random = ~site,
    data = d,
    backend = "inla",
    verbose = FALSE
  )))
}

.c6_serialized_mb <- function(x) length(serialize(x, NULL)) / 2^20

.c6_object_mb <- function(x) as.numeric(utils::object.size(x)) / 2^20

test_that("DEF-204: a saved INLA fit is no larger than the fit in memory", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )

  fit <- .c6_fit(2000L)

  # object.size() counts no environment: the margin is what environments add.
  expect_lt(.c6_serialized_mb(fit), 1.25 * .c6_object_mb(fit) + 5)
  expect_lt(.c6_serialized_mb(fit$extras$formula), 5)
  expect_lt(.c6_serialized_mb(fit$inla$all.hyper), 1)
})

test_that("vcov() factorises one configuration's precision and reuses its ordering", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )

  fit <- .c6_fit(400L)
  cs <- fit$inla$misc$configs
  expect_gt(length(cs$config), 1L)

  # Reference: a fresh factorisation of every configuration's precision.
  fixed_names <- rownames(fit$inla$summary.fixed)
  offset <- sum(cs$contents$length) - cs$n
  idx <- cs$contents$start[match(fixed_names, cs$contents$tag)] - offset
  unit <- Matrix::sparseMatrix(
    i = idx,
    j = seq_along(idx),
    x = 1,
    dims = c(cs$n, length(idx))
  )
  w <- .fb_inla_config_weights(cs)
  second <- matrix(0, length(idx), length(idx))
  first <- numeric(length(idx))
  for (k in seq_along(cs$config)) {
    chol_q <- Matrix::Cholesky(
      Matrix::forceSymmetric(cs$config[[k]]$Q, uplo = "U"),
      LDL = FALSE
    )
    cov_f <- as.matrix(Matrix::solve(chol_q, unit))[idx, , drop = FALSE]
    mu <- cs$config[[k]]$improved.mean[idx]
    second <- second + w[k] * (cov_f + tcrossprod(mu))
    first <- first + w[k] * mu
  }
  reference <- second - tcrossprod(first)
  reference <- (reference + t(reference)) / 2
  reference <- reference *
    tcrossprod(.fb_inla_draw_map(fit, fixed_names)$b)
  v <- vcov(fit)
  expect_identical(unname(v), unname(reference))

  n_chol <- 0L
  cholesky <- Matrix::Cholesky
  testthat::local_mocked_bindings(
    Cholesky = function(...) {
      n_chol <<- n_chol + 1L
      cholesky(...)
    },
    .package = "Matrix"
  )
  invisible(vcov(fit))
  expect_identical(n_chol, 1L)

  n_chol <- 0L
  set.seed(1L)
  invisible(.fb_inla_config_sample(fit, 200L))
  expect_identical(n_chol, 1L)
})

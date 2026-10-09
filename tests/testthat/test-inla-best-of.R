# Best of k INLA fits: the fit in this session plus k - 1 fits in fresh R
# processes, keeping the highest marginal likelihood (R/inla_best_of.R).

# A stand-in for an INLA fit: a marginal likelihood, a mode status and a
# narrow marginal for each precision.
.bo_fake <- function(mlik, prec = c("the Gaussian observations" = 4,
                                    site = 1),
                     status = 0L) {
  marg <- function(p) {
    x <- seq(0.99 * p, 1.01 * p, length.out = 201L)
    cbind(x = x, y = stats::dnorm(x, p, 0.001 * p))
  }
  nm <- paste("Precision for", names(prec))
  list(
    mlik = matrix(c(mlik, mlik), 2L, 1L),
    mode = list(mode.status = status),
    summary.hyperpar = data.frame(mean = unname(prec), row.names = nm),
    marginals.hyperpar = stats::setNames(lapply(unname(prec), marg), nm)
  )
}

# Replaces the fresh-process fits with `fits`, in order, and counts calls.
.bo_mock_isolated <- function(fits, env = parent.frame()) {
  calls <- new.env()
  calls$n <- 0L
  calls$args <- list()
  testthat::local_mocked_bindings(
    .fb_inla_isolated_fit = function(args) {
      calls$n <- calls$n + 1L
      calls$args[[calls$n]] <- args
      fits[[calls$n]]
    },
    .fb_isolation_callr_available = function() TRUE,
    .env = env
  )
  calls
}

.bo_data <- function() {
  set.seed(3)
  d <- data.frame(
    site = factor(rep(1:6, each = 12)),
    trt = factor(rep(c("a", "b", "c"), 24)),
    x = stats::rnorm(72)
  )
  d$y <- 1 + c(0, 0.5, -0.3)[d$trt] + 0.4 * d$x +
    stats::rnorm(6)[d$site] + stats::rnorm(72)
  d
}

.bo_fit <- function() {
  suppressWarnings(suppressMessages(flexybayes(
    y ~ trt + x, random = ~site, data = .bo_data(), backend = "inla",
    verbose = FALSE
  )))
}

test_that("the option is 3 by default and refuses anything but a count", {
  withr::local_options(flexyBayes.inla_best_of = NULL)
  expect_identical(.fb_inla_best_of_k(), 3L)
  withr::local_options(flexyBayes.inla_best_of = 1)
  expect_identical(.fb_inla_best_of_k(), 1L)
  withr::local_options(flexyBayes.inla_best_of = 5L)
  expect_identical(.fb_inla_best_of_k(), 5L)
  for (bad in list(0L, -1, 2.5, "3", NA, NA_integer_, c(1L, 2L), Inf, TRUE)) {
    withr::local_options(flexyBayes.inla_best_of = bad)
    expect_error(
      .fb_inla_best_of_k(),
      class = "flexybayes_refusal_inla_best_of_invalid"
    )
  }
})

test_that("a bad option value refuses before INLA is called", {
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_best_of = 0L
  )
  testthat::local_mocked_bindings(
    .inla_call = function(args) stop("INLA must not be called")
  )
  expect_error(
    flexybayes(y ~ trt, random = ~site, data = .bo_data(), backend = "inla",
               verbose = FALSE),
    class = "flexybayes_refusal_inla_best_of_invalid"
  )
})

test_that("the fit with the highest marginal likelihood is kept, the first of equals", {
  skip_if_not_installed("INLA")
  calls <- .bo_mock_isolated(list(.bo_fake(-9), .bo_fake(-9)))
  args <- list(formula = y ~ 1, marker = "the call")
  out <- .fb_inla_best_of(.bo_fake(-10), args, k = 3L)
  expect_identical(calls$n, 2L)
  expect_identical(calls$args, list(args, args))
  expect_identical(out$record$chosen, 2L)
  expect_identical(out$fit$mlik[1L, 1L], -9)
  expect_equal(out$record$mlik, c(-10, -9, -9))
  expect_equal(out$record$mlik_range, 1)

  calls <- .bo_mock_isolated(list(.bo_fake(-9.5), .bo_fake(-12)))
  out <- .fb_inla_best_of(.bo_fake(-9), args, k = 3L)
  expect_identical(out$record$chosen, 1L)
  expect_identical(out$fit$mlik[1L, 1L], -9)
})

test_that("the record carries each fit's status and each SD's spread", {
  skip_if_not_installed("INLA")
  .bo_mock_isolated(list(
    .bo_fake(-8, c("the Gaussian observations" = 4, site = 1 / 1.21),
             status = 1L),
    .bo_fake(-9, c("the Gaussian observations" = 4, site = 1))
  ))
  rec <- .fb_inla_best_of(.bo_fake(-10), list(), k = 3L)$record
  expect_identical(rec$k, 3L)
  expect_identical(rec$mode_status, c(0L, 1L, 0L))
  expect_identical(dim(rec$sd), c(3L, 2L))
  expect_identical(
    colnames(rec$sd),
    c("the Gaussian observations", "site")
  )
  # SD = 1 / sqrt(precision): 1 against 1.1 for site, 0.5 throughout for
  # the residual.
  expect_equal(unname(rec$sd[, "site"]), c(1, 1.1, 1), tolerance = 1e-4)
  expect_equal(
    unname(rec$sd_rel_diff),
    c(0, 0.1),
    tolerance = 1e-4
  )
  line <- .fb_best_of_line(rec)
  expect_match(line, "fit 2 of 3 kept (highest marginal likelihood)", fixed = TRUE)
  expect_match(
    line,
    "the SD for site differed by 10.0% across the 3 fits",
    fixed = TRUE
  )
})

test_that("no line is printed when every SD agrees within 5%", {
  rec <- list(k = 3L, chosen = 1L, sd_rel_diff = c(a = 0.01, b = 0.05))
  expect_null(.fb_best_of_line(rec))
  expect_null(.fb_best_of_line(NULL))
})

test_that("a failed fit in a separate process stops the fit with its message", {
  skip_if_not_installed("INLA")
  testthat::local_mocked_bindings(
    .fb_inla_isolated_fit = function(args) {
      stop(.fb_refusal_condition(
        reason_code = "isolated_process_died",
        message = "the process was ended by signal 9 (SIGKILL)"
      ))
    },
    .fb_isolation_callr_available = function() TRUE
  )
  err <- tryCatch(
    .fb_inla_best_of(.bo_fake(-10), list(), k = 3L),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_isolated_process_died")
  expect_match(
    conditionMessage(err),
    "INLA fit 2 of 3, run in a separate R process, failed: ",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "signal 9 (SIGKILL)", fixed = TRUE)

  n <- 0L
  testthat::local_mocked_bindings(
    .fb_inla_isolated_fit = function(args) {
      n <<- n + 1L
      if (n == 2L) stop("the inla-program exited with an error")
      .bo_fake(-1)
    }
  )
  expect_error(
    .fb_inla_best_of(.bo_fake(-10), list(), k = 3L),
    "INLA fit 3 of 3, run in a separate R process, failed: the inla-program",
    fixed = TRUE
  )
})

test_that("k above 1 without callr refuses by name", {
  testthat::local_mocked_bindings(
    .fb_isolation_callr_available = function() FALSE
  )
  withr::local_options(flexyBayes.inla_best_of = 2L)
  expect_error(
    .fb_inla_best_of_k(),
    class = "flexybayes_refusal_isolation_requires_callr"
  )
  withr::local_options(flexyBayes.inla_best_of = 1L)
  expect_identical(.fb_inla_best_of_k(), 1L)
})

test_that("k = 1 is the single fit, with nothing added to the record", {
  testthat::local_mocked_bindings(
    .fb_inla_isolated_fit = function(args) stop("no further fit at k = 1")
  )
  fake <- .bo_fake(-10)
  out <- .fb_inla_best_of(fake, list(), k = 1L)
  expect_identical(out$fit, fake)
  expect_null(out$record)

  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1",
    flexyBayes.inla_best_of = 1L
  )
  fit <- .bo_fit()
  rec <- fit$extras$reproducibility
  expect_identical(
    names(rec),
    c(
      "r_version", "flexybayes_version", "engine", "engine_versions",
      "threads", "seed"
    )
  )
  expect_false(any(grepl("Best of", utils::capture.output(print(fit)))))
})

test_that("flexybayes() keeps the better fit from a separate process and says so", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1",
    flexyBayes.inla_best_of = 3L
  )
  calls <- new.env()
  calls$n <- 0L
  testthat::local_mocked_bindings(
    .fb_inla_isolated_fit = function(args) {
      calls$n <- calls$n + 1L
      fit <- .inla_call(args)
      if (calls$n == 1L) {
        # 1 log unit better, with the site SD 10% larger.
        fit$mlik[, 1L] <- fit$mlik[, 1L] + 1
        nm <- "Precision for site"
        m <- fit$marginals.hyperpar[[nm]]
        m[, 1L] <- m[, 1L] / 1.21
        m[, 2L] <- m[, 2L] * 1.21
        fit$marginals.hyperpar[[nm]] <- m
      } else {
        fit$mlik[, 1L] <- fit$mlik[, 1L] - 1
      }
      fit
    },
    .fb_isolation_callr_available = function() TRUE
  )
  fit <- .bo_fit()
  expect_identical(calls$n, 2L)
  rec <- fit$extras$reproducibility$best_of
  expect_identical(rec$chosen, 2L)
  expect_equal(rec$mlik[[2L]] - rec$mlik[[1L]], 1, tolerance = 1e-3)
  expect_equal(rec$mlik_range, 2, tolerance = 1e-3)
  # The marginal likelihood recorded is the one the returned fit reports.
  expect_equal(rec$mlik[[2L]], unname(fit$inla$mlik[1L, 1L]), tolerance = 1e-8)
  expect_equal(unname(rec$sd_rel_diff[["site"]]), 0.1, tolerance = 1e-3)
  printed <- utils::capture.output(print(fit))
  expect_true(any(grepl(
    "Best of  : fit 2 of 3 kept (highest marginal likelihood)",
    printed,
    fixed = TRUE
  )))
  summarised <- utils::capture.output(print(summary(fit)))
  expect_true(any(grepl("the SD for site differed by", summarised, fixed = TRUE)))
})

test_that("best of 3 runs two fits in fresh R processes on a package fixture", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("callr")
  skip_if_not_long_tier()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1",
    flexyBayes.inla_best_of = 3L
  )
  fit <- .bo_fit()
  rec <- fit$extras$reproducibility$best_of
  expect_identical(rec$k, 3L)
  expect_length(rec$mlik, 3L)
  expect_true(all(is.finite(rec$mlik)))
  expect_identical(rec$chosen, which.max(rec$mlik))
  expect_equal(
    rec$mlik[[rec$chosen]],
    unname(fit$inla$mlik[1L, 1L]),
    tolerance = 1e-8
  )
  expect_identical(rec$mode_status, rep(0L, 3L))
  expect_identical(nrow(rec$sd), 3L)
  expect_true(all(is.finite(rec$sd_rel_diff)))
})

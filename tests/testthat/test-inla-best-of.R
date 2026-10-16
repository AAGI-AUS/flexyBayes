# Best of k INLA fits: the fit in this session plus k - 1 fits in fresh R
# processes, keeping the highest marginal likelihood (R/inla_best_of.R).

# A stand-in for an INLA fit: a marginal likelihood, a mode status and a
# narrow marginal for each precision.
.bo_fake <- function(
  mlik,
  prec = c("the Gaussian observations" = 4, site = 1),
  status = 0L
) {
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

# Replaces the fresh-process fits with `fits`, in the order they are started,
# and counts calls. Fit j finishes `delays[[j]]` seconds after it starts. A
# condition in `fits` is raised when that fit is collected.
.bo_mock_isolated <- function(fits, delays = 0, env = parent.frame()) {
  calls <- new.env()
  calls$n <- 0L
  calls$args <- list()
  calls$running <- 0L
  calls$max_running <- 0L
  calls$finished <- integer(0)
  calls$stopped <- integer(0)
  delays <- rep_len(delays, length(fits))
  testthat::local_mocked_bindings(
    .fb_inla_isolated_start = function(args) {
      calls$n <- calls$n + 1L
      calls$args[[calls$n]] <- args
      calls$running <- calls$running + 1L
      calls$max_running <- max(calls$max_running, calls$running)
      list(i = calls$n, done_at = Sys.time() + delays[[calls$n]])
    },
    .fb_inla_isolated_running = function(handle) Sys.time() < handle$done_at,
    .fb_inla_isolated_collect = function(handle) {
      calls$running <- calls$running - 1L
      calls$finished <- c(calls$finished, handle$i)
      fit <- fits[[handle$i]]
      if (inherits(fit, "condition")) {
        stop(fit)
      }
      list(ok = TRUE, value = fit, seed = get0(".Random.seed", globalenv()))
    },
    .fb_inla_isolated_stop = function(handle) {
      calls$running <- calls$running - 1L
      calls$stopped <- c(calls$stopped, handle$i)
      invisible(NULL)
    },
    .fb_isolation_callr_available = function() TRUE,
    .env = env
  )
  calls
}

# Free memory and the memory of one fit, in GB.
.bo_mock_memory <- function(free_gb = 64, fit_gb = 1.5, env = parent.frame()) {
  testthat::local_mocked_bindings(
    .fb_free_memory_bytes = function() free_gb * 1024^3,
    .fb_inla_fit_bytes = function(fb, data) fit_gb * 1024^3,
    .env = env
  )
}

# The record without the wall times, which differ between runs.
.bo_strip <- function(rec) {
  rec$wall_s <- NULL
  rec
}

.bo_data <- function() {
  set.seed(3)
  d <- data.frame(
    site = factor(rep(1:6, each = 12)),
    trt = factor(rep(c("a", "b", "c"), 24)),
    x = stats::rnorm(72)
  )
  d$y <- 1 +
    c(0, 0.5, -0.3)[d$trt] +
    0.4 * d$x +
    stats::rnorm(6)[d$site] +
    stats::rnorm(72)
  d
}

.bo_fit <- function() {
  suppressWarnings(suppressMessages(flexybayes(
    y ~ trt + x,
    random = ~site,
    data = .bo_data(),
    backend = "inla",
    verbose = FALSE
  )))
}

test_that("the option is 3 by default and refuses anything but a count", {
  withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = ""))
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

test_that("the default is 1 while R CMD check limits cores", {
  withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = "TRUE"))
  withr::local_options(flexyBayes.inla_best_of = NULL)
  expect_identical(.fb_inla_best_of_k(), 1L)
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
    flexybayes(
      y ~ trt,
      random = ~site,
      data = .bo_data(),
      backend = "inla",
      verbose = FALSE
    ),
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
    .bo_fake(
      -8,
      c("the Gaussian observations" = 4, site = 1 / 1.21),
      status = 1L
    ),
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
  expect_match(
    line,
    "fit 2 of 3 kept (highest marginal likelihood)",
    fixed = TRUE
  )
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
  died <- .fb_refusal_condition(
    reason_code = "isolated_process_died",
    message = "the process was ended by signal 9 (SIGKILL)"
  )
  .bo_mock_isolated(list(died, died))
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

  .bo_mock_isolated(list(
    .bo_fake(-1),
    simpleError("the inla-program exited with an error")
  ))
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
    .fb_inla_isolated_start = function(args) stop("no further fit at k = 1"),
    .fb_free_memory_bytes = function() stop("no memory probe at k = 1"),
    .fb_inla_fit_bytes = function(fb, data) stop("no memory estimate at k = 1")
  )
  withr::local_options(flexyBayes.inla_best_of_workers = 0L)
  expect_null(.fb_inla_best_of_launch(list(), 1L, NULL, NULL))
  fake <- .bo_fake(-10)
  out <- .fb_inla_best_of(fake, list(), k = 1L)
  expect_identical(out$fit, fake)
  expect_null(out$record)
  out <- .fb_inla_best_of(fake, list(), k = 1L, pool = NULL)
  expect_identical(out, list(fit = fake, record = NULL))

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
      "r_version",
      "flexybayes_version",
      "engine",
      "engine_versions",
      "threads",
      "seed"
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
    .fb_inla_isolated_start = function(args) {
      calls$n <- calls$n + 1L
      list(i = calls$n, args = args)
    },
    .fb_inla_isolated_running = function(handle) FALSE,
    .fb_inla_isolated_collect = function(handle) {
      fit <- .inla_call(handle$args)
      if (handle$i == 1L) {
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
      list(ok = TRUE, value = fit, seed = get0(".Random.seed", globalenv()))
    },
    .fb_isolation_callr_available = function() TRUE
  )
  .bo_mock_memory()
  fit <- .bo_fit()
  expect_identical(calls$n, 2L)
  rec <- fit$extras$reproducibility$best_of
  expect_identical(rec$chosen, 2L)
  expect_true(rec$parallel)
  expect_length(rec$wall_s, 3L)
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
  expect_true(any(grepl(
    "the SD for site differed by",
    summarised,
    fixed = TRUE
  )))
})

test_that("the fits started before the session fit run beside it", {
  skip_if_not_installed("INLA")
  .bo_mock_memory()
  calls <- .bo_mock_isolated(list(.bo_fake(-9), .bo_fake(-12)))
  args <- list(formula = y ~ 1, marker = "the call")
  pool <- .fb_inla_best_of_launch(args, 3L, NULL, NULL)
  # Both further fits are running before the session fit starts.
  expect_identical(calls$n, 2L)
  expect_identical(calls$args, list(args, args))
  out <- .fb_inla_best_of(.bo_fake(-10), args, 3L, pool)
  rec <- out$record
  expect_identical(rec$chosen, 2L)
  expect_equal(rec$mlik, c(-10, -9, -12))
  expect_identical(rec$workers, 3L)
  expect_true(rec$parallel)
  expect_length(rec$wall_s, 3L)
  expect_true(all(is.finite(rec$wall_s) & rec$wall_s >= 0))
})

test_that("the fit kept and the record do not depend on which fit finishes first", {
  skip_if_not_installed("INLA")
  .bo_mock_memory()
  # Fits 3 and 4 tie on the highest marginal likelihood, and fit 3 is kept.
  fits <- list(.bo_fake(-9), .bo_fake(-8, status = 1L), .bo_fake(-8))
  run <- function(delays) {
    calls <- .bo_mock_isolated(fits, delays = delays)
    pool <- .fb_inla_best_of_launch(list(), 4L, NULL, NULL)
    out <- .fb_inla_best_of(.bo_fake(-10), list(), 4L, pool)
    list(out = out, finished = calls$finished)
  }
  a <- run(c(0.6, 0.3, 0.05))
  b <- run(c(0.05, 0.3, 0.6))
  # Children are numbered by start: child 1 is fit 2.
  expect_identical(a$finished, c(3L, 2L, 1L))
  expect_identical(b$finished, c(1L, 2L, 3L))
  for (x in list(a, b)) {
    expect_identical(x$out$record$chosen, 3L)
    expect_identical(x$out$fit, fits[[2L]])
    expect_equal(x$out$record$mlik, c(-10, -9, -8, -8))
    expect_identical(x$out$record$mode_status, c(0L, 0L, 1L, 0L))
  }
  expect_identical(.bo_strip(a$out$record), .bo_strip(b$out$record))
})

test_that("option flexyBayes.inla_best_of_workers caps the fits running at once", {
  skip_if_not_installed("INLA")
  .bo_mock_memory()
  fits <- list(.bo_fake(-9), .bo_fake(-8), .bo_fake(-7))

  withr::local_options(flexyBayes.inla_best_of_workers = 2L)
  calls <- .bo_mock_isolated(fits, delays = 0.1)
  pool <- .fb_inla_best_of_launch(list(), 4L, NULL, NULL)
  # The session fit is one of the two.
  expect_identical(calls$n, 1L)
  rec <- .fb_inla_best_of(.bo_fake(-10), list(), 4L, pool)$record
  expect_identical(calls$n, 3L)
  expect_identical(calls$max_running, 2L)
  expect_identical(rec$workers, 2L)
  expect_true(rec$parallel)
  expect_identical(rec$chosen, 4L)

  withr::local_options(flexyBayes.inla_best_of_workers = 1L)
  calls <- .bo_mock_isolated(fits)
  pool <- .fb_inla_best_of_launch(list(), 4L, NULL, NULL)
  expect_identical(calls$n, 0L)
  rec <- .fb_inla_best_of(.bo_fake(-10), list(), 4L, pool)$record
  expect_identical(calls$max_running, 1L)
  expect_identical(rec$workers, 1L)
  expect_false(rec$parallel)

  withr::local_options(flexyBayes.inla_best_of_workers = 10L)
  expect_identical(.fb_inla_best_of_workers(3L, NULL, NULL), 3L)
  withr::local_options(flexyBayes.inla_best_of_workers = NULL)
  expect_identical(.fb_inla_best_of_workers(3L, NULL, NULL), 3L)
  for (bad in list(0L, -1, 2.5, "2", NA, c(1L, 2L), Inf)) {
    withr::local_options(flexyBayes.inla_best_of_workers = bad)
    expect_error(
      .fb_inla_best_of_workers(3L, NULL, NULL),
      "option flexyBayes.inla_best_of_workers must be a whole number"
    )
  }
})

test_that("the memory guard lowers the fits at once, and says why", {
  .bo_mock_memory(free_gb = 12, fit_gb = 8)
  expect_message(
    w <- .fb_inla_best_of_workers(3L, NULL, NULL),
    paste0(
      "the 3 INLA fits run one after another: 12.0 GB of memory is free ",
      "and one fit needs about 8.0 GB"
    ),
    fixed = TRUE
  )
  expect_identical(w, 1L)

  .bo_mock_memory(free_gb = 20, fit_gb = 8)
  expect_message(
    w <- .fb_inla_best_of_workers(3L, NULL, NULL),
    "2 of the 3 INLA fits run at once: 20.0 GB of memory is free",
    fixed = TRUE
  )
  expect_identical(w, 2L)

  .bo_mock_memory(free_gb = 64, fit_gb = 8)
  expect_silent(w <- .fb_inla_best_of_workers(3L, NULL, NULL))
  expect_identical(w, 3L)

  .bo_mock_memory(free_gb = NA, fit_gb = 8)
  expect_message(
    w <- .fb_inla_best_of_workers(3L, NULL, NULL),
    "the 3 INLA fits run one after another: the free memory could not be read",
    fixed = TRUE
  )
  expect_identical(w, 1L)

  .bo_mock_memory(free_gb = 64, fit_gb = NA)
  expect_message(
    w <- .fb_inla_best_of_workers(3L, NULL, NULL),
    "the memory of one fit could not be estimated",
    fixed = TRUE
  )
  expect_identical(w, 1L)

  # No probe when the option asks for one fit at a time.
  testthat::local_mocked_bindings(
    .fb_free_memory_bytes = function() stop("no probe"),
    .fb_inla_fit_bytes = function(fb, data) stop("no estimate")
  )
  withr::local_options(flexyBayes.inla_best_of_workers = 1L)
  expect_silent(w <- .fb_inla_best_of_workers(3L, NULL, NULL))
  expect_identical(w, 1L)
})

test_that("a sequential fallback starts no fit before the session fit", {
  skip_if_not_installed("INLA")
  .bo_mock_memory(free_gb = 2, fit_gb = 1.5)
  calls <- .bo_mock_isolated(list(.bo_fake(-9), .bo_fake(-8)))
  expect_message(
    pool <- .fb_inla_best_of_launch(list(), 3L, NULL, NULL),
    "run one after another"
  )
  expect_identical(calls$n, 0L)
  rec <- .fb_inla_best_of(.bo_fake(-10), list(), 3L, pool)$record
  expect_identical(calls$max_running, 1L)
  expect_false(rec$parallel)
  expect_identical(rec$chosen, 3L)
})

test_that("in parallel the lowest-numbered failure is reported and the rest stopped", {
  skip_if_not_installed("INLA")
  .bo_mock_memory()
  # Fit 3 fails first, fit 2 fails later: fit 2 is named.
  calls <- .bo_mock_isolated(
    list(simpleError("late"), simpleError("early"), .bo_fake(-1)),
    delays = c(0.3, 0.05, 5)
  )
  pool <- .fb_inla_best_of_launch(list(), 4L, NULL, NULL)
  expect_error(
    .fb_inla_best_of(.bo_fake(-10), list(), 4L, pool),
    "INLA fit 2 of 4, run in a separate R process, failed: late",
    fixed = TRUE
  )
  # Fit 4, still running, was stopped.
  .fb_inla_best_of_stop(pool)
  expect_identical(calls$stopped, 3L)
})

test_that("the memory of one fit is estimated from the design-memory preflight", {
  fb <- fb_from_brms(y ~ trt + x + (1 | site), data = .bo_data())
  est <- .fb_inla_fit_bytes(fb, .bo_data())
  pf <- .fb_preflight(fb, .fb_dataset(.bo_data()))
  expect_equal(
    est,
    .FB_BEST_OF_BASE_BYTES +
      .FB_BEST_OF_ESTIMATE_FACTOR * as.numeric(pf$memory_estimate)
  )
  expect_true(is.na(.fb_inla_fit_bytes(NULL, NULL)))
})

test_that("free memory is read from vm_stat on macOS", {
  skip_on_os(c("windows", "linux", "solaris"))
  free <- .fb_free_memory_bytes()
  expect_true(is.finite(free) && free > 0)
  expect_lt(free, .fb_probe_system_ram())
})

test_that("best of 3 runs two fits in fresh R processes, one after another, on a package fixture", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("callr")
  skip_if_not_long_tier()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1",
    flexyBayes.inla_best_of = 3L,
    flexyBayes.inla_best_of_workers = 1L
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
  expect_false(rec$parallel)
})

test_that("best of 3 runs its fits at once in fresh R processes on a package fixture", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("callr")
  skip_if_not_long_tier()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1",
    flexyBayes.inla_best_of = 3L,
    flexyBayes.inla_best_of_workers = 3L
  )
  # The fixture is small. The guard is tested above.
  testthat::local_mocked_bindings(.fb_free_memory_bytes = function() {
    64 * 1024^3
  })
  t0 <- proc.time()[["elapsed"]]
  fit <- .bo_fit()
  elapsed <- proc.time()[["elapsed"]] - t0
  rec <- fit$extras$reproducibility$best_of
  expect_true(rec$parallel)
  expect_identical(rec$workers, 3L)
  expect_true(all(is.finite(rec$mlik)))
  expect_identical(rec$chosen, which.max(rec$mlik))
  expect_equal(
    rec$mlik[[rec$chosen]],
    unname(fit$inla$mlik[1L, 1L]),
    tolerance = 1e-8
  )
  expect_identical(rec$mode_status, rep(0L, 3L))
  # The fits overlapped: together they took longer than the whole call.
  expect_gt(sum(rec$wall_s), elapsed)
})

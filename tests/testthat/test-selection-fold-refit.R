# =============================================================================
# What an INLA fold refit inside fb_select() computes (R/selection_fold_refit.R).
#
# A fold refit is read for its variance components, its fitted values or its
# posterior draws on the held-out rows. It asks INLA for that readout only and
# runs no field test. On one fold of a 2,880-row fixture with a per-trial AR1
# field, the lean and full refits took the same 391 function evaluations to
# the same marginal likelihood and the same variance components. Through four
# folds of fb_select(covariance = "base_model") the fold variance components
# agreed to 2.8e-3 relative (2.2e-3 between two full runs) and the scores to
# 1.0e-5 (2.5e-5 between two full runs), with the same selected support.
# =============================================================================

fbfr_quiet <- function(env = parent.frame()) {
  withr::local_options(
    list(
      flexyBayes.inla_threads = "1:1",
      flexyBayes.silence_default_prior_note = TRUE,
      flexyBayes.silence_high_missingness_warning = TRUE,
      flexyBayes.silence_spatial_collapse_warning = TRUE
    ),
    .local_envir = env
  )
}

fbfr_site_data <- function(seed = 20261001L, n_site = 12L, per = 10L) {
  .fb_sel_with_seed(seed, {
    n <- n_site * per
    d <- data.frame(
      site = factor(rep(sprintf("S%02d", seq_len(n_site)), each = per))
    )
    for (v in c("x", "z1", "z2", "z3")) {
      d[[v]] <- stats::rnorm(n)
    }
    u <- stats::rnorm(n_site, sd = 0.8)
    d$y <- 2 + 0.5 * d$x + 0.7 * d$z1 + u[as.integer(d$site)] +
      stats::rnorm(n, sd = 0.5)
    d
  })
}

# Every INLA call a block makes, with its arguments.
fbfr_capture <- function(env = parent.frame()) {
  seen <- new.env()
  seen$calls <- list()
  real <- .inla_call
  testthat::local_mocked_bindings(
    .inla_call = function(args) {
      seen$calls[[length(seen$calls) + 1L]] <- args
      real(args)
    },
    .env = env
  )
  seen
}

fbfr_rel <- function(a, b) {
  max(abs(a - b)) / max(abs(b))
}


test_that("a fold refit asks INLA for its readout only", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbfr_quiet()
  d <- fbfr_site_data()
  base <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                     verbose = FALSE)
  fd <- base$data
  fd$y[fd$site %in% c("S01", "S02", "S03")] <- NA
  seen <- fbfr_capture()

  full <- .fb_sel_bayes_refit(base, fd, y ~ x, "inla", na_action = "augment")
  expect_identical(
    seen$calls[[1L]]$control.compute,
    list(config = TRUE, return.marginals = TRUE, dic = TRUE, waic = TRUE)
  )
  expect_null(seen$calls[[1L]]$control.mode)

  lean <- list()
  for (readout in c("varcomp", "fitted", "draws")) {
    seen$calls <- list()
    lean[[readout]] <- .fb_sel_bayes_refit(
      base, fd, y ~ x, "inla", na_action = "augment", readout = readout
    )
    expect_identical(
      seen$calls[[1L]]$control.compute,
      list(
        config = identical(readout, "draws"),
        return.marginals = FALSE,
        dic = FALSE,
        waic = FALSE
      ),
      info = readout
    )
    expect_null(seen$calls[[1L]]$control.mode)
    expect_null(lean[[readout]]$inla$marginals.fixed)
    expect_false(is.list(lean[[readout]]$inla$dic))
  }
  # The plan ends with the refit.
  expect_null(.fb_fold_refit_active())

  # Each readout agrees with the full refit's.
  vc_full <- .fb_sel_summary(full)$varcomp
  vc_lean <- .fb_sel_summary(lean$varcomp)$varcomp
  expect_identical(vc_lean$component, vc_full$component)
  expect_lt(fbfr_rel(vc_lean$estimate, vc_full$estimate), 5e-3)
  expect_lt(fbfr_rel(vc_lean$std.error, vc_full$std.error), 5e-3)
  expect_lt(fbfr_rel(stats::fitted(lean$fitted), stats::fitted(full)), 5e-3)
  assess <- which(is.na(fd$y))
  draws <- .fb_sel_posterior_predict_fold(
    lean$draws, assess, fd, prediction = "known_cluster", n_draws = 200L,
    seed = 1L
  )
  expect_identical(dim(draws$mu), c(200L, length(assess)))
})


test_that("a fold refit runs no field test, and a full fit still does", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbfr_quiet()
  # No spatial signal, so the full fit's field test refits without the field.
  d <- .fb_sel_with_seed(5L, {
    grid <- expand.grid(col = seq_len(6L), row = seq_len(8L))
    out <- data.frame(row = factor(grid$row), col = factor(grid$col))
    out$y <- 3 + stats::rnorm(nrow(out), sd = 0.6)
    out
  })
  base <- flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = d,
                     backend = "inla", verbose = FALSE)
  fd <- base$data
  fd$y[1:8] <- NA
  seen <- fbfr_capture()

  # The field test's fit is the one call that asks for the marginal
  # likelihood alone.
  field_test_calls <- function() {
    sum(vapply(
      seen$calls,
      function(a) identical(a$control.compute, list(mlik = TRUE)),
      logical(1L)
    ))
  }

  full <- suppressWarnings(
    .fb_sel_bayes_refit(base, fd, y ~ 1, "inla", na_action = "augment")
  )
  expect_true(isTRUE(full$extras$field_test$refit))
  expect_identical(field_test_calls(), 1L)

  seen$calls <- list()
  lean <- .fb_sel_bayes_refit(
    base, fd, y ~ 1, "inla", na_action = "augment", readout = "varcomp"
  )
  expect_null(lean$extras$field_test)
  expect_identical(field_test_calls(), 0L)
})


test_that("fb_select() with lean fold refits selects and scores as with full ones", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  skip_if_not_long_tier()
  fbfr_quiet()
  d <- fbfr_site_data()
  d$w <- rep(c(1, 0.6), length.out = nrow(d))
  base <- flexybayes(y ~ x, random = ~site, data = d, weights = d$w,
                     backend = "inla", verbose = FALSE)
  run <- function() {
    fb_select(
      base,
      ~ z1 + z2 + z3,
      penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
      workflow = "penalized_optimizer",
      resamples = fb_resamples(unit = ~site, v = 4L, weights_by = ~site),
      covariance = "base_model",
      seed = 1L
    )
  }
  lean <- run()
  full <- testthat::with_mocked_bindings(
    run(),
    .fb_fold_refit_plan = function(readout) NULL
  )
  expect_identical(lean$selection$chosen_terms, full$selection$chosen_terms)
  # The splits are exact; each fold's weights come from an engine refit.
  for (f in seq_along(lean$folds)) {
    a <- lean$folds[[f]]
    b <- full$folds[[f]]
    expect_identical(a$train_rows, b$train_rows)
    expect_identical(a$assess_rows, b$assess_rows)
    expect_lt(fbfr_rel(a$weights, b$weights), 1e-4)
  }
  expect_lt(max(abs(lean$scores$score - full$scores$score)), 1e-3)
  for (f in seq_along(lean$certificate$covariance$folds)) {
    a <- lean$certificate$covariance$folds[[f]]
    b <- full$certificate$covariance$folds[[f]]
    expect_lt(fbfr_rel(c(a$residual, a$terms), c(b$residual, b$terms)), 0.02)
  }
})

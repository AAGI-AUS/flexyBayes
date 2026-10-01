# =============================================================================
# fb_select(keep_fold_data = TRUE) and the three functions that read the kept
# fold data: fb_select_fold_data(), fb_select_path() and fb_select_refit().
#
# Oracles.
#   - The run itself: the refit of the empty and the full column set equals
#     the benchmark's protected and full models, and the path at the run's
#     alpha equals the run's supports.
#   - A second run at another alpha: the path re-solved at that alpha from
#     the first run's fold data equals the second run's supports.
#   - The fold's whitening written out with the stored factor.
# The GLS against a written-out V is in test-known-covariance-two-stage.R.
# =============================================================================

fbsa_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE,
    .local_envir = envir
  )
  if (requireNamespace("INLA", quietly = TRUE)) {
    old <- INLA::inla.getOption("num.threads")
    INLA::inla.setOption(num.threads = "1:1")
    withr::defer(INLA::inla.setOption(num.threads = old), envir = envir)
  }
}

fbsa_data <- function(seed = 3L) {
  .fb_sel_with_seed(seed, {
    d <- data.frame(
      site = factor(rep(sprintf("s%02d", 1:8), each = 6L)),
      x1 = stats::rnorm(48L),
      x2 = stats::rnorm(48L),
      x3 = stats::rnorm(48L),
      x4 = stats::rnorm(48L)
    )
    d$y <- 1 + 0.6 * d$x1 - 0.3 * d$x2 + stats::rnorm(8L)[d$site] +
      stats::rnorm(48L)
    d
  })
}

fbsa_supports <- function(sel) {
  pkg <- as.data.frame(sel$benchmark$glmnet$support$by_lambda)
  pkg <- pkg[pkg$arm == "benchmark", ]
  pkg$support <- vapply(
    strsplit(pkg$terms, ", ", fixed = TRUE),
    function(t) paste(sort(t[nzchar(t)]), collapse = ";"),
    character(1L)
  )
  pkg[order(pkg$fold, pkg$lambda_index), c("fold", "lambda_index", "s",
                                            "support")]
}

fbsa_net <- function(alpha) {
  fb_penalty(term(everything()) ~ elastic_net(lambda = 1, alpha = alpha))
}


test_that("fold data is kept on request only, and its absence refuses by name", {
  refused <- "flexybayes_refusal_selection_fold_data_not_kept"
  expect_error(fb_select(NULL, ~x1, keep_fold_data = TRUE), class = refused)
  expect_error(fb_select(NULL, ~x1, keep_fold_data = NA),
               "`keep_fold_data` must be TRUE or FALSE")
  bare <- structure(list(), class = "fb_select")
  expect_error(fb_select_fold_data(bare), class = refused)
  expect_error(fb_select_path(bare), class = refused)
  expect_error(fb_select_refit(bare, "x1"), class = refused)
})

test_that("the refit of the empty and full sets is the benchmark's protected and full model", {
  skip_if_not_installed("INLA")
  fbsa_quiet()
  d <- fbsa_data()
  base <- flexybayes(y ~ 1, random = ~site, data = d, backend = "inla",
                     verbose = FALSE)
  sel <- fb_select(base, ~ x1 + x2 + x3 + x4,
                   workflow = "unpenalized_benchmark",
                   resamples = fb_resamples(unit = ~site, v = 4L),
                   seed = 1L, keep_fold_data = TRUE)
  clusters <- as.data.frame(sel$clusters)
  for (set in list(
    list(model = "protected_only", columns = character(0)),
    list(model = "full", columns = paste0("x", 1:4))
  )) {
    got <- fb_select_refit(sel, set$columns)$scores
    want <- clusters[clusters$model == set$model, ]
    want <- want[match(paste(got$fold, got$cluster),
                       paste(want$fold, want$cluster)), ]
    expect_equal(got$mse, want$mse, tolerance = 1e-12)
    expect_equal(got$nlpd, want$nlpd, tolerance = 1e-12)
  }

  # The fold data: rows as the run's folds, whitening by the stored factor.
  folds <- fb_select_fold_data(sel)
  expect_identical(names(folds), vapply(sel$folds, function(f) f$label, ""))
  for (split in sel$folds) {
    fold <- folds[[split$label]]
    expect_identical(fold$train_rows, split$train_rows)
    expect_identical(fold$assess_rows, split$assess_rows)
    expect_identical(fold$response, d$y)
    expect_equal(
      colMeans(fold$x_candidate[fold$train_rows, ]),
      c(x1 = 0, x2 = 0, x3 = 0, x4 = 0),
      tolerance = 1e-12
    )
  }
  expect_named(fb_select_fold_data(sel, folds = "2"), "2")
  expect_error(fb_select_fold_data(sel, folds = "9"), "does not have: 9")
  expect_error(fb_select_refit(sel, "z9"), "not candidates of the run: z9")
  expect_error(fb_select_path(sel),
               class = "flexybayes_refusal_selection_penalty_spec_malformed")
})

test_that("the path at the run's alpha is the run's, and at another alpha is that run's", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  fbsa_quiet()
  d <- fbsa_data()
  base <- flexybayes(y ~ 1, random = ~site, data = d, backend = "inla",
                     verbose = FALSE)
  run <- function(penalty) {
    fb_select(base, ~ x1 + x2 + x3 + x4, penalty = penalty,
              workflow = "penalized_optimizer",
              decision = fb_decision(method = "full", metric = "rmse",
                                     n_lambda = 12L),
              resamples = fb_resamples(unit = ~site, v = 4L), seed = 1L,
              keep_fold_data = TRUE)
  }
  at_one <- run(fbsa_net(1))
  at_low <- run(fbsa_net(0.2))
  columns <- c("fold", "lambda_index", "s", "support")

  same <- fb_select_path(at_one)
  expect_identical(attr(same, "alpha"), 1)
  rownames(same) <- NULL
  expect_identical(same[, columns],
                   `rownames<-`(fbsa_supports(at_one), NULL))

  moved <- fb_select_path(at_one, alpha = 0.2)
  expect_identical(attr(moved, "alpha"), 0.2)
  expect_identical(moved[, columns],
                   `rownames<-`(fbsa_supports(at_low), NULL))
  expect_false(identical(moved$support, same$support))

  # Another grid and a subset of folds.
  short <- fb_select_path(at_one, n_lambda = 5L, folds = c("1", "3"))
  expect_identical(nrow(short), 10L)
  expect_setequal(short$fold, c("1", "3"))

  expect_error(fb_select_path(at_one, alpha = 0), "single number in \\(0, 1\\]")
  laplace <- run(fb_penalty(term(everything()) ~ laplace(scale = 0.3)))
  expect_error(
    fb_select_path(laplace, alpha = 0.5),
    class = "flexybayes_refusal_selection_penalty_spec_malformed"
  )
  expect_s3_class(fb_select_path(laplace), "data.frame")
})

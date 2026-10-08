# predict(allow_new_levels = TRUE) against ASReml 4.2's plug-in prediction of
# trial, cell and locality means at a locality the fit did not see. The
# ASReml side is computed by data-raw/asreml_predict_new_level.R and stored in
# fixtures/asreml_predict_new_level.rds: the REML variances, the fixed-part
# cell predictions with their covariance (predict(vcov = TRUE)), and the
# variance of each set mean under the REML variances, computed there as a
# dense covariance and by the sum-of-shares formula of the comparison study.

.pna_cache <- new.env(parent = emptyenv())

.pna_setup <- function() {
  if (is.null(.pna_cache$fit)) {
    fx <- readRDS(testthat::test_path(
      "fixtures",
      "asreml_predict_new_level.rds"
    ))
    train <- droplevels(fx$data[fx$data$loc != fx$held_out, ])
    test <- fx$data[fx$data$loc == fx$held_out, ]
    rownames(test) <- NULL
    withr::local_options(flexyBayes.inla_threads = "1:1")
    fit <- suppressMessages(flexybayes(
      stats::as.formula(fx$call$fixed),
      random = stats::as.formula(fx$call$random),
      data = train,
      backend = "inla",
      verbose = FALSE
    ))
    .pna_cache$fx <- fx
    .pna_cache$test <- test
    .pna_cache$fit <- fit
  }
  .pna_cache
}

test_that("at ASReml's REML variances the set variance is ASReml's plug-in", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  s <- .pna_setup()
  fx <- s$fx
  expect_identical(fx$stamp$asreml, "4.2.0.482")
  # One "draw" holding ASReml's variances, as standard deviations, and its
  # field correlations.
  at_reml <- as.data.frame(
    as.list(c(sqrt(fx$variance), fx$correlation)),
    check.names = FALSE
  )
  for (nm in names(fx$sets)) {
    nd <- s$test[fx$sets[[nm]], ]
    plan <- .fb_predict_new_plan(s$fit, nd, positions = TRUE)
    one_set <- list(id = rep(1L, nrow(nd)), size = nrow(nd))
    v <- .fb_predict_new_moments(
      plan,
      one_set,
      matrix(0, nrow(nd), 1L),
      at_reml
    )$variance
    expected <- fx$prediction$var_random[fx$prediction$set == nm]
    expect_equal(drop(v), expected, tolerance = 1e-10, label = nm)
  }
})

test_that("the posterior predictive agrees with ASReml's plug-in prediction", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  s <- .pna_setup()
  fx <- s$fx
  set.seed(20261008L)
  by_trial <- predict(
    s$fit,
    newdata = s$test,
    allow_new_levels = TRUE,
    mean_by = ~TRIAL,
    n_draws = 4000L
  )
  set.seed(20261009L)
  by_cell <- predict(
    s$fit,
    newdata = s$test,
    allow_new_levels = TRUE,
    mean_by = ~ TRIAL + VARIETY + NRATE,
    n_draws = 4000L
  )
  set.seed(20261010L)
  by_loc <- predict(
    s$fit,
    newdata = s$test,
    allow_new_levels = TRUE,
    mean_by = ~loc,
    n_draws = 4000L
  )
  cell_key <- paste(by_cell$TRIAL, by_cell$VARIETY, by_cell$NRATE, sep = "_")
  ours <- rbind(
    by_trial[c("estimate", "sd")],
    by_cell[
      match(c("T31_V1_lo", "T31_V6_hi", "T32_V3_hi"), cell_key),
      c("estimate", "sd")
    ],
    by_loc[c("estimate", "sd")]
  )
  theirs <- fx$prediction
  expect_identical(
    theirs$set,
    c(
      "trial_T31",
      "trial_T32",
      "cell_T31_V1_lo",
      "cell_T31_V6_hi",
      "cell_T32_V3_hi",
      "locality_L16"
    )
  )

  # Mean: both estimate the same fixed-effect function, the posterior under a
  # vague prior and ASReml by GLS at the REML variances, so they differ by a
  # small fraction of its standard error. A quarter of it is allowed; the
  # Monte Carlo error of 4000 draws is under a tenth.
  se_fixed <- sqrt(theirs$var_fixed)
  expect_lt(max(abs(ours$estimate - theirs$mean) / se_fixed), 0.25)

  # Variance: the posterior predictive variance averages each component's
  # posterior over the draws, where ASReml plugs in its REML estimate, so the
  # ratio of the two lies between the smallest and the largest ratio of a
  # component's posterior mean variance to its REML value. Those ratios are
  # read from summary()$varcomp (posterior mean of a variance = estimate^2 +
  # std.error^2), independently of predict(). The allowance of 0.05 covers
  # the parts that do not follow this rule: the fixed-effect spread (4.5% of
  # the plug-in variance of a trial mean), the field's correlations, and the
  # Monte Carlo error of 4000 draws (about 1%).
  vc <- summary(s$fit)$varcomp
  vc <- vc[match(names(fx$variance), vc$component), ]
  ratio <- (vc$estimate^2 + vc$std.error^2) / fx$variance
  observed <- ours$sd^2 / theirs$sd^2
  expect_true(all(observed > min(ratio) - 0.05), label = "variance ratio floor")
  expect_true(
    all(observed < max(ratio) + 0.05),
    label = "variance ratio ceiling"
  )
})

# predict(allow_new_levels = TRUE): the posterior predictive at levels the fit
# did not see. The correspondence to ASReml's plug-in prediction is in
# test-predict-new-levels-asreml.R; this file holds the contract.

.pnl_cache <- new.env(parent = emptyenv())

.pnl_fixture <- function() {
  readRDS(testthat::test_path("fixtures", "asreml_predict_new_level.rds"))
}

.pnl_split <- function() {
  fx <- .pnl_fixture()
  train <- droplevels(fx$data[fx$data$loc != fx$held_out, ])
  test <- fx$data[fx$data$loc == fx$held_out, ]
  rownames(test) <- NULL
  list(fx = fx, train = train, test = test)
}

.pnl_fit <- function() {
  if (is.null(.pnl_cache$fit)) {
    s <- .pnl_split()
    withr::local_options(flexyBayes.inla_threads = "1:1")
    .pnl_cache$fit <- suppressMessages(flexybayes(
      stats::as.formula(s$fx$call$fixed),
      random = stats::as.formula(s$fx$call$random),
      data = s$train,
      backend = "inla",
      verbose = FALSE
    ))
  }
  .pnl_cache$fit
}

.pnl_small_fit <- function(name, ...) {
  if (is.null(.pnl_cache[[name]])) {
    withr::local_options(flexyBayes.inla_threads = "1:1")
    .pnl_cache[[name]] <- suppressMessages(flexybayes(
      ...,
      backend = "inla",
      verbose = FALSE
    ))
  }
  .pnl_cache[[name]]
}


test_that("a new locality gets the posterior predictive table", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .pnl_fit()
  test <- .pnl_split()$test
  set.seed(1)
  p <- predict(fit, newdata = test, allow_new_levels = TRUE, mean_by = ~TRIAL)
  expect_s3_class(p, "fb_predict_new")
  expect_named(p, c("TRIAL", "n_rows", "estimate", "sd", "lower", "upper"))
  expect_equal(as.character(p$TRIAL), c("T31", "T32"))
  expect_equal(p$n_rows, c(24L, 24L))
  expect_true(all(p$lower < p$estimate & p$estimate < p$upper))
  draws <- attr(p, "draws")
  expect_equal(dim(draws), c(1000L, 2L))
  # The joint draws are a sample from the mixture the summary describes.
  expect_lt(max(abs(colMeans(draws) - p$estimate) / (p$sd / sqrt(1000))), 4)
  expect_lt(max(abs(apply(draws, 2L, stats::sd) / p$sd - 1)), 0.1)
  # The two trials share a locality, so their draws are correlated.
  expect_gt(stats::cor(draws[, 1L], draws[, 2L]), 0.2)
  expect_output(
    print(p),
    "Posterior predictive: mean over the rows of each TRIAL"
  )
})

test_that("the posterior predictive repeats under set.seed()", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .pnl_fit()
  test <- .pnl_split()$test
  set.seed(7)
  a <- predict(
    fit,
    newdata = test,
    allow_new_levels = TRUE,
    mean_by = ~TRIAL,
    n_draws = 200L
  )
  set.seed(7)
  b <- predict(
    fit,
    newdata = test,
    allow_new_levels = TRUE,
    mean_by = ~TRIAL,
    n_draws = 200L
  )
  expect_s3_class(a, "fb_predict_new")
  expect_identical(a, b)
})

test_that("mean_by takes a formula or column names, and NULL predicts rows", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .pnl_fit()
  test <- .pnl_split()$test
  set.seed(2)
  a <- predict(
    fit,
    newdata = test,
    allow_new_levels = TRUE,
    mean_by = ~ TRIAL + VARIETY,
    n_draws = 200L
  )
  set.seed(2)
  b <- predict(
    fit,
    newdata = test,
    allow_new_levels = TRUE,
    mean_by = c("TRIAL", "VARIETY"),
    n_draws = 200L
  )
  expect_identical(a, b)
  expect_equal(nrow(a), 12L)
  expect_true(all(a$n_rows == 4L))
  rows <- predict(
    fit,
    newdata = test[1:3, ],
    allow_new_levels = TRUE,
    n_draws = 200L
  )
  expect_named(rows, c("row", "n_rows", "estimate", "sd", "lower", "upper"))
  expect_equal(rows$row, 1:3)
})

test_that("levels the fit saw enter at their posterior draws", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .pnl_fit()
  train <- .fb_fit_data(fit)
  # Every level seen: the predictive mean of a fitted row is the posterior
  # mean of its linear predictor, which INLA reports for the fitted rows.
  set.seed(3)
  p <- predict(
    fit,
    newdata = train[1:60, ],
    allow_new_levels = TRUE,
    n_draws = 2000L
  )
  fv <- fit$inla$summary.fitted.values
  z <- (p$estimate - fv$mean[1:60]) / (fv$sd[1:60] / sqrt(2000))
  expect_lt(max(abs(z)), 4.5)
  # A row's predictive SD adds the residual to the SD of its linear predictor.
  sigma2 <- mean(fb_variance_draws(fit, 4000L, seed = 1L)$sigma^2)
  expect_lt(max(abs(p$sd^2 - fv$sd[1:60]^2 - sigma2)), 0.01)
})

test_that("a new trial at a seen locality takes that locality's effect", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .pnl_fit()
  test <- .pnl_split()$test
  at_seen <- test
  at_seen$loc <- factor("L01", levels = levels(fit$data$loc))
  set.seed(4)
  new_loc <- predict(
    fit,
    newdata = test,
    allow_new_levels = TRUE,
    mean_by = ~TRIAL,
    n_draws = 4000L
  )
  set.seed(4)
  seen_loc <- predict(
    fit,
    newdata = at_seen,
    allow_new_levels = TRUE,
    mean_by = ~TRIAL,
    n_draws = 4000L
  )
  effect <- fit$inla$summary.random$loc
  effect <- effect[effect$ID == "L01", ]
  expect_lt(
    abs(seen_loc$estimate[1L] - new_loc$estimate[1L] - effect$mean),
    4 * effect$sd / sqrt(4000)
  )
  # The locality variance leaves the prediction, and the fixed part is
  # replaced by the fixed part plus the locality's effect, drawn jointly.
  sd_loc2 <- mean(fb_variance_draws(fit, 4000L, seed = 1L)$sd_loc^2)
  set.seed(9)
  post <- fb_as_draws_simple(fit, n_samples = 4000L)
  X <- stats::model.matrix(
    ~ VARIETY * NRATE,
    data = test[test$TRIAL == "T31", ]
  )
  fixed_names <- rownames(fit$inla$summary.fixed)
  beta <- vapply(
    paste0(fixed_names, ":1"),
    function(k) post[[k]],
    numeric(4000L)
  )
  fixed_part <- drop(beta %*% colMeans(X)[fixed_names])
  u <- post[[paste0("loc:", which(fit$inla$summary.random$loc$ID == "L01"))]]
  expected <- sd_loc2 + stats::var(fixed_part) - stats::var(fixed_part + u)
  expect_equal(new_loc$sd[1L]^2 - seen_loc$sd[1L]^2, expected, tolerance = 0.1)
})

test_that("without allow_new_levels the INLA prediction is population-level", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .pnl_fit()
  test <- .pnl_split()$test
  X <- stats::model.matrix(~ VARIETY * NRATE, data = test)
  expect_equal(
    predict(fit, newdata = test),
    drop(X %*% coef(fit)[colnames(X)]),
    ignore_attr = TRUE
  )
})

test_that("the new-level path refuses arguments it does not combine with", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .pnl_fit()
  test <- .pnl_split()$test
  expect_error(predict(fit, newdata = test, mean_by = ~TRIAL), "mean_by")
  expect_error(
    predict(fit, newdata = test, allow_new_levels = TRUE, se.fit = TRUE),
    "se.fit = TRUE"
  )
  expect_error(
    predict(fit, newdata = test, allow_new_levels = TRUE, classify = "VARIETY"),
    "classify"
  )
  expect_error(predict(fit, allow_new_levels = TRUE), "newdata")
  expect_error(
    predict(fit, newdata = test, allow_new_levels = TRUE, mean_by = ~site),
    "site"
  )
})

test_that("a new level with no fitted variance is refused by name", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .pnl_fit()
  test <- .pnl_split()$test
  new_design <- test
  levels(new_design$design) <- c(levels(new_design$design), "lattice")
  new_design$design[] <- "lattice"
  expect_error(
    predict(fit, newdata = new_design, allow_new_levels = TRUE),
    class = "flexybayes_refusal"
  )
  cnd <- tryCatch(
    predict(fit, newdata = new_design, allow_new_levels = TRUE),
    error = function(e) e
  )
  expect_identical(cnd$reason_code, "predict_new_levels_no_variance")
  expect_match(conditionMessage(cnd), "diag(design):mainplot", fixed = TRUE)

  off_lattice <- test
  off_lattice$RANGE <- factor(
    as.character(off_lattice$RANGE),
    levels = c(levels(off_lattice$RANGE), "9")
  )
  off_lattice$RANGE[1L] <- "9"
  cnd <- tryCatch(
    predict(fit, newdata = off_lattice, allow_new_levels = TRUE),
    error = function(e) e
  )
  expect_identical(cnd$reason_code, "predict_new_levels_no_variance")
  expect_match(conditionMessage(cnd), "off it: 9", fixed = TRUE)

  s <- .pnl_split()
  dsum_fit <- .pnl_small_fit(
    "dsum",
    y ~ VARIETY,
    random = ~TRIAL,
    residual = ~ dsum(~ units | design),
    data = s$train[s$train$design == "rcb" | s$train$loc %in% c("L01", "L02"), ]
  )
  new_section <- s$test
  levels(new_section$design) <- c(levels(new_section$design), "lattice")
  new_section$design[] <- "lattice"
  cnd <- tryCatch(
    predict(dsum_fit, newdata = new_section, allow_new_levels = TRUE),
    error = function(e) e
  )
  expect_identical(cnd$reason_code, "predict_new_levels_no_variance")
  expect_match(conditionMessage(cnd), "dsum(~ units | design)", fixed = TRUE)
  # A section the fit saw is drawn from its own residual variance.
  set.seed(5)
  ok <- predict(
    dsum_fit,
    newdata = s$test,
    allow_new_levels = TRUE,
    mean_by = ~design,
    n_draws = 200L
  )
  expect_equal(sort(as.character(ok$design)), c("rcb", "split"))
})

test_that("a structure the new-level draw is not built for is refused", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  s <- .pnl_split()
  d <- s$train
  d$x <- as.numeric(as.character(d$RANGE))
  slope_fit <- .pnl_small_fit("slope", y ~ x + (x || TRIAL), data = d)
  nd <- s$test
  nd$x <- as.numeric(as.character(nd$RANGE))
  cnd <- tryCatch(
    predict(slope_fit, newdata = nd, allow_new_levels = TRUE),
    error = function(e) e
  )
  expect_s3_class(cnd, "flexybayes_refusal")
  expect_identical(
    cnd$reason_code,
    "predict_new_levels_structure_not_supported"
  )
  expect_match(
    conditionMessage(cnd),
    "simple_slope_uncor(TRIAL, x)",
    fixed = TRUE
  )
})

test_that("a response that is not Gaussian on the identity link is refused", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  s <- .pnl_split()
  d <- s$train
  set.seed(6)
  d$count <- stats::rpois(nrow(d), 3)
  pois_fit <- .pnl_small_fit(
    "poisson",
    count ~ VARIETY,
    random = ~TRIAL,
    family = "poisson",
    data = d
  )
  cnd <- tryCatch(
    predict(pois_fit, newdata = s$test, allow_new_levels = TRUE),
    error = function(e) e
  )
  expect_s3_class(cnd, "flexybayes_refusal")
  expect_identical(cnd$reason_code, "predict_new_levels_family_not_supported")
})

test_that("brms predicts a new locality from its own group-level SDs", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  skip_if_not_long_tier()
  s <- .pnl_split()
  fit <- suppressMessages(suppressWarnings(flexybayes(
    y ~ VARIETY,
    random = ~ loc + TRIAL,
    data = s$train,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  set.seed(8)
  p <- predict(fit, newdata = s$test, allow_new_levels = TRUE, mean_by = ~TRIAL)
  expect_s3_class(p, "fb_predict_new")
  expect_identical(attr(p, "engine"), "brms")
  expect_equal(dim(attr(p, "draws")), c(500L, 2L))
  # The predictive variance of a trial mean is the locality, trial and
  # residual variances (the last over 24 plots), plus the spread of the
  # fixed part, read from brms's own draws.
  post <- as.data.frame(posterior::as_draws_df(fit$brms))
  X <- stats::model.matrix(~VARIETY, data = s$test)
  b <- as.matrix(post[, paste0(
    "b_",
    sub("\\(Intercept\\)", "Intercept", colnames(X))
  )])
  fixed_part <- b %*% colMeans(X)
  expected <- mean(post$sd_loc__Intercept^2) +
    mean(post$sd_TRIAL__Intercept^2) +
    mean(post$sigma^2) / 24 +
    stats::var(drop(fixed_part))
  expect_equal(p$sd^2, rep(expected, 2L), tolerance = 0.15)
  expect_equal(
    p$estimate,
    rep(mean(fixed_part), 2L),
    tolerance = 0.1 * sqrt(expected)
  )
})

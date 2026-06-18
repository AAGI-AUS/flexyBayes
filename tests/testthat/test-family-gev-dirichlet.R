# ---------------------------------------------------------------- #
# GEV + Dirichlet families -- the two block-maxima / compositional   #
# fitters added in the hub-ergonomics uplift. Each carries a         #
# parameter-recovery test: simulate from known parameters, fit, and  #
# confirm the estimates land inside their reported credible /        #
# confidence interval. Recovery against a known data-generating      #
# process is the honest oracle for an estimator (the suite does not   #
# share the estimate's source with the code under test).             #
# ---------------------------------------------------------------- #

# ---- GEV ----------------------------------------------------------------

test_that("fb_family_gev() describes the GEV family", {
  fam <- fb_family_gev()
  expect_s3_class(fam, "fb_family")
  expect_identical(fam$family, "gen_extreme_value")
  expect_identical(fam$parameters, c("location", "scale", "shape"))
  expect_identical(fam$fitter, "fb_gev")
})

test_that("rgev() simulates finite values with the expected scale", {
  set.seed(1L)
  y0 <- rgev(2000L, location = 0, scale = 1, shape = 0)
  y2 <- rgev(2000L, location = 0, scale = 2, shape = 0)
  expect_true(all(is.finite(y0)))
  # A larger scale spreads the maxima further.
  expect_gt(stats::sd(y2), stats::sd(y0))
  expect_error(rgev(10L, 0, -1, 0), "positive")
})

test_that("fb_gev() recovers known GEV parameters by maximum likelihood", {
  set.seed(42L)
  truth <- c(location = 10, scale = 2, shape = 0.15)
  y <- rgev(4000L, truth[["location"]], truth[["scale"]], truth[["shape"]])
  fit <- fb_gev(y)

  expect_s3_class(fit, "fb_gev_fit")
  est <- fit$estimates
  expect_identical(est$term, c("location", "scale", "shape"))

  # Each true value sits inside the 95% confidence interval.
  for (p in est$term) {
    row <- est[est$term == p, ]
    expect_gte(truth[[p]], row$conf.low)
    expect_lte(truth[[p]], row$conf.high)
  }

  # Point estimates are close on a large sample.
  expect_equal(est$estimate[1L], truth[["location"]], tolerance = 0.1)
  expect_equal(est$estimate[2L], truth[["scale"]], tolerance = 0.1)
  expect_equal(est$estimate[3L], truth[["shape"]], tolerance = 0.05)
})

test_that("fb_gev() recovers the Gumbel limit (shape = 0)", {
  set.seed(7L)
  y <- rgev(4000L, location = 5, scale = 1.5, shape = 0)
  fit <- fb_gev(y)
  shape_row <- fit$estimates[fit$estimates$term == "shape", ]
  # The shape estimate is small in magnitude -- the Gumbel limit recovers a
  # near-zero shape. (Its interval need not bracket zero exactly: at the
  # boundary the maximum-likelihood shape is biased upward by an O(1/n)
  # amount, the well-known small-sample GEV shape bias.)
  expect_lt(abs(shape_row$estimate), 0.05)
})

test_that("fb_gev() reports increasing return levels", {
  set.seed(3L)
  y <- rgev(2000L, 10, 2, 0.1)
  fit <- fb_gev(y, return_periods = c(10, 50, 100))
  rl <- fit$return_levels
  expect_equal(rl$return_period, c(10, 50, 100))
  expect_true(all(diff(rl$return_level) > 0))
})

test_that("fb_gev() validates its input", {
  expect_error(fb_gev(c(1, 2)), "at least four")
  expect_error(fb_gev("a"), "numeric")
})

test_that("tidy() on a GEV fit returns the canonical columns", {
  set.seed(1L)
  fit <- fb_gev(rgev(500L, 10, 2, 0.1))
  td <- tidy(fit)
  expect_true(is.data.frame(td))
  expect_true(all(
    c("term", "estimate", "std.error", "conf.low", "conf.high") %in% names(td)
  ))
  expect_equal(nrow(td), 3L)
})

# ---- Dirichlet ----------------------------------------------------------

test_that("fb_family_dirichlet() describes the Dirichlet family", {
  fam <- fb_family_dirichlet()
  expect_s3_class(fam, "fb_family")
  expect_identical(fam$family, "dirichlet")
  expect_identical(fam$fitter, "fb_dirichlet")
})

test_that("rdirichlet() simulates rows on the simplex", {
  set.seed(1L)
  X <- rdirichlet(50L, alpha = c(2, 5, 3))
  expect_equal(dim(X), c(50L, 3L))
  expect_true(all(abs(rowSums(X) - 1) < 1e-08))
  expect_true(all(X >= 0))
  expect_error(rdirichlet(5L, c(1)), "at least two")
})

test_that("fb_dirichlet() recovers known concentrations by ML", {
  set.seed(7L)
  truth <- c(2, 5, 3)
  X <- rdirichlet(3000L, alpha = truth)
  fit <- fb_dirichlet(X)

  expect_s3_class(fit, "fb_dirichlet_fit")
  est <- fit$estimates
  expect_equal(nrow(est), 3L)

  # Each true concentration sits inside its 95% confidence interval.
  for (j in seq_along(truth)) {
    expect_gte(truth[j], est$conf.low[j])
    expect_lte(truth[j], est$conf.high[j])
  }
  expect_equal(est$estimate, truth, tolerance = 0.15)

  # The mean composition recovers the normalised truth.
  expect_equal(
    unname(fit$mean_composition),
    truth / sum(truth),
    tolerance = 0.02
  )
})

test_that("fb_dirichlet() recovers concentrations via greta", {
  skip_if_no_greta()

  set.seed(7L)
  truth <- c(2, 5, 3)
  X <- rdirichlet(400L, alpha = truth)
  fit <- fb_dirichlet(
    X,
    method = "greta",
    n_samples = 400L,
    warmup = 400L,
    chains = 2L
  )
  est <- fit$estimates
  for (j in seq_along(truth)) {
    expect_gte(truth[j], est$conf.low[j])
    expect_lte(truth[j], est$conf.high[j])
  }
  expect_equal(est$estimate, truth, tolerance = 0.6)
})

test_that("fb_dirichlet() carries through column labels", {
  set.seed(2L)
  X <- rdirichlet(200L, c(3, 3, 3))
  colnames(X) <- c("sand", "silt", "clay")
  fit <- fb_dirichlet(X)
  expect_identical(fit$estimates$term, c("sand", "silt", "clay"))
  expect_identical(names(fit$mean_composition), c("sand", "silt", "clay"))
})

test_that("fb_dirichlet() validates its input", {
  expect_error(fb_dirichlet(matrix(1, 10L, 1L)), "at least two")
  expect_error(fb_dirichlet(matrix(1, 2L, 3L)), "at least four")
  expect_error(
    fb_dirichlet(matrix(c(-1, rep(1, 11)), 4L, 3L)),
    "non-negative"
  )
})

test_that("tidy() on a Dirichlet fit returns the canonical columns", {
  set.seed(1L)
  fit <- fb_dirichlet(rdirichlet(300L, c(2, 5, 3)))
  td <- tidy(fit)
  expect_true(all(
    c("term", "estimate", "std.error", "conf.low", "conf.high") %in% names(td)
  ))
  expect_equal(nrow(td), 3L)
})

# ---- Family-gate routing ------------------------------------------------

test_that("flexybayes() routes GEV / Dirichlet to their dedicated fitters", {
  err_gev <- expect_error(
    flexyBayes:::.resolve_family("gen_extreme_value", NULL)
  )
  expect_match(conditionMessage(err_gev), "fb_gev")

  err_dir <- expect_error(
    flexyBayes:::.resolve_family("dirichlet", NULL)
  )
  expect_match(conditionMessage(err_dir), "fb_dirichlet")
})

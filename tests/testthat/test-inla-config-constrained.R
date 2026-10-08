# An INLA fit of random = ~ site + spl(x) whose stored precision, at one
# hyperparameter configuration, is singular to working precision: the rw2
# smooth and the intercept share a constant that only the smooth's
# sum-to-zero constraint fixes. The oracle is INLA::inla.posterior.sample()
# on the same fit, 100,000 draws, stored with the fit by
# data-raw/inla_spline_constrained.R. Moments are compared within four Monte
# Carlo standard errors.

.cc_ref <- function() {
  readRDS(testthat::test_path("fixtures", "inla_spline_constrained.rds"))
}

test_that("DEF-205: the stored precision does not factorise as it stands", {
  ref <- .cc_ref()
  k <- ref$failing_configs[[1L]]
  q <- Matrix::forceSymmetric(ref$fit$inla$misc$configs$config[[k]]$Q,
                              uplo = "U")
  expect_error(
    suppressWarnings(Matrix::Cholesky(q, LDL = FALSE)),
    "not positive"
  )
})

test_that("DEF-205: vcov() matches INLA's own posterior draws", {
  skip_if_not_installed("INLA")
  ref <- .cc_ref()
  fit <- ref$fit
  n <- ref$oracle$n_draws
  scale <- fit$extras$response_standardisation$scale
  v <- vcov(fit) / scale^2
  sd_inla <- sqrt(diag(ref$oracle$fixed_cov))
  # An SD estimated from n normal draws has standard error about sd / sqrt(2n).
  gap <- (sqrt(diag(v)) - sd_inla) / (sd_inla / sqrt(2 * n))
  expect_true(all(abs(gap) <= 4), label = paste(round(gap, 2), collapse = ", "))
  cor_fb <- stats::cov2cor(v)
  cor_inla <- stats::cov2cor(ref$oracle$fixed_cov)
  off <- upper.tri(cor_fb)
  expect_lt(max(abs(cor_fb[off] - cor_inla[off])), 4 / sqrt(n))
  # INLA's own marginal SDs, from its integration over the configurations.
  expect_equal(
    unname(sqrt(diag(vcov(fit)))),
    fit$inla$summary.fixed$sd,
    tolerance = 1e-3
  )
})

test_that("DEF-205: posterior draws meet the constraint and match INLA's", {
  skip_if_not_installed("INLA")
  ref <- .cc_ref()
  fit <- ref$fit
  n <- 4000L
  draws <- withr::with_seed(5L, .fb_inla_config_sample(fit, n))
  std <- fit$extras$response_standardisation
  spline <- (draws$latent[grep("^x:", rownames(draws$latent)), ]) / std$scale
  expect_lt(max(abs(colSums(spline))), 1e-8)
  fixed_rows <- paste0(rownames(fit$inla$summary.fixed), ":1")
  fixed <- draws$latent[fixed_rows, ]
  fixed[1L, ] <- fixed[1L, ] - std$centre
  fixed <- fixed / std$scale
  sd_inla <- sqrt(diag(ref$oracle$fixed_cov))
  mean_se <- sd_inla * sqrt(1 / n + 1 / ref$oracle$n_draws)
  gap_mean <- (rowMeans(fixed) - ref$oracle$fixed_mean) / mean_se
  expect_true(all(abs(gap_mean) <= 4),
              label = paste(round(gap_mean, 2), collapse = ", "))
  sd_se <- ref$oracle$spline_sd * sqrt(1 / (2 * n) + 1 / (2 * ref$oracle$n_draws))
  gap_sd <- (apply(spline, 1L, stats::sd) - ref$oracle$spline_sd) / sd_se
  expect_lt(max(abs(gap_sd)), 4.5)
})

test_that("DEF-205: the user-facing summaries of the fit are computed", {
  skip_if_not_installed("INLA")
  ref <- .cc_ref()
  fit <- ref$fit
  se <- predict(fit, se.fit = TRUE)$se.fit
  expect_true(all(is.finite(se)) && all(se > 0))
  cls <- as.data.frame(predict(fit, classify = "trt"))
  expect_identical(nrow(cls), 3L)
  expect_true(all(is.finite(cls$std.error)))
  simple <- withr::with_seed(5L, fb_as_draws_simple(fit))
  expect_true(all(is.finite(unlist(simple))))
})

test_that("a precision no constraint repairs is refused by name", {
  ref <- .cc_ref()
  fit <- ref$fit
  fit$inla$misc$configs$constr <- NULL
  err <- tryCatch(vcov(fit), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_inla_config_precision_not_positive_definite")
  expect_match(conditionMessage(err), "configuration", fixed = TRUE)
})

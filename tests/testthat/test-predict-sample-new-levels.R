# test-predict-sample-new-levels.R -- Stage 3B file-output sister
# (ADR 0023 §"sample" activation, v0.3.5). Covers the active
# `allow_new_levels = "sample"` branch: success on mixed
# known+unknown levels, caller-seed-driven non-determinism,
# variance scaling with the RE posterior tau, backward-compat with
# legacy fits, the no-draws fit refusal, and confirmation that the
# v0.3.4 deferred-stop is gone.

suppressPackageStartupMessages({
  library(testthat)
})

old_opts <- options(
  flexyBayes.silence_default_prior_note = TRUE,
  flexyBayes.silence_uniform_inla_approx = TRUE,
  flexyBayes.silence_auto_fallback_note = TRUE,
  flexyBayes.silence_auto_inla_missing_note = TRUE
)
on.exit(options(old_opts), add = TRUE)


mk_predict_fit <- function(seed = 2026L, N = 60L, J = 4L, prior_vc_sd = 1) {
  set.seed(seed)
  dat <- data.frame(
    x = rnorm(N),
    g = factor(sample(letters[seq_len(J)], N, replace = TRUE)),
    y = rnorm(N)
  )
  fit <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = dat,
    backend = "greta",
    n_samples = 50L,
    warmup = 30L,
    chains = 1L,
    prior_fixed_sd = 10,
    prior_vc_sd = prior_vc_sd,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  list(fit = fit, dat = dat)
}


# ---------------------------------------------------------------- #
# Active: success on mixed known+unknown levels, no warning         #
# ---------------------------------------------------------------- #

test_that("sample activation: known+unknown returns numeric of length nrow(newdata) without warning", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:5, ]
  nd$g <- as.character(nd$g)
  nd$g[c(2L, 4L)] <- c("ZZ", "WW")
  set.seed(42L)
  p <- expect_silent(
    predict(fx$fit, newdata = nd, allow_new_levels = "sample")
  )
  expect_length(p, 5L)
  expect_true(all(is.finite(p)))
})


# ---------------------------------------------------------------- #
# Caller-seed-driven non-determinism on unknown rows                #
# ---------------------------------------------------------------- #

test_that("sample activation: different seeds give different predictions on unknown rows", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:5, ]
  nd$g <- as.character(nd$g)
  nd$g[2L] <- "ZZ"
  nd$g[4L] <- "WW"
  set.seed(101L)
  p_a <- predict(fx$fit, newdata = nd, allow_new_levels = "sample")
  set.seed(202L)
  p_b <- predict(fx$fit, newdata = nd, allow_new_levels = "sample")
  # Unknown rows (positions 2 and 4) must differ across seeds.
  expect_false(isTRUE(all.equal(p_a[2L], p_b[2L])))
  expect_false(isTRUE(all.equal(p_a[4L], p_b[4L])))
})

test_that("sample activation: known rows are seed-independent", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:5, ]
  nd$g <- as.character(nd$g)
  nd$g[2L] <- "ZZ"
  set.seed(101L)
  p_a <- predict(fx$fit, newdata = nd, allow_new_levels = "sample")
  set.seed(202L)
  p_b <- predict(fx$fit, newdata = nd, allow_new_levels = "sample")
  # Known rows (positions 1, 3, 4, 5) must be identical across
  # seeds -- only the unknown-row RE sample consumes the RNG.
  known_idx <- c(1L, 3L, 4L, 5L)
  expect_identical(p_a[known_idx], p_b[known_idx])
})


# ---------------------------------------------------------------- #
# Variance scaling with posterior tau (file-output interval)        #
# ---------------------------------------------------------------- #

test_that("sample activation: file-output interval on unknown row reflects RE variance", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit(prior_vc_sd = 1)
  nd <- fx$dat[1:5, ]
  nd$g <- as.character(nd$g)
  nd$g[2L] <- "ZZ"

  f_sample <- tempfile(fileext = ".rds")
  f_pop <- tempfile(fileext = ".rds")
  on.exit(unlink(c(f_sample, f_pop)), add = TRUE)

  set.seed(101L)
  suppressWarnings(predict(
    fx$fit,
    newdata = nd,
    output_file = f_pop,
    format = "rds",
    allow_new_levels = "population"
  ))
  set.seed(101L)
  predict(
    fx$fit,
    newdata = nd,
    output_file = f_sample,
    format = "rds",
    allow_new_levels = "sample"
  )

  out_pop <- readRDS(f_pop)
  out_sample <- readRDS(f_sample)

  # Unknown row (position 2). Under "population", the row is NA
  # (model.matrix drops the unknown level). Under "sample", the row
  # has a finite prediction with a posterior interval reflecting
  # the sampled-RE variance. The interval width is the diagnostic
  # of RE-variance contribution.
  sample_width_unknown <- out_sample$upper[2L] - out_sample$lower[2L]
  sample_width_known <- out_sample$upper[1L] - out_sample$lower[1L]
  expect_true(is.finite(sample_width_unknown))
  # Sampled-RE row's interval is no narrower than the known-row
  # interval (RE contribution adds variance, never subtracts).
  expect_gte(sample_width_unknown, sample_width_known)
})


# ---------------------------------------------------------------- #
# Backward-compat: legacy fit (no extras$fb_dataset)                #
# ---------------------------------------------------------------- #

test_that("sample + legacy fit silently skips dictionary resolution + sampling", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  legacy_fit <- fx$fit
  legacy_fit$extras$fb_dataset <- NULL
  nd <- fx$dat[1:5, ]
  nd$g <- as.character(nd$g)
  nd$g[1L] <- "ZZ"
  # No warning, no stop -- legacy path bypasses both the dictionary
  # resolver and the sample branch.
  p <- expect_silent(
    predict(legacy_fit, newdata = nd, allow_new_levels = "sample")
  )
  expect_length(p, nrow(nd))
})


# ---------------------------------------------------------------- #
# Deferred-stop removal: v0.3.4 reserved label is gone              #
# ---------------------------------------------------------------- #

test_that("sample activation: v0.3.4 deferred-stop refusal is no longer raised", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:5, ]
  nd$g <- as.character(nd$g)
  nd$g[1L] <- "ZZ"
  err <- tryCatch(
    predict(fx$fit, newdata = nd, allow_new_levels = "sample"),
    error = function(e) e
  )
  # If err is NOT a condition object, predict succeeded -- the
  # deferred-stop is gone.
  if (inherits(err, "condition")) {
    expect_failure(
      expect_match(conditionMessage(err), "reserved at v0.3\\.4")
    )
    expect_failure(
      expect_match(conditionMessage(err), "deferred to v0.3\\.5")
    )
  } else {
    succeed("sample mode ran without raising the v0.3.4 deferred-stop")
  }
})


# ---------------------------------------------------------------- #
# Edge: no-draws fit raises structured refusal                      #
# ---------------------------------------------------------------- #

test_that("sample activation: fit without $greta$draws raises structured refusal", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  no_draws_fit <- fx$fit
  no_draws_fit$greta$draws <- NULL
  nd <- fx$dat[1:5, ]
  nd$g <- as.character(nd$g)
  nd$g[1L] <- "ZZ"
  err <- expect_error(
    predict(no_draws_fit, newdata = nd, allow_new_levels = "sample"),
    "posterior draws on \\$greta\\$draws"
  )
  expect_match(conditionMessage(err), "INLA draws adapter is")
})

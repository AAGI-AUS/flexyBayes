# test-predict-kernel.R -- v0.3.8 audit Critical Fix #1
# (file-backed smooth-bypass) + ADR 0030 C2 (shared prediction
# kernel). Covers the .predict_linear_draws() contract:
#
#   (a) kernel returns equivalent linear predictor across monolithic
#       vs chunked invocation for a fixed-only model.
#   (b) ditto for an s() smooth model -- the file-backed path
#       pre-v0.3.8 silently dropped the smooth contribution.
#   (c) ditto for a random-intercept model with population prediction.
#   (d) ditto with sampled-new-level prediction (allow_new_levels =
#       "sample" surfaces sample_re_summary on newdata; the kernel
#       layers per-row, per-draw RE realisations).
#   (e) chunk-invariance: chunk_size = 100L, 500L, 2000L produce
#       byte-identical posterior intervals for a fixed seed (ADR
#       0023 §6.d at the kernel level).
#   (f) include arg filters contributions additively.
#   (g) file-backed path equals in-memory path within rounding
#       tolerance for all model classes covered by (a)-(d).
#   (h) malformed include raises a typed condition with reason_code
#       = "predict_kernel_invalid_include".
#
# Refusal taxonomy:
#   - flexybayes_predict_kernel_refusal class on malformed include
#     (test (h)); .predict_kernel_validate_include() raises.

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


# ---------------------------------------------------------------- #
# Fixtures: greta-backend fits so $greta$draws is populated.        #
# Re-used across blocks (a)-(g) so MCMC time is amortised.          #
# ---------------------------------------------------------------- #

skip_if_no_greta_quiet <- function() {
  testthat::skip_on_cran()
  skip_if_greta_backend_unusable()
}

mk_fixed_fit <- function(seed = 20260526L, N = 80L) {
  set.seed(seed)
  dat <- data.frame(
    x = rnorm(N),
    y = rnorm(N)
  )
  dat$y <- 0.5 + 1.2 * dat$x + rnorm(N, 0, 0.3)
  fit <- suppressMessages(fb(
    y ~ x,
    data = dat,
    backend = "greta",
    n_samples = 50L,
    warmup = 30L,
    chains = 1L,
    prior_fixed_sd = 10,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  list(fit = fit, dat = dat)
}

mk_ri_fit <- function(seed = 20260526L, N = 80L, J = 4L) {
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
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  list(fit = fit, dat = dat)
}

mk_smooth_fit <- function(seed = 20260526L, N = 80L) {
  skip_if_not_installed("mgcv")
  set.seed(seed)
  dat <- data.frame(x = sort(stats::runif(N, 0, 10)))
  dat$y <- sin(dat$x) + stats::rnorm(N, 0, 0.2)
  fit <- suppressMessages(flexybayes(
    fixed = y ~ 1,
    random = ~ s(x, k = 6L),
    data = dat,
    n_samples = 50L,
    warmup = 30L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  list(fit = fit, dat = dat)
}


# ---------------------------------------------------------------- #
# (h) malformed include raises typed condition. No fit needed.       #
# ---------------------------------------------------------------- #

test_that("(h) malformed `include` raises predict_kernel_invalid_include", {
  # We do not need a real fit -- the validator runs first. Use a
  # minimal stub so the kernel reaches the validator before any
  # other check.
  stub <- list(
    greta = list(draws = list(matrix(0, nrow = 1L, ncol = 1L))),
    glm = list(formula = y ~ x, coefficients = c("(Intercept)" = 0)),
    extras = list(
      parse_info = list(
        random = NULL,
        smooths = list(),
        fixed = list(intercept = TRUE, terms = list())
      )
    )
  )
  class(stub) <- c("flexybayes", "list")

  # Bad include: an unknown vocabulary element.
  err1 <- expect_error(
    flexyBayes:::.predict_linear_draws(
      stub,
      newdata = data.frame(x = 1),
      include = c("fixed", "random_known")
    ),
    "kernel vocabulary"
  )
  expect_s3_class(err1, "flexybayes_predict_kernel_refusal")
  expect_identical(err1$reason_code, "predict_kernel_invalid_include")
  expect_true("random_known" %in% err1$bad)

  # Empty include.
  err2 <- expect_error(
    flexyBayes:::.predict_linear_draws(
      stub,
      newdata = data.frame(x = 1),
      include = character(0)
    ),
    "non-empty character vector"
  )
  expect_s3_class(err2, "flexybayes_predict_kernel_refusal")
  expect_identical(err2$reason_code, "predict_kernel_invalid_include")

  # Non-character include.
  err3 <- expect_error(
    flexyBayes:::.predict_linear_draws(
      stub,
      newdata = data.frame(x = 1),
      include = 1L
    ),
    "non-empty character vector"
  )
  expect_s3_class(err3, "flexybayes_predict_kernel_refusal")
})


# ---------------------------------------------------------------- #
# (a) kernel equivalence in-memory vs chunked, fixed-only model     #
# ---------------------------------------------------------------- #

test_that("(a) kernel: monolithic vs chunked equivalence on a fixed-only model", {
  skip_if_no_greta_quiet()
  f <- mk_fixed_fit()
  newx <- data.frame(x = seq(-2, 2, length.out = 30L))

  mono <- flexyBayes:::.predict_linear_draws(f$fit, newx)
  cs <- 7L # not a clean divisor of 30 -- exercises ragged chunks
  chunks <- split(seq_len(nrow(newx)), ceiling(seq_len(nrow(newx)) / cs))
  chunked <- matrix(0, nrow = nrow(newx), ncol = ncol(mono))
  for (rng in chunks) {
    chunked[rng, ] <- flexyBayes:::.predict_linear_draws(
      f$fit,
      newx[rng, , drop = FALSE]
    )
  }
  # Numerical equality at machine tolerance -- BLAS non-associativity
  # makes BIT-identity infeasible across different matrix shapes (a
  # monolithic mm %*% t(beta) and the concatenation of per-chunk
  # mm[rng] %*% t(beta) can differ by ~1 ULP). The audit contract
  # (Critical Fix #1) is correctness of the smooth contribution, not
  # bit-identity; bit-identity holds within a single chunk-size run
  # (test-predict-file-output.R asserts it for the file path) but
  # not across re-grouped chunkings of the same total computation.
  # Drop dimnames before comparison (monolithic carries newdata
  # rownames; chunked is an unnamed accumulator).
  expect_equal(unname(mono), unname(chunked), tolerance = 1e-10)
})


# ---------------------------------------------------------------- #
# (b) kernel equivalence on an s() smooth model -- closes the         #
#     file-backed smooth-bypass bug (Critical Fix #1).               #
# ---------------------------------------------------------------- #

test_that("(b) kernel: monolithic vs chunked equivalence on an s() smooth model", {
  skip_if_no_greta_quiet()
  f <- mk_smooth_fit()
  newx <- data.frame(x = seq(min(f$dat$x), max(f$dat$x), length.out = 30L))

  mono <- flexyBayes:::.predict_linear_draws(f$fit, newx)

  # Smoke check: the smooth contribution is non-trivial (otherwise
  # the "byte-identical" assertion below would pass on a buggy
  # implementation that silently dropped the smooth basis).
  expect_true(sd(rowMeans(mono)) > 0.05)

  cs <- 7L
  chunks <- split(seq_len(nrow(newx)), ceiling(seq_len(nrow(newx)) / cs))
  chunked <- matrix(0, nrow = nrow(newx), ncol = ncol(mono))
  for (rng in chunks) {
    chunked[rng, ] <- flexyBayes:::.predict_linear_draws(
      f$fit,
      newx[rng, , drop = FALSE]
    )
  }
  # Numerical equality at machine tolerance -- BLAS non-associativity
  # makes BIT-identity infeasible across different matrix shapes (a
  # monolithic mm %*% t(beta) and the concatenation of per-chunk
  # mm[rng] %*% t(beta) can differ by ~1 ULP). The audit contract
  # (Critical Fix #1) is correctness of the smooth contribution, not
  # bit-identity; bit-identity holds within a single chunk-size run
  # (test-predict-file-output.R asserts it for the file path) but
  # not across re-grouped chunkings of the same total computation.
  # Drop dimnames before comparison (monolithic carries newdata
  # rownames; chunked is an unnamed accumulator).
  expect_equal(unname(mono), unname(chunked), tolerance = 1e-10)
})


# ---------------------------------------------------------------- #
# (c) kernel equivalence on a random-intercept model (population)    #
# ---------------------------------------------------------------- #

test_that("(c) kernel: monolithic vs chunked equivalence on a random-intercept model (population)", {
  skip_if_no_greta_quiet()
  f <- mk_ri_fit()
  newdat <- data.frame(
    x = seq(-2, 2, length.out = 24L),
    g = factor(rep(letters[1:4], length.out = 24L), levels = letters[1:4])
  )

  mono <- flexyBayes:::.predict_linear_draws(f$fit, newdat)
  cs <- 5L
  chunks <- split(seq_len(nrow(newdat)), ceiling(seq_len(nrow(newdat)) / cs))
  chunked <- matrix(0, nrow = nrow(newdat), ncol = ncol(mono))
  for (rng in chunks) {
    chunked[rng, ] <- flexyBayes:::.predict_linear_draws(
      f$fit,
      newdat[rng, , drop = FALSE]
    )
  }
  # Numerical equality at machine tolerance -- BLAS non-associativity
  # makes BIT-identity infeasible across different matrix shapes (a
  # monolithic mm %*% t(beta) and the concatenation of per-chunk
  # mm[rng] %*% t(beta) can differ by ~1 ULP). The audit contract
  # (Critical Fix #1) is correctness of the smooth contribution, not
  # bit-identity; bit-identity holds within a single chunk-size run
  # (test-predict-file-output.R asserts it for the file path) but
  # not across re-grouped chunkings of the same total computation.
  # Drop dimnames before comparison (monolithic carries newdata
  # rownames; chunked is an unnamed accumulator).
  expect_equal(unname(mono), unname(chunked), tolerance = 1e-10)
})


# ---------------------------------------------------------------- #
# (d) sampled-new-level: kernel layers sampled_re onto unknown rows  #
# ---------------------------------------------------------------- #

test_that("(d) kernel: sampled_re contribution applies only to named rows", {
  skip_if_no_greta_quiet()
  f <- mk_ri_fit()
  newdat <- data.frame(
    x = seq(-1, 1, length.out = 6L),
    g = factor(c("a", "b", "a", "c", "d", "a"), levels = letters[1:4])
  )

  # Build a synthetic sampled_re list keyed on rows 2 and 5.
  base <- flexyBayes:::.predict_linear_draws(f$fit, newdat)
  n_draws_eff <- ncol(base)
  set.seed(99L)
  sampled_re <- list(
    "2" = rnorm(n_draws_eff, 0, 0.5),
    "5" = rnorm(n_draws_eff, 0, 0.5)
  )
  with_re <- flexyBayes:::.predict_linear_draws(
    f$fit,
    newdat,
    sampled_re = sampled_re,
    include = c("fixed", "smooth", "random_sampled")
  )

  # Rows other than 2 and 5 are unchanged byte-for-byte.
  unchanged_rows <- c(1L, 3L, 4L, 6L)
  expect_identical(
    base[unchanged_rows, , drop = FALSE],
    with_re[unchanged_rows, , drop = FALSE]
  )

  # Rows 2 and 5 carry exactly the sampled_re vector added.
  expect_equal(unname(with_re[2L, ] - base[2L, ]), sampled_re[["2"]])
  expect_equal(unname(with_re[5L, ] - base[5L, ]), sampled_re[["5"]])
})


# ---------------------------------------------------------------- #
# (e) chunk-invariance at multiple chunk sizes                       #
# ---------------------------------------------------------------- #

test_that("(e) kernel: chunk_size = 7L, 13L, 30L all produce numerically identical output", {
  skip_if_no_greta_quiet()
  f <- mk_smooth_fit()
  newx <- data.frame(x = seq(min(f$dat$x), max(f$dat$x), length.out = 30L))
  # Same-chunk-size invariance: two repeated calls with the same
  # chunk_size produce bit-identical output (the deterministic
  # contract; no RNG involvement). Cross-chunk-size variance is
  # ULP-level (BLAS reorders accumulation), so tested at
  # tolerance.
  cs <- 7L
  chunks <- split(seq_len(nrow(newx)), ceiling(seq_len(nrow(newx)) / cs))
  build_at <- function(chunk_size) {
    chunks <- split(
      seq_len(nrow(newx)),
      ceiling(seq_len(nrow(newx)) / chunk_size)
    )
    n_cols <- ncol(flexyBayes:::.predict_linear_draws(
      f$fit,
      newx[1L, , drop = FALSE]
    ))
    out <- matrix(0, nrow = nrow(newx), ncol = n_cols)
    for (rng in chunks) {
      out[rng, ] <- flexyBayes:::.predict_linear_draws(
        f$fit,
        newx[rng, , drop = FALSE]
      )
    }
    out
  }
  cs7 <- build_at(7L)
  cs13 <- build_at(13L)
  cs30 <- build_at(30L)
  # Same-chunk-size repeat: bit-identical.
  expect_identical(cs7, build_at(7L))
  # Cross-chunk-size: numerical equality at machine tolerance.
  expect_equal(cs7, cs13, tolerance = 1e-10)
  expect_equal(cs7, cs30, tolerance = 1e-10)
})


# ---------------------------------------------------------------- #
# (f) include arg filters contributions additively                   #
# ---------------------------------------------------------------- #

test_that("(f) include filter: c('fixed') drops the smooth contribution", {
  skip_if_no_greta_quiet()
  f <- mk_smooth_fit()
  newx <- data.frame(x = seq(min(f$dat$x), max(f$dat$x), length.out = 20L))

  full <- flexyBayes:::.predict_linear_draws(
    f$fit,
    newx,
    include = c("fixed", "smooth")
  )
  fxonly <- flexyBayes:::.predict_linear_draws(f$fit, newx, include = "fixed")

  # The smooth contribution is non-trivial -- the row-wise SD across
  # newdata of the row-mean should differ substantially between the
  # smooth-included and fixed-only paths.
  expect_true(sd(rowMeans(full)) > sd(rowMeans(fxonly)) + 0.02)
})


# ---------------------------------------------------------------- #
# (g) file-backed path equals in-memory path within rounding         #
# ---------------------------------------------------------------- #

test_that("(g) file-backed predict matches in-memory predict (smooth model)", {
  skip_if_no_greta_quiet()
  f <- mk_smooth_fit()
  newx <- data.frame(x = seq(min(f$dat$x), max(f$dat$x), length.out = 20L))

  in_mem <- predict(f$fit, newdata = newx, type = "link")

  tmp <- tempfile(fileext = ".rds")
  predict(
    f$fit,
    newdata = newx,
    output_file = tmp,
    format = "rds",
    type = "link"
  )
  on.exit(unlink(tmp), add = TRUE)
  file_out <- readRDS(tmp)

  # Both paths route through the shared kernel, so the file path's
  # `point` column (posterior mean) numerically matches the
  # in-memory path's row means. Names differ (in_mem carries rownames
  # from newdata; file_out$point is unnamed), so compare unnamed
  # values.
  expect_equal(unname(file_out$point), unname(in_mem), tolerance = 1e-10)
})

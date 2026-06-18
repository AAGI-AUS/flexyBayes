# test-predict-file-output-equivalence.R -- v0.4.0 section-11 trust
# gate. The file-backed prediction path (ADR 0023) must reproduce the
# in-memory prediction kernel, not just for the fixed-only path. Two
# equivalences are checked per term class:
#
#   (1) the file payload's `point` column equals the in-memory
#       predict() point vector (the file path computes the same point
#       prediction as the in-memory kernel); and
#   (2) the chunked file path is byte-equivalent (to floating-point
#       tolerance; BLAS non-associativity at the chunk boundary admits
#       ULP-level slack) to the unchunked single-pass file path, so
#       chunking -- the file-path-specific machinery -- does not perturb
#       the posterior intervals.
#
# Four term classes (ADR 0023 + v040-plan section 11):
#   (i)   fixed-only                  y ~ x
#   (ii)  smooth                      y ~ s(x)
#   (iii) random-intercept (pop.)     y ~ x + (1 | g)
#   (iv)  sampled new level           y ~ x + (1 | g), allow_new_levels
#
# greta-backed (the per-draw file path requires $greta$draws); the file
# runs only when greta is available.

suppressPackageStartupMessages(library(testthat))

# Silence the fit-time notes, scoped to the CALLING test_that frame so
# the option never leaks into other test files. (A session-wide
# teardown_env() scope would silence the default-prior note that
# test-smooth.R asserts fires.)
.fbe_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = envir
  )
}

# (1) in-memory point == file point; (2) chunked file == unchunked file.
.expect_file_equiv <- function(fit, newdata, ...) {
  mem <- predict(fit, newdata = newdata, ...) # point vector

  f_whole <- tempfile(fileext = ".rds")
  f_chunk <- tempfile(fileext = ".rds")
  on.exit(unlink(c(f_whole, f_chunk)), add = TRUE)

  # Unchunked single pass (the in-memory kernel computed in one go).
  predict(fit, newdata = newdata, output_file = f_whole, format = "rds", ...)
  # Chunked file path (chunk_size below nrow forces multi-chunk).
  predict(
    fit,
    newdata = newdata,
    output_file = f_chunk,
    format = "rds",
    chunk_size = max(1L, nrow(newdata) %/% 3L),
    ...
  )

  whole <- readRDS(f_whole)
  chunk <- readRDS(f_chunk)

  expect_true(all(c("point", "lower", "upper") %in% names(whole)))
  expect_equal(nrow(whole), nrow(newdata))

  # (1) file point matches the in-memory predict() point vector.
  expect_equal(whole$point, as.numeric(mem))
  # (2) chunk-invariance of the full intervals.
  expect_equal(chunk$point, whole$point)
  expect_equal(chunk$lower, whole$lower)
  expect_equal(chunk$upper, whole$upper)
}

# ---------------------------------------------------------------- #
# (i) Fixed-only                                                    #
# ---------------------------------------------------------------- #

test_that("file output matches the in-memory kernel for a fixed-only model", {
  skip_if_greta_backend_unusable()
  .fbe_quiet()
  set.seed(101L)
  dat <- data.frame(x = rnorm(60L), y = rnorm(60L))
  fit <- suppressMessages(fb(
    y ~ x,
    data = dat,
    backend = "greta",
    n_samples = 50L,
    warmup = 30L,
    chains = 1L,
    prior_fixed_sd = 10,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  .expect_file_equiv(fit, dat[1:21, , drop = FALSE])
})

# ---------------------------------------------------------------- #
# (ii) Smooth s(x)                                                  #
# ---------------------------------------------------------------- #

test_that("file output matches the in-memory kernel for a smooth s(x) model", {
  skip_if_greta_backend_unusable()
  .fbe_quiet()
  skip_if_not_installed("mgcv")
  set.seed(102L)
  dat <- data.frame(x = sort(rnorm(60L)))
  dat$y <- sin(dat$x) + rnorm(60L, sd = 0.3)
  # Smooths are ingested through the asreml surface (s() in `random`),
  # not the brms surface -- fb_brms() does not yet support s().
  fit <- suppressMessages(flexybayes(
    y ~ 1,
    random = ~ s(x),
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
  .expect_file_equiv(fit, dat[1:21, , drop = FALSE])
})

# ---------------------------------------------------------------- #
# (iii) Random-intercept, population-level prediction               #
# ---------------------------------------------------------------- #

test_that("file output matches the in-memory kernel for a population-level RI prediction", {
  skip_if_greta_backend_unusable()
  .fbe_quiet()
  set.seed(103L)
  dat <- data.frame(
    x = rnorm(80L),
    g = factor(sample(letters[1:5], 80L, replace = TRUE)),
    y = rnorm(80L)
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
  .expect_file_equiv(
    fit,
    dat[1:21, , drop = FALSE],
    allow_new_levels = "population"
  )
})

# ---------------------------------------------------------------- #
# (iv) Sampled new level (allow_new_levels = "sample")              #
# ---------------------------------------------------------------- #
#
# Sample mode draws the unknown level's effect per draw from the RNG, so
# the file path's chunked iteration re-samples per chunk (ADR 0023 §6.d:
# chunked sample is supported only via output_file and is not byte-equal
# to the unchunked pass). The honest claim is reproducibility under a
# fixed seed: the unchunked file output is identical across runs, and
# the intervals are well-formed (finite, ordered).

test_that("the sampled-new-level file path is reproducible under a fixed seed and well-formed", {
  skip_if_greta_backend_unusable()
  .fbe_quiet()
  set.seed(104L)
  dat <- data.frame(
    x = rnorm(80L),
    g = factor(sample(letters[1:5], 80L, replace = TRUE)),
    y = rnorm(80L)
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
  nd <- data.frame(
    x = rnorm(15L),
    g = factor(rep("z", 15L), levels = c(letters[1:5], "z"))
  )

  f1 <- tempfile(fileext = ".rds")
  f2 <- tempfile(fileext = ".rds")
  on.exit(unlink(c(f1, f2)), add = TRUE)
  set.seed(99L)
  predict(
    fit,
    newdata = nd,
    allow_new_levels = "sample",
    output_file = f1,
    format = "rds"
  )
  set.seed(99L)
  predict(
    fit,
    newdata = nd,
    allow_new_levels = "sample",
    output_file = f2,
    format = "rds"
  )
  o1 <- readRDS(f1)
  o2 <- readRDS(f2)

  expect_equal(o1$point, o2$point) # reproducible under the seed
  expect_equal(o1$lower, o2$lower)
  expect_equal(o1$upper, o2$upper)
  expect_true(all(is.finite(o1$point)))
  expect_true(all(o1$lower <= o1$point & o1$point <= o1$upper))
})

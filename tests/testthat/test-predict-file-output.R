# test-predict-file-output.R -- Stage 3B file-output (ADR 0023
# Decision 4). Covers .resolve_format() pure-path table, the
# fst-not-installed structured refusal, csv / rds / fst round-trip
# with factor-type preservation, chunked vs unchunked equivalence
# at floating-point precision, format = "auto" end-to-end, and the
# interop branch. fst availability is flipped via
# testthat::local_mocked_bindings() to keep the matrix coverage
# independent of the host's fst install state.

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


# Test fixture: small greta-fit (identical shape to
# test-predict-chunked.R fixture). Greta-backend so $greta$draws
# is populated, which the file-output per-draw path requires.
mk_predict_fit <- function(seed = 2026L, N = 60L, J = 4L) {
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


# ---------------------------------------------------------------- #
# ADR §6.a: .resolve_format() pure-path table                       #
# ---------------------------------------------------------------- #

test_that(".resolve_format(): format='auto' + N<1e6 returns 'rds'", {
  expect_identical(
    flexyBayes:::.resolve_format("auto", 1e5, FALSE, TRUE),
    "rds"
  )
  expect_identical(
    flexyBayes:::.resolve_format("auto", 1e5, FALSE, FALSE),
    "rds"
  )
})

test_that(".resolve_format(): format='auto' + N>=1e6 + fst returns 'fst'", {
  expect_identical(
    flexyBayes:::.resolve_format("auto", 1e6, FALSE, TRUE),
    "fst"
  )
  expect_identical(
    flexyBayes:::.resolve_format("auto", 1e7, FALSE, TRUE),
    "fst"
  )
})

test_that(".resolve_format(): format='auto' + N>=1e6 + !fst returns 'rds'", {
  expect_identical(
    flexyBayes:::.resolve_format("auto", 1e6, FALSE, FALSE),
    "rds"
  )
  expect_identical(
    flexyBayes:::.resolve_format("auto", 1e7, FALSE, FALSE),
    "rds"
  )
})

test_that(".resolve_format(): format='auto' + interop returns 'csv'", {
  expect_identical(
    flexyBayes:::.resolve_format("auto", 1e6, TRUE, TRUE),
    "csv"
  )
  expect_identical(
    flexyBayes:::.resolve_format("auto", 1e7, TRUE, TRUE),
    "csv"
  )
})

test_that(".resolve_format(): explicit format pins identity for csv/rds", {
  expect_identical(
    flexyBayes:::.resolve_format("csv", 1e7, FALSE, FALSE),
    "csv"
  )
  expect_identical(
    flexyBayes:::.resolve_format("rds", 1e7, FALSE, FALSE),
    "rds"
  )
})


# ---------------------------------------------------------------- #
# ADR §6.c: structured refusal on format='fst' + !fst               #
# ---------------------------------------------------------------- #

test_that(".resolve_format(): format='fst' + fst_available=FALSE raises structured refusal", {
  err <- expect_error(
    flexyBayes:::.resolve_format("fst", 1e6, FALSE, FALSE),
    "'fst' package is not installed"
  )
  msg <- conditionMessage(err)
  # Contains the install command verbatim per ADR 0023 §4 message.
  expect_match(msg, "install\\.packages\\(\"fst\"\\)")
  # Contains the rds override hint per ADR 0023 §4 message.
  expect_match(msg, "format = \"rds\"")
})

test_that(".resolve_format(): format='fst' + fst_available=TRUE returns 'fst'", {
  expect_identical(
    flexyBayes:::.resolve_format("fst", 1e6, FALSE, TRUE),
    "fst"
  )
})


# ---------------------------------------------------------------- #
# ADR §6.b: csv / rds / fst round-trip with factor preservation     #
# ---------------------------------------------------------------- #

test_that("rds round-trip: predict file output preserves factor type", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  res <- predict(
    fx$fit,
    newdata = fx$dat[1:20, ],
    output_file = f,
    format = "rds"
  )
  expect_identical(res, f) # invisible(path)
  out <- readRDS(f)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("point", "lower", "upper") %in% names(out)))
  expect_identical(nrow(out), 20L)
  expect_true(is.factor(out$g)) # rds preserves factor type
})

test_that("csv round-trip: predict file output reads factor as character", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  predict(fx$fit, newdata = fx$dat[1:20, ], output_file = f, format = "csv")
  out <- data.table::fread(f)
  expect_true(all(c("point", "lower", "upper") %in% names(out)))
  expect_identical(nrow(out), 20L)
  # csv loses factor type by design -- column reads back as character.
  expect_true(is.character(out$g))
})

test_that("fst round-trip: predict file output preserves factor type (when fst installed)", {
  skip_if_greta_backend_unusable()
  testthat::skip_if_not_installed("fst")
  fx <- mk_predict_fit()
  f <- tempfile(fileext = ".fst")
  on.exit(unlink(f), add = TRUE)
  predict(fx$fit, newdata = fx$dat[1:20, ], output_file = f, format = "fst")
  out <- fst::read_fst(f, as.data.table = TRUE)
  expect_true(all(c("point", "lower", "upper") %in% names(out)))
  expect_identical(nrow(out), 20L)
  expect_true(is.factor(out$g))
})


# ---------------------------------------------------------------- #
# ADR §6.d: chunk_size invariance on file output                    #
# ---------------------------------------------------------------- #

test_that("chunk_size invariance: chunked vs unchunked file output equivalent at fp precision", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd_big <- fx$dat[rep(seq_len(nrow(fx$dat)), 2L), ]
  f1 <- tempfile(fileext = ".rds")
  f2 <- tempfile(fileext = ".rds")
  on.exit(unlink(c(f1, f2)), add = TRUE)
  predict(fx$fit, newdata = nd_big, output_file = f1, format = "rds")
  predict(
    fx$fit,
    newdata = nd_big,
    output_file = f2,
    format = "rds",
    chunk_size = 25L
  )
  o1 <- readRDS(f1)
  o2 <- readRDS(f2)
  # ADR 0023 §6.d says "bitwise-identical" -- in practice BLAS
  # non-associativity at the matrix-size boundary introduces ~ULP-
  # level differences. The test commits to floating-point equivalence
  # at testthat::expect_equal default tolerance (~1.5e-8).
  expect_equal(o1$point, o2$point)
  expect_equal(o1$lower, o2$lower)
  expect_equal(o1$upper, o2$upper)
})


# ---------------------------------------------------------------- #
# ADR §6.a (end-to-end): format='auto' resolution                   #
# ---------------------------------------------------------------- #

test_that("format='auto' end-to-end: small N writes .rds payload (not fst)", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  f <- tempfile(fileext = ".out")
  on.exit(unlink(f), add = TRUE)
  predict(fx$fit, newdata = fx$dat[1:20, ], output_file = f, format = "auto")
  # The file is rds-formatted (small N -> rds branch); readRDS must
  # succeed.
  out <- readRDS(f)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("point", "lower", "upper") %in% names(out)))
})

test_that("format='auto' + mocked !fst at large-N reports 'rds'", {
  # Pure-path test: mock .fst_available to FALSE; .resolve_format
  # at N = 1e6 must fall through to rds.
  testthat::local_mocked_bindings(
    .fst_available = function() FALSE,
    .package = "flexyBayes"
  )
  expect_identical(
    flexyBayes:::.resolve_format(
      "auto",
      1e6,
      FALSE,
      flexyBayes:::.fst_available()
    ),
    "rds"
  )
})

test_that("format='auto' + mocked fst at large-N reports 'fst'", {
  testthat::local_mocked_bindings(
    .fst_available = function() TRUE,
    .package = "flexyBayes"
  )
  expect_identical(
    flexyBayes:::.resolve_format(
      "auto",
      1e6,
      FALSE,
      flexyBayes:::.fst_available()
    ),
    "fst"
  )
})


# ---------------------------------------------------------------- #
# interop branch                                                    #
# ---------------------------------------------------------------- #

test_that("interop=TRUE + format='auto' writes csv (predict body level)", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  predict(
    fx$fit,
    newdata = fx$dat[1:20, ],
    output_file = f,
    format = "auto",
    interop = TRUE
  )
  # Must be readable as csv.
  out <- data.table::fread(f)
  expect_true(all(c("point", "lower", "upper") %in% names(out)))
  expect_true(is.character(out$g)) # csv -> character readback
})

test_that("interop=TRUE + explicit format='rds' wins (auto-only rule)", {
  # interop only flips the auto branch -- explicit format pins
  # identity.
  expect_identical(
    flexyBayes:::.resolve_format("rds", 1e6, TRUE, TRUE),
    "rds"
  )
  expect_identical(
    flexyBayes:::.resolve_format("fst", 1e6, TRUE, TRUE),
    "fst"
  )
})


# ---------------------------------------------------------------- #
# Overwrite refusal + signature edges                               #
# ---------------------------------------------------------------- #

test_that("overwrite refusal: existing output_file raises structured stop", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  predict(fx$fit, newdata = fx$dat[1:5, ], output_file = f, format = "rds")
  err <- expect_error(
    predict(fx$fit, newdata = fx$dat[1:5, ], output_file = f, format = "rds"),
    "already exists"
  )
  expect_match(conditionMessage(err), "Refusing to overwrite silently")
  expect_match(conditionMessage(err), "unlink")
})

test_that("output_file + newdata=NULL raises structured refusal", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  err <- expect_error(
    predict(fx$fit, newdata = NULL, output_file = f),
    "`output_file` is only supported with `newdata`"
  )
  expect_false(file.exists(f)) # no file written on refusal
})

test_that("output_file with se.fit=TRUE writes a 4th column", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  predict(
    fx$fit,
    newdata = fx$dat[1:10, ],
    output_file = f,
    format = "rds",
    se.fit = TRUE
  )
  out <- readRDS(f)
  expect_true(all(
    c("point", "lower", "upper", "se.fit") %in%
      names(out)
  ))
  expect_true(all(is.finite(out$se.fit)))
  expect_true(all(out$se.fit >= 0))
})

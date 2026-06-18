# test-predict-chunked.R -- Stage 3B core (ADR 0023 Decisions 1, 2, 3).
#
# Covers the trimmed v0.3.4 scope: dictionary-backed factor handling
# (Decision 1), `allow_new_levels` policy (Decision 2 with "population"
# + "refuse" active and "sample" reserved), and chunked iteration
# (Decision 3). File-backed output formats + fst Suggests (Decision 4)
# are deferred to v0.3.5; the corresponding test blocks land at that
# release.

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


# Test fixture: small greta-fit on a 4-level grouping factor + one
# continuous covariate. Greta is intentional -- the dispatch attaches
# extras$fb_dataset on every emit path; greta keeps the test wall-time
# under 30 s on the dev laptop and exercises the full class chain.
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
# Decision 1: dictionary-backed factor handling                     #
# ---------------------------------------------------------------- #

test_that("Decision 1: fit carries extras$fb_dataset with dictionaries", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  expect_false(is.null(fx$fit$extras$fb_dataset))
  expect_s3_class(fx$fit$extras$fb_dataset, "fb_dataset")
  expect_true("g" %in% names(fx$fit$extras$fb_dataset$dictionaries))
  expect_setequal(fx$fit$extras$fb_dataset$dictionaries$g, letters[1:4])
  expect_identical(fx$fit$extras$fb_dataset$n_rows, 60L)
})

test_that("Decision 1: extras$fb_dataset is metadata-only (no $data)", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  expect_identical(fx$fit$extras$fb_dataset$origin, "metadata-only")
  expect_null(fx$fit$extras$fb_dataset$data)
})

test_that("Decision 1: known-level newdata round-trips factor codes", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:10, ]
  # Pass character (not factor) -- resolver should coerce to factor
  # with the fit-time level set.
  nd$g <- as.character(nd$g)
  p <- predict(fx$fit, newdata = nd)
  expect_length(p, nrow(nd))
  expect_true(all(is.finite(p)))
})


# ---------------------------------------------------------------- #
# Decision 2: allow_new_levels policy                               #
# ---------------------------------------------------------------- #

test_that("Decision 2: allow_new_levels match.arg rejects unknown values", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:5, ]
  expect_error(
    predict(fx$fit, newdata = nd, allow_new_levels = "ignore"),
    "should be one of"
  )
})

test_that("Decision 2: allow_new_levels='population' default succeeds on known levels", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:5, ]
  # Default allow_new_levels -- should NOT warn when no unknown levels
  p <- expect_silent(predict(fx$fit, newdata = nd))
  expect_length(p, 5L)
})

test_that("Decision 2: allow_new_levels='population' warns + sets unknown rows to NA factor", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:5, ]
  nd$g <- as.character(nd$g)
  nd$g[2] <- "ZZ"
  nd$g[4] <- "WW"
  expect_warning(
    p <- predict(fx$fit, newdata = nd, allow_new_levels = "population"),
    "unknown factor level"
  )
  expect_length(p, nrow(nd))
})

test_that("Decision 2: allow_new_levels='refuse' raises structured stop on unknown level", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:5, ]
  nd$g <- as.character(nd$g)
  nd$g[2] <- "ZZ"
  err <- expect_error(
    predict(fx$fit, newdata = nd, allow_new_levels = "refuse"),
    "level\\(s\\) not present in the fit-time dictionary"
  )
  expect_match(conditionMessage(err), "column `g`")
  expect_match(conditionMessage(err), "ZZ")
})

test_that("Decision 2: allow_new_levels='refuse' lists multiple unknown levels", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:6, ]
  nd$g <- as.character(nd$g)
  nd$g[1] <- "ZZ"
  nd$g[3] <- "YY"
  nd$g[5] <- "XX"
  err <- expect_error(
    predict(fx$fit, newdata = nd, allow_new_levels = "refuse"),
    "3 unknown level\\(s\\)"
  )
  expect_match(conditionMessage(err), "XX")
  expect_match(conditionMessage(err), "YY")
  expect_match(conditionMessage(err), "ZZ")
})

test_that("Decision 2: allow_new_levels='population' is the default", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:5, ]
  p_default <- predict(fx$fit, newdata = nd)
  p_pop <- predict(fx$fit, newdata = nd, allow_new_levels = "population")
  expect_identical(p_default, p_pop)
})


# ---------------------------------------------------------------- #
# Decision 3: chunked iteration                                     #
# ---------------------------------------------------------------- #

test_that("Decision 3: chunk_size=NULL preserves single-pass behaviour", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:20, ]
  p_default <- predict(fx$fit, newdata = nd)
  p_null <- predict(fx$fit, newdata = nd, chunk_size = NULL)
  expect_identical(p_default, p_null)
})

test_that("Decision 3: chunk_size > nrow(newdata) bypasses chunked path", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd <- fx$dat[1:10, ]
  p_unchunked <- predict(fx$fit, newdata = nd)
  p_above <- predict(fx$fit, newdata = nd, chunk_size = 1000L)
  expect_identical(p_unchunked, p_above)
})

test_that("Decision 3: chunked vs unchunked produce identical output", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  # Replicate the data to get 180 rows, well above the chunk_size
  nd_big <- fx$dat[rep(seq_len(nrow(fx$dat)), 3), ]
  p_full <- predict(fx$fit, newdata = nd_big)
  p_chunk <- predict(fx$fit, newdata = nd_big, chunk_size = 50L)
  expect_equal(p_full, p_chunk)
})

test_that("Decision 3: chunk_size + se.fit=TRUE assembles fit + se.fit", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  nd_big <- fx$dat[rep(seq_len(nrow(fx$dat)), 3), ]
  out_full <- predict(fx$fit, newdata = nd_big, se.fit = TRUE)
  out_chunk <- predict(
    fx$fit,
    newdata = nd_big,
    se.fit = TRUE,
    chunk_size = 50L
  )
  expect_named(out_chunk, c("fit", "se.fit"))
  expect_length(out_chunk$fit, nrow(nd_big))
  expect_length(out_chunk$se.fit, nrow(nd_big))
  # Single-pass se.fit carries rowSums() names; chunked path
  # concatenates per-chunk vectors without names. Values are equal;
  # strip names for the comparison.
  expect_equal(unname(out_full$fit), unname(out_chunk$fit))
  expect_equal(unname(out_full$se.fit), unname(out_chunk$se.fit))
})


# ---------------------------------------------------------------- #
# Backward-compat: legacy fits without extras$fb_dataset             #
# ---------------------------------------------------------------- #

test_that("Backward-compat: legacy fit (no extras$fb_dataset) skips dictionary resolution", {
  skip_if_greta_backend_unusable()
  fx <- mk_predict_fit()
  # Simulate a legacy fit by stripping the new slot
  legacy_fit <- fx$fit
  legacy_fit$extras$fb_dataset <- NULL
  nd <- fx$dat[1:5, ]
  # No warning even with an unknown level -- the dictionary path is
  # skipped entirely for legacy fits.
  nd$g <- as.character(nd$g)
  nd$g[1] <- "ZZ"
  expect_warning(
    p <- predict(legacy_fit, newdata = nd, allow_new_levels = "population"),
    NA # asserts NO warning fires
  )
  # The legacy path produces whatever the v0.3.3 behaviour produced;
  # the test commits only to "no crash, no spurious dictionary
  # warning". The actual numeric output is the v0.3.3 reference.
  expect_length(p, nrow(nd))
})

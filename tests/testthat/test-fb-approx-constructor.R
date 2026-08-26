# Tests for fb_approx() -- the approximation-scheme constructor noun
# (v0.4.0 Wave 1 Phase 1D).

suppressPackageStartupMessages(library(testthat))


# ---------------------------------------------------------------- #
# Construction + shape                                              #
# ---------------------------------------------------------------- #

test_that("fb_approx() returns a classed list with scheme + kwargs", {
  a <- fb_approx("low_rank_smooth", rank = 5L)
  expect_s3_class(a, "fb_approx")
  expect_true(is_fb_approx(a))
  expect_identical(class(a), c("fb_approx", "list"))
  # scheme + tuning are plain list elements so the parse path reads
  # spec$scheme / spec$rank directly.
  expect_identical(a$scheme, "low_rank_smooth")
  expect_identical(a$rank, 5L)
})

test_that("fb_approx() carries a literature-referenced bias-bound promise", {
  a <- fb_approx("low_rank_smooth", rank = 5L)
  bb <- attr(a, "bias_bound_promise")
  expect_true(is.character(bb) && nzchar(bb))
  expect_match(bb, "Frobenius")
  expect_match(bb, "Wood")
  # Provenance (internal design records) must not leak onto the object.
  expect_false(any(grepl("ADR", c(bb, unlist(a)))))
})


# ---------------------------------------------------------------- #
# Validation                                                        #
# ---------------------------------------------------------------- #

test_that("fb_approx() refuses an unregistered scheme with a catchable class", {
  err <- tryCatch(fb_approx("not_a_scheme"), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_approximation_scheme_unknown")
  expect_match(conditionMessage(err), "low_rank_smooth")
})

test_that("fb_approx() rejects a non-string scheme", {
  expect_error(fb_approx(1L), "non-empty single string")
  expect_error(fb_approx(character(0)), "non-empty single string")
})


# ---------------------------------------------------------------- #
# print                                                             #
# ---------------------------------------------------------------- #

test_that("print.fb_approx() shows scheme, kwargs, and bias bound", {
  out <- utils::capture.output(print(fb_approx("low_rank_smooth", rank = 5L)))
  expect_true(any(grepl("<fb_approx> scheme = \"low_rank_smooth\"", out)))
  expect_true(any(grepl("rank = 5", out)))
  expect_true(any(grepl("bias bound", out)))
})


# ---------------------------------------------------------------- #
# Consumption: s(x, representation = fb_approx(...))                 #
# ---------------------------------------------------------------- #

test_that("fb_approx() is accepted as a smooth representation spec", {
  skip_if_not_installed("mgcv")
  skip_on_cran()
  skip_on_ci()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(1L)
  d <- data.frame(x = runif(80))
  d$y <- sin(2 * pi * d$x) + rnorm(80, 0, 0.3)
  # plan = TRUE parses + resolves the approximation spec without fitting
  # (and without needing an active engine -- the low_rank_smooth scheme
  # has no consumer on either active engine, so this only checks that
  # the fb_approx() form and the equivalent list form both build a
  # plan; neither backend is asked to fit it). Smooths live in `random`,
  # not `fixed` -- flexyBayes refuses a smoother written in the fixed
  # part with fixed_smoother_not_supported.
  p_obj <- suppressMessages(flexybayes(
    fixed = y ~ 1,
    random = ~ s(x, representation = fb_approx("low_rank_smooth", rank = 5L)),
    data = d,
    plan = TRUE,
    verbose = FALSE
  ))
  p_lst <- suppressMessages(flexybayes(
    fixed = y ~ 1,
    random = ~ s(
      x,
      representation = list(scheme = "low_rank_smooth", rank = 5L)
    ),
    data = d,
    plan = TRUE,
    verbose = FALSE
  ))
  expect_s3_class(p_obj, "fb_plan")
  expect_s3_class(p_lst, "fb_plan")
})

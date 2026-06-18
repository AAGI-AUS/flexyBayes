# Tests for validate_approximation() --- ADR 0030 C5 + ADR 0027
# (v0.4.0 Wave 1 Phase 1B). The sixth exported verb. Dispatch is on
# the fit's registered approximation scheme via the .approximation_
# registry dispatch table; an exact fit refuses. The validation-result
# object is <fb_approximation_validation>. These tests use synthetic
# fits carrying the emit-time fit$extras$parse_info$approx slot so the
# verb is exercised without a live MCMC fit (the end-to-end fit is
# covered, stress-gated, in test-emit-smooth-low-rank.R).

# .mk_lr_fit() --- a synthetic flexybayes fit carrying one or more
# low_rank_smooth approximations with the given captures.
.mk_lr_fit <- function(captures = c(x = 0.9891)) {
  approx <- lapply(seq_along(captures), function(i) {
    list(
      scheme = "low_rank_smooth",
      rank = 4L,
      k = 9L,
      frobenius_capture = unname(captures[[i]]),
      V_K = matrix(0, 9L, 4L),
      singular_values = rep(1, 9L)
    )
  })
  names(approx) <- names(captures)
  structure(
    list(
      exactness = "approximate_low_rank_smooth",
      extras = list(parse_info = list(approx = approx))
    ),
    class = c("flexybayes", "list")
  )
}

# ---------------------------------------------------------------- #
# (a) exported generic + result object shape                        #
# ---------------------------------------------------------------- #

test_that("validate_approximation() returns a structured result", {
  expect_true(is.function(validate_approximation))
  v <- validate_approximation(.mk_lr_fit(c(x = 0.9995)))
  expect_s3_class(v, "fb_approximation_validation")
  expect_identical(v$scheme, "low_rank_smooth")
  expect_true(v$pass)
  expect_identical(v$threshold, 0.99)
  expect_named(v$per_smooth, "x")
  r <- v$per_smooth$x
  expect_setequal(
    names(r),
    c(
      "smooth",
      "scheme",
      "rank",
      "k",
      "frobenius_capture",
      "bias_bound",
      "threshold",
      "pass"
    )
  )
  expect_equal(r$bias_bound, 1 - r$frobenius_capture)
})

# ---------------------------------------------------------------- #
# (b) threshold behaviour                                           #
# ---------------------------------------------------------------- #

test_that("the pass verdict tracks the threshold", {
  fit <- .mk_lr_fit(c(x = 0.9891))
  expect_false(validate_approximation(fit)$pass) # < 0.99
  expect_false(validate_approximation(fit, threshold = 0.999)$pass)
  expect_true(validate_approximation(fit, threshold = 0.90)$pass)
})

test_that("a fit passes only when every smooth clears the threshold", {
  fit <- .mk_lr_fit(c(x = 0.999, z = 0.95))
  v <- validate_approximation(fit, threshold = 0.99)
  expect_false(v$pass) # z fails
  expect_true(v$per_smooth$x$pass)
  expect_false(v$per_smooth$z$pass)
  expect_true(validate_approximation(fit, threshold = 0.90)$pass)
})

# ---------------------------------------------------------------- #
# (c) refusal on an exact fit                                       #
# ---------------------------------------------------------------- #

test_that("validate_approximation() refuses an exact fit", {
  exact <- structure(
    list(exactness = "exact", extras = list(parse_info = list())),
    class = c("flexybayes", "list")
  )
  err <- tryCatch(
    validate_approximation(exact),
    flexybayes_approximation_scheme_unknown = function(e) e
  )
  expect_s3_class(err, "flexybayes_approximation_scheme_unknown")
  expect_identical(err$reason_code, "approximation_scheme_unknown")
  expect_match(conditionMessage(err), "no recognised approximation")
})

# ---------------------------------------------------------------- #
# (d) scheme resolution helper                                      #
# ---------------------------------------------------------------- #

test_that(".fit_approximation_scheme() resolves the scheme", {
  expect_identical(
    flexyBayes:::.fit_approximation_scheme(.mk_lr_fit()),
    "low_rank_smooth"
  )
  exact <- structure(
    list(extras = list(parse_info = list())),
    class = c("flexybayes", "list")
  )
  expect_null(flexyBayes:::.fit_approximation_scheme(exact))
})

# ---------------------------------------------------------------- #
# (e) print method snapshot                                         #
# ---------------------------------------------------------------- #

test_that("print.fb_approximation_validation snapshot", {
  v <- validate_approximation(.mk_lr_fit(c(x = 0.95)))
  expect_snapshot(print(v))
})

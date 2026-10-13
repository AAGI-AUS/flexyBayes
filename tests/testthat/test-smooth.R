# mgcv-style univariate smooth ingest s(x) (ADR 0004 D3).
#
# s() is parsed, and te(), ti() and t2() are not. The basis is built at parse
# time with mgcv::smoothCon() into a `smooth_mgcv` random term. Neither
# engine fits that term (INLA's smooth route is spl(), a native rw2), so
# `flexybayes(random = ~ s(x))` refuses on both backends before any code
# is emitted.

skip_if_no_mgcv <- function() skip_if_not_installed("mgcv")

mk_smooth_data <- function() {
  set.seed(2026L)
  n <- 50L
  x <- sort(runif(n, 0, 10))
  y <- sin(x) + rnorm(n, sd = 0.2)
  data.frame(x = x, y = y, g = factor(rep(1:5, length.out = n)))
}

# ---------------------------------------------------------------- #
# Parse-time basis construction                                    #
# ---------------------------------------------------------------- #

test_that(".enrich() builds the smoothCon basis for s(x)", {
  skip_if_no_mgcv()
  d <- mk_smooth_data()
  fb <- flexyBayes:::fb_from_asreml(
    fixed = y ~ 1,
    random = ~ s(x),
    data = d
  )
  expect_length(fb$random_terms, 1L)
  rt <- fb$random_terms[[1]]
  expect_identical(rt$type, "smooth_mgcv")
  expect_identical(rt$var, "x")
  expect_true(is.matrix(rt$X))
  expect_identical(nrow(rt$X), nrow(d))
  expect_true(rt$k >= 1L)
  expect_true(is.character(rt$smooth_label))
})

test_that("s(x, k = 6) honours user-supplied k", {
  skip_if_no_mgcv()
  d <- mk_smooth_data()
  fb <- flexyBayes:::fb_from_asreml(
    fixed = y ~ 1,
    random = ~ s(x, k = 6),
    data = d
  )
  rt <- fb$random_terms[[1]]
  expect_identical(rt$type, "smooth_mgcv")
  # absorb.cons = TRUE drops one column to fold the sum-to-zero
  # constraint, so k = 6 gives 5 effective columns.
  expect_identical(rt$k, 5L)
})

test_that(".enrich() errors politely when s() variable is missing from data", {
  skip_if_no_mgcv()
  d <- mk_smooth_data()
  expect_error(
    flexyBayes:::fb_from_asreml(
      fixed = y ~ 1,
      random = ~ s(missing_var),
      data = d
    ),
    regexp = "not found in data"
  )
})

test_that(".brms_term_refusal_message() distinguishes spl() (INLA route) from s() (no active route)", {
  # brms lowers neither smoother. spl() reaches INLA as a second-order
  # random walk on an ASReml-style basis. s() (smooth_mgcv) reaches
  # neither engine, so its message must not offer INLA.
  spl_msg <- flexyBayes:::.brms_term_refusal_message(
    list(type = "spline", var = "x")
  )
  expect_match(spl_msg, "spl(x)", fixed = TRUE)
  expect_match(spl_msg, "second-order random walk", fixed = TRUE)
  expect_match(spl_msg, 'backend = "inla"', fixed = TRUE)

  smooth_msg <- flexyBayes:::.brms_term_refusal_message(
    list(type = "smooth_mgcv", var = "x")
  )
  # Recommends spl(x), not s(x), on INLA.
  expect_match(smooth_msg, "spl(x)", fixed = TRUE)
  expect_false(grepl('s(x) with backend = "inla"', smooth_msg, fixed = TRUE))
  # brms fits s() in its own formula (test-engine-claims.R), so the message
  # says flexyBayes fits it on neither engine, not that no engine does.
  expect_match(smooth_msg, "flexyBayes fits an s(x) smooth on neither engine",
    fixed = TRUE)
  expect_false(grepl("no active engine", smooth_msg, fixed = TRUE))
})

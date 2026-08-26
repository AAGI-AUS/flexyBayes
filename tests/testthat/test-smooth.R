# Tests for mgcv-style univariate smooth ingest s(x) (ADR 0004 D3).
#
# v0.1 supports s() only; te() / ti() / t2() defer to v0.2. The basis is
# constructed at parse time via mgcv::smoothCon() into a `smooth_mgcv`
# random-term type, tested below. 0.9.3: the only engine that ever emitted
# this term type was withdrawn entirely (see NEWS.md) -- neither active
# engine represents it (brms has no lowering for it; INLA's own smooth
# route uses a native rw2 term instead, a different random-term type), so
# the code-generation tests this file used to carry are gone. A live
# `flexybayes(random = ~ s(x))` call now refuses on both `backend = "brms"`
# and `backend = "inla"` before any code is emitted.

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
  # constraint; k = 6 -> 5 effective columns.
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
  # 0.9.3: brms has no lowering for either smoother, but the two
  # messages must not claim the same alternative. spl() genuinely
  # reaches INLA as a second-order random walk; s() (smooth_mgcv) does
  # not -- INLA's rw2 route is spl(), a different term type built on
  # an ASReml-style basis, not a lowering of the mgcv smooth basis
  # s() constructs. Before this fix the smooth_mgcv message wrongly
  # claimed s(x) itself would run on INLA as a second-order random
  # walk, which it does not: neither active engine fits it.
  spl_msg <- flexyBayes:::.brms_term_refusal_message(
    list(type = "spline", var = "x")
  )
  expect_match(spl_msg, "spl(x)", fixed = TRUE)
  expect_match(spl_msg, "second-order random walk", fixed = TRUE)
  expect_match(spl_msg, 'backend = "inla"', fixed = TRUE)

  smooth_msg <- flexyBayes:::.brms_term_refusal_message(
    list(type = "smooth_mgcv", var = "x")
  )
  # Recommends spl(x), not s(x), on INLA -- and does not repeat the
  # false claim that s(x) itself is what runs as a random walk there.
  expect_match(smooth_msg, "spl(x)", fixed = TRUE)
  expect_false(grepl('s(x) with backend = "inla"', smooth_msg, fixed = TRUE))
  expect_match(smooth_msg, "no active engine", fixed = TRUE)
})

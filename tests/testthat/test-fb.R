# Tests for fb() -- alias for flexybayes() under ADR 0004 D1.
#
# fb is a literal alias for flexybayes(); see R/fb.R. The original
# brms-format ingest path is deferred to v0.2 as fb_brms(); the
# internal helper fb_from_brms() in R/fb_from_brms.R remains
# unexported and continues to be tested in test-fb-from-brms.R for
# v0.2 work continuity.

mk_fb_data <- function() {
  set.seed(2026L)
  n <- 30L
  data.frame(
    y = rnorm(n),
    x = rnorm(n),
    g = factor(rep(1:5, length.out = n))
  )
}

# ---------------------------------------------------------------- #
# Alias relationship                                                #
# ---------------------------------------------------------------- #

test_that("fb is identical to flexybayes (ADR 0004 D1)", {
  expect_identical(fb, flexybayes)
})

test_that("fb() and flexybayes() byte-identical on the same call", {
  d <- mk_fb_data()

  via_fb <- fb(
    fixed = y ~ x,
    random = ~g,
    data = d,
    verbose = FALSE,
    return_code = TRUE
  )
  via_flexy <- flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    verbose = FALSE,
    return_code = TRUE
  )
  expect_identical(via_fb, via_flexy)
})

# ---------------------------------------------------------------- #
# Code generation (auto routes return_code = TRUE to brms/Stan --   #
# the v0.1 default fitting path was the native engine withdrawn     #
# entirely in 0.9.3, see NEWS.md; INLA has no return_code text form)#
# ---------------------------------------------------------------- #

test_that("fb() generates brms/Stan code for fixed-only model", {
  d <- mk_fb_data()
  code <- fb(
    fixed = y ~ x,
    data = d,
    verbose = FALSE,
    return_code = TRUE
  )
  expect_type(code, "character")
  expect_true(nzchar(code))
  expect_match(code, "sigma")
})

test_that("fb() generates RE code for asreml-style random intercept", {
  d <- mk_fb_data()
  code <- fb(
    fixed = y ~ x,
    random = ~g,
    data = d,
    verbose = FALSE,
    return_code = TRUE
  )
  expect_type(code, "character")
  expect_match(code, "sd_1")
  expect_match(code, "z_1")
})

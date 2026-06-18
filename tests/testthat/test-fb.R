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
  skip_if_no_greta()
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
# Code generation via greta backend (the v0.1 default and only      #
# fitting path; INLA emit and brms passthrough deferred to v0.2)    #
# ---------------------------------------------------------------- #

test_that("fb() generates greta code for fixed-only model", {
  skip_if_no_greta()
  d <- mk_fb_data()
  code <- fb(
    fixed = y ~ x,
    data = d,
    verbose = FALSE,
    return_code = TRUE
  )
  expect_type(code, "character")
  expect_true(nzchar(code))
  expect_match(code, "normal\\(")
})

test_that("fb() generates RE code for asreml-style random intercept", {
  skip_if_no_greta()
  d <- mk_fb_data()
  code <- fb(
    fixed = y ~ x,
    random = ~g,
    data = d,
    verbose = FALSE,
    return_code = TRUE
  )
  expect_type(code, "character")
  expect_match(code, "sigma_g")
  expect_match(code, "g_raw")
})

# DEF-200: a non-Gaussian INLA fit with a three-way random term reads the
# term's grouping variables by exact name.

test_that("DEF-200: covariate standardisation reads a:b:cc without partial matching", {
  set.seed(3)
  d <- expand.grid(a = factor(1:3), b = factor(1:2), cc = factor(1:2), rep = 1:4)
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rpois(nrow(d), 3)
  fb <- .fb_from_asreml_impl(
    fixed = y ~ x, random = ~ a:b:cc, data = d, family = "poisson"
  )
  withr::local_options(warnPartialMatchDollar = TRUE)
  expect_no_warning(protected <- .fb_inla_covariate_protected(fb))
  expect_setequal(protected, c("y", "a", "b", "cc"))
  expect_identical(.fb_inla_covariate_columns(fb, d), "x")
})

# DEF-227, found by the 0.16.0 release gate's long tier. A stub stands in
# for the brms fit, so the test runs without the long tier.

test_that("DEF-227: canonical draws of a brms fit do not pass n_samples", {
  seen <- NULL
  testthat::local_mocked_bindings(
    fb_as_draws_simple = function(fit, ...) {
      seen <<- names(list(...))
      list(b_Intercept = c(1, 2))
    }
  )
  fit <- structure(list(extras = list()), class = c("flexybayes_brms", "flexybayes"))
  flexyBayes:::.fb_canonical_draws(fit, n_samples = 20L, drop = FALSE)
  expect_false("n_samples" %in% seen)
  fit_inla <- structure(list(extras = list()), class = c("flexybayes_inla", "flexybayes"))
  flexyBayes:::.fb_canonical_draws(fit_inla, n_samples = 20L, drop = FALSE)
  expect_true("n_samples" %in% seen)
})

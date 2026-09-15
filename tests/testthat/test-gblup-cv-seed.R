# fb_gblup_cv(seed = ) fixes its folds with set.seed(). It used to leave the
# caller's random-number stream reset to that seed, so every draw the caller
# made afterwards depended on a seed they had given to one function.

.gblup_seed_inputs <- function() {
  set.seed(11L)
  n <- 40L
  m <- 200L
  z <- matrix(stats::rbinom(n * m, 2L, 0.3), n, m)
  zc <- scale(z, scale = FALSE)
  g <- tcrossprod(zc) / m + diag(n) * 1e-3
  y <- as.vector(t(chol(g)) %*% stats::rnorm(n)) + stats::rnorm(n)
  list(y = y, g = g)
}

test_that("fb_gblup_cv(seed = ) leaves the caller's random-number state as it found it", {
  inputs <- .gblup_seed_inputs()
  set.seed(2026L)
  before <- get(".Random.seed", envir = globalenv())

  cv <- fb_gblup_cv(inputs$y, inputs$g, folds = 4L, seed = 1L)

  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("the same seed still gives the same folds", {
  inputs <- .gblup_seed_inputs()
  first <- fb_gblup_cv(inputs$y, inputs$g, folds = 4L, seed = 1L)
  set.seed(5L)
  second <- fb_gblup_cv(inputs$y, inputs$g, folds = 4L, seed = 1L)
  expect_identical(first$predicted, second$predicted)
})

test_that("an unstarted random-number stream is left unstarted", {
  inputs <- .gblup_seed_inputs()
  withr::local_preserve_seed()
  rm(".Random.seed", envir = globalenv())

  cv <- fb_gblup_cv(inputs$y, inputs$g, folds = 4L, seed = 1L)

  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

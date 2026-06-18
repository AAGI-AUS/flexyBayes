# End-to-end gretaR backend: fit a real model through the out-of-process worker
# and confirm the C1 backend contract (canonical-named posterior::draws_array,
# correct recovery, draws feed triangulate()). Skipped unless gretaR is
# resolvable -- installed at the version floor, or a dev source tree pointed to
# by options(flexyBayes.gretaR_home = ...). It is skipped on CRAN and in CI
# where gretaR (torch) is not installed; it runs locally and is the functional
# proof of the activated backend.

.gretaR_resolvable <- function() {
  if (nzchar(getOption("flexyBayes.gretaR_home", ""))) {
    return(TRUE)
  }
  if (!nzchar(system.file(package = "gretaR"))) {
    return(FALSE)
  }
  v <- tryCatch(utils::packageVersion("gretaR"), error = function(e) NULL)
  !is.null(v) && v >= "0.4.0"
}

test_that("backend = 'gretaR' fits a Gaussian GLM with canonical draws", {
  skip_on_cran()
  skip_if_not(.gretaR_resolvable(),
              "gretaR not installed/resolvable; set flexyBayes.gretaR_home")

  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(1)
  n <- 120L
  x <- rnorm(n)
  y <- 1.5 - 0.8 * x + rnorm(n, 0, 0.7)
  d <- data.frame(y = y, x = x)

  fit <- flexybayes(y ~ x, data = d, family = "gaussian", backend = "gretaR",
                    n_samples = 400L, warmup = 400L, chains = 2L,
                    verbose = FALSE)

  expect_s3_class(fit, "flexybayes_gretaR")
  expect_s3_class(fit$draws, "draws_array")
  # Canonical names flow straight from model_from_arrays(names=) -- no relabel.
  expect_setequal(posterior::variables(fit$draws), c("(Intercept)", "x", "sigma"))

  s <- posterior::summarise_draws(fit$draws, mean = mean)
  m <- stats::setNames(s$mean, s$variable)
  expect_lt(abs(m[["(Intercept)"]] - 1.5), 0.3)
  expect_lt(abs(m[["x"]] - (-0.8)), 0.3)

  # fb_as_draws_simple contract: a named list of numeric vectors, one per param.
  ds <- fb_as_draws_simple(fit)
  expect_type(ds, "list")
  expect_setequal(names(ds), c("(Intercept)", "x", "sigma"))
  expect_true(all(vapply(ds, is.numeric, logical(1L))))
})

test_that("a gretaR fit is a working triangulate() member", {
  skip_on_cran()
  skip_if_not(.gretaR_resolvable(),
              "gretaR not installed/resolvable; set flexyBayes.gretaR_home")

  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(2)
  n <- 100L
  x <- rnorm(n)
  y <- 1.0 + 0.5 * x + rnorm(n, 0, 0.6)
  d <- data.frame(y = y, x = x)

  f1 <- flexybayes(y ~ x, data = d, family = "gaussian", backend = "gretaR",
                   n_samples = 300L, warmup = 300L, chains = 2L, verbose = FALSE)
  f2 <- flexybayes(y ~ x, data = d, family = "gaussian", backend = "gretaR",
                   n_samples = 300L, warmup = 300L, chains = 2L, verbose = FALSE)

  tr <- triangulate(f1, f2)
  expect_s3_class(tr, "triangulate_result")
})

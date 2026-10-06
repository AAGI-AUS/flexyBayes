# DEF-199: under prior_vc_sd, INLA and brms put the same prior on the
# standard deviation of a three-way random term.

.ilpm_data <- function() {
  set.seed(3)
  d <- expand.grid(a = factor(1:3), b = factor(1:2), cc = factor(1:2), rep = 1:4)
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(12)[interaction(d$a, d$b, d$cc)] +
    stats::rnorm(nrow(d))
  d
}

test_that("DEF-199: INLA and brms give a:b:cc the same prior under prior_vc_sd", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  d <- .ilpm_data()
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ x, random = ~ a:b + a:b:cc, data = d, prior_vc_sd = 2,
    backend = "inla", verbose = FALSE
  )))
  applied <- prior_summary(fit)$legacy_vc_applied

  fb <- .fb_from_asreml_impl(
    fixed = y ~ x, random = ~ a:b + a:b:cc, data = d, family = "gaussian"
  )
  specs <- .brms_legacy_specs(fb, prior_fixed_sd = 100, prior_vc_sd = 2)
  brms_sd <- Filter(function(s) identical(s$class, "sd"), specs)
  brms_prior <- stats::setNames(
    vapply(brms_sd, `[[`, "", "string"),
    paste0("sd_", vapply(brms_sd, `[[`, "", "group"))
  )

  expect_identical(brms_prior[["sd_a:b:cc"]], "lognormal(0, 2)")
  expect_identical(applied[["sd_a:b:cc"]], brms_prior[["sd_a:b:cc"]])
  expect_identical(applied[["sd_a:b"]], brms_prior[["sd_a:b"]])
  expect_length(prior_summary(fit)$filled_priors, 0L)
})

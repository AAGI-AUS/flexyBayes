# DEF-198: on brms, prior_vc_sd reaches an interaction random effect
# (a:b, a:b:cc) as lognormal(0, prior_vc_sd) on its standard deviation.

.lpm_data <- function() {
  set.seed(3)
  d <- expand.grid(a = factor(1:3), b = factor(1:2), cc = factor(1:2), rep = 1:4)
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(12)[interaction(d$a, d$b, d$cc)] + stats::rnorm(nrow(d))
  d
}

.lpm_sd_rows <- function(random, d) {
  fb <- .fb_from_asreml_impl(fixed = y ~ x, random = random, data = d, family = "gaussian")
  specs <- .brms_legacy_specs(fb, prior_fixed_sd = 100, prior_vc_sd = 2)
  Filter(function(s) identical(s$class, "sd"), specs)
}

test_that("DEF-198: a three-way grouping term gets the sd prior", {
  rows <- .lpm_sd_rows(~ a:b:cc, .lpm_data())
  expect_length(rows, 1L)
  expect_identical(rows[[1L]]$group, "a:b:cc")
  expect_identical(rows[[1L]]$string, "lognormal(0, 2)")
  expect_identical(rows[[1L]]$input, "scalar:prior_vc_sd")
})

test_that("DEF-198: a two-way grouping term gets the sd prior", {
  rows <- .lpm_sd_rows(~ a + a:b, .lpm_data())
  expect_identical(vapply(rows, `[[`, "", "group"), c("a", "a:b"))
})

test_that("DEF-198: brms accepts every prior row for the formula it is given", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- .lpm_data()
  fb <- .fb_from_asreml_impl(fixed = y ~ x, random = ~ a:b:cc, data = d, family = "gaussian")
  pr <- .priors_to_brms(NULL, fb, prior_fixed_sd = 100, prior_vc_sd = 2)
  vp <- suppressMessages(brms::validate_prior(pr, y ~ x + (1 | a:b:cc), data = d))
  row <- vp[vp$class == "sd" & vp$group == "a:b:cc" & vp$coef == "", , drop = FALSE]
  expect_identical(as.character(row$prior), "lognormal(0, 2)")
})

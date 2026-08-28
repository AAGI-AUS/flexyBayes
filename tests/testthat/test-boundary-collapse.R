# A variance component pinned at the boundary must say so, on any term
# type, not only on the spatial field.
#
# These are unit tests on constructed fit objects rather than live fits.
# That is deliberate: the trigger is a degenerate INLA mode, and three
# runs of the identical call on besag.met C1 returned an sd_gen upper
# bound of 0.00396, 0.0176 and 11.75. A live fixture would be a flaky
# test asserting a non-reproducible mode. The logic under test is a
# reading of the canonical variance-component table, so the table is
# what the test supplies.

make_vc <- function(components, q97.5, sigma_median = 10, estimate = NULL) {
  vc <- data.frame(
    component = components,
    estimate = if (is.null(estimate)) q97.5 / 2 else estimate,
    sd = rep(1, length(components)),
    q2.5 = rep(0, length(components)),
    q97.5 = q97.5,
    stringsAsFactors = FALSE
  )
  attr(vc, "posterior_median") <- c(sigma = sigma_median)
  vc
}

make_fit <- function(vc) {
  structure(
    list(extras = list(variance_comps = vc)),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}

test_that("a boundary-pinned iid component is reported", {
  fit <- make_fit(make_vc(
    c("sigma", "sd_gen"),
    q97.5 = c(17, 0.004),
    sigma_median = 15.6
  ))
  reasons <- .fb_boundary_collapse_reasons(fit)
  expect_length(reasons, 1L)
  expect_match(reasons, "sd_gen", fixed = TRUE)
  expect_match(reasons, "0.004", fixed = TRUE)
  expect_warning(.fb_warn_boundary_collapse(fit), "at the boundary")
})

test_that("a well-identified component is not reported", {
  fit <- make_fit(make_vc(
    c("sigma", "sd_gen"),
    q97.5 = c(3.15, 1.37),
    sigma_median = 3.15
  ))
  expect_length(.fb_boundary_collapse_reasons(fit), 0L)
  expect_silent(.fb_warn_boundary_collapse(fit))
})

test_that("sd_spatial is left to the spatial detector, not double-reported", {
  fit <- make_fit(make_vc(
    c("sigma", "sd_spatial"),
    q97.5 = c(15.5, 0.015),
    sigma_median = 15.5
  ))
  expect_length(.fb_boundary_collapse_reasons(fit), 0L)
})

test_that("several collapsed components are all named", {
  fit <- make_fit(make_vc(
    c("sigma", "sd_gen", "sd_block"),
    q97.5 = c(12, 0.001, 0.002),
    sigma_median = 12
  ))
  reasons <- .fb_boundary_collapse_reasons(fit)
  expect_length(reasons, 2L)
  expect_true(any(grepl("sd_gen", reasons, fixed = TRUE)))
  expect_true(any(grepl("sd_block", reasons, fixed = TRUE)))
})

test_that("the residual scale falls back to the sigma row", {
  vc <- make_vc(c("sigma", "sd_gen"), q97.5 = c(20, 0.01), sigma_median = 15)
  attr(vc, "posterior_median") <- NULL
  vc$estimate <- c(15, 0.005)
  expect_equal(.fb_residual_scale(vc), 15)
  expect_length(.fb_boundary_collapse_reasons(make_fit(vc)), 1L)
})

test_that("an absent or unusable table is not an error", {
  expect_length(.fb_boundary_collapse_reasons(make_fit(NULL)), 0L)
  expect_length(
    .fb_boundary_collapse_reasons(make_fit(data.frame())),
    0L
  )
  # a table with no residual scale to compare against
  vc <- make_vc(c("sd_gen"), q97.5 = 0.01, sigma_median = 10)
  attr(vc, "posterior_median") <- NULL
  vc$estimate <- 0.005
  expect_length(.fb_boundary_collapse_reasons(make_fit(vc)), 0L)
})

test_that("the warning can be silenced by option", {
  fit <- make_fit(make_vc(
    c("sigma", "sd_gen"),
    q97.5 = c(17, 0.004),
    sigma_median = 15.6
  ))
  withr::with_options(
    list(flexyBayes.silence_boundary_collapse_warning = TRUE),
    expect_silent(.fb_warn_boundary_collapse(fit))
  )
})

# The calibration itself, pinned against the measured sweep rather than
# against a number typed in here. The threshold started at 0.05 of the
# residual SD and fired on seven suite fixtures whose components are null
# by construction (`mk_inla_data()` is `y = rnorm(30)` with no group
# effect). `inst/validation/boundary_calibration.csv` holds 112 INLA fits
# with the group SD set to exactly zero, produced by
# `tools/calibrate_boundary_collapse.R`; the floor read off it is the
# number the threshold has to sit below. Reading it here rather than
# restating it means a wrong threshold cannot pass by also being written
# into the test.

.calib_path <- function() {
  local <- testthat::test_path("..", "..", "inst", "validation",
                               "boundary_calibration.csv")
  if (file.exists(local)) {
    return(local)
  }
  installed <- system.file("validation", "boundary_calibration.csv",
                           package = "flexyBayes")
  if (nzchar(installed)) installed else NA_character_
}

.null_floor <- function() {
  p <- .calib_path()
  if (is.na(p)) {
    return(NA_real_)
  }
  d <- utils::read.csv(p, stringsAsFactors = FALSE)
  d <- d[d$arm == "null" & is.finite(d$ratio), , drop = FALSE]
  if (!nrow(d)) NA_real_ else min(d$ratio)
}

test_that("the threshold sits below the measured null floor", {
  floor_ratio <- .null_floor()
  skip_if(is.na(floor_ratio), "calibration record not in this tree")
  expect_lt(.FB_COLLAPSE_FRACTION, floor_ratio)
  # and not so far below it that a real collapse escapes: the degenerate
  # mode measured on besag.met C1 is 0.00396 / 15.6
  expect_gt(.FB_COLLAPSE_FRACTION, 0.00396 / 15.6)
})

test_that("a null component at the measured floor does not warn", {
  floor_ratio <- .null_floor()
  skip_if(is.na(floor_ratio), "calibration record not in this tree")
  sigma <- 1.088
  fit <- make_fit(make_vc(
    c("sigma", "sd_gen"),
    q97.5 = c(sigma, floor_ratio * sigma),
    sigma_median = sigma
  ))
  expect_length(.fb_boundary_collapse_reasons(fit), 0L)
  expect_silent(.fb_warn_boundary_collapse(fit))
})

test_that("the measured degenerate mode still warns", {
  fit <- make_fit(make_vc(
    c("sigma", "sd_gen"),
    q97.5 = c(15.6, 0.00396),
    sigma_median = 15.6
  ))
  expect_length(.fb_boundary_collapse_reasons(fit), 1L)
  expect_warning(.fb_warn_boundary_collapse(fit), "at the boundary")
})

test_that("the calibration record is the shape the threshold argues from", {
  p <- .calib_path()
  skip_if(is.na(p), "calibration record not in this tree")
  d <- utils::read.csv(p, stringsAsFactors = FALSE)
  expect_true(all(c("arm", "n", "k", "seed", "ratio") %in% names(d)))
  expect_gte(sum(d$arm == "null"), 100L)
  # the positive control: a real component is orders of magnitude higher
  re <- d[d$arm == "real" & is.finite(d$ratio), ]
  expect_gt(stats::median(re$ratio), 10 * .null_floor())
})

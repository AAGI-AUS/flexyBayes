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


# =============================================================================
# Signed components are not scale components (audit 2026-09-11, F18).
#
# The display flag used a denylist: everything that was not a correlation and
# not `sigma` was treated as a variance component and tested against
# `upper < 0.01 * sigma`. Factor-analytic loadings are signed, so a negative
# loading -- a crossover interaction, which is a real and interpretable result
# -- has a negative upper credible bound, satisfies that inequality for
# arithmetic reasons, and was printed as `collapsed`. The audit found it on
# the shipped FA(1) and FA(2) summaries.
#
# The fit-time detector above has always used the allowlist, so the two
# surfaces disagreed. These tests pin the agreement in both directions: the
# loading must go unflagged, and a genuinely collapsed specific SD must still
# be flagged by both.
# =============================================================================

test_that("a negative factor-analytic loading is not marked collapsed", {
  components <- c("sigma", "sd_fa1", "fa1_env1", "fa1_env4", "sd_gen__env1")
  # fa1_env4 is the crossover loading: the whole interval sits below zero.
  upper <- c(17, 2.1, 1.8, -0.42, 1.2)
  medians <- c(sigma = 15.6)
  notes <- .fb_varcomp_notes(components, upper, medians)
  expect_equal(notes[[4]], "")
  expect_true(all(notes == ""))
})

test_that("a genuinely collapsed specific SD is still flagged", {
  # The falsifier. Without this the fix could simply have stopped flagging
  # anything and passed the test above.
  components <- c("sigma", "fa1_env4", "sd_gen__env3")
  upper <- c(17, -0.42, 0.004)
  medians <- c(sigma = 15.6)
  notes <- .fb_varcomp_notes(components, upper, medians)
  expect_equal(notes[[2]], "")
  expect_equal(notes[[3]], "collapsed")
})

test_that("the table and the fit-time diagnostic agree on the same components", {
  components <- c("sigma", "fa1_env4", "fa2_env5", "sd_gen__env3", "sd_fa1")
  upper <- c(17, -0.42, -1.03, 0.004, 2.1)
  vc <- make_vc(components, q97.5 = upper, sigma_median = 15.6)
  notes <- .fb_varcomp_notes(components, upper, c(sigma = 15.6))
  flagged_table <- components[notes == "collapsed"]
  flagged_fit <- .fb_boundary_collapse_reasons(make_fit(vc))
  expect_equal(flagged_table, "sd_gen__env3")
  expect_length(flagged_fit, 1L)
  expect_match(flagged_fit, "sd_gen__env3", fixed = TRUE)
})

test_that("component domains are classified from the canonical name", {
  expect_equal(.fb_component_domain("sigma"), "positive")
  expect_equal(.fb_component_domain("sd_gen"), "positive")
  expect_equal(.fb_component_domain("sd_spatial"), "positive")
  expect_equal(.fb_component_domain("sd_fa1"), "positive")
  expect_equal(.fb_component_domain("cor_g__Intercept__x"), "unit_interval")
  expect_equal(.fb_component_domain("rho_row"), "unit_interval")
  expect_equal(.fb_component_domain("fa1_env4"), "signed")
  expect_equal(.fb_component_domain("fa12_env4"), "signed")
  # A component the canonical namer did not produce fails safe: it is not
  # admitted to a test written for standard deviations.
  expect_equal(.fb_component_domain("shape"), "unknown")
  expect_equal(.fb_component_domain("shape"), "unknown")
})

test_that("an unknown-domain component is never flagged", {
  notes <- .fb_varcomp_notes(
    c("sigma", "shape"),
    c(17, 0.0001),
    c(sigma = 15.6)
  )
  expect_true(all(notes == ""))
})


# =============================================================================
# The Representation: line names a conditional translation (audit F06).
#
# `Representation: exact (model, not inference)` means the emit carries no
# structural approximation. A reader takes it for more: "this is the model
# ASReml would fit". For a factor-analytic term the two come apart -- the
# loadings are anchored to one at a level taken from factor-level order, which
# spans the factor-analytic family only where the leading loading block
# supports the anchoring. The audit found FA(1) and FA(2) summaries printing
# the bare `exact` badge on exactly that lowering.
# =============================================================================

make_fit_with_terms <- function(types) {
  structure(
    list(
      exactness = "exact",
      extras = list(
        parse_info = list(random = lapply(types, function(t) list(type = t)))
      )
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}

test_that("a factor-analytic fit does not report a bare exact representation", {
  fit <- make_fit_with_terms(c("fa_gxe"))
  label <- .repr_label_for_fit(fit, list())
  expect_false(identical(label, "exact (model, not inference)"))
  expect_match(label, "conditional", fixed = TRUE)
  expect_match(label, "TRANSLATION.md", fixed = TRUE)
})

test_that("a bare fa() term is treated the same as the interaction form", {
  expect_match(
    .repr_label_for_fit(make_fit_with_terms(c("fa")), list()),
    "conditional",
    fixed = TRUE
  )
})

test_that("a fit with no conditional term keeps the plain exact label", {
  # The falsifier: without this, always returning the caveat would pass above.
  for (types in list(c("ide"), c("nested", "at"), c("us_gxe"), character(0))) {
    expect_equal(
      .repr_label_for_fit(make_fit_with_terms(types), list()),
      "exact (model, not inference)"
    )
  }
})

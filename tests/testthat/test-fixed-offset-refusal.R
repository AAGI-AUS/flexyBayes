# =============================================================================
# offset() in the fixed formula (2026-09-13).
#
# terms() files an offset apart from the term labels, and the fixed-formula
# parser read only the labels, so `Reaction ~ Days + offset(Days)` fitted
# `Reaction ~ Days` on INLA and on brms with no message: a Days coefficient of
# 10.47 where lme4 with the offset returns 9.47. The refusal sits where the
# fixed formula is ingested, so every route that ingests one is exercised.
# =============================================================================

.off_sleepstudy <- function() {
  e <- new.env()
  utils::data("sleepstudy", package = "lme4", envir = e)
  e$sleepstudy
}

# A design the aggregated and streamed routes accept: every fixed term a
# factor, so without the refusal those routes would compress and fit rather
# than refuse for a reason of their own (a continuous fixed term forces one
# cell per row, which they refuse as continuous_cell_key_data_dependent).
.off_cells <- function(n = 120L, seed = 11L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(rep(sprintf("f%d", seq_len(4L)), length.out = n)),
    g = factor(rep(sprintf("g%02d", seq_len(12L)), each = 10L))
  )
  d$z <- stats::rnorm(n)
  d$y <- stats::rnorm(n, 10) + d$z
  d
}

.off_expect_refusal <- function(expr, offsets = "offset(Days)") {
  err <- tryCatch(expr, error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_fixed_offset_unsupported")
  # A route that fits instead of refusing returns no condition to read, and
  # each route is to fail on its own rather than end the test at the first.
  if (inherits(err, "condition")) {
    for (o in offsets) {
      expect_match(conditionMessage(err), o, fixed = TRUE)
    }
  }
}

test_that("every route that ingests a fixed formula refuses an offset", {
  skip_if_not_installed("lme4")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_sleepstudy()

  .off_expect_refusal(fb_from_asreml(
    fixed = Reaction ~ Days + offset(Days), random = ~ Subject, data = d
  ))
  .off_expect_refusal(
    fb_from_brms(Reaction ~ Days + offset(Days) + (1 | Subject), data = d)
  )
  .off_expect_refusal(
    fb_plan(Reaction ~ Days + offset(Days), random = ~ Subject, data = d)
  )
  for (backend in c("auto", "inla", "brms")) {
    .off_expect_refusal(flexybayes(
      Reaction ~ Days + offset(Days),
      random = ~ Subject,
      data = d,
      backend = backend,
      verbose = FALSE
    ))
  }

  cells <- .off_cells()
  .off_expect_refusal(
    flexybayes(
      y ~ f + offset(z),
      random = ~ g,
      data = cells,
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    ),
    offsets = "offset(z)"
  )
  .off_expect_refusal(
    suppressMessages(flexybayes_stream(
      y ~ f + offset(z),
      random = ~ g,
      source = cells,
      backend = "inla",
      fit = FALSE,
      verbose = FALSE
    )),
    offsets = "offset(z)"
  )
})

test_that("the refusal names every offset the formula carries", {
  skip_if_not_installed("lme4")
  d <- .off_sleepstudy()
  .off_expect_refusal(
    fb_from_asreml(
      fixed = Reaction ~ offset(Days) + Days + offset(log(Days + 1)),
      random = ~ Subject,
      data = d
    ),
    offsets = c("offset(Days)", "offset(log(Days + 1))")
  )
})

test_that("a column named offset is a covariate, not an offset", {
  skip_if_not_installed("lme4")
  d <- .off_sleepstudy()
  d$offset <- d$Days
  fb <- fb_from_asreml(fixed = Reaction ~ offset, random = ~ Subject, data = d)
  expect_identical(fb$fixed[[1L]]$var, "offset")
  expect_identical(fb$fixed[[1L]]$type, "continuous")
})

test_that("the response minus the offset fits the offset model", {
  # The spelling the refusal names, checked against lme4 with the offset.
  # With and without the offset the Days coefficients differ by exactly 1.
  # INLA's posterior mean and lme4's estimate on this balanced design
  # differed by 2.06e-6 in each of five identical calls (2026-09-13), so the
  # tolerance is 1e-4.
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_sleepstudy()
  fit <- flexybayes(
    Reaction - Days ~ Days,
    random = ~ Subject,
    data = d,
    backend = "inla",
    verbose = FALSE
  )
  ref <- lme4::fixef(
    lme4::lmer(Reaction ~ Days + offset(Days) + (1 | Subject), data = d)
  )
  expect_lt(abs(coef(fit)[["Days"]] - ref[["Days"]]), 1e-4)
})

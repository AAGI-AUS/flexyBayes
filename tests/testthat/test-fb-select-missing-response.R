# =============================================================================
# fb_select() on a base whose data carry a missing response.
#
# The unpenalised benchmark, the penalised optimiser and the comparison
# `benchmark = TRUE` adds all estimate a fold's covariance by REML from the
# response of every training row, and score every held-out row. A missing
# response stopped inside stats::optim() with the untyped "L-BFGS-B needs
# finite values of 'fn'". Each route now refuses by name before it runs, and
# `plan = TRUE` returns the same refusal.
# =============================================================================

fbmr_data <- function(missing_rows = c(4L, 19L)) {
  set.seed(20260919L)
  J <- 12L
  m <- 6L
  d <- data.frame(site = factor(rep(sprintf("S%02d", seq_len(J)), each = m)))
  d$clay <- rep(stats::rnorm(J), each = m)
  d$rain <- rep(stats::rnorm(J), each = m)
  d$plot_z <- stats::rnorm(nrow(d))
  d$yield <- 3 + 0.8 * d$clay + 0.2 * d$plot_z +
    rep(stats::rnorm(J, sd = 0.6), each = m) + stats::rnorm(nrow(d), sd = 0.7)
  d$yield[missing_rows] <- NA
  d
}

fbmr_base <- function(d) {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  flexybayes(yield ~ plot_z, random = ~ site, data = d, backend = "inla",
             verbose = FALSE)
}

fbmr_folds <- function() fb_resamples(unit = ~ site, v = 3L)

fbmr_class <- "flexybayes_refusal_selection_missing_response_not_carried"


test_that("the unpenalised benchmark refuses a missing response by name", {
  base <- fbmr_base(fbmr_data())
  err <- expect_error(
    fb_select(base, candidates = ~ clay + rain,
              workflow = "unpenalized_benchmark",
              resamples = fbmr_folds(), seed = 1L),
    class = fbmr_class
  )
  expect_identical(err$rows, c(4L, 19L))
  expect_match(conditionMessage(err), "(rows 4, 19)", fixed = TRUE)
  expect_match(conditionMessage(err), "workflow = \"unpenalized_benchmark\"",
               fixed = TRUE)
})


test_that("the penalised optimiser refuses it, and so does its plan", {
  skip_if_not_installed("glmnet")
  base <- fbmr_base(fbmr_data())
  pen <- fb_penalty(term(everything()) ~ ridge(scale = 1))
  for (plan in c(FALSE, TRUE)) {
    expect_error(
      fb_select(base, candidates = ~ clay + rain, penalty = pen,
                workflow = "penalized_optimizer",
                resamples = fbmr_folds(), seed = 1L, plan = plan),
      class = fbmr_class
    )
  }
})


test_that("benchmark = TRUE on the Bayesian route refuses before any fold refit", {
  base <- fbmr_base(fbmr_data())
  for (plan in c(FALSE, TRUE)) {
    err <- expect_error(
      fb_select(base, candidates = ~ clay + rain, benchmark = TRUE,
                resamples = fbmr_folds(), seed = 1L, plan = plan),
      class = fbmr_class
    )
  }
  expect_match(conditionMessage(err), "`benchmark = TRUE`", fixed = TRUE)
})


test_that("benchmark = TRUE with a shrinkage prior refuses before the reference fit", {
  # The budget stops the call before anything is sampled if the refusal
  # does not come first, so this block never fits the reference model.
  skip_if_not_installed("brms")
  skip_if_not_installed("projpred")
  base <- fbmr_base(fbmr_data())
  pen <- fb_penalty(term(everything()) ~ ridge(scale = 0.5))
  for (plan in c(FALSE, TRUE)) {
    expect_error(
      fb_select(base, candidates = ~ clay + rain, penalty = pen,
                benchmark = TRUE, resamples = fbmr_folds(), seed = 1L,
                budget = 0.001, plan = plan),
      class = fbmr_class
    )
  }
})


test_that("a base with every response observed is not refused", {
  base <- fbmr_base(fbmr_data(missing_rows = integer(0)))
  sel <- fb_select(base, candidates = ~ clay + rain,
                   workflow = "unpenalized_benchmark",
                   resamples = fbmr_folds(), seed = 1L)
  expect_s3_class(sel, "fb_select")
})


test_that("the refusal is registered", {
  expect_true(
    "selection_missing_response_not_carried" %in% fb_refusals()$reason_code
  )
})

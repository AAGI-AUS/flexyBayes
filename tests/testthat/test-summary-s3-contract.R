# =============================================================================
# summary() obeys the S3 contract: compute and return, do not print.
#
# summary.flexybayes(), summary.flexybayes_aggregated() and
# summary.flexybayes_inla() return visibly and print nothing.
# print.summary.flexybayes() does the printing, so `s <- summary(fit)` is
# silent and `print(summary(fit))` prints once.
#
# capture.output() auto-prints a visible value, so
# `capture.output(summary(fit))` still captures. expect_output() does not,
# so those sites print explicitly.
# =============================================================================

fb_s3_fixture <- function() {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  d <- agridat::john.alpha
  d$gen <- factor(d$gen)
  d$rep <- factor(d$rep)
  d$block <- factor(d$block)
  # `rep` has three levels, so sd_rep is weakly identified.
  withr::local_options(flexyBayes.silence_weak_identification_warning = TRUE)
  flexybayes(
    fixed = yield ~ gen,
    random = ~ rep + rep:block,
    data = d,
    backend = "inla",
    verbose = FALSE,
    prior = .fb_test_john_alpha_prior(d)
  )
}


test_that("assigning summary() prints nothing", {
  fit <- fb_s3_fixture()
  out <- utils::capture.output(s <- summary(fit))
  expect_identical(length(out), 0L)
  expect_s3_class(s, "summary.flexybayes")
})


test_that("summary() returns visibly, so the console still prints it", {
  fit <- fb_s3_fixture()
  # withVisible is the direct check. An invisible return would mean a bare
  # summary(fit) at the prompt printed nothing, which would be a regression
  # in the other direction.
  expect_true(withVisible(summary(fit))$visible)
  expect_gt(length(utils::capture.output(summary(fit))), 10L)
})


test_that("print(summary()) prints once, not twice", {
  fit <- fb_s3_fixture()
  once <- length(utils::capture.output(print(summary(fit))))
  auto <- length(utils::capture.output(summary(fit)))
  expect_identical(once, auto)
})


test_that("the summary object carries what a programmatic consumer needs", {
  fit <- fb_s3_fixture()
  s <- summary(fit)
  expect_true(all(
    c("fixed", "varcomp", "converge", "n_observed", "engine") %in% names(s)
  ))
  expect_gt(nrow(s$varcomp), 0L)
})

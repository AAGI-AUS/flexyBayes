# =============================================================================
# summary() obeys the S3 contract: compute and return, do not print.
#
# summary.flexybayes(), summary.flexybayes_aggregated() and
# summary.flexybayes_inla() used to call print() unconditionally and return
# invisibly. Three consequences, all of them wrong and none of them visible
# to a test that only checked the returned value:
#
#   * `s <- summary(fit)` emitted 67 lines of model summary;
#   * `print(summary(fit))` emitted 134, printing the whole thing twice;
#   * a programmatic consumer had no way to read the object quietly. The
#     selection substrate reads it once per fold, so a ten-fold run would
#     have written several hundred lines of summary into the console.
#
# print.summary.flexybayes() already existed. The printing was simply on the
# wrong side of the contract.
#
# The fix is behaviour-preserving where it matters, which is why the blast
# radius was two lines rather than the fifteen a grep suggests:
# capture.output() auto-prints a visible value, so every
# `capture.output(summary(fit))` site still captures. testthat's
# expect_output() does not, so those two sites now print explicitly.
# =============================================================================

fb_s3_fixture <- function() {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  d <- agridat::john.alpha
  d$gen <- factor(d$gen)
  d$rep <- factor(d$rep)
  d$block <- factor(d$block)
  flexybayes(
    fixed = yield ~ gen, random = ~ rep + rep:block,
    data = d, backend = "inla", verbose = FALSE,
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

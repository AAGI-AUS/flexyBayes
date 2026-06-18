# Fit-time convergence warning (.fb_warn_poor_convergence). Tested on
# synthetic fit objects so the behaviour is deterministic and fast: the
# helper's contract is a pure function of the attached convergence slot.

mk_fit <- function(psrf_pt = NULL, n_eff = NULL) {
  gelman <- if (is.null(psrf_pt)) {
    list(psrf = NULL)
  } else {
    list(psrf = matrix(
      c(psrf_pt, rep(1.0, length(psrf_pt))),
      ncol = 2L,
      dimnames = list(NULL, c("Point est.", "Upper C.I."))
    ))
  }
  fit <- list(extras = list(convergence = list(gelman = gelman, n_eff = n_eff)))
  class(fit) <- "flexybayes"
  fit
}

test_that("warns when any Rhat exceeds the threshold", {
  withr::local_options(flexyBayes.silence_convergence_warning = FALSE)
  expect_warning(
    flexyBayes:::.fb_warn_poor_convergence(mk_fit(c(1.01, 1.43))),
    "may not have converged"
  )
})

test_that("warns when effective sample size is below the floor", {
  withr::local_options(flexyBayes.silence_convergence_warning = FALSE)
  expect_warning(
    flexyBayes:::.fb_warn_poor_convergence(mk_fit(c(1.0), n_eff = c(5, 8))),
    "effective sample size"
  )
})

test_that("is quiet on a well-mixed fit", {
  withr::local_options(flexyBayes.silence_convergence_warning = FALSE)
  expect_no_warning(
    flexyBayes:::.fb_warn_poor_convergence(mk_fit(c(1.01, 1.02), n_eff = c(900, 950)))
  )
})

test_that("is a no-op for a deterministic (INLA / Laplace) fit with no psrf", {
  withr::local_options(flexyBayes.silence_convergence_warning = FALSE)
  expect_no_warning(flexyBayes:::.fb_warn_poor_convergence(mk_fit(NULL)))
})

test_that("is a no-op on non-fit objects (code / plan)", {
  withr::local_options(flexyBayes.silence_convergence_warning = FALSE)
  expect_no_warning(flexyBayes:::.fb_warn_poor_convergence(list(a = 1)))
})

test_that("respects the silence option", {
  withr::local_options(flexyBayes.silence_convergence_warning = TRUE)
  expect_no_warning(flexyBayes:::.fb_warn_poor_convergence(mk_fit(c(1.9))))
})

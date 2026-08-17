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

# The structured-term note. The August 2026 MET probe found the `us`
# branch of this warning diagnosing the wrong thing: it blamed
# non-identified factor-analytic loadings (brms parameterises us() by
# standard deviations and a Cholesky correlation, which is identified)
# and sent the reader to fb_structured_cov(), which abstains for us()
# terms. Reading it at face value would turn a genuine mixing failure
# into a labelling artefact.

mk_struct_fit <- function(type) {
  fit <- mk_fit(c(1.01, 1.43))
  fit$extras$parse_info <- list(random = list(list(type = type)))
  fit
}

test_that("the us() note names the confounding, not a rotation artefact", {
  withr::local_options(flexyBayes.silence_convergence_warning = FALSE)
  w <- tryCatch(
    flexyBayes:::.fb_warn_poor_convergence(mk_struct_fit("us_gxe")),
    warning = function(w) w
  )
  msg <- conditionMessage(w)
  expect_match(msg, "one observation per cell", fixed = TRUE)
  expect_match(msg, "confounded with the diagonal", fixed = TRUE)
  expect_match(msg, "replication within cell", fixed = TRUE)
  # It must not send the reader anywhere that abstains for us() terms.
  expect_false(grepl("fb_structured_cov", msg, fixed = TRUE))
  expect_false(grepl("rotation/sign", msg, fixed = TRUE))
})

test_that("the fa() note keeps its own, different diagnosis", {
  withr::local_options(flexyBayes.silence_convergence_warning = FALSE)
  w <- tryCatch(
    flexyBayes:::.fb_warn_poor_convergence(mk_struct_fit("fa_gxe")),
    warning = function(w) w
  )
  msg <- conditionMessage(w)
  expect_match(msg, "rotation/sign", fixed = TRUE)
  expect_match(msg, "fb_structured_cov()", fixed = TRUE)
  expect_false(grepl("one observation per cell", msg, fixed = TRUE))
})

test_that("an unstructured fit with no convergence problem stays quiet", {
  withr::local_options(flexyBayes.silence_convergence_warning = FALSE)
  fit <- mk_fit(c(1.00, 1.01), n_eff = c(900, 1200))
  fit$extras$parse_info <- list(random = list(list(type = "us_gxe")))
  expect_no_warning(flexyBayes:::.fb_warn_poor_convergence(fit))
})

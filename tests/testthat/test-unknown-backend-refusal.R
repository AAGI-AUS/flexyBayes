# =============================================================================
# Any backend name outside the two registered engines is unknown to every
# fitting entry point. `.check_known_backend_name()` refuses anything not
# in the active set, so made-up names stand in for any unregistered one.
# A new fitting route that skips the shared check fails here.
# =============================================================================

test_that("naming an unregistered backend raises unknown_backend at every known fitting entry point", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(1L)
  d <- data.frame(y = stats::rnorm(20L), x = stats::rnorm(20L))

  # Every unregistered name raises the same reason code, and both active
  # engines are named in the message.
  for (nm in c("not_a_registered_backend", "still_not_one")) {
    err <- tryCatch(
      suppressMessages(flexybayes(y ~ x, data = d, backend = nm)),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_unknown_backend_refusal")
    expect_identical(err$reason_code, "unknown_backend", label = nm)
    expect_match(conditionMessage(err), '"inla"', fixed = TRUE, label = nm)
    expect_match(conditionMessage(err), '"brms"', fixed = TRUE, label = nm)
  }
})

test_that("fb_dirichlet() rejects any method beyond the sole ml estimator", {
  # "ml" is the only method, so any other request fails match.arg()
  # before any fitting code runs.
  m <- matrix(stats::rgamma(60L, 2), 20L, 3L)
  m <- m / rowSums(m)
  err <- tryCatch(
    fb_dirichlet(m, method = "not_a_real_method"),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_false(inherits(err, "flexybayes_refusal"))
  expect_match(conditionMessage(err), "ml", fixed = TRUE)
})

test_that("the maximum-likelihood Dirichlet fit is unaffected", {
  set.seed(2L)
  m <- matrix(stats::rgamma(90L, 2), 30L, 3L)
  m <- m / rowSums(m)
  fit <- fb_dirichlet(m)
  expect_s3_class(fit, "fb_dirichlet_fit")
  expect_identical(fit$method, "ml")
  expect_equal(sum(fit$mean_composition), 1, tolerance = 1e-8)
})

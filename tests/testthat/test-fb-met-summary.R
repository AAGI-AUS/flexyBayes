# Breeder MET summary (G2, collapsed at 0.9.3). fb_met_summary() computed
# overall performance, stability, G x E BLUPs, and the environment
# correlation matrix from the *realised* factor-analytic effects the
# withdrawn native engine's fa(env, k):gen emit monitored (see NEWS.md,
# 0.9.3). Neither active engine emits an fa() term at all -- INLA and brms
# both refuse it as a structured-covariance term outside their
# representable set -- so the function now abstains unconditionally with
# `met_summary_not_available` regardless of what `fit` carries. These tests
# were rewritten from a computational-correctness suite (pinned arithmetic
# on a hand-built fit, a real MCMC-fitted end-to-end check) to a refusal
# suite, because there is no longer a code path that computes anything to
# check.

test_that("fb_met_summary() refuses unconditionally on a well-formed flexybayes fit", {
  fake <- structure(
    list(extras = list(parse_info = list(random = list(
      list(type = "fa_gxe", inner = "gen", outer = "env", k = 2L)
    )))),
    class = c("flexybayes", "list")
  )
  err <- tryCatch(fb_met_summary(fake), error = function(e) e)
  expect_s3_class(err, "error")
  expect_identical(err$reason_code, "met_summary_not_available")
  expect_match(conditionMessage(err), "not available", fixed = TRUE)
  expect_match(conditionMessage(err), "summary\\(\\)")
})

test_that("fb_met_summary() refuses unconditionally even without a factor-analytic term", {
  # Before 0.9.3 the absence of a fa() term raised a *different* refusal
  # ("no factor-analytic term"); now the function abstains before it ever
  # inspects `parse_info$random`, so the reason code is the same regardless
  # of what the fit carries.
  bad <- structure(
    list(extras = list(parse_info = list(random = list()))),
    class = c("flexybayes", "list")
  )
  err <- tryCatch(fb_met_summary(bad), error = function(e) e)
  expect_identical(err$reason_code, "met_summary_not_available")
})

test_that("fb_met_summary() refuses a non-flexybayes object with the class-check message", {
  expect_error(fb_met_summary(list()), "flexybayes")
  expect_error(fb_met_summary(1:5), "flexybayes")
})

test_that("fb_met_summary() refuses identically on an INLA-shaped and a brms-shaped fit", {
  # The refusal is unconditional on the class, not just the bare
  # `flexybayes` parent -- confirms an inla/brms fit does not accidentally
  # reach a different, class-specific method.
  for (cls in c("flexybayes_inla", "flexybayes_brms")) {
    fit <- structure(list(extras = list()), class = c(cls, "flexybayes", "list"))
    err <- tryCatch(fb_met_summary(fit), error = function(e) e)
    expect_identical(err$reason_code, "met_summary_not_available", label = cls)
  }
})

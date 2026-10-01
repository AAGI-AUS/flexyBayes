# brms, lme4 and glmmTMB export ranef() from nlme, and brms exports
# prior_summary() from rstantools. Attached after flexyBayes, those generics
# mask the package's own, so a bare call dispatches through them.

# Evaluates `expr` where only base R is visible, as in a user session, so a
# package method is found through its S3 registration alone.
.gm_eval <- function(expr, fit) {
  eval(expr, envir = list2env(list(fit = fit), parent = baseenv()))
}

.gm_fits <- function() {
  list(
    inla = structure(list(), class = c("flexybayes_inla", "flexybayes")),
    brms = structure(list(), class = c("flexybayes_brms", "flexybayes"))
  )
}

test_that("prior_summary() through rstantools' generic reaches the package method", {
  skip_if_not_installed("rstantools")
  for (fit in .gm_fits()) {
    masked <- .gm_eval(quote(rstantools::prior_summary(fit)), fit)
    expect_s3_class(masked, "prior_summary_flexybayes")
    expect_identical(masked, prior_summary(fit))
  }
})

test_that("ranef() through nlme's generic reaches the package method", {
  skip_if_not_installed("nlme")
  for (fit in .gm_fits()) {
    expect_identical(.gm_eval(quote(nlme::ranef(fit)), fit), ranef(fit))
  }
})
